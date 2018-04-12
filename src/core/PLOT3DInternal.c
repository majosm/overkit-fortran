// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

typedef enum {
  P3D_LITTLE_ENDIAN = 1,
  P3D_BIG_ENDIAN = 2
} p3d_endian;

typedef enum {
  P3D_STANDARD = 1,
  P3D_EXTENDED = 2
} p3d_format;

typedef long long p3d_offset;

static void ReadRecordWrapper(p3d_endian Endian, p3d_format Format, long long *Value, FILE *In);
static void WriteRecordWrapper(p3d_endian Endian, p3d_format Format, long long Value, FILE *Out);

static size_t fread_endian(p3d_endian Endian, void *Data, size_t ElementSize,
  size_t NumElements, FILE *In);
static size_t fwrite_endian(p3d_endian Endian, void *Data, size_t ElementSize,
  size_t NumElements, FILE *Out);

static void SwapEndian(void *Data, size_t ElementSize, size_t NumElements, void *SwappedData);

// Report the machine's endianness
p3d_endian P3DInternalMachineEndian() {

  unsigned char EndianTest[2] = {1, 0};

  if(*(short *)EndianTest == 1) {
    return P3D_LITTLE_ENDIAN;
  } else {
    return P3D_BIG_ENDIAN;
  }

}

// Figure out a PLOT3D grid file's format
void P3DInternalDetectGridFormat(char *FilePath, p3d_endian *Endian, p3d_format *Format,
  int *NumDims, int *NumGrids, int *WithIBlank, int *Error) {

  FILE *GridFile;
  char ErrorString[256];
  int i, m;
  int *NumPointsAll;
  int *NumPoints;
  size_t GridSize;

  *Error = 0;

  GridFile = fopen(FilePath, "rb");
  if (!GridFile) {
    sprintf(ErrorString, "%s %s %s", "ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  unsigned char InitialBytes[sizeof(long long)];
  fread(InitialBytes, sizeof(unsigned char), sizeof(long long), GridFile);

  // If little endian, the first byte will be the size of the record containing the grid count
  if (InitialBytes[0] != 0) {
    *Endian = P3D_LITTLE_ENDIAN;
  } else {
    *Endian = P3D_BIG_ENDIAN;
  }

  int RecordSizeInt;
  long long RecordSizeLongLong;
  memcpy(&RecordSizeInt, InitialBytes, sizeof(int));
  memcpy(&RecordSizeLongLong, InitialBytes, sizeof(long long));
  if (*Endian != P3DInternalMachineEndian()) {
    int RecordSizeIntCopy = RecordSizeInt;
    SwapEndian(&RecordSizeIntCopy, sizeof(int), 1, &RecordSizeInt);
    long long RecordSizeLongLongCopy = RecordSizeLongLong;
    SwapEndian(&RecordSizeLongLongCopy, sizeof(long long), 1, &RecordSizeLongLong);
  }

  if (RecordSizeLongLong == sizeof(int)) {
    *Format = P3D_EXTENDED;
  } else if (RecordSizeInt == sizeof(int)) {
    *Format = P3D_STANDARD;
  } else {
    sprintf(ErrorString, "ERROR: Failed to detect PLOT3D grid file format.\n");
    perror(ErrorString);
    *Error = 1;
    goto close_file;
  }

  rewind(GridFile);

  int RecordWrapperSize = *Format == P3D_STANDARD ? sizeof(int) : sizeof(long long);
  long long RecordWrapper;

  // Read the number of grids
  fseek(GridFile, RecordWrapperSize, SEEK_CUR);
  fread_endian(*Endian, NumGrids, sizeof(int), 1, GridFile);
  fseek(GridFile, RecordWrapperSize, SEEK_CUR);

  // The next record contains the number of points over all grids;
  // Use the first record wrapper to figure out the grid dimension
  ReadRecordWrapper(*Endian, *Format, &RecordWrapper, GridFile);
  *NumDims = RecordWrapper/(sizeof(int)*(*NumGrids));

  // Read the number of points over all grids
  NumPointsAll = (int *)malloc(sizeof(int)*(*NumGrids)*3);
  for (m = 0; m < *NumGrids; ++m) {
    NumPoints = NumPointsAll + 3*m;
    fread_endian(*Endian, NumPoints, sizeof(int), *NumDims, GridFile);
    for (i = *NumDims; i < 3; ++i) {
      NumPoints[i] = 1;
    }
  }
  fseek(GridFile, RecordWrapperSize, SEEK_CUR);

  // The next record contains grid 1
  // Use the first record wrapper to figure out whether the grid has IBlank or not
  ReadRecordWrapper(*Endian, *Format, &RecordWrapper, GridFile);
  GridSize = (size_t)NumPointsAll[0]*(size_t)NumPointsAll[1]*(size_t)NumPointsAll[2];
  *WithIBlank = RecordWrapper > sizeof(double)*(*NumDims)*GridSize;

  free(NumPointsAll);

  close_file:
    fclose(GridFile);

}

// Read the size of a grid in a PLOT3D grid file
void P3DInternalGetGridSize(char *FilePath, p3d_endian Endian, p3d_format Format, int NumDims,
  int GridID, int *NumPoints, int *Error) {

  FILE *GridFile;
  char ErrorString[256];
  int i;
  size_t Offset;

  *Error = 0;

  GridFile = fopen(FilePath, "rb");
  if (!GridFile) {
    sprintf(ErrorString, "%s %s %s", "ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  int RecordWrapperSize = Format == P3D_STANDARD ? sizeof(int) : sizeof(long long);

  // Skip past NumGrids record and part of NumPoints record that contains prior grids' sizes
  Offset = 0;
  Offset += RecordWrapperSize;
  Offset += sizeof(int);
  Offset += RecordWrapperSize;
  Offset += RecordWrapperSize;
  Offset += sizeof(int) * GridID * NumDims;
  fseek(GridFile, Offset, SEEK_SET);

  // Read the number of points for the current grid
  fread_endian(Endian, NumPoints, sizeof(int), NumDims, GridFile);
  for (i = NumDims; i < 3; ++i) {
    NumPoints[i] = 1;
  }

  fclose(GridFile);

}

// Report the offset of a grid record in a PLOT3D grid file
p3d_offset P3DInternalGetGridOffset(p3d_format Format, int NumDims, int NumGrids,
  int WithIBlank, int *NumPointsAll, int GridID) {

  p3d_offset Offset;
  int m;
  int *NumPoints;
  size_t GridSize;

  int RecordWrapperSize = Format == P3D_STANDARD ? sizeof(int) : sizeof(long long);

  Offset = 0;

  // Number of grids record
  Offset += RecordWrapperSize;
  Offset += sizeof(int);
  Offset += RecordWrapperSize;

  // Number of points over all grids record
  Offset += RecordWrapperSize;
  Offset += sizeof(int)*NumDims*NumGrids;
  Offset += RecordWrapperSize;

  // Previous grid records
  for (m = 0; m < GridID; ++m) {
    NumPoints = NumPointsAll + 3*m;
    GridSize = (size_t)NumPoints[0]*(size_t)NumPoints[1]*(size_t)NumPoints[2];
    Offset += RecordWrapperSize;
    Offset += sizeof(double)*NumDims*GridSize;
    if (WithIBlank) {
      Offset += sizeof(int)*GridSize;
    }
    Offset += RecordWrapperSize;
  }

  return Offset;

}

// Create a PLOT3D grid file, write the header, and pad the rest with zeros
void P3DInternalCreateGridFile(char *FilePath, p3d_endian Endian, p3d_format Format,
  int NumDims, int NumGrids, int WithIBlank, int *NumPointsAll, int *Error) {

  FILE *GridFile = NULL;
  char ErrorString[256];
  int m;
  long long RecordSize;
  int *NumPoints;
  size_t GridSize;
  const char Zeros[4096] = { 0 };
  size_t NumBytes;
  int WriteSize;

  *Error = 0;

  GridFile = fopen(FilePath, "wb");
  if (!GridFile) {
    sprintf(ErrorString,"%s %s %s","ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  int RecordWrapperSize = Format == P3D_STANDARD ? sizeof(int) : sizeof(long long);

  // Write the number of grids
  RecordSize = sizeof(int);
  WriteRecordWrapper(Endian, Format, RecordSize, GridFile);
  fwrite_endian(Endian, &NumGrids, sizeof(int), 1, GridFile);
  WriteRecordWrapper(Endian, Format, RecordSize, GridFile);

  // Write the number of points over all grids
  RecordSize = sizeof(int)*NumDims*NumGrids;
  WriteRecordWrapper(Endian, Format, RecordSize, GridFile);
  for (m = 0; m < NumGrids; ++m) {
    fwrite_endian(Endian, NumPointsAll+3*m, sizeof(int), NumDims, GridFile);
  }
  WriteRecordWrapper(Endian, Format, RecordSize, GridFile);

  // Pad the rest of the file with zeros
  for (m = 0; m < NumGrids; ++m) {

    NumPoints = NumPointsAll + 3*m;
    GridSize = (size_t)NumPoints[0]*(size_t)NumPoints[1]*(size_t)NumPoints[2];

    NumBytes = 0;
    NumBytes += RecordWrapperSize;
    NumBytes += sizeof(double)*NumDims*GridSize;
    if (WithIBlank) NumBytes += sizeof(int)*GridSize;
    NumBytes += RecordWrapperSize;

    while (NumBytes > 0) {
      WriteSize = min(4096, NumBytes);
      fwrite(Zeros, 1, WriteSize, GridFile);
      NumBytes -= WriteSize;
    }

  }

  fclose(GridFile);

}

// Read a PLOT3D grid from a grid file
void P3DInternalReadSingleGrid(char *FilePath, p3d_endian Endian, p3d_format Format,
  int NumDims, int WithIBlank, int *NumPoints, p3d_offset Offset, double *X, double *Y,
  double *Z, int *IBlank, int *Error) {

  FILE *GridFile = NULL;
  char ErrorString[256];
  size_t GridSize;

  *Error = 0;

  GridFile = fopen(FilePath, "rb+");
  if (!GridFile) {
    sprintf(ErrorString, "%s %s %s", "ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  int RecordWrapperSize = Format == P3D_STANDARD ? sizeof(int) : sizeof(long long);

  // Jump to grid location in file
  p3d_offset FileLoc = 0;
  while (FileLoc < Offset) {
    long int SeekAmount = min(LONG_MAX, Offset-FileLoc);
    fseek(GridFile, SeekAmount, SEEK_CUR);
    FileLoc += SeekAmount;
  }

  // Read the grid data
  GridSize = (size_t)NumPoints[0]*(size_t)NumPoints[1]*(size_t)NumPoints[2];
  fseek(GridFile, RecordWrapperSize, SEEK_CUR);
  fread_endian(Endian, X, sizeof(double), GridSize, GridFile);
  fread_endian(Endian, Y, sizeof(double), GridSize, GridFile);
  if (NumDims == 3) {
    fread_endian(Endian, Z, sizeof(double), GridSize, GridFile);
  }
  if (WithIBlank) {
    fread_endian(Endian, IBlank, sizeof(int), GridSize, GridFile);
  }

  fclose(GridFile);

}

// Can't easily pass null pointers from Fortran, so call these wrappers instead
void P3DInternalReadSingleGrid2D(char *FilePath, p3d_endian Endian, p3d_format Format,
  int *NumPoints, p3d_offset Offset, double *X, double *Y, int *Error) {
  P3DInternalReadSingleGrid(FilePath, Endian, Format, 2, 0, NumPoints, Offset, X, Y, NULL, NULL,
    Error);
}
void P3DInternalReadSingleGrid2DWithIBlank(char *FilePath, p3d_endian Endian, p3d_format Format,
  int *NumPoints, p3d_offset Offset, double *X, double *Y, int *IBlank, int *Error) {
  P3DInternalReadSingleGrid(FilePath, Endian, Format, 2, 1, NumPoints, Offset, X, Y, NULL, IBlank,
    Error);
}
void P3DInternalReadSingleGrid3D(char *FilePath, p3d_endian Endian, p3d_format Format,
  int *NumPoints, p3d_offset Offset, double *X, double *Y, double *Z, int *Error) {
  P3DInternalReadSingleGrid(FilePath, Endian, Format, 3, 0, NumPoints, Offset, X, Y, Z, NULL,
    Error);
}
void P3DInternalReadSingleGrid3DWithIBlank(char *FilePath, p3d_endian Endian, p3d_format Format,
  int *NumPoints, p3d_offset Offset, double *X, double *Y, double *Z, int *IBlank,
  int *Error) {
  P3DInternalReadSingleGrid(FilePath, Endian, Format, 3, 1, NumPoints, Offset, X, Y, Z, IBlank,
    Error);
}

// Write a PLOT3D grid to a grid file
void P3DInternalWriteSingleGrid(char *FilePath, p3d_endian Endian, p3d_format Format,
  int NumDims, int WithIBlank, int *NumPoints, p3d_offset Offset, double *X, double *Y,
  double *Z, int *IBlank, int *Error) {

  FILE *GridFile = NULL;
  char ErrorString[256];
  size_t GridSize;
  long long RecordSize;

  *Error = 0;

  GridFile = fopen(FilePath, "rb+");
  if (!GridFile) {
    sprintf(ErrorString, "%s %s %s", "ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  // Jump to grid location in file
  p3d_offset FileLoc = 0;
  while (FileLoc < Offset) {
    long int SeekAmount = min(LONG_MAX, Offset-FileLoc);
    fseek(GridFile, SeekAmount, SEEK_CUR);
    FileLoc += SeekAmount;
  }

  // Write the grid data
  GridSize = (size_t)NumPoints[0]*(size_t)NumPoints[1]*(size_t)NumPoints[2];
  RecordSize = sizeof(double)*NumDims*GridSize;
  if (WithIBlank == 1) RecordSize += sizeof(int)*GridSize;
  WriteRecordWrapper(Endian, Format, RecordSize, GridFile);
  fwrite_endian(Endian, X, sizeof(double), GridSize, GridFile);
  fwrite_endian(Endian, Y, sizeof(double), GridSize, GridFile);
  if (NumDims == 3) {
    fwrite_endian(Endian, Z, sizeof(double), GridSize, GridFile);
  }
  if (WithIBlank == 1) {
    fwrite_endian(Endian, IBlank, sizeof(int), GridSize, GridFile);
  }
  WriteRecordWrapper(Endian, Format, RecordSize, GridFile);

  fclose(GridFile);

}

// Can't easily pass null pointers from Fortran, so call these wrappers instead
void P3DInternalWriteSingleGrid2D(char *FilePath, p3d_endian Endian, p3d_format Format,
  int *NumPoints, p3d_offset Offset, double *X, double *Y, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, Endian, Format, 2, 0, NumPoints, Offset, X, Y, NULL, NULL,
    Error);
}
void P3DInternalWriteSingleGrid2DWithIBlank(char *FilePath, p3d_endian Endian,
  p3d_format Format, int *NumPoints, p3d_offset Offset, double *X, double *Y, int *IBlank,
  int *Error) {
  P3DInternalWriteSingleGrid(FilePath, Endian, Format, 2, 1, NumPoints, Offset, X, Y, NULL, IBlank,
    Error);
}
void P3DInternalWriteSingleGrid3D(char *FilePath, p3d_endian Endian, p3d_format Format,
  int *NumPoints, p3d_offset Offset, double *X, double *Y, double *Z, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, Endian, Format, 3, 0, NumPoints, Offset, X, Y, Z, NULL,
    Error);
}
void P3DInternalWriteSingleGrid3DWithIBlank(char *FilePath, p3d_endian Endian,
  p3d_format Format, int *NumPoints, p3d_offset Offset, double *X, double *Y, double *Z,
  int *IBlank, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, Endian, Format, 3, 1, NumPoints, Offset, X, Y, Z, IBlank,
    Error);
}

static void ReadRecordWrapper(p3d_endian Endian, p3d_format Format, long long *Value, FILE *In) {

  if (Format == P3D_STANDARD) {
    int ValueInt;
    fread_endian(Endian, &ValueInt, sizeof(int), 1, In);
    *Value = ValueInt;
  } else {
    fread_endian(Endian, Value, sizeof(long long), 1, In);
  }

}

static void WriteRecordWrapper(p3d_endian Endian, p3d_format Format, long long Value, FILE *Out) {

  if (Format == P3D_STANDARD) {
    int ValueInt = (int)Value;
    fwrite_endian(Endian, &ValueInt, sizeof(int), 1, Out);
  } else {
    fwrite_endian(Endian, &Value, sizeof(long long), 1, Out);
  }

}

// Read data, possibly with an endianness that is different from the machine's
static size_t fread_endian(p3d_endian Endian, void *Data, size_t ElementSize,
  size_t NumElements, FILE *In) {

  int DifferentEndian;
  void *DataToRead;
  size_t Result;

  DifferentEndian = Endian != P3DInternalMachineEndian();

  if (DifferentEndian) {
    DataToRead = malloc(ElementSize * NumElements);
  } else {
    DataToRead = Data;
  }

  Result = fread(DataToRead, ElementSize, NumElements, In);

  if (DifferentEndian) {
    SwapEndian(DataToRead, ElementSize, NumElements, Data);
    free(DataToRead);
  }

  return Result;

}

// Write data, possibly with an endianness that is different from the machine's
static size_t fwrite_endian(p3d_endian Endian, void *Data, size_t ElementSize,
  size_t NumElements, FILE *Out) {

  int DifferentEndian;
  void *DataToWrite;
  size_t Result;

  DifferentEndian = Endian != P3DInternalMachineEndian();

  if (DifferentEndian) {
    DataToWrite = malloc(ElementSize * NumElements);
    SwapEndian(Data, ElementSize, NumElements, DataToWrite);
  } else {
    DataToWrite = Data;
  }

  Result = fwrite(DataToWrite, ElementSize, NumElements, Out);

  if (DifferentEndian) {
    free(DataToWrite);
  }

  return Result;

}

static void SwapEndian(void *Data, size_t ElementSize, size_t NumElements, void *SwappedData) {

  int i, j;
  char *A, *B;

  for (i = 0; i < NumElements; ++i) {
    A = (char *)Data + i*ElementSize;
    B = (char *)SwappedData + i*ElementSize;
    for (j = 0; j < ElementSize; ++j) {
      B[j] = A[ElementSize-j-1];
    }
  }

}

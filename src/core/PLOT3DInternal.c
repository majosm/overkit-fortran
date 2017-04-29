// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <stdio.h>
#include <stdlib.h>

typedef enum {
  P3D_LITTLE_ENDIAN = 1,
  P3D_BIG_ENDIAN = 2
} p3d_endian_type;

typedef int p3d_size_type;
typedef long long int p3d_offset_type;

static p3d_size_type fread_endian(p3d_endian_type Endian, void *Data, p3d_size_type ElementSize,
  p3d_size_type NumElements, FILE *In);
static p3d_size_type fwrite_endian(p3d_endian_type Endian, void *Data, p3d_size_type ElementSize, 
  p3d_size_type NumElements, FILE *Out);

// Report the machine's endianness
p3d_endian_type P3DInternalMachineEndian() {

  unsigned char EndianTest[2] = {1, 0};

  if(*(short *)EndianTest == 1) {
    return P3D_LITTLE_ENDIAN;
  } else {
    return P3D_BIG_ENDIAN;
  }

}

// Figure out a PLOT3D grid file's format
void P3DInternalDetectGridFormat(char *FilePath, int *NumDims, int *NumGrids, int *WithIBlank,
  p3d_endian_type *Endian, int *Error) {

  FILE *GridFile;
  char ErrorString[256];
  int i, m;
  p3d_size_type Record;
  int *NumPointsAll;
  int *NumPoints;
  p3d_size_type GridSize;

  *Error = 0;

  GridFile = fopen(FilePath, "rb");
  if (!GridFile) {
    sprintf(ErrorString, "%s %s %s", "ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  // Determine the endianness by checking the first record indicator, which should
  // have a value of sizeof(int)
  *Endian = P3D_LITTLE_ENDIAN;
  fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
  if (Record != sizeof(int)) *Endian = P3D_BIG_ENDIAN;

  rewind(GridFile);

  // Read the number of grids
  fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
  fread_endian(*Endian, NumGrids, sizeof(int), 1, GridFile);
  fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

  if (*NumGrids != 0) {

    // The next record contains the number of points over all grids;
    // Use the first record indicator to figure out the grid dimension
    fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
    *NumDims = Record/(sizeof(int)*(*NumGrids));

    // Read the number of points over all grids
    NumPointsAll = (int *)malloc(sizeof(int)*(*NumGrids)*3);
    for (m = 0; m < *NumGrids; ++m) {
      NumPoints = NumPointsAll + 3*m;
      fread_endian(*Endian, NumPoints, sizeof(int), *NumDims, GridFile);
      for (i = *NumDims; i < 3; ++i) {
        NumPoints[i] = 1;
      }
    }
    fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

    // The next record contains grid 1
    // Use the first record indicator to figure out whether the grid has IBlank or not
    fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
    GridSize = NumPointsAll[0]*NumPointsAll[1]*NumPointsAll[2];
    *WithIBlank = Record > sizeof(double)*(*NumDims)*GridSize;

    free(NumPointsAll);

  } else {

    // Degenerate case, empty grid file -- just pick some values for the remaining
    // format parameters
    *NumDims = 3;
    *WithIBlank = 0;
    *Endian = P3DInternalMachineEndian();

  }

  fclose(GridFile);

}

// Read the size of a grid in a PLOT3D grid file
void P3DInternalGetGridSize(char *FilePath, int NumDims, p3d_endian_type Endian, int GridID,
  int *NumPoints, int *Error) {

  FILE *GridFile;
  char ErrorString[256];
  int i;
  p3d_offset_type Offset;

  *Error = 0;

  GridFile = fopen(FilePath, "rb");
  if (!GridFile) {
    sprintf(ErrorString, "%s %s %s", "ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  // Skip past NumGrids record and part of NumPoints record that contains prior grids' sizes
  Offset = 2*sizeof(p3d_size_type) + sizeof(int);
  Offset += sizeof(p3d_size_type) + sizeof(int) * GridID * NumDims;
  fseek(GridFile, (long int)Offset, SEEK_SET);

  // Read the number of points for the current grid
  fread_endian(Endian, NumPoints, sizeof(int), NumDims, GridFile);
  for (i = NumDims; i < 3; ++i) {
    NumPoints[i] = 1;
  }

  fclose(GridFile);

}

// Report the offset of a grid record in a PLOT3D grid file
p3d_offset_type P3DInternalGetGridOffset(int NumDims, int NumGrids, int *NumPointsAll,
  int WithIBlank, int GridID) {

  p3d_offset_type Offset;
  int m;
  int *NumPoints;
  p3d_size_type GridSize;

  Offset = 0;

  // Number of grids record
  Offset += sizeof(int) + 2*sizeof(p3d_size_type);

  // Number of points over all grids record
  Offset += sizeof(int)*NumDims*NumGrids + 2*sizeof(p3d_size_type);

  // Previous grid records
  for (m = 0; m < GridID; ++m) {
    NumPoints = NumPointsAll + 3*m;
    GridSize = NumPoints[0]*NumPoints[1]*NumPoints[2];
    Offset += sizeof(double)*NumDims*GridSize;
    if (WithIBlank) {
      Offset += sizeof(int)*GridSize;
    }
    Offset += 2*sizeof(p3d_size_type);
  }

  return Offset;

}

// Create a PLOT3D grid file, write the header, and pad the rest with zeros
void P3DInternalCreateGridFile(char *FilePath, int NumDims, int NumGrids, int *NumPointsAll,
  int WithIBlank, int Endian, int *Error) {

  FILE *GridFile = NULL;
  char ErrorString[256];
  p3d_size_type Record;
  int m;
  int *NumPoints;
  p3d_size_type GridSize;
  const char Zeros[4096] = { 0 };
  p3d_size_type NumBytes;
  int WriteSize;

  *Error = 0;

  GridFile = fopen(FilePath, "wb");
  if (!GridFile) {
    sprintf(ErrorString,"%s %s %s","ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  // Write the number of grids
  Record = sizeof(int);
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
  fwrite_endian(Endian, &NumGrids, sizeof(int), 1, GridFile);
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

  // Write the number of points over all grids
  Record = sizeof(int)*NumDims*NumGrids;
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
  for (m = 0; m < NumGrids; ++m) {
    fwrite_endian(Endian, NumPointsAll+3*m, sizeof(int), NumDims, GridFile);
  }
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

  // Pad the rest of the file with zeros
  for (m = 0; m < NumGrids; ++m) {

    NumPoints = NumPointsAll + 3*m;
    GridSize = NumPoints[0]*NumPoints[1]*NumPoints[2];
    NumBytes = sizeof(double)*NumDims*GridSize;
    if (WithIBlank) {
      NumBytes += sizeof(int)*GridSize;
    }
    NumBytes += 2*sizeof(p3d_size_type);

    while (NumBytes > 0) {
      WriteSize = NumBytes > 4096 ? 4096 : NumBytes;
      fwrite(Zeros, 1, WriteSize, GridFile);
      NumBytes -= WriteSize;
    }

  }

  fclose(GridFile);

}

// Read a PLOT3D grid from a grid file
void P3DInternalReadSingleGrid(char *FilePath, int NumDims, int *NumPoints, int WithIBlank,
  p3d_endian_type Endian, p3d_offset_type Offset, double *X, double *Y, double *Z, int *IBlank,
  int *Error) {

  FILE *GridFile = NULL;
  char ErrorString[256];
  p3d_size_type GridSize;

  *Error = 0;

  GridFile = fopen(FilePath, "rb+");
  if (!GridFile) {
    sprintf(ErrorString, "%s %s %s", "ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  // Jump to grid location in file
  // TODO: Figure out how to deal with offsets greater than long int max size
  fseek(GridFile, (long int)(Offset + sizeof(p3d_size_type)), SEEK_SET);

  // Read the grid data
  GridSize = NumPoints[0]*NumPoints[1]*NumPoints[2];
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
void P3DInternalReadSingleGrid2D(char *FilePath, int *NumPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, int *Error) {
  P3DInternalReadSingleGrid(FilePath, 2, NumPoints, 0, Endian, Offset, X, Y, NULL, NULL, Error);
}
void P3DInternalReadSingleGrid2DWithIBlank(char *FilePath, int *NumPoints,
  p3d_endian_type Endian, p3d_offset_type Offset, double *X, double *Y, int *IBlank, int *Error) {
  P3DInternalReadSingleGrid(FilePath, 2, NumPoints, 1, Endian, Offset, X, Y, NULL, IBlank, Error);
}
void P3DInternalReadSingleGrid3D(char *FilePath, int *NumPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, double *Z, int *Error) {
  P3DInternalReadSingleGrid(FilePath, 3, NumPoints, 0, Endian, Offset, X, Y, Z, NULL, Error);
}
void P3DInternalReadSingleGrid3DWithIBlank(char *FilePath, int *NumPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, double *Z, int *IBlank, int *Error) {
  P3DInternalReadSingleGrid(FilePath, 3, NumPoints, 1, Endian, Offset, X, Y, Z, IBlank, Error);
}

// Write a PLOT3D grid to a grid file
void P3DInternalWriteSingleGrid(char *FilePath, int NumDims, int *NumPoints, int WithIBlank,
  p3d_endian_type Endian, p3d_offset_type Offset, double *X, double *Y, double *Z, int *IBlank,
  int *Error) {

  FILE *GridFile = NULL;
  char ErrorString[256];
  p3d_size_type GridSize;
  p3d_size_type Record;

  *Error = 0;

  GridFile = fopen(FilePath, "rb+");
  if (!GridFile) {
    sprintf(ErrorString, "%s %s %s", "ERROR: Unable to open file ", FilePath, ".\n");
    perror(ErrorString);
    *Error = 1;
    return;
  }

  // Jump to grid location in file
  // TODO: Figure out how to deal with offsets greater than long int max size
  fseek(GridFile, (long int)Offset, SEEK_SET);

  // Write the grid data
  GridSize = NumPoints[0]*NumPoints[1]*NumPoints[2];
  Record = sizeof(double)*NumDims*GridSize + sizeof(int)*GridSize;
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
  fwrite_endian(Endian, X, sizeof(double), GridSize, GridFile);
  fwrite_endian(Endian, Y, sizeof(double), GridSize, GridFile);
  if (NumDims == 3) {
    fwrite_endian(Endian, Z, sizeof(double), GridSize, GridFile);
  }
  if (WithIBlank == 1) {
    fwrite_endian(Endian, IBlank, sizeof(int), GridSize, GridFile);
  }
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

  fclose(GridFile);

}

// Can't easily pass null pointers from Fortran, so call these wrappers instead
void P3DInternalWriteSingleGrid2D(char *FilePath, int *NumPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, 2, NumPoints, 0, Endian, Offset, X, Y, NULL, NULL, Error);
}
void P3DInternalWriteSingleGrid2DWithIBlank(char *FilePath, int *NumPoints,
  p3d_endian_type Endian, p3d_offset_type Offset, double *X, double *Y, int *IBlank, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, 2, NumPoints, 1, Endian, Offset, X, Y, NULL, IBlank, Error);
}
void P3DInternalWriteSingleGrid3D(char *FilePath, int *NumPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, double *Z, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, 3, NumPoints, 0, Endian, Offset, X, Y, Z, NULL, Error);
}
void P3DInternalWriteSingleGrid3DWithIBlank(char *FilePath, int *NumPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, double *Z, int *IBlank, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, 3, NumPoints, 1, Endian, Offset, X, Y, Z, IBlank, Error);
}

static void SwapEndian(void *Data, p3d_size_type ElementSize, p3d_size_type NumElements,
  void *SwappedData) {

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

// Read data, possibly with an endianness that is different from the machine's
static p3d_size_type fread_endian(p3d_endian_type Endian, void *Data, p3d_size_type ElementSize,
  p3d_size_type NumElements, FILE *In) {

  int DifferentEndian;
  void *DataToRead;
  p3d_size_type Result;

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
static p3d_size_type fwrite_endian(p3d_endian_type Endian, void *Data, p3d_size_type ElementSize, 
  p3d_size_type NumElements, FILE *Out) {

  int DifferentEndian;
  void *DataToWrite;
  p3d_size_type Result;

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

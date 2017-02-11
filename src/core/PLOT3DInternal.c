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
  p3d_size_type nElements, FILE *In);
static p3d_size_type fwrite_endian(p3d_endian_type Endian, void *Data, p3d_size_type ElementSize, 
  p3d_size_type nElements, FILE *Out);

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
void P3DInternalDetectGridFormat(char *FilePath, int *nDims, int *nGrids, int *WithIBlank,
  p3d_endian_type *Endian, int *Error) {

  FILE *GridFile;
  char ErrorString[256];
  int i, m;
  p3d_size_type Record;
  int *nPointsAll;
  int *nPoints;
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
  fread_endian(*Endian, nGrids, sizeof(int), 1, GridFile);
  fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

  if (*nGrids != 0) {

    // The next record contains the number of points over all grids;
    // Use the first record indicator to figure out the grid dimension
    fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
    *nDims = Record/(sizeof(int)*(*nGrids));

    // Read the number of points over all grids
    nPointsAll = (int *)malloc(sizeof(int)*(*nGrids)*3);
    for (m = 0; m < *nGrids; ++m) {
      nPoints = nPointsAll + 3*m;
      fread_endian(*Endian, nPoints, sizeof(int), *nDims, GridFile);
      for (i = *nDims; i < 3; ++i) {
        nPoints[i] = 1;
      }
    }
    fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

    // The next record contains grid 1
    // Use the first record indicator to figure out whether the grid has IBlank or not
    fread_endian(*Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
    GridSize = nPointsAll[0]*nPointsAll[1]*nPointsAll[2];
    *WithIBlank = Record > sizeof(double)*(*nDims)*GridSize;

    free(nPointsAll);

  } else {

    // Degenerate case, empty grid file -- just pick some values for the remaining
    // format parameters
    *nDims = 3;
    *WithIBlank = 0;
    *Endian = P3DInternalMachineEndian();

  }

  fclose(GridFile);

}

// Read the size of a grid in a PLOT3D grid file
void P3DInternalGetGridSize(char *FilePath, int nDims, p3d_endian_type Endian, int GridID,
  int *nPoints, int *Error) {

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

  // Skip past nGrids record and part of nPoints record that contains prior grids' sizes
  Offset = 2*sizeof(p3d_size_type) + sizeof(int);
  Offset += sizeof(p3d_size_type) + sizeof(int) * GridID * nDims;
  fseek(GridFile, (long int)Offset, SEEK_SET);

  // Read the number of points for the current grid
  fread_endian(Endian, nPoints, sizeof(int), nDims, GridFile);
  for (i = nDims; i < 3; ++i) {
    nPoints[i] = 1;
  }

  fclose(GridFile);

}

// Report the offset of a grid record in a PLOT3D grid file
p3d_offset_type P3DInternalGetGridOffset(int nDims, int nGrids, int *nPointsAll, int WithIBlank,
  int GridID) {

  p3d_offset_type Offset;
  int m;
  int *nPoints;
  p3d_size_type GridSize;

  Offset = 0;

  // Number of grids record
  Offset += sizeof(int) + 2*sizeof(p3d_size_type);

  // Number of points over all grids record
  Offset += sizeof(int)*nDims*nGrids + 2*sizeof(p3d_size_type);

  // Previous grid records
  for (m = 0; m < GridID; ++m) {
    nPoints = nPointsAll + 3*m;
    GridSize = nPoints[0]*nPoints[1]*nPoints[2];
    Offset += sizeof(double)*nDims*GridSize;
    if (WithIBlank) {
      Offset += sizeof(int)*GridSize;
    }
    Offset += 2*sizeof(p3d_size_type);
  }

  return Offset;

}

// Create a PLOT3D grid file, write the header, and pad the rest with zeros
void P3DInternalCreateGridFile(char *FilePath, int nDims, int nGrids, int *nPointsAll,
  int WithIBlank, int Endian, int *Error) {

  FILE *GridFile = NULL;
  char ErrorString[256];
  p3d_size_type Record;
  int m;
  int *nPoints;
  p3d_size_type GridSize;
  const char Zeros[4096] = { 0 };
  p3d_size_type nBytes;
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
  fwrite_endian(Endian, &nGrids, sizeof(int), 1, GridFile);
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

  // Write the number of points over all grids
  Record = sizeof(int)*nDims*nGrids;
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
  for (m = 0; m < nGrids; ++m) {
    fwrite_endian(Endian, nPointsAll+3*m, sizeof(int), nDims, GridFile);
  }
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

  // Pad the rest of the file with zeros
  for (m = 0; m < nGrids; ++m) {

    nPoints = nPointsAll + 3*m;
    GridSize = nPoints[0]*nPoints[1]*nPoints[2];
    nBytes = sizeof(double)*nDims*GridSize;
    if (WithIBlank) {
      nBytes += sizeof(int)*GridSize;
    }
    nBytes += 2*sizeof(p3d_size_type);

    while (nBytes > 0) {
      WriteSize = nBytes > 4096 ? 4096 : nBytes;
      fwrite(Zeros, 1, WriteSize, GridFile);
      nBytes -= WriteSize;
    }

  }

  fclose(GridFile);

}

// Read a PLOT3D grid from a grid file
void P3DInternalReadSingleGrid(char *FilePath, int nDims, int *nPoints, int WithIBlank,
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
  GridSize = nPoints[0]*nPoints[1]*nPoints[2];
  fread_endian(Endian, X, sizeof(double), GridSize, GridFile);
  fread_endian(Endian, Y, sizeof(double), GridSize, GridFile);
  if (nDims == 3) {
    fread_endian(Endian, Z, sizeof(double), GridSize, GridFile);
  }
  if (WithIBlank) {
    fread_endian(Endian, IBlank, sizeof(int), GridSize, GridFile);
  }

  fclose(GridFile);

}

// Can't easily pass null pointers from Fortran, so call these wrappers instead
void P3DInternalReadSingleGrid2D(char *FilePath, int *nPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, int *Error) {
  P3DInternalReadSingleGrid(FilePath, 2, nPoints, 0, Endian, Offset, X, Y, NULL, NULL, Error);
}
void P3DInternalReadSingleGrid2DWithIBlank(char *FilePath, int *nPoints,
  p3d_endian_type Endian, p3d_offset_type Offset, double *X, double *Y, int *IBlank, int *Error) {
  P3DInternalReadSingleGrid(FilePath, 2, nPoints, 1, Endian, Offset, X, Y, NULL, IBlank, Error);
}
void P3DInternalReadSingleGrid3D(char *FilePath, int *nPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, double *Z, int *Error) {
  P3DInternalReadSingleGrid(FilePath, 3, nPoints, 0, Endian, Offset, X, Y, Z, NULL, Error);
}
void P3DInternalReadSingleGrid3DWithIBlank(char *FilePath, int *nPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, double *Z, int *IBlank, int *Error) {
  P3DInternalReadSingleGrid(FilePath, 3, nPoints, 1, Endian, Offset, X, Y, Z, IBlank, Error);
}

// Write a PLOT3D grid to a grid file
void P3DInternalWriteSingleGrid(char *FilePath, int nDims, int *nPoints, int WithIBlank,
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
  GridSize = nPoints[0]*nPoints[1]*nPoints[2];
  Record = sizeof(double)*nDims*GridSize + sizeof(int)*GridSize;
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);
  fwrite_endian(Endian, X, sizeof(double), GridSize, GridFile);
  fwrite_endian(Endian, Y, sizeof(double), GridSize, GridFile);
  if (nDims == 3) {
    fwrite_endian(Endian, Z, sizeof(double), GridSize, GridFile);
  }
  if (WithIBlank == 1) {
    fwrite_endian(Endian, IBlank, sizeof(int), GridSize, GridFile);
  }
  fwrite_endian(Endian, &Record, sizeof(p3d_size_type), 1, GridFile);

  fclose(GridFile);

}

// Can't easily pass null pointers from Fortran, so call these wrappers instead
void P3DInternalWriteSingleGrid2D(char *FilePath, int *nPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, 2, nPoints, 0, Endian, Offset, X, Y, NULL, NULL, Error);
}
void P3DInternalWriteSingleGrid2DWithIBlank(char *FilePath, int *nPoints,
  p3d_endian_type Endian, p3d_offset_type Offset, double *X, double *Y, int *IBlank, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, 2, nPoints, 1, Endian, Offset, X, Y, NULL, IBlank, Error);
}
void P3DInternalWriteSingleGrid3D(char *FilePath, int *nPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, double *Z, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, 3, nPoints, 0, Endian, Offset, X, Y, Z, NULL, Error);
}
void P3DInternalWriteSingleGrid3DWithIBlank(char *FilePath, int *nPoints, p3d_endian_type Endian,
  p3d_offset_type Offset, double *X, double *Y, double *Z, int *IBlank, int *Error) {
  P3DInternalWriteSingleGrid(FilePath, 3, nPoints, 1, Endian, Offset, X, Y, Z, IBlank, Error);
}

static void SwapEndian(void *Data, p3d_size_type ElementSize, p3d_size_type nElements,
  void *SwappedData) {

  int i, j;
  char *A, *B;

  for (i = 0; i < nElements; ++i) {
    A = (char *)Data + i*ElementSize;
    B = (char *)SwappedData + i*ElementSize;
    for (j = 0; j < ElementSize; ++j) {
      B[j] = A[ElementSize-j-1];
    }
  }

}

// Read data, possibly with an endianness that is different from the machine's
static p3d_size_type fread_endian(p3d_endian_type Endian, void *Data, p3d_size_type ElementSize,
  p3d_size_type nElements, FILE *In) {

  int DifferentEndian;
  void *DataToRead;
  p3d_size_type Result;

  DifferentEndian = Endian != P3DInternalMachineEndian();

  if (DifferentEndian) {
    DataToRead = malloc(ElementSize * nElements);
  } else {
    DataToRead = Data;
  }

  Result = fread(DataToRead, ElementSize, nElements, In);

  if (DifferentEndian) {
    SwapEndian(DataToRead, ElementSize, nElements, Data);
    free(DataToRead);
  }

  return Result;

}

// Write data, possibly with an endianness that is different from the machine's
static p3d_size_type fwrite_endian(p3d_endian_type Endian, void *Data, p3d_size_type ElementSize, 
  p3d_size_type nElements, FILE *Out) {

  int DifferentEndian;
  void *DataToWrite;
  p3d_size_type Result;

  DifferentEndian = Endian != P3DInternalMachineEndian();

  if (DifferentEndian) {
    DataToWrite = malloc(ElementSize * nElements);
    SwapEndian(Data, ElementSize, nElements, DataToWrite);
  } else {
    DataToWrite = Data;
  }

  Result = fwrite(DataToWrite, ElementSize, nElements, Out);

  if (DifferentEndian) {
    free(DataToWrite);
  }

  return Result;

}

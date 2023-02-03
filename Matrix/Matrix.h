// File: $Id: Matrix.h,v 1.1 2016/10/27 16:49:54 csci243 Exp csci243 $
//
// Description:  Interface to the Matrix module.
//
// Author: wrc 
//

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdio.h>
#include <stdbool.h>

// // // // // // // // // // // // // // // // //
//
// Type declarations
//
// // // // // // // // // // // // // // // // //

/// The opaque Matrix type.

typedef struct matrix_st *Matrix;

/// Status and Error codes

typedef enum {
    Success, BadRowNumber, BadColNumber
} Status;

// // // // // // // // // // // // // // // // //
//
// Function prototypes
//
// // // // // // // // // // // // // // // // //

/// Create a Matrix with the specified dimension.
///
/// If the matrix is square, it is initialized to the identity
/// matrix; otherwise, it is initialized to the zero matrix.
///
/// @return a pointer to a Matrix on success, or NULL if the creation fails

Matrix mat_create( size_t rows, size_t cols );


/// Destroy the supplied Matrix, releasing any dynamic memory it used.
///
/// @param mat the Matrix to be destroyed

void mat_destroy( Matrix mat );


/// Initialize a Matrix from an array of values.
///
/// Assumes that the data parameter contains enough values
/// to completely initialize the matrix, with the values
/// held in row-major order (i.e., for M rows and N columns,
/// the first N values initialize the first row, the next N
/// values initialize the second row, etc.).
///
/// @param mat  the matrix to be initialized
/// @param data the array to be copied into the matrix

void mat_init( Matrix mat, const float data[] );


/// Duplicate a Matrix.
///
/// @param mat the matrix to be duplicated
///
/// @return a pointer to the new Matrix on success, or NULL if any
///         error occurred (such as an allocation failure)

Matrix mat_duplicate( const Matrix mat );


/// Compare two matrices for equality.
///
/// Two matrices are equal if they have the same dimensions,
/// and for all rows and columns, the corresponding entries
/// in the matrices are equal.
///
/// @param m1 the first operand Matrix
/// @param m2 the second operand Matrix
///
/// @return true if the matrices are equal, else false

bool mat_equals( const Matrix m1, const Matrix m2 );


/// Multiply a matrix by a scalar.
///
/// @param mat  the Matrix to be modified
/// @param data the scalar multiplier

void mat_scalar_mult( Matrix mat, float data );


/// Multiply two matrices, returning a new result matrix.
///
/// If the matrices are compatible, creates a new Matrix of the proper
/// size and populates it with the matrix product of the two source
/// matrices.
///
/// @param m1 the first matrix operand
/// @param m2 the second matrix operand
///
/// @return the product matrix on success, or NULL if an error occurred

Matrix mat_mult( const Matrix m1, const Matrix m2 );


/// Retrieve the value from one cell in a Matrix.
///
/// @param mat  the matrix to examine
/// @param data a pointer to a float into which the value is to be placed
/// @param row  the row number of the desired cell (1..M)
/// @param col  the column number of the desired cell
///
/// @return a status code indicating the result of the operation

Status mat_get_cell( const Matrix mat, float *data, size_t row, size_t col );


/// Retrieve the values from one row in a Matrix.
///
/// Assumes that the data parameter is an array large enough
/// to hold the contents of the supplied Matrix.  Copies the
/// data from the Matrix into the data array in row-major
/// storage order (i.e., the first row in the matrix is saved
/// into the first four elements of the data array, the second
/// row into the second four elements, etc.).
///
/// @param mat  the matrix to update
/// @param data an array into which the row values are to be placed
/// @param row  the row number of the desired cell (1..M)
///
/// @return a status code indicating the result of the operation

Status mat_get_row( const Matrix mat, float data[], size_t row );


/// Assign a value to one cell in a Matrix.
///
/// @param mat  the matrix to update
/// @param data the value to assign
/// @param row  the row number of the desired cell (1..M)
/// @param col  the column number of the desired cell
///
/// @return a status code indicating the result of the operation

Status mat_set_cell( Matrix mat, float data, size_t row, size_t col );


/// Assign a set of values to one row in a Matrix.
///
/// Assumes that the data parameter contains enough values
/// to completely initialize the specified row of the Matrix.
///
/// @param mat  the matrix to update
/// @param data an array containing values to be placed into the row
/// @param row  the row number of the desired cell (1..M)
///
/// @return a status code indicating the result of the operation

Status mat_set_row( Matrix mat, const float data[], size_t row );

/// Transpose a matrix
///
/// @param mat the Matrix to be transposed
///
/// @return a pointer to a new Matrix on success, or NULL

Matrix mat_transpose( const Matrix mat );


/// Print a matrix
///
/// Prints the matrix in a standard format:
///
///    M rows, N columns:
///       v11   v12   ...   v1N
///       v21   v22   ...   v2N
///        ...
///       vM1   vM2   ...   vMN
///
/// Each "vIJ" is printed using an 8-character field with
/// three digits to the right of the decimal point.
///
/// @param stream where to print the Matrix
/// @param mat    the Matrix to be printed

void mat_print( const Matrix mat, FILE *stream );

#endif

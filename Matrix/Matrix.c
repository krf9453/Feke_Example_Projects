/*
 * Matrix.c program implements all of the functions from Matrix.h
 * @name Kelly Feke
 * @login krf9453
 */

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <strings.h>


#include "Matrix.h"

/// the Matrix struct holds the number of rows in matrix,
/// the number of colunms in matrix, and the pointer where the 
/// values of the matrix are held
struct matrix_st {
    size_t rows;    /// num rows in matrix 
    size_t cols;    /// num cols in matrix
    float *contents;/// contents of matrix
};

/*
 * Create a Matrix with the specified dimension.
 * If the matrix is square, it is initialized to the identity
 * matrix; otherwise, it is initialized to the zero matrix.
 *
 * @return a pointer to a Matrix on success, or NULL if the creation fails
*/
Matrix mat_create ( size_t rows, size_t cols ) {
    Matrix new; /// new Matrix being created
    new = (Matrix) malloc( sizeof(struct matrix_st));
    if (new != 0) {
        
        new -> contents = malloc (rows * cols * sizeof(float));

        /// if a square matrix
        if (rows == cols ){
            for (size_t i = 0; i < rows; i++){
                for (size_t j = 0; j < cols; j++){
                    //make identity matrix
                    if ( i == j ){
                        new -> contents[i * cols + j] = 1;
                    } else {
                        new -> contents[i * cols + j] = 0;
                }
            }
        }
        }
        /// not a square matrix
        else {
            for (size_t i = 0; i < rows; i++){
                for (size_t j = 0; j < cols; j++){
                    new -> contents[i * cols + j] = 0;
                }
            }
        }
        new -> rows = rows;
        new -> cols = cols;
    }
    return (new);
}

/*
 * Destroy the supplied Matrix, releasing any dynamic memory it used.
 * 
 * @param mat the Matrix to be destroyed
 */
void mat_destroy( Matrix mat ) {
    assert( mat != 0 );

    if ( mat -> contents != 0 ) {
        free( mat  -> contents );
    }

    free(mat);
}




/*
 * Initialize a Matrix from an array of values.
 * Assumes that the data parameter contains enough values
 * to completely initialize the matrix, with the values
 * held in row-major order (i.e., for M rows and N columns,
 * the first N values initialize the first row, the next N
 * values initialize the second row, etc.).
 *
 * @param mat  the matrix to be initialized
 * @param data the array to be copied into the matrix
 */
void mat_init (Matrix mat, const float data[]) {
    
    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            mat -> contents[i * mat->cols + j] = data[i * mat->cols + j];

        }
    }
}

/*
 * Duplicate a Matrix.
 *
 * @param mat the matrix to be duplicated
 *
 * @return a pointer to the new Matrix on success, or NULL if any 
 * error occurred (such as an allocation failure)
 */
Matrix mat_duplicate ( const Matrix mat) {
    Matrix new = mat_create(mat->rows, mat->cols); /// duped matrix

    if (new != 0) {
        new->cols = mat->cols;
        new->rows = mat->rows;
        for (size_t i = 0; i < mat->rows; i++) {
            for (size_t j = 0; j < mat->cols; j++) {
                new -> contents[i * mat->cols + j] = mat -> contents[i * mat->cols + j];
            }
        }
        return new;
    } else {
        return NULL;
    }
}



/*
 * Compare two matrices for equality.
 *
 * Two matrices are equal if they have the same dimensions,
 * and for all rows and columns, the corresponding entries
 * in the matrices are equal.
 *
 * @param m1 the first operand Matrix
 * @param m2 the second operand Matrix
 * @return true if the matrices are equal, else false
 */
bool mat_equals( const Matrix m1, const Matrix m2 ){
    if (m1->rows != m2->rows || m1->cols != m2->cols) {
        return false;
    }
    for (size_t i = 0; i < m1->rows; i++) {
        for (size_t j = 0; j < m1->cols; j++) {
            if (m1 -> contents[i * m1->cols + j] != m2 -> contents[i * m2->cols + j]) {
                return false;
            }
        }
    }
    return true;
}

/*
 * Multiply a matrix by a scalar.
 *
 * @param mat  the Matrix to be modified
 * @param data the scalar multiplier
 */
void mat_scalar_mult( Matrix mat, float data ) {
    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            mat -> contents[i * mat->cols + j] = data * mat -> contents[i * mat->cols + j];
        }
    }
}

/*
 * Multiply two matrices, returning a new result matrix.
 * If the matrices are compatible, creates a new Matrix of the proper
 * size and populates it with the matrix product of the two source
 * matrices.
 *
 * @param m1 the first matrix operand
 * @param m2 the second matrix operand
 * @return the product matrix on success, or NULL if an error occurred
 */
Matrix mat_mult( const Matrix m1, const Matrix m2 ) {
    if (m1->cols != m2->rows) {
        return NULL;
    }

    Matrix new = mat_create(m1->rows, m2->cols); /// new product matrix
    float sum = 0;

    /// k iterate through cols of m2
    for (size_t k = 0; k < m2->cols; k++) {
        /// i iterate through rows of m1
        for (size_t i = 0; i < m1->rows; i++) {
            /// j iterate through rows of m2 (= cols of m1)
            for (size_t j = 0; j < m2->rows; j++) {
                sum += (m1 -> contents[i * m1->cols + j]) * (m2 -> contents[j * m2->cols + k]);
            }   
            new -> contents[i * new->cols + k] = sum;
            sum = 0;

        }
    }
    return new;
}

/*
 * Retrieve the value from one cell in a Matrix.
 * 
 * @param mat  the matrix to examine
 * @param data a pointer to a float into which the value is to be placed
 * @param row  the row number of the desired cell (1..M)
 * @param col  the column number of the desired cell
 * @return a status code indicating the result of the operation
 */
Status mat_get_cell( const Matrix mat, float *data, size_t row, size_t col) { 
    if (row > mat->rows) {
        return BadRowNumber;
    }
    if (col > mat->cols) {
        return BadColNumber;
    }
    float cell;
    cell = mat -> contents[(row-1)*mat->cols +(col-1)];
    *data = cell;
    return Success;

}

/*
 * Retrieve the values from one row in a Matrix.
 *
 * Assumes that the data parameter is an array large enough
 * to hold the contents of the supplied Matrix.  Copies the
 * data from the Matrix into the data array in row-major
 * storage order (i.e., the first row in the matrix is saved
 * into the first four elements of the data array, the second
 * row into the second four elements, etc.).
 *
 * @param mat  the matrix to update
 * @param data an array into which the row values are to be placed
 * @param row  the row number of the desired cell (1..M)
 * @return a status code indicating the result of the operation
 */
Status mat_get_row( const Matrix mat, float data[], size_t row) {

    if (row > mat->rows) {
        return BadRowNumber;
    }

    for (size_t j = 0; j < mat->cols; j++) {
        data[j] = mat -> contents[(row-1) * mat->cols + j];
    }

    return Success;
}

/*
 * Assign a value to one cell in a Matrix.
 *
 * @param mat  the matrix to update
 * @param data the value to assign
 * @param row  the row number of the desired cell (1..M)
 * @param col  the column number of the desired cell
 * @return a status code indicating the result of the operation
 */
Status mat_set_cell( Matrix mat, float data, size_t row, size_t col) {
    if (row > mat->rows) {
        return BadRowNumber;
    }
    if (col > mat->cols) {
        return BadColNumber;
    }

    mat -> contents[(row-1) * mat->cols + (col-1)] = data;

    return Success;
}

/*
 * Assign a set of values to one row in a Matrix.
 * Assumes that the data parameter contains enough values
 * to completely initialize the specified row of the Matrix.
 * @param mat  the matrix to update
 * @param data an array containing values to be placed into the row
 * @param row  the row number of the desired cell (1..M)
 * @return a status code indicating the result of the operation
 */
Status mat_set_row( Matrix mat, const float data[], size_t row) {

    if (row > mat->rows) {
        return BadRowNumber;
    }

    for (size_t j = 0; j < mat->cols; j++) {
        mat -> contents[(row-1) * mat->cols + j] = data[j];
    }

    return Success;
}


/*
 * Transpose a matrix
 *
 * @param mat the Matrix to be transposed
 * @return a pointer to a new Matrix on success, or NULL
 */
Matrix mat_transpose( const Matrix mat) {

    Matrix new = mat_create(mat->cols, mat->rows);

    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            new -> contents[j * new->cols + i] = mat -> contents[i * mat->cols + j];
        }
    }
    return new;
}

/*
 * Print a matrix
 * Prints the matrix in a standard format:
 *
 *    M rows, N columns:
 *       v11   v12   ...   v1N
 *       v21   v22   ...   v2N
 *        ...
 *       vM1   vM2   ...   vMN
 *
 * Each "vIJ" is printed using an 8-character field with
 * three digits to the right of the decimal point.
 *
 * @param stream where to print the Matrix
 * @param mat    the Matrix to be printed


 */
void mat_print( const Matrix mat, FILE *stream) {


    fprintf(stream, "%zu rows, %zu columns:\n", mat->rows, mat->cols);

    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            fprintf(stream, "%*.3f", 8, mat -> contents[i * mat->cols + j]);
        }
        fprintf(stream, "\n");
    }
}


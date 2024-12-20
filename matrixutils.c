#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixutils.h"
#include "symnmf.h"

/** Multiplies two matrices and returns the resulting matrix. */
double** multiply_matrix(double **matrix1, int rows_matrix1, int cols_matrix1, double **matrix2, int rows_matrix2, int cols_matrix2)
{

    double **result;
    int i, j, k;
    if(cols_matrix1 != rows_matrix2)
    {
        throw_error();
    }
    result = matrix_memory_calloc(rows_matrix1, cols_matrix2);
    for (i = 0; i < rows_matrix1; i++){
        for (j = 0; j < cols_matrix2; j++){
            result[i][j] = 0;
            for (k = 0; k < cols_matrix1; k++){
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}

/** Transposes a matrix and returns the resulting matrix. */
double** transpose_matrix(double **matrix, int rows, int cols)
{
    double **result;
    int i, j;
    result = matrix_memory_calloc(cols, rows);
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}

/** Copies the contents of the source matrix to the target matrix. */
void copy_matrix(double **source_matrix, double **target_matrix, int rows, int cols){
    int i, j;
    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            target_matrix[i][j] = source_matrix[i][j];
        }
    }
}

/** Computes the Frobenius norm of a matrix. */
double frobenius_norm(double **matrix, int rows, int cols)
{
    double sum = 0;
    int i, j;
    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            sum += pow(matrix[i][j], 2);
        }
    }
    return sum;
}

/** Computes the difference between two matrices and returns the resulting matrix. */
double** difference_matrix(double **A, double **B, int rows, int cols)
{
    double **diff;
    int i, j = 0;
    diff = matrix_memory_calloc(rows, cols);
    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            diff[i][j] = A[i][j] - B[i][j];
        }
    }
    return diff;
}

/** Prints the matrix with values formatted to 4 decimal places, separated by commas. */
void print_output(double **matrix, int N, int M)
{
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < M; j++){
            printf("%.4f", matrix[i][j]);
            if(j != M - 1){
                printf(",");
            }
        }
        printf("\n");
    }
}

/** Allocates memory for a matrix using calloc and returns the allocated matrix. */
double** matrix_memory_calloc(int rows, int cols){
    double **result;
    int i, j;

    result = (double**)calloc(rows, sizeof(double*));
    if (result == NULL){
        throw_error();
    }

    for (i = 0; i < rows; i++){
        result[i] = (double*)calloc(cols, sizeof(double));
        if (result[i] == NULL){
            for (j = 0; j < i; j++){
                free(result[j]);
            }
            free(result);
            throw_error();
        }
    }
    return result;
}

/** Frees the memory allocated for a matrix */
void matrix_memory_free(double **matrix, int rows){
    int i;
    if (matrix == NULL) {
        return;
    }
    for (i = 0; i < rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}

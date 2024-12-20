# ifndef MATRIXUTILS_H
# define MATRIXUTILS_H

double** multiply_matrix(double **matrix1, int rows_matrix1, int cols_matrix1, double **matrix2, int rows_matrix2, int cols_matrix2);
double** transpose_matrix(double **matrix, int rows, int cols);
void copy_matrix(double **source_matrix, double **target_matrix, int rows, int cols);
double frobenius_norm(double **matrix, int rows, int cols);
double** difference_matrix(double **A, double **B, int rows, int cols);
void print_output(double **matrix, int N, int M);
double** matrix_memory_calloc(int rows, int cols);
void matrix_memory_free(double **matrix, int rows);

# endif


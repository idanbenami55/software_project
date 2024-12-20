#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"
#include "matrixutils.h"
#define BETA 0.5
#define SYM "sym"
#define DDG "ddg"
#define NORM "norm"

/** Prints an error message and exits the program. */
void throw_error(void)
{
    fprintf(stderr, "An Error Has Occurred\n");
    exit(1);
}

/** Calculates the squared Euclidean distance between two vectors. */
double euclidean_distance(double *a, double *b, int d)
{
    double sum = 0, diff;
    int i;
    for (i = 0; i < d; i++){
        diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sum;
}

/** Computes and returns the symmetric similarity matrix. */
double** sym(double **points, int N, int d)
{
    double **symatrix, distance;
    int i, j;
    symatrix = matrix_memory_calloc(N, N);
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i == j){
                symatrix[i][j] = 0.0;
            }
            else{
                distance = euclidean_distance(points[i], points[j], d);
                symatrix[i][j] = exp(-distance / 2.0);
            }
        }
    }
    return symatrix;
}

/** Computes and returns the diagonal degree matrix. */
double** ddg(double **points, int N, int d)
{
    double **diag, **symatrix, sum;
    int i, j;
    symatrix = sym(points, N, d);
    if (symatrix == NULL) {
        throw_error();
    }
    diag = matrix_memory_calloc(N, N);
    for (i = 0; i < N; i++){
        sum = 0;
        for (j = 0; j < N; j++){
            sum += symatrix[i][j];
        }
        diag[i][i] = sum;
    }
    matrix_memory_free(symatrix, N);
    return diag;
}

/** Computes and returns diag^power for a given diagonal matrix. */
double** diag_power(double **diag, double power, int N)
{
    double **diag_pow;
    int i;
    diag_pow = matrix_memory_calloc(N, N);
    for (i = 0; i < N; i++){
        if (diag[i][i] == 0.0 && power < 0.0){
            matrix_memory_free(diag_pow, N);
            return NULL;
        } else {
            diag_pow[i][i] = pow(diag[i][i], power);
        }
    }
    return diag_pow;
}

/** Computes and returns the normalized similarity matrix using the formula D^(-0.5) * A * D^(-0.5). */
double** norm(double **points, int N, int d)
{
    double **symatrix, **diag, **norm, **squared_diag, **temp;
    symatrix = sym(points, N, d);
    if (symatrix == NULL) {
        throw_error();
    }
    diag = ddg(points, N, d);
    if (diag == NULL) {
        matrix_memory_free(symatrix, N);
        throw_error();
    }
    squared_diag = diag_power(diag, -0.5, N);
    if (squared_diag == NULL) {
        matrix_memory_free(symatrix, N);
        matrix_memory_free(diag, N);
        throw_error();
    }
    temp = multiply_matrix(squared_diag, N, N, symatrix, N, N);
    if (temp == NULL) {
        matrix_memory_free(symatrix, N);
        matrix_memory_free(diag, N);
        matrix_memory_free(squared_diag, N);
        throw_error();
    }
    norm = multiply_matrix(temp, N, N, squared_diag, N, N);
    if (norm == NULL) {
        matrix_memory_free(symatrix, N);
        matrix_memory_free(diag, N);
        matrix_memory_free(squared_diag, N);
        matrix_memory_free(temp, N);
        throw_error();
    }
    matrix_memory_free(symatrix, N);
    matrix_memory_free(diag, N);
    matrix_memory_free(temp, N);
    matrix_memory_free(squared_diag, N);
    return norm;
}

/** Updates matrix H using the Symmetric NMF multiplicative update rule. */
void update_h(double **next_matrix, double **cur_matrix, double **norm_matrix, int N, int k) 
{
    double **wh_matrix, **transpose_cur_matrix, **temp_matrix, **hhth_matrix;
    int i, j;
    wh_matrix = multiply_matrix(norm_matrix, N, N, cur_matrix, N, k);
    if (wh_matrix == NULL) {
        throw_error();
    }
    transpose_cur_matrix = transpose_matrix(cur_matrix, N, k);
    if (transpose_cur_matrix == NULL) {
        matrix_memory_free(wh_matrix, N);
        throw_error();
    }
    temp_matrix = multiply_matrix(cur_matrix, N, k, transpose_cur_matrix, k, N);
    if (temp_matrix == NULL) {
        matrix_memory_free(wh_matrix, N);
        matrix_memory_free(transpose_cur_matrix, k);
        throw_error();
    }
    hhth_matrix = multiply_matrix(temp_matrix, N, N, cur_matrix, N, k);
    if (hhth_matrix == NULL) {
        matrix_memory_free(wh_matrix, N);
        matrix_memory_free(transpose_cur_matrix, k);
        matrix_memory_free(temp_matrix, N);
        throw_error();
    }
    for (i = 0; i < N; i++){
        for (j = 0; j < k; j++){
            next_matrix[i][j] = cur_matrix[i][j] * (1 - BETA * (1 - wh_matrix[i][j] / hhth_matrix[i][j]));
        }
    }
    matrix_memory_free(wh_matrix, N);
    matrix_memory_free(transpose_cur_matrix, k);
    matrix_memory_free(temp_matrix, N);
    matrix_memory_free(hhth_matrix, N);
}

/** Optimizes matrix H using iterative updates until convergence or max iterations. */
double** optimize_h(double **init_matrix, double **norm_matrix, int N, int k, int max_iter, double epsilon)
{
    double **cur_matrix, **next_matrix, **diff;
    int iter_counter = 0;
    next_matrix = matrix_memory_calloc(N, k);
    if (next_matrix == NULL) {
        throw_error();
    }
    cur_matrix = matrix_memory_calloc(N, k);
    if (cur_matrix == NULL) {
        matrix_memory_free(next_matrix, N);
        throw_error();
    }
    copy_matrix(init_matrix, cur_matrix, N, k);
    update_h(next_matrix, init_matrix, norm_matrix, N, k);
    iter_counter += 1;
    diff = difference_matrix(next_matrix, cur_matrix, N, k);
    if (diff == NULL) {
        matrix_memory_free(next_matrix, N);
        matrix_memory_free(cur_matrix, N);
        throw_error();
    }
    while (frobenius_norm(diff, N, k) >= epsilon && iter_counter < max_iter)
    {
        copy_matrix(next_matrix, cur_matrix, N, k);
        matrix_memory_free(diff, N);
        update_h(next_matrix, cur_matrix, norm_matrix, N, k);
        diff = difference_matrix(next_matrix, cur_matrix, N, k);
        if (diff == NULL) {
            matrix_memory_free(next_matrix, N);
            matrix_memory_free(cur_matrix, N);
            throw_error();
        }
        iter_counter += 1;
    }
    matrix_memory_free(diff, N);
    matrix_memory_free(cur_matrix, N);
    return next_matrix;
}

/** Counts the number of rows in a file. */
int num_of_rows(FILE *file) {
    int row = 0;
    int c;
    while ((c = fgetc(file)) != EOF) {
        if (c == '\n') {
            row++;
        }
    }
    if (ftell(file) > 0 && c != '\n') {
        row++;
    }
    fseek(file, 0, SEEK_SET);
    return row - 1;
}

/** Counts the number of columns in a file based on the first line. */
int num_of_columns(FILE *input_file) {
    int c;
    int col = 1; 
    while ((c = fgetc(input_file)) != '\n' && c != EOF) {
        if (c == ',') {
            col++;
        }
    }
    fseek(input_file, 0, SEEK_SET);
    return col;
}

/** Initializes a matrix A with values read from a file. */
int init_matrix(double **A, FILE *input_file, int n, int m){
    int i = 0;
    int j = 0;
    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
            if (fscanf(input_file, "%lf,", &A[i][j]) != 1){
                return 1;
            }
        }
    }
    return 0;
}

int main(int argc, char **argv)
{
    char *goal, *file_name;
    int dim = 0, points_number = 0;
    double **points, **res = NULL;
    FILE *file;
    if (argc != 3) { throw_error(); }
    goal = argv[1];
    file_name = argv[2];
    file = fopen(file_name,"r");
    if(file == NULL) {throw_error();}
    points_number = num_of_rows(file);
    dim = num_of_columns(file);
    points = matrix_memory_calloc(points_number, dim);
    if (init_matrix(points, file, points_number, dim) == 1){
        matrix_memory_free(points, points_number);
        fclose(file);
        throw_error();
    }
    fclose(file);
    if (strcmp(goal, SYM) == 0) {
        res = sym(points, points_number, dim);
    }
    else if(strcmp(goal, DDG) == 0) {
        res = ddg(points, points_number, dim);
    }
    else if(strcmp(goal, NORM) == 0) {
        res = norm(points, points_number, dim);
    }
    else {
        throw_error();
    }
    print_output(res, points_number, points_number);
    matrix_memory_free(points, points_number);
    matrix_memory_free(res, points_number);
    return 0;
}


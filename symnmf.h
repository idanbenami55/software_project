# ifndef SYMNMF_H
# define SYMNMF_H

void throw_error(void);
double euclidean_distance(double *a, double *b, int d);
double** sym(double** points, int N, int d);
double** ddg(double **points, int N, int d);
double** diag_power(double **diag, double power, int N);
double** norm(double **points, int N, int d);
void update_h(double **next_matrix, double **cur_matrix, double **norm_matrix, int N, int k);
double** optimize_h(double **init_matrix, double **norm_matrix, int N, int k, int max_iter, double epsilon);
double** read_points(char *file_name, int *points_number, int *dim);

#endif


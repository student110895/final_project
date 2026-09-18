#ifndef SYMNMF_H
#define SYMNMF_H

double **allocate_matrix(int rows, int cols);
void free_matrix(double **matrix, int rows);

double **sym(double **points, int n, int d);
double **ddg(double **points, int n, int d);
double **norm(double **points, int n, int d);

double **symnmf(double **H, double **W, int n, int k);

#endif

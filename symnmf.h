#ifndef SYMNMF_H
#define SYMNMF_H

/* --- Memory Management --- */

/* Allocates contiguous memory for a 2D array of doubles */
double **allocate_matrix(int rows, int cols);

/* Safely frees a dynamically allocated 2D array */
void free_matrix(double **matrix, int rows);


/* --- Mathematical Matrix Operations --- */

/* Calculates the Similarity Matrix (A) using exponential distances */
double **sym(double **points, int n, int d);

/* Generates the Degree Matrix (D) based on the Similarity Matrix */
double **ddg(double **points, int n, int d);

/* Calculates the Normalized Similarity Matrix (W) */
double **norm(double **points, int n, int d);


/* --- Core SymNMF Algorithm --- */

/* Computes the optimized H matrix until convergence or max iterations */
double **symnmf(double **H, double **W, int n, int k);

#endif
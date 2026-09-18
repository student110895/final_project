#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

#define ERROR_MSG "An Error Has Occurred"
#define BETA 0.5
#define EPSILON 0.0001
#define MAX_ITER 300
#define DENOM_EPSILON 0.000001

/* Allocates contiguous memory for a 2D array of doubles */
double **allocate_matrix(int rows, int cols) {
    double **matrix;
    int i, j;

    matrix = (double **)malloc(rows * sizeof(double *));
    if (!matrix) return NULL;

    for (i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(cols * sizeof(double));
        if (!matrix[i]) {
            for (j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

/* Safely frees a dynamically allocated 2D array */
void free_matrix(double **matrix, int rows) {
    int i;
    if (!matrix) return;
    for (i = 0; i < rows; i++) free(matrix[i]);
    free(matrix);
}

/* Creates an independent deep copy of a given matrix */
static double **copy_matrix(double **matrix, int rows, int cols) {
    double **copy;
    int i, j;

    copy = allocate_matrix(rows, cols);
    if (!copy) return NULL;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            copy[i][j] = matrix[i][j];
        }
    }
    return copy;
}

/* Parses the input file to count the number of data points (n) and dimensions (d) */
static int get_dimensions(FILE *file, int *n, int *d) {
    int c, commas = 0, lines = 0, first_line = 1, last = '\n';

    while ((c = fgetc(file)) != EOF) {
        if (first_line && c == ',') commas++;
        if (c == '\n' && last != '\n') { lines++; first_line = 0; }
        last = c;
    }
    if (last != '\n') lines++;
    
    *n = lines;
    *d = commas + 1;
    rewind(file);
    return lines > 0;
}

/* Reads comma-separated float values from the file into an n x d matrix */
static double **read_points(FILE *file, int n, int d) {
    double **points;
    int i, j, separator;

    points = allocate_matrix(n, d);
    if (!points) return NULL;

    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            if (fscanf(file, "%lf", &points[i][j]) != 1) {
                free_matrix(points, n); return NULL;
            }
            if (j < d - 1) {
                separator = fgetc(file);
                if (separator != ',') { free_matrix(points, n); return NULL; }
            }
        }
    }
    return points;
}

/* Calculates the squared Euclidean distance between two vectors */
static double squared_distance(double *x, double *y, int d) {
    double sum = 0.0, diff;
    int i;
    for (i = 0; i < d; i++) {
        diff = x[i] - y[i];
        sum += diff * diff;
    }
    return sum;
}

/* Calculates the Similarity Matrix (A) using exponential distances */
double **sym(double **points, int n, int d) {
    double **A;
    int i, j;

    A = allocate_matrix(n, n);
    if (!A) return NULL;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) A[i][j] = 0.0;
            else A[i][j] = exp(-squared_distance(points[i], points[j], d) / 2.0);
        }
    }
    return A;
}

/* Calculates the Degree Matrix (D) as a diagonal matrix of row sums */
static double **degree_from_similarity(double **A, int n) {
    double **D, sum;
    int i, j;

    D = allocate_matrix(n, n);
    if (!D) return NULL;

    for (i = 0; i < n; i++) {
        sum = 0.0;
        for (j = 0; j < n; j++) {
            D[i][j] = 0.0;
            sum += A[i][j];
        }
        D[i][i] = sum;
    }
    return D;
}

/* Generates the Similarity Matrix and returns its corresponding Degree Matrix */
double **ddg(double **points, int n, int d) {
    double **A, **D;
    
    A = sym(points, n, d);
    if (!A) return NULL;
    
    D = degree_from_similarity(A, n);
    free_matrix(A, n);
    return D;
}

/* Calculates the Normalized Similarity Matrix (W) */
double **norm(double **points, int n, int d) {
    double **A, **D, **W;
    int i, j;

    A = sym(points, n, d);
    if (!A) return NULL;
    
    D = degree_from_similarity(A, n);
    if (!D) { free_matrix(A, n); return NULL; }
    
    W = allocate_matrix(n, n);
    if (!W) { free_matrix(A, n); free_matrix(D, n); return NULL; }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (D[i][i] == 0.0 || D[j][j] == 0.0) W[i][j] = 0.0;
            else W[i][j] = A[i][j] / (sqrt(D[i][i]) * sqrt(D[j][j]));
        }
    }
    free_matrix(A, n); free_matrix(D, n);
    return W;
}

/* Performs standard matrix multiplication for matrices of compliant dimensions */
static double **multiply(double **A, double **B, int rows, int common, int cols) {
    double **res, sum;
    int i, j, k;

    res = allocate_matrix(rows, cols);
    if (!res) return NULL;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            sum = 0.0;
            for (k = 0; k < common; k++) sum += A[i][k] * B[k][j];
            res[i][j] = sum;
        }
    }
    return res;
}

/* Calculates H^T * H without allocating the transpose matrix in memory */
static double **calculate_ht_h(double **H, int n, int k) {
    double **res, sum;
    int i, j, row;

    res = allocate_matrix(k, k);
    if (!res) return NULL;

    for (i = 0; i < k; i++) {
        for (j = 0; j < k; j++) {
            sum = 0.0;
            for (row = 0; row < n; row++) sum += H[row][i] * H[row][j];
            res[i][j] = sum;
        }
    }
    return res;
}

/* Calculates the squared Frobenius norm of the difference between two matrices */
static double frobenius_difference(double **A, double **B, int rows, int cols) {
    double sum = 0.0, diff;
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            diff = A[i][j] - B[i][j];
            sum += diff * diff;
        }
    }
    return sum;
}

/* Computes the next iteration matrix for H in the SymNMF algorithm */
static double **update_h(double **H, double **W, int n, int k) {
    double **WH, **HtH, **denom, **next;
    int i, j;

    WH = multiply(W, H, n, n, k);
    HtH = calculate_ht_h(H, n, k);
    if (!WH || !HtH) { free_matrix(WH, n); free_matrix(HtH, k); return NULL; }

    denom = multiply(H, HtH, n, k, k);
    next = allocate_matrix(n, k);
    if (!denom || !next) {
        free_matrix(WH, n); free_matrix(HtH, k);
        free_matrix(denom, n); free_matrix(next, n);
        return NULL;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            next[i][j] = H[i][j] * ((1.0 - BETA) + BETA * WH[i][j] / (denom[i][j] + DENOM_EPSILON));
        }
    }
    free_matrix(WH, n); free_matrix(HtH, k); free_matrix(denom, n);
    return next;
}

/* Core SymNMF optimization algorithm looping until convergence or MAX_ITER */
double **symnmf(double **H, double **W, int n, int k) {
    double **curr, **next, diff;
    int iter;

    curr = copy_matrix(H, n, k);
    if (!curr) return NULL;

    for (iter = 0; iter < MAX_ITER; iter++) {
        next = update_h(curr, W, n, k);
        if (!next) { free_matrix(curr, n); return NULL; }

        diff = frobenius_difference(next, curr, n, k);
        free_matrix(curr, n);
        curr = next;

        if (diff < EPSILON) break;
    }
    return curr;
}

/* Prints a given matrix to stdout with exactly 4 decimal places per element */
static void print_matrix(double **matrix, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", matrix[i][j]);
            if (j < cols - 1) printf(",");
        }
        printf("\n");
    }
}

/* Validates that the requested goal matches one of the three C operations */
static int valid_goal(char *goal) {
    return strcmp(goal, "sym") == 0 || strcmp(goal, "ddg") == 0 || strcmp(goal, "norm") == 0;
}

/* Routes execution to the appropriate mathematical function based on user input */
static double **run_goal(char *goal, double **points, int n, int d) {
    if (strcmp(goal, "sym") == 0) return sym(points, n, d);
    if (strcmp(goal, "ddg") == 0) return ddg(points, n, d);
    return norm(points, n, d);
}

/* Entry point for standalone C execution and argument validation */
int main(int argc, char **argv) {
    FILE *file;
    double **points, **result;
    int n, d;

    if (argc != 3 || !valid_goal(argv[1])) { printf("%s\n", ERROR_MSG); return 1; }
    
    file = fopen(argv[2], "r");
    if (!file) { printf("%s\n", ERROR_MSG); return 1; }
    
    if (!get_dimensions(file, &n, &d)) { fclose(file); printf("%s\n", ERROR_MSG); return 1; }
    
    points = read_points(file, n, d);
    fclose(file);
    if (!points) { printf("%s\n", ERROR_MSG); return 1; }

    result = run_goal(argv[1], points, n, d);
    if (!result) { free_matrix(points, n); printf("%s\n", ERROR_MSG); return 1; }

    print_matrix(result, n, n);
    free_matrix(result, n);
    free_matrix(points, n);
    return 0;
}
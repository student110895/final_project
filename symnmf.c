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

double **allocate_matrix(int rows, int cols)
{
    double **matrix;
    int i;
    int j;

    matrix = (double **)malloc(rows * sizeof(double *));
    if (matrix == NULL)
        return NULL;

    for (i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(cols * sizeof(double));

        if (matrix[i] == NULL) {
            for (j = 0; j < i; j++)
                free(matrix[j]);

            free(matrix);
            return NULL;
        }
    }

    return matrix;
}


void free_matrix(double **matrix, int rows)
{
    int i;

    if (matrix == NULL)
        return;

    for (i = 0; i < rows; i++)
        free(matrix[i]);

    free(matrix);
}


static double **copy_matrix(double **matrix, int rows, int cols)
{
    double **copy;
    int i;
    int j;

    copy = allocate_matrix(rows, cols);
    if (copy == NULL)
        return NULL;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++)
            copy[i][j] = matrix[i][j];
    }

    return copy;
}


static int get_dimensions(FILE *file, int *n, int *d)
{
    int c;
    int commas = 0;
    int lines = 0;
    int first_line = 1;
    int last = '\n';

    while ((c = fgetc(file)) != EOF) {
        if (first_line && c == ',')
            commas++;

        /* Only increment line count if the line actually contains data */
        if (c == '\n' && last != '\n') {
            lines++;
            first_line = 0;
        }

        last = c;
    }

    /* Account for a file that doesn't end with a newline */
    if (last != '\n')
        lines++;

    *n = lines;
    *d = commas + 1;
    rewind(file);

    return lines > 0;
}


static double **read_points(FILE *file, int n, int d)
{
    double **points;
    int i;
    int j;
    int separator;

    points = allocate_matrix(n, d);
    if (points == NULL)
        return NULL;

    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            if (fscanf(file, "%lf", &points[i][j]) != 1) {
                free_matrix(points, n);
                return NULL;
            }

            if (j < d - 1) {
                separator = fgetc(file);

                if (separator != ',') {
                    free_matrix(points, n);
                    return NULL;
                }
            }
        }
    }

    return points;
}


static double squared_distance(double *x, double *y, int d)
{
    double sum;
    double diff;
    int i;

    sum = 0.0;

    for (i = 0; i < d; i++) {
        diff = x[i] - y[i];
        sum += diff * diff;
    }

    return sum;
}


double **sym(double **points, int n, int d)
{
    double **A;
    int i;
    int j;

    A = allocate_matrix(n, n);
    if (A == NULL)
        return NULL;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j)
                A[i][j] = 0.0;
            else
                A[i][j] =
                    exp(-squared_distance(points[i], points[j], d) / 2.0);
        }
    }

    return A;
}


static double **degree_from_similarity(double **A, int n)
{
    double **D;
    double sum;
    int i;
    int j;

    D = allocate_matrix(n, n);
    if (D == NULL)
        return NULL;

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


double **ddg(double **points, int n, int d)
{
    double **A;
    double **D;

    A = sym(points, n, d);
    if (A == NULL)
        return NULL;

    D = degree_from_similarity(A, n);

    free_matrix(A, n);

    return D;
}


double **norm(double **points, int n, int d)
{
    double **A;
    double **D;
    double **W;
    int i;
    int j;

    A = sym(points, n, d);
    if (A == NULL)
        return NULL;

    D = degree_from_similarity(A, n);
    if (D == NULL) {
        free_matrix(A, n);
        return NULL;
    }

    W = allocate_matrix(n, n);
    if (W == NULL) {
        free_matrix(A, n);
        free_matrix(D, n);
        return NULL;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (D[i][i] == 0.0 || D[j][j] == 0.0)
                W[i][j] = 0.0;
            else
                W[i][j] =
                    A[i][j] /
                    (sqrt(D[i][i]) * sqrt(D[j][j]));
        }
    }

    free_matrix(A, n);
    free_matrix(D, n);

    return W;
}


static double **multiply(double **A, double **B,
                         int rows, int common, int cols)
{
    double **result;
    double sum;
    int i;
    int j;
    int k;

    result = allocate_matrix(rows, cols);
    if (result == NULL)
        return NULL;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            sum = 0.0;

            for (k = 0; k < common; k++)
                sum += A[i][k] * B[k][j];

            result[i][j] = sum;
        }
    }

    return result;
}


static double **calculate_ht_h(double **H, int n, int k)
{
    double **result;
    double sum;
    int i;
    int j;
    int row;

    result = allocate_matrix(k, k);
    if (result == NULL)
        return NULL;

    for (i = 0; i < k; i++) {
        for (j = 0; j < k; j++) {
            sum = 0.0;

            for (row = 0; row < n; row++)
                sum += H[row][i] * H[row][j];

            result[i][j] = sum;
        }
    }

    return result;
}


static double frobenius_difference(double **A, double **B,
                                   int rows, int cols)
{
    double sum;
    double diff;
    int i;
    int j;

    sum = 0.0;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            diff = A[i][j] - B[i][j];
            sum += diff * diff;
        }
    }

    return sum;
}


static double **update_h(double **H, double **W, int n, int k)
{
    double **WH;
    double **HtH;
    double **denominator;
    double **next;
    int i;
    int j;

    WH = multiply(W, H, n, n, k);
    HtH = calculate_ht_h(H, n, k);

    if (WH == NULL || HtH == NULL) {
        free_matrix(WH, n);
        free_matrix(HtH, k);
        return NULL;
    }

    denominator = multiply(H, HtH, n, k, k);
    next = allocate_matrix(n, k);

    if (denominator == NULL || next == NULL) {
        free_matrix(WH, n);
        free_matrix(HtH, k);
        free_matrix(denominator, n);
        free_matrix(next, n);
        return NULL;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            next[i][j] =
                H[i][j] *
                ((1.0 - BETA) +
                BETA * WH[i][j] /
                (denominator[i][j] + DENOM_EPSILON));
        }
    }

    free_matrix(WH, n);
    free_matrix(HtH, k);
    free_matrix(denominator, n);

    return next;
}


double **symnmf(double **H, double **W, int n, int k)
{
    double **current;
    double **next;
    double difference;
    int iteration;

    current = copy_matrix(H, n, k);
    if (current == NULL)
        return NULL;

    for (iteration = 0; iteration < MAX_ITER; iteration++) {
        next = update_h(current, W, n, k);

        if (next == NULL) {
            free_matrix(current, n);
            return NULL;
        }

        difference =
            frobenius_difference(next, current, n, k);

        free_matrix(current, n);
        current = next;

        if (difference < EPSILON)
            break;
    }

    return current;
}


static void print_matrix(double **matrix, int rows, int cols)
{
    int i;
    int j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", matrix[i][j]);

            if (j < cols - 1)
                printf(",");
        }

        printf("\n");
    }
}


static int valid_goal(char *goal)
{
    return strcmp(goal, "sym") == 0 ||
           strcmp(goal, "ddg") == 0 ||
           strcmp(goal, "norm") == 0;
}


static double **run_goal(char *goal, double **points, int n, int d)
{
    if (strcmp(goal, "sym") == 0)
        return sym(points, n, d);

    if (strcmp(goal, "ddg") == 0)
        return ddg(points, n, d);

    return norm(points, n, d);
}


int main(int argc, char **argv)
{
    FILE *file;
    double **points;
    double **result;
    int n;
    int d;

    if (argc != 3 || !valid_goal(argv[1])) {
        printf("%s\n", ERROR_MSG);
        return 1;
    }

    file = fopen(argv[2], "r");
    if (file == NULL) {
        printf("%s\n", ERROR_MSG);
        return 1;
    }

    if (!get_dimensions(file, &n, &d)) {
        fclose(file);
        printf("%s\n", ERROR_MSG);
        return 1;
    }

    points = read_points(file, n, d);
    fclose(file);

    if (points == NULL) {
        printf("%s\n", ERROR_MSG);
        return 1;
    }

    result = run_goal(argv[1], points, n, d);

    if (result == NULL) {
        free_matrix(points, n);
        printf("%s\n", ERROR_MSG);
        return 1;
    }

    print_matrix(result, n, n);

    free_matrix(result, n);
    free_matrix(points, n);

    return 0;
}

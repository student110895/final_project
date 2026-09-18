#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

#define ERROR_MSG "An Error Has Occurred"


double **allocate_matrix(int rows, int cols)
{
    double **mat;
    int i;

    mat = (double **)malloc(rows * sizeof(double *));
    if (mat == NULL)
        return NULL;

    for (i = 0; i < rows; i++) {
        mat[i] = (double *)malloc(cols * sizeof(double));
        if (mat[i] == NULL)
            return NULL;
    }

    return mat;
}


void free_matrix(double **mat, int rows)
{
    int i;

    for (i = 0; i < rows; i++)
        free(mat[i]);

    free(mat);
}


void get_dimensions(FILE *file, int *n, int *d)
{
    int c;
    int commas = 0;
    int lines = 0;
    int first_line = 1;
    int last = '\n';

    while ((c = fgetc(file)) != EOF) {
        if (first_line && c == ',')
            commas++;

        if (c == '\n') {
            lines++;
            first_line = 0;
        }

        last = c;
    }

    if (last != '\n')
        lines++;

    *n = lines;
    *d = commas + 1;
    rewind(file);
}


double **read_points(FILE *file, int n, int d)
{
    double **points;
    int i;
    int j;

    points = allocate_matrix(n, d);

    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            fscanf(file, "%lf", &points[i][j]);

            if (j < d - 1)
                fgetc(file);
        }
    }

    return points;
}


double squared_distance(double *x, double *y, int d)
{
    double sum = 0.0;
    double diff;
    int i;

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


double **ddg(double **points, int n, int d)
{
    double **A;
    double **D;
    double sum;
    int i;
    int j;

    A = sym(points, n, d);
    D = allocate_matrix(n, n);

    for (i = 0; i < n; i++) {
        sum = 0.0;

        for (j = 0; j < n; j++)
            sum += A[i][j];

        for (j = 0; j < n; j++)
            D[i][j] = 0.0;

        D[i][i] = sum;
    }

    free_matrix(A, n);

    return D;
}


void print_matrix(double **mat, int rows, int cols)
{
    int i;
    int j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", mat[i][j]);

            if (j < cols - 1)
                printf(",");
        }

        printf("\n");
    }
}


int main(int argc, char **argv)
{
    FILE *file;
    double **points;
    double **result;
    int n;
    int d;

    if (argc != 3) {
        printf("%s\n", ERROR_MSG);
        return 1;
    }

    if (strcmp(argv[1], "sym") != 0 &&
        strcmp(argv[1], "ddg") != 0) {
        printf("%s\n", ERROR_MSG);
        return 1;
    }

    file = fopen(argv[2], "r");

    if (file == NULL) {
        printf("%s\n", ERROR_MSG);
        return 1;
    }

    get_dimensions(file, &n, &d);
    points = read_points(file, n, d);
    fclose(file);

    if (strcmp(argv[1], "sym") == 0)
        result = sym(points, n, d);
    else
        result = ddg(points, n, d);

    print_matrix(result, n, n);

    free_matrix(result, n);
    free_matrix(points, n);

    return 0;
}

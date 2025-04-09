#ifndef HEADER_H
#define HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#define FLOAT 1
#define DOUBLE 2
#define ZERO 0
#define RANDOM 3
#define SID 520489042
#define A_TIMES_B 4
#define B_TIMES_C 5
#define L1_CACHE 65536

typedef struct {
    int thread_id;
    int start_row;
    int end_row;
    int start_col;
    int end_col;
    int common_dim;
    int common_dim_limit;
    float **M_1;
    float **M_2;
    float **M_3;
} ThreadDataFloat;

typedef struct {
    int thread_id;
    int start_row;
    int end_row;
    int start_col;
    int end_col;
    int common_dim;
    int common_dim_limit;
    double **M_1;
    double **M_2;
    double **M_3;
} ThreadDataDouble;

//Compares final output matrices, of type double
void compare_matrices_double(double **M_1, double **M_2, int rows, int cols);

//Compares final output matrices, of type float
void compare_matrices_float(float **M_1, float **M_2, int rows, int cols);

//Executes the multiplication of sub-matrices, of type double, in a 2D blocking algorithm, performed by a thread in parallel
void* multiply_block_double_par(void* arg);

//Executes the multiplication of sub-matrices, of type float, in a 2D blocking algorithm, performed by a thread in parallel
void* multiply_block_float_par(void* arg);

//Prints a matrix, of type double, to an output file for reference
void print_matrix_to_file_double(double **M, int rows, int cols, const char *filename);

//Prints a matrix, of type float, to an output file for reference
void print_matrix_to_file_float(float **M, int rows, int cols, const char *filename);

//Handles the multiplication of two matrices, of type double, in parallel. This oversees the 2D blocking orchestration and thread management
void multiply_matrix_double_par(double **M_1, double **M_2, double **M_3, int rows, int common_dim, int cols, int num_threads, int block_size);

//Handles the multiplication of two matrices, of type float, in parallel. This oversees the 2D blocking orchestration and thread management
void multiply_matrix_float_par(float **M_1, float **M_2, float **M_3, int rows, int common_dim, int cols, int num_threads, int block_size);

//Handles the multiplication of two matrices, of type double, in sequential order
void multiply_matrix_double_seq(double **M_1, double **M_2, double **M_3, int rows, int common_dim, int cols);

//Handles the multiplication of two matrices, of type float, in sequential order
void multiply_matrix_float_seq(float **M_1, float **M_2, float **M_3, int rows, int common_dim, int cols);

//Initialises a matrix, of type double, with either randomly generated values or 0
double ** setup_matrix_double(double **M, int rows, int cols, int array_type, int input_type, const char *filename);

//Initialises a matrix, of type float, with either randomly generated values or 0
float ** setup_matrix_float(float **M, int rows, int cols, int array_type, int input_type, const char *filename);

//Confirms if the input provided is a valid non-zero integer
bool is_valid_int(const char *buff);

#endif // HEADER_H

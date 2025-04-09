//main must generate random values based on SID
//matrices must be determined by user input, (dimensions and type)
//matrix multiplication must be manually verified as correct
//wall time must be outputted
//must intelligently determine if it is faster to multiply a*(b*c) or (a*b)*c with matrix T
//a * b * c = d

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
#include <sys/stat.h>

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

void write_double(){

}

void write_float(){
    
}

void compare_files_double(const char *filename_1, const char *filename_2) {
    FILE *file_1 = fopen(filename_1, "r");
    if (file_1 == NULL) {
        perror("Failed to open file 1");
        exit(EXIT_FAILURE);
    }

    FILE *file_2 = fopen(filename_2, "r");
    if (file_2 == NULL) {
        perror("Failed to open file 2");
        fclose(file_1);
        exit(EXIT_FAILURE);
    }

    struct stat st_1, st_2;
    if (stat(filename_1, &st_1) != 0) {
        perror("stat");
        fclose(file_1);
        fclose(file_2);
        exit(EXIT_FAILURE);
    }

    if (stat(filename_2, &st_2) != 0) {
        perror("stat");
        fclose(file_1);
        fclose(file_2);
        exit(EXIT_FAILURE);
    }

    if (st_1.st_size != st_2.st_size) {
        printf("Files are not identical (different sizes) File size of sequential is %lld, file size of parallel is %lld\n", st_1.st_size, st_2.st_size);
        fclose(file_1);
        fclose(file_2);
        return;
    }

    char *buffer_1 = (char *)malloc(st_1.st_size);
    char *buffer_2 = (char *)malloc(st_2.st_size);

    if (buffer_1 == NULL || buffer_2 == NULL) {
        perror("Failed to allocate memory");
        fclose(file_1);
        fclose(file_2);
        free(buffer_1);
        free(buffer_2);
        exit(EXIT_FAILURE);
    }

    fread(buffer_1, 1, st_1.st_size, file_1);
    fread(buffer_2, 1, st_2.st_size, file_2);

    bool identical = true;

    double epsilon = 1.0E-10;

    for (size_t i = 0; i< st_1.st_size / sizeof(double); i++){
        float val_1 = ((double *)buffer_1)[i];
        float val_2 = ((double *)buffer_2)[i];
        if (fabs(val_1 - val_2) > epsilon){
            identical = false;
            printf("Files are not identical at element %zu\n", i);
            break;
        }
    }

    if (identical){
        printf("Files are identical\n");
    } else {
        printf("Files are not identical\n");
    }

    fclose(file_1);
    fclose(file_2);
    free(buffer_1);
    free(buffer_2);
}

void compare_files_float(const char *filename_1, const char *filename_2) {
    printf("Testing\n");
    FILE *file_1 = fopen(filename_1, "r");
    if (file_1 == NULL) {
        perror("Failed to open file 1");
        exit(EXIT_FAILURE);
    }

    FILE *file_2 = fopen(filename_2, "r");
    if (file_2 == NULL) {
        perror("Failed to open file 2");
        fclose(file_1);
        exit(EXIT_FAILURE);
    }

    struct stat st_1, st_2;
    if (stat(filename_1, &st_1) != 0) {
        perror("stat");
        fclose(file_1);
        fclose(file_2);
        exit(EXIT_FAILURE);
    }

    if (stat(filename_2, &st_2) != 0) {
        perror("stat");
        fclose(file_1);
        fclose(file_2);
        exit(EXIT_FAILURE);
    }

    if (st_1.st_size != st_2.st_size) {
        printf("Files are not identical (different sizes) File size of sequential is %lld, file size of parallel is %lld\n", st_1.st_size, st_2.st_size);
        fclose(file_1);
        fclose(file_2);
        return;
    }

    char *buffer_1 = (char *)malloc(st_1.st_size);
    char *buffer_2 = (char *)malloc(st_2.st_size);

    if (buffer_1 == NULL || buffer_2 == NULL) {
        perror("Failed to allocate memory");
        fclose(file_1);
        fclose(file_2);
        free(buffer_1);
        free(buffer_2);
        exit(EXIT_FAILURE);
    }

    fread(buffer_1, 1, st_1.st_size, file_1);
    fread(buffer_2, 1, st_2.st_size, file_2);

    bool identical = true;
    float epsilon = 1.0E-7;

    for (size_t i = 0; i< st_1.st_size / sizeof(float); i++){
        float val_1 = ((float *)buffer_1)[i];
        float val_2 = ((float *)buffer_2)[i];
        if (fabs(val_1 - val_2) > epsilon){
            identical = false;
            printf("Files are not identical at element %zu\n", i);
            break;
        }
    }

    if (identical == true){
        printf("Files are identical\n");
    } else {
        printf("Files are not identical\n");
    }

    fclose(file_1);
    fclose(file_2);
    free(buffer_1);
    free(buffer_2);
}

void multiply_matrix_double_seq(double **M_1, double **M_2, double **M_3, int rows, int common_dim, int cols){
    for (int i = 0; i < rows; i++){
        for (int k = 0; k < common_dim; k++){
            for (int j = 0; j < cols; j++){
                M_3[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
}

void multiply_matrix_float_seq(float **M_1, float **M_2, float **M_3, int rows, int common_dim, int cols){
    for (int i = 0; i < rows; i++){
        for (int k = 0; k < common_dim; k++){
            for (int j = 0; j < cols; j++){
                M_3[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
}

void* multiply_block_float_par(void* arg) {
    ThreadDataFloat* data = (ThreadDataFloat*)arg;
    for (int i = data->start_row; i < data->end_row; i++) {
        for (int j = data->start_col; j < data->end_col; j+=4) {
            for (int k = 0; k < data->common_dim_limit; k+=4){
                data->M_3[i][j] += data->M_1[i][k] * data->M_2[k][j];
                data->M_3[i][j+1] += data->M_1[i][k] * data->M_2[k][j+1];
                data->M_3[i][j+2] += data->M_1[i][k] * data->M_2[k][j+2];
                data->M_3[i][j+3] += data->M_1[i][k] * data->M_2[k][j+3];

                data->M_3[i][j] += data->M_1[i][k+1] * data->M_2[k+1][j];
                data->M_3[i][j+1] += data->M_1[i][k+1] * data->M_2[k+1][j+1];
                data->M_3[i][j+2] += data->M_1[i][k+1] * data->M_2[k+1][j+2];
                data->M_3[i][j+3] += data->M_1[i][k+1] * data->M_2[k+1][j+3];

                data->M_3[i][j] += data->M_1[i][k+2] * data->M_2[k+2][j];
                data->M_3[i][j+1] += data->M_1[i][k+2] * data->M_2[k+2][j+1];
                data->M_3[i][j+2] += data->M_1[i][k+2] * data->M_2[k+2][j+2];
                data->M_3[i][j+3] += data->M_1[i][k+2] * data->M_2[k+2][j+3];

                data->M_3[i][j] += data->M_1[i][k+3] * data->M_2[k+3][j];
                data->M_3[i][j+1] += data->M_1[i][k+3] * data->M_2[k+3][j+1];
                data->M_3[i][j+2] += data->M_1[i][k+3] * data->M_2[k+3][j+2];
                data->M_3[i][j+3] += data->M_1[i][k+3] * data->M_2[k+3][j+3];
            }
            for (int k = data->common_dim_limit; k < data->common_dim; k++){
                data->M_3[i][j] += data->M_1[i][k] * data->M_2[k][j];
                data->M_3[i][j+1] += data->M_1[i][k] * data->M_2[k][j+1];
                data->M_3[i][j+2] += data->M_1[i][k] * data->M_2[k][j+2];
                data->M_3[i][j+3] += data->M_1[i][k] * data->M_2[k][j+3];
            }
            // for (int k = 0; k < data->common_dim_limit; k+=4){
            //     sum0 += data->M_1[i][k] * data->M_2[k][j] +
            //             data->M_1[i][k+1] * data->M_2[k+1][j] +
            //             data->M_1[i][k+2] * data->M_2[k+2][j] +
            //             data->M_1[i][k+3] * data->M_2[k+3][j];

            //     sum1 += data->M_1[i][k] * data->M_2[k][j+1] +
            //            data->M_1[i][k+1] * data->M_2[k+1][j+1] +
            //            data->M_1[i][k+2] * data->M_2[k+2][j+1] +
            //            data->M_1[i][k+3] * data->M_2[k+3][j+1];

            //     sum2 += data->M_1[i][k] * data->M_2[k][j+2] +
            //            data->M_1[i][k+1] * data->M_2[k+1][j+2] +
            //            data->M_1[i][k+2] * data->M_2[k+2][j+2] +
            //            data->M_1[i][k+3] * data->M_2[k+3][j+2];

            //     sum3 += data->M_1[i][k] * data->M_2[k][j+3] +
            //            data->M_1[i][k+1] * data->M_2[k+1][j+3] +
            //            data->M_1[i][k+2] * data->M_2[k+2][j+3] +
            //            data->M_1[i][k+3] * data->M_2[k+3][j+3];
            // }

            // // Handle remaining k iterations
            // for (int k = data->common_dim_limit; k < data->common_dim; k++) {
            //     sum0 += data->M_1[i][k] * data->M_2[k][j];
            //     sum1 += data->M_1[i][k] * data->M_2[k][j+1];
            //     sum2 += data->M_1[i][k] * data->M_2[k][j+2];
            //     sum3 += data->M_1[i][k] * data->M_2[k][j+3];
            // }
        }
    }
    return NULL;
}

void* multiply_block_double_par(void* arg) {
    ThreadDataDouble* data = (ThreadDataDouble*)arg;
    for (int i = data->start_row; i < data->end_row; i++) {
        for (int j = data->start_col; j < data->end_col; j+=4) {
            for (int k = 0; k < data->common_dim_limit; k+=4){
                data->M_3[i][j] += data->M_1[i][k] * data->M_2[k][j];
                data->M_3[i][j+1] += data->M_1[i][k] * data->M_2[k][j+1];
                data->M_3[i][j+2] += data->M_1[i][k] * data->M_2[k][j+2];
                data->M_3[i][j+3] += data->M_1[i][k] * data->M_2[k][j+3];

                data->M_3[i][j] += data->M_1[i][k+1] * data->M_2[k+1][j];
                data->M_3[i][j+1] += data->M_1[i][k+1] * data->M_2[k+1][j+1];
                data->M_3[i][j+2] += data->M_1[i][k+1] * data->M_2[k+1][j+2];
                data->M_3[i][j+3] += data->M_1[i][k+1] * data->M_2[k+1][j+3];

                data->M_3[i][j] += data->M_1[i][k+2] * data->M_2[k+2][j];
                data->M_3[i][j+1] += data->M_1[i][k+2] * data->M_2[k+2][j+1];
                data->M_3[i][j+2] += data->M_1[i][k+2] * data->M_2[k+2][j+2];
                data->M_3[i][j+3] += data->M_1[i][k+2] * data->M_2[k+2][j+3];

                data->M_3[i][j] += data->M_1[i][k+3] * data->M_2[k+3][j];
                data->M_3[i][j+1] += data->M_1[i][k+3] * data->M_2[k+3][j+1];
                data->M_3[i][j+2] += data->M_1[i][k+3] * data->M_2[k+3][j+2];
                data->M_3[i][j+3] += data->M_1[i][k+3] * data->M_2[k+3][j+3];
            }
            for (int k = data->common_dim_limit; k < data->common_dim; k++){
                data->M_3[i][j] += data->M_1[i][k] * data->M_2[k][j];
                data->M_3[i][j+1] += data->M_1[i][k] * data->M_2[k][j+1];
                data->M_3[i][j+2] += data->M_1[i][k] * data->M_2[k][j+2];
                data->M_3[i][j+3] += data->M_1[i][k] * data->M_2[k][j+3];
            }    
            // for (int k = 0; k < data->common_dim_limit; k+=4){
            //     sum0 += data->M_1[i][k] * data->M_2[k][j] +
            //             data->M_1[i][k+1] * data->M_2[k+1][j] +
            //             data->M_1[i][k+2] * data->M_2[k+2][j] +
            //             data->M_1[i][k+3] * data->M_2[k+3][j];

            //     sum1 += data->M_1[i][k] * data->M_2[k][j+1] +
            //            data->M_1[i][k+1] * data->M_2[k+1][j+1] +
            //            data->M_1[i][k+2] * data->M_2[k+2][j+1] +
            //            data->M_1[i][k+3] * data->M_2[k+3][j+1];

            //     sum2 += data->M_1[i][k] * data->M_2[k][j+2] +
            //            data->M_1[i][k+1] * data->M_2[k+1][j+2] +
            //            data->M_1[i][k+2] * data->M_2[k+2][j+2] +
            //            data->M_1[i][k+3] * data->M_2[k+3][j+2];

            //     sum3 += data->M_1[i][k] * data->M_2[k][j+3] +
            //            data->M_1[i][k+1] * data->M_2[k+1][j+3] +
            //            data->M_1[i][k+2] * data->M_2[k+2][j+3] +
            //            data->M_1[i][k+3] * data->M_2[k+3][j+3];
            // }

            // // Handle remaining k iterations
            // for (int k = data->common_dim_limit; k < data->common_dim; k++) {
            //     sum0 += data->M_1[i][k] * data->M_2[k][j];
            //     sum1 += data->M_1[i][k] * data->M_2[k][j+1];
            //     sum2 += data->M_1[i][k] * data->M_2[k][j+2];
            //     sum3 += data->M_1[i][k] * data->M_2[k][j+3];
            // }

            // data->M_3[i][j] = sum0;
            // data->M_3[i][j+1] = sum1;
            // data->M_3[i][j+2] = sum2;
            // data->M_3[i][j+3] = sum3;
        }
    }
    return NULL;
}


void print_matrix_to_file_double(double **M, int rows, int cols, const char *filename){
    printf("Printing matrix of dimensions %d x %d to file %s\n", rows, cols, filename);
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file");
        exit(1);
    }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            fprintf(file, "%.12lf ", M[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void print_matrix_to_file_float(float **M, int rows, int cols, const char *filename){
    printf("Printing matrix of dimensions %d x %d to file %s\n", rows, cols, filename);
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file");
        exit(1);
    }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            fprintf(file, "%.6f ", M[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void multiply_matrix_double_par(double **M_1, double **M_2, double **M_3, int rows, int common_dim, int cols, int num_threads, int block_size){
    pthread_t threads[num_threads];
    ThreadDataDouble thread_data[num_threads];

    // Partition the matrix into blocks and assign to worker threads
    int block_rows = floor(rows / block_size);
    int block_cols = floor(cols / block_size);
    int thread_index = 0;
    int common_dim_limit = floor(common_dim/4) * 4;

    bool firstIteration = true;


    for (int i = 0; i < block_rows; i++) {
        for (int j = 0; j < block_cols; j++) {
            if (thread_index == num_threads){
                if (firstIteration){
                    firstIteration = false;
                }
                thread_index = 0;
            }
            if (!firstIteration){
                //waiting for it to join before continuing
                pthread_join(threads[thread_index], NULL);
            }
            thread_data[thread_index].thread_id = thread_index;
            thread_data[thread_index].start_row = i * block_size;
            thread_data[thread_index].end_row = (i + 1) * block_size;
            thread_data[thread_index].start_col = j * block_size;
            thread_data[thread_index].end_col = (j + 1) * block_size;
            thread_data[thread_index].common_dim = common_dim;
            thread_data[thread_index].common_dim_limit = common_dim_limit;
            thread_data[thread_index].M_1 = M_1;
            thread_data[thread_index].M_2 = M_2;
            thread_data[thread_index].M_3 = M_3;
            pthread_create(&threads[thread_index], NULL, multiply_block_double_par, &thread_data[thread_index]);
            // pthread_join(threads[thread_index], NULL);
            thread_index++;
        }
    }

    // Wait for all threads to complete
    for (int i = 0; i < thread_index; i++) {
        pthread_join(threads[i], NULL);
    }

    // Process the leftover rows and columns
    for (int i = block_rows * block_size; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < common_dim; k++){
                M_3[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
    for (int j = block_cols * block_size; j < cols; j++) {
        for (int i = 0; i < block_rows * block_size; i++) {
            for (int k = 0; k < common_dim; k++){
                M_3[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
}

void multiply_matrix_float_par(float **M_1, float **M_2, float **M_3, int rows, int common_dim, int cols, int num_threads, int block_size){
    pthread_t threads[num_threads];
    ThreadDataFloat thread_data[num_threads];

    // Partition the matrix into blocks and assign to worker threads
    int block_rows = floor(rows / block_size);
    int block_cols = floor(cols / block_size);
    int thread_index = 0;
    int common_dim_limit = floor(common_dim/4) * 4;

    bool firstIteration = true;

    for (int i = 0; i < block_rows; i++) {
        for (int j = 0; j < block_cols; j++) {
            if (thread_index == num_threads){
                if (firstIteration){
                    firstIteration = false;
                }
                thread_index = 0;
            }
            if (!firstIteration){
                pthread_join(threads[thread_index], NULL);
            }
            thread_data[thread_index].thread_id = thread_index;
            thread_data[thread_index].start_row = i * block_size;
            thread_data[thread_index].end_row = (i + 1) * block_size;
            thread_data[thread_index].start_col = j * block_size;
            thread_data[thread_index].end_col = (j + 1) * block_size;
            thread_data[thread_index].common_dim = common_dim;
            thread_data[thread_index].common_dim_limit = common_dim_limit;
            thread_data[thread_index].M_1 = M_1;
            thread_data[thread_index].M_2 = M_2;
            thread_data[thread_index].M_3 = M_3;
            pthread_create(&threads[thread_index], NULL, multiply_block_float_par, &thread_data[thread_index]);
            // pthread_join(threads[thread_index], NULL);
            thread_index++;
        }
    }

    // Wait for all worker threads to complete
    for (int i = 0; i < thread_index; i++) {
        pthread_join(threads[i], NULL);
    }

    // Process the leftover rows and columns
    for (int i = block_rows * block_size; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < common_dim; k++){
                M_3[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
    for (int j = block_cols * block_size; j < cols; j++) {
        for (int i = 0; i < block_rows * block_size; i++) {
            for (int k = 0; k < common_dim; k++){
                M_3[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
}

double ** setup_matrix_double(double **M, int rows, int cols, int array_type, int input_type, const char *filename){
    M = (double **)malloc(rows * sizeof(double *));

    // allocate memory to columns
    for (int i = 0; i < rows; i++){
        M[i] = (double *)malloc(cols * sizeof(double));
    }

    if (input_type == ZERO){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                M[i][j] = 0.0;
            }
        }  
    } else if (input_type == RANDOM){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                M[i][j] = (double) rand() / RAND_MAX;
            }
        }
    }
    print_matrix_to_file_double(M, rows, cols, filename);
        
    return M;
}

float ** setup_matrix_float(float **M, int rows, int cols, int array_type, int input_type, const char *filename){
    M = (float **)malloc(rows * sizeof(float *));

    // allocate memory to columns
    for (int i = 0; i < rows; i++){
        M[i] = (float *)malloc(cols * sizeof(float));
    }

    if (input_type == ZERO){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                M[i][j] = 0.0f;
            }
        }  
    } else if (input_type == RANDOM){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                M[i][j] = (float) rand() / RAND_MAX;
            }
        }
    }

    print_matrix_to_file_float(M, rows, cols, filename);
        
    return M;
}

bool is_valid_int(const char *buff) {
    char *end;
    bool is_valid_int = false;
    errno = 0;
 
    const long int_to_test = strtol(buff, &end, 10);
 
    if (end == buff) {
      (void) fprintf(stderr, "%s: not a decimal number\n", buff);
    } else if ('\0' != *end) {
      (void) fprintf(stderr, "%s: extra characters at end of input: %s\n", buff, end);
    } else if ((LONG_MIN == int_to_test || LONG_MAX == int_to_test) && ERANGE == errno) {
      (void) fprintf(stderr, "%s out of range of type long\n", buff);
    } else if (int_to_test > INT_MAX) {
      (void) fprintf(stderr, "%ld greater than INT_MAX\n", int_to_test);
    } else if (int_to_test < 0) {
      (void) fprintf(stderr, "%ld less than 0\n", int_to_test);
    } else {
      is_valid_int = true;
    }
    return is_valid_int;
}

int main(int argc, char *argv[]){
    if (argc != 7) {
        printf("Usage: %s <m> <k> <l> <n> <type> <num_threads>\n", argv[0]);
        return 1;
    }
    int m, k, l, n, num_threads = 0;

    if (is_valid_int(argv[1]) && is_valid_int(argv[2]) && is_valid_int(argv[3]) && is_valid_int(argv[4])){
        m = atoi(argv[1]); 
        k = atoi(argv[2]); 
        l = atoi(argv[3]); 
        n = atoi(argv[4]); 
    } else {
        (void) fprintf(stderr, "Exiting due to invalid command-line input\n");
        exit(EXIT_FAILURE);
    }

    char *type = argv[5];
    int array_type;

    if (!isalpha(type[0])){
        (void) fprintf(stderr, "Exiting due to invalid command-line input\n");
        exit(EXIT_FAILURE);
    } else if (!(strcmp(type, "double") == 0 || strcmp(type, "float") == 0)){
        (void) fprintf(stderr, "Exiting due to invalid command-line input\n");
        exit(EXIT_FAILURE);
    } else if (strcmp(type, "double") == 0){
        array_type = DOUBLE;
    } else if (strcmp(type, "float") == 0){
        array_type = FLOAT;
    }

    if (is_valid_int(argv[6])){
        if (atoi(argv[6]) > 64){
            (void) fprintf(stderr, "Exiting due to invalid command-line input\n");
            exit(EXIT_FAILURE);
        } else {
            num_threads = atoi(argv[6]); 
        }
    }

    int t_rows, t_cols;
    int t_type;
    //determining optimal T by comparing the number of scalar multiplications which will occur
    long long first_case = (long long)m*k*l+(long long)m*l*n;
    long long second_case = (long long)k*l*n+(long long)m*k*n;
    if (first_case < second_case){
        // in this case, (A * B) * C is more efficient
        t_rows = m;
        t_cols = l;
        t_type = A_TIMES_B;
    } else {
        // otherwise, A * (B * C) is more efficient
        t_rows = k;
        t_cols = n;
        t_type = B_TIMES_C;
    }
    srand(SID);

    if (array_type == FLOAT){

        int tempBlockSize = floor(sqrt(L1_CACHE/(3 * sizeof(float))));
        int blockSize = floor(tempBlockSize/4) * 4;

        float **A = setup_matrix_float(A, m, k, array_type, RANDOM, "A");
        float **B = setup_matrix_float(B, k, l, array_type, RANDOM, "B");
        float **C = setup_matrix_float(C, l, n, array_type, RANDOM, "C");
        float **D_seq = setup_matrix_float(D_seq, m, n, array_type, ZERO, "D_seq");
        float **T_seq = setup_matrix_float(T_seq, t_rows, t_cols, array_type, ZERO, "T_seq");
        float **D_par = setup_matrix_float(D_par, m, n, array_type, ZERO, "D_par");
        float **T_par = setup_matrix_float(T_par, t_rows, t_cols, array_type, ZERO, "T_par");

        struct timeval start_seq, end_seq;
        gettimeofday(&start_seq, NULL);
        if (t_type == A_TIMES_B){
            multiply_matrix_float_seq(A, B, T_seq, m, k, l);
            multiply_matrix_float_seq(T_seq, C, D_seq, t_rows, t_cols, n);
        } else if (t_type == B_TIMES_C){
            multiply_matrix_float_seq(B, C, T_seq, k, l, n);
            multiply_matrix_float_seq(A, T_seq, D_seq, m, t_rows, t_cols);
        }
        gettimeofday(&end_seq, NULL);
        long seconds_seq = end_seq.tv_sec - start_seq.tv_sec;
        long microseconds_seq = end_seq.tv_usec - start_seq.tv_usec;
        double elapsed_seq = seconds_seq + microseconds_seq * 1e-6;
        printf("Sequential Elapsed time: %.6f seconds\n", elapsed_seq);

        printf("Printing sequential matrices to file\n");
        print_matrix_to_file_float(T_seq, t_rows, t_cols, "T_seq");
        print_matrix_to_file_float(D_seq, m, n, "D_seq");

        struct timeval start_par, end_par;
        gettimeofday(&start_par, NULL);
        if (t_type == A_TIMES_B){
            multiply_matrix_float_par(A, B, T_par, m, k, l, num_threads, blockSize);
            multiply_matrix_float_par(T_par, C, D_par, t_rows, t_cols, n, num_threads, blockSize);
        } else if (t_type == B_TIMES_C){
            multiply_matrix_float_par(B, C, T_par, k, l, n, num_threads, blockSize);
            multiply_matrix_float_par(A, T_par, D_par, m, t_rows, t_cols, num_threads, blockSize);
        }
        gettimeofday(&end_par, NULL);
        long seconds_par = end_par.tv_sec - start_par.tv_sec;
        long microseconds_par = end_par.tv_usec - start_par.tv_usec;
        double elapsed_par = seconds_par + microseconds_par * 1e-6;
        printf("Parallel Elapsed time: %.6f seconds\n", elapsed_par);

        printf("Printing parallel matrices to file\n");
        printf("T_par\n");
        print_matrix_to_file_float(T_par, t_rows, t_cols, "T_par");
        printf("D_par\n");
        print_matrix_to_file_float(D_par, m, n, "D_par");

        bool identical = true;
        float epsilon = 1.0E-7;
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                if (fabs(D_par[i][j] - D_seq[i][j]) > epsilon){
                    printf("Mistmatch on index [%d][%d]\n", i, j);
                    identical = false;
                    break;
                }
            }
        }
        if (identical){
            printf("Outputs are identical\n");
        } else {
            printf("Outputs are not identical\n");
        }

        // printf("Comparing T\n");
        // compare_files_float("T_seq", "T_par");
        // printf("Comparing D\n");
        // compare_files_float("D_seq", "D_par");

        double speedup = elapsed_seq - elapsed_par;
        printf("Speedup: %.6f seconds\n", speedup);

        // Free allocated memory
        for (int i = 0; i < m; i++) {
            free(A[i]);
            free(D_par[i]);
            free(D_seq[i]);
        }
        for (int i = 0; i < k; i++) {
            free(B[i]);
        }
        for (int i = 0; i < l; i++) {
            free(C[i]);
        }
        for (int i = 0; i < t_rows; i++){
            free(T_par[i]);
            free(T_seq[i]);
        }
        free(A);
        free(B);
        free(C);
        free(D_seq);
        free(T_seq);
        free(D_par);
        free(T_par);

    } else if (array_type == DOUBLE){

        int tempBlockSize = floor(sqrt(L1_CACHE/(3 * sizeof(double))));
        int blockSize = floor(tempBlockSize/4) * 4;

        double **A = setup_matrix_double(A, m, k, array_type, RANDOM, "A");
        double **B = setup_matrix_double(B, k, l, array_type, RANDOM, "B");
        double **C = setup_matrix_double(C, l, n, array_type, RANDOM, "C");
        double **D_seq = setup_matrix_double(D_seq, m, n, array_type, ZERO, "D_seq");
        double **T_seq = setup_matrix_double(T_seq, t_rows, t_cols, array_type, ZERO, "T_seq");
        double **D_par = setup_matrix_double(D_par, m, n, array_type, ZERO, "D_par");
        double **T_par = setup_matrix_double(T_par, t_rows, t_cols, array_type, ZERO, "T_par");

        struct timeval start_seq, end_seq;
        gettimeofday(&start_seq, NULL);
        if (t_type == A_TIMES_B){
            multiply_matrix_double_seq(A, B, T_seq, m, k, l);
            multiply_matrix_double_seq(T_seq, C, D_seq, t_rows, t_cols, n);
        } else if (t_type == B_TIMES_C){
            multiply_matrix_double_seq(B, C, T_seq, k, l, n);
            multiply_matrix_double_seq(A, T_seq, D_seq, m, t_rows, t_cols);
        }
        gettimeofday(&end_seq, NULL);
        long seconds_seq = end_seq.tv_sec - start_seq.tv_sec;
        long microseconds_seq = end_seq.tv_usec - start_seq.tv_usec;
        double elapsed_seq = seconds_seq + microseconds_seq * 1e-6;
        printf("Sequential Elapsed time: %.6f seconds\n", elapsed_seq);

        printf("Printing sequential matrices to file\n");
        print_matrix_to_file_double(T_seq, t_rows, t_cols, "T_seq");
        print_matrix_to_file_double(D_seq, m, n, "D_seq");

        struct timeval start_par, end_par;
        gettimeofday(&start_par, NULL);
        if (t_type == A_TIMES_B){
            multiply_matrix_double_par(A, B, T_par, m, k, l, num_threads, blockSize);
            multiply_matrix_double_par(T_par, C, D_par, t_rows, t_cols, n, num_threads, blockSize);
        } else if (t_type == B_TIMES_C){
            multiply_matrix_double_par(B, C, T_par, k, l, n, num_threads, blockSize);
            multiply_matrix_double_par(A, T_par, D_par, m, t_rows, t_cols, num_threads, blockSize);
        }
        gettimeofday(&end_par, NULL);
        long seconds_par = end_par.tv_sec - start_par.tv_sec;
        long microseconds_par = end_par.tv_usec - start_par.tv_usec;
        double elapsed_par = seconds_par + microseconds_par * 1e-6;
        printf("Parallel Elapsed time: %.6f seconds\n", elapsed_par);

        printf("Printing parallel matrices to file\n");
        printf("T_par\n");
        print_matrix_to_file_double(T_par, t_rows, t_cols, "T_par");
        printf("D_par\n");
        print_matrix_to_file_double(D_par, m, n, "D_par");

        bool identical = true;
        double epsilon = 1.0E-10;
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                if (fabs(D_par[i][j] - D_seq[i][j]) > epsilon){
                    printf("Mistmatch on index [%d][%d]\n", i, j);
                    identical = false;
                    break;
                }
            }
        }
        if (identical){
            printf("Outputs are identical\n");
        } else {
            printf("Outputs are not identical\n");
        }
        double speedup = elapsed_seq - elapsed_par;
        printf("Speedup: %.6f seconds\n", speedup);

        // Free allocated memory
        for (int i = 0; i < m; i++) {
            free(A[i]);
            free(D_seq[i]);
            free(D_par[i]);
        }
        for (int i = 0; i < k; i++) {
            free(B[i]);
        }
        for (int i = 0; i < l; i++) {
            free(C[i]);
        }
        for (int i = 0; i < t_rows; i++){
            free(T_seq[i]);
            free(T_par[i]);
        }
        free(A);
        free(B);
        free(C);
        free(D_seq);
        free(T_seq);
        free(D_par);
        free(T_par);
    }
}
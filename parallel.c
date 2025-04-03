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

#define FLOAT 1
#define DOUBLE 2
#define ZERO 0
#define RANDOM 3
#define SID 520489042
#define A_TIMES_B 4
#define B_TIMES_C 5
#define L1_CACHE 64
#define LINE_SIZE 128

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

void* multiply_block_float(void* arg) {
    ThreadDataFloat* data = (ThreadDataFloat*)arg;
    printf("Thread %d - Multiply block for start_row: %d end_row: %d\nstart_col: %d end_col: %d\n", data->thread_id, data->start_row, data->end_row, data->start_col, data->end_col);
    // for (int i = data->start_row; i < data->end_row; i++) {
    //     for (int j = data->start_col; j < data->end_col; j++) {
    //         for (int k = 0; k < data->common_dim; k++){
    //             data->M_3[i][j] += data->M_1[i][k] * data->M_2[k][j];
    //             // printf("M_3[%d][%d] += %f * %f]\n", i, j, data->M_1[i][k], data->M_2[k][j]);
    //         }
    //     }
    // }
    for (int i = data->start_row; i < data->end_row; i++) {
        for (int j = data->start_col; j < data->end_col; j+=4) {
            int j_1 = j;
            int j_2 = j+1;
            int j_3 = j+2;
            int j_4 = j+3;
            // for (int k = 0; k < data->common_dim_limit; k+=4){
            //     int k_1 = k;
            //     int k_2 = k+1;
            //     int k_3 = k+2;
            //     int k_4 = k+3;

            //     data->M_3[i][j_1] += data->M_1[i][k_1] * data->M_2[k_1][j_1];
            //     data->M_3[i][j_1] += data->M_1[i][k_2] * data->M_2[k_2][j_1];
            //     data->M_3[i][j_1] += data->M_1[i][k_3] * data->M_2[k_3][j_1];
            //     data->M_3[i][j_1] += data->M_1[i][k_4] * data->M_2[k_4][j_1];

            //     data->M_3[i][j_2] += data->M_1[i][k_1] * data->M_2[k_1][j_2];
            //     data->M_3[i][j_2] += data->M_1[i][k_2] * data->M_2[k_2][j_2];
            //     data->M_3[i][j_2] += data->M_1[i][k_3] * data->M_2[k_3][j_2];
            //     data->M_3[i][j_2] += data->M_1[i][k_4] * data->M_2[k_4][j_2];

            //     data->M_3[i][j_3] += data->M_1[i][k_1] * data->M_2[k_1][j_3];
            //     data->M_3[i][j_3] += data->M_1[i][k_2] * data->M_2[k_2][j_3];
            //     data->M_3[i][j_3] += data->M_1[i][k_3] * data->M_2[k_3][j_3];
            //     data->M_3[i][j_3] += data->M_1[i][k_4] * data->M_2[k_4][j_3];

            //     data->M_3[i][j_4] += data->M_1[i][k_1] * data->M_2[k_1][j_4];
            //     data->M_3[i][j_4] += data->M_1[i][k_2] * data->M_2[k_2][j_4];
            //     data->M_3[i][j_4] += data->M_1[i][k_3] * data->M_2[k_3][j_4];
            //     data->M_3[i][j_4] += data->M_1[i][k_4] * data->M_2[k_4][j_4];
            //     // data->M_3[i][j] += data->M_1[i][k] * data->M_2[k][j];
            //     // printf("M_3[%d][%d] += %f * %f]\n", i, j, data->M_1[i][k], data->M_2[k][j]);
            // }
            for (int k = 0; k < data->common_dim; k++){
                data->M_3[i][j_1] += data->M_1[i][k] * data->M_2[k][j_1];
                data->M_3[i][j_2] += data->M_1[i][k] * data->M_2[k][j_2];
                data->M_3[i][j_3] += data->M_1[i][k] * data->M_2[k][j_3];
                data->M_3[i][j_4] += data->M_1[i][k] * data->M_2[k][j_4];
            }
        }
    }
    return NULL;
}

void* multiply_block_double(void* arg) {
    ThreadDataDouble* data = (ThreadDataDouble*)arg;
    printf("Thread %d - Multiply block for start_row: %d end_row: %d\nstart_col: %d end_col: %d\n", data->thread_id, data->start_row, data->end_row, data->start_col, data->end_col);
    for (int i = data->start_row; i < data->end_row; i++) {
        for (int j = data->start_col; j < data->end_col; j+=4) {
            int j_1 = j;
            int j_2 = j+1;
            int j_3 = j+2;
            int j_4 = j+3;
            // for (int k = 0; k < data->common_dim_limit; k+=4){
            //     int k_1 = k;
            //     int k_2 = k+1;
            //     int k_3 = k+2;
            //     int k_4 = k+3;

            //     data->M_3[i][j_1] += data->M_1[i][k_1] * data->M_2[k_1][j_1];
            //     data->M_3[i][j_1] += data->M_1[i][k_2] * data->M_2[k_2][j_1];
            //     data->M_3[i][j_1] += data->M_1[i][k_3] * data->M_2[k_3][j_1];
            //     data->M_3[i][j_1] += data->M_1[i][k_4] * data->M_2[k_4][j_1];

            //     data->M_3[i][j_2] += data->M_1[i][k_1] * data->M_2[k_1][j_2];
            //     data->M_3[i][j_2] += data->M_1[i][k_2] * data->M_2[k_2][j_2];
            //     data->M_3[i][j_2] += data->M_1[i][k_3] * data->M_2[k_3][j_2];
            //     data->M_3[i][j_2] += data->M_1[i][k_4] * data->M_2[k_4][j_2];

            //     data->M_3[i][j_3] += data->M_1[i][k_1] * data->M_2[k_1][j_3];
            //     data->M_3[i][j_3] += data->M_1[i][k_2] * data->M_2[k_2][j_3];
            //     data->M_3[i][j_3] += data->M_1[i][k_3] * data->M_2[k_3][j_3];
            //     data->M_3[i][j_3] += data->M_1[i][k_4] * data->M_2[k_4][j_3];

            //     data->M_3[i][j_4] += data->M_1[i][k_1] * data->M_2[k_1][j_4];
            //     data->M_3[i][j_4] += data->M_1[i][k_2] * data->M_2[k_2][j_4];
            //     data->M_3[i][j_4] += data->M_1[i][k_3] * data->M_2[k_3][j_4];
            //     data->M_3[i][j_4] += data->M_1[i][k_4] * data->M_2[k_4][j_4];
            //     // data->M_3[i][j] += data->M_1[i][k] * data->M_2[k][j];
            //     // printf("M_3[%d][%d] += %f * %f]\n", i, j, data->M_1[i][k], data->M_2[k][j]);
            // }
            for (int k = 0; k < data->common_dim; k++){
                data->M_3[i][j_1] += data->M_1[i][k] * data->M_2[k][j_1];
                data->M_3[i][j_2] += data->M_1[i][k] * data->M_2[k][j_2];
                data->M_3[i][j_3] += data->M_1[i][k] * data->M_2[k][j_3];
                data->M_3[i][j_4] += data->M_1[i][k] * data->M_2[k][j_4];
            }
        }
    }
    return NULL;
}

void print_matrix_to_file_double(double **M, int rows, int cols, const char *filename){
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file");
        exit(1);
    }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            fprintf(file, "%lf ", M[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void print_matrix_to_file_float(float **M, int rows, int cols, const char *filename){
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file");
        exit(1);
    }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            fprintf(file, "%lf ", M[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void multiply_matrix_double(double **M_1, double **M_2, double **M_3, int rows, int common_dim, int cols, const char *filename, int num_threads, int block_size){
    pthread_t threads[num_threads];
    ThreadDataDouble thread_data[num_threads];

    // Partition the matrix into blocks and assign to worker threads
    int block_rows = floor(rows / block_size);
    int block_cols = floor(cols / block_size);
    int thread_index = 0;
    int common_dim_limit = floor(common_dim/4);

    // struct timeval start, end;
    // gettimeofday(&start, NULL);

    for (int i = 0; i < block_rows; i++) {
        for (int j = 0; j < block_cols; j++) {
            if (thread_index == num_threads){
                thread_index = 0;
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
            pthread_create(&threads[thread_index], NULL, multiply_block_float, &thread_data[thread_index]);
            thread_index++;
        }
    }

    // Wait for all worker threads to complete
    for (int i = 0; i < thread_index; i++) {
        pthread_join(threads[i], NULL);
    }

    // Process the leftover rows and columns
    printf("Cleaning up\n");
    printf("Multiply block for start_row: %d end_row: %d\nstart_col: %d end_col: %d\n", block_rows * block_size, rows-1, 0, cols-1);
    for (int i = block_rows * block_size; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < common_dim; k++){
                M_3[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
    printf("Multiply block for start_row: %d end_row: %d\nstart_col: %d end_col: %d\n", 0, rows-1, block_cols * block_size, cols-1);
    for (int j = block_cols * block_size; j < cols; j++) {
        for (int i = 0; i < block_rows * block_size; i++) {
            for (int k = 0; k < common_dim; k++){
                M_3[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
}

void multiply_matrix_float(float **M_1, float **M_2, float **M_3, int rows, int common_dim, int cols, const char *filename, int num_threads, int block_size){
    pthread_t threads[num_threads];
    ThreadDataFloat thread_data[num_threads];

    // Partition the matrix into blocks and assign to worker threads
    int block_rows = floor(rows / block_size);
    int block_cols = floor(cols / block_size);
    int thread_index = 0;
    int common_dim_limit = floor(common_dim/4);

    // struct timeval start, end;
    // gettimeofday(&start, NULL);

    for (int i = 0; i < block_rows; i++) {
        for (int j = 0; j < block_cols; j++) {
            if (thread_index == num_threads){
                thread_index = 0;
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
            pthread_create(&threads[thread_index], NULL, multiply_block_float, &thread_data[thread_index]);
            thread_index++;
        }
    }

    // Wait for all worker threads to complete
    for (int i = 0; i < thread_index; i++) {
        pthread_join(threads[i], NULL);
    }

    // Process the leftover rows and columns
    printf("Cleaning up\n");
    printf("Multiply block for start_row: %d end_row: %d\nstart_col: %d end_col: %d\n", block_rows * block_size, rows-1, 0, cols-1);
    for (int i = block_rows * block_size; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < common_dim; k++){
                M_3[i][j] += M_1[i][k] * M_2[k][j];
            }
        }
    }
    printf("Multiply block for start_row: %d end_row: %d\nstart_col: %d end_col: %d\n", 0, rows-1, block_cols * block_size, cols-1);
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
                M[i][j] = 0;
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
                M[i][j] = 0;
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

    printf("m is %d\n", m);
    printf("k is %d\n", k);
    printf("l is %d\n", l);
    printf("n is %d\n", n);

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
        printf("A TIMES B\n");
    } else {
        // otherwise, A * (B * C) is more efficient
        t_rows = k;
        t_cols = n;
        t_type = B_TIMES_C;
        printf("B times C\n");
    }
    srand(SID);

    int blockSize = floor(sqrt(L1_CACHE/sizeof(float)));
    printf("blockSize: %d\n", blockSize);
    int roundedBlockSize;
    if (blockSize % LINE_SIZE < LINE_SIZE/2){
        //round down
        roundedBlockSize = fmax(floor(blockSize/LINE_SIZE), 1) * LINE_SIZE;
    } else {
        //round up
        roundedBlockSize = ceil(blockSize/LINE_SIZE) * LINE_SIZE;
    }
    printf("roundedBlockSize: %d\n", roundedBlockSize);

    if (array_type == FLOAT){
        printf("type is float\n");
        printf("num_threads is %d\n", num_threads);

        float **A = setup_matrix_float(A, m, k, array_type, RANDOM, "A");
        float **B = setup_matrix_float(B, k, l, array_type, RANDOM, "B");
        float **C = setup_matrix_float(C, l, n, array_type, RANDOM, "C");
        float **D = setup_matrix_float(D, m, n, array_type, ZERO, "D");
        float **T = setup_matrix_float(T, t_rows, t_cols, array_type, ZERO, "T");


        struct timeval start, end;
        gettimeofday(&start, NULL);
        if (t_type == A_TIMES_B){
            printf("Multiplying matrix A * B = T, where A = (%d x %d), B = (%d x %d), T = (%d x %d)\n",
            m, k, k, l, t_rows, t_cols);
            multiply_matrix_float(A, B, T, m, k, l, "T", num_threads, roundedBlockSize);
            printf("Multiplying matrix T * C = D, where T = (%d x %d), C = (%d x %d), D = (%d x %d)\n",
                t_rows, t_cols, l, n, m, n);
            multiply_matrix_float(T, C, D, t_rows, t_cols, n, "D", num_threads, roundedBlockSize);
        } else if (t_type == B_TIMES_C){
            printf("Multiplying matrix B * C = T, where B = (%d x %d), C = (%d x %d), T = (%d x %d)\n",
                k, l, l, n, t_rows, t_cols);
            multiply_matrix_float(B, C, T, k, l, n, "T", num_threads, roundedBlockSize);
            printf("Multiplying matrix A * T = D, where A = (%d x %d), T = (%d x %d), D = (%d x %d)\n",
                m, k, t_rows, t_cols, m, n);
            multiply_matrix_float(A, T, D, m, t_rows, t_cols, "D", num_threads, roundedBlockSize);
        }
        gettimeofday(&end, NULL);
        long seconds = end.tv_sec - start.tv_sec;
        long microseconds = end.tv_usec - start.tv_usec;
        double elapsed = seconds + microseconds * 1e-6;
        printf("Elapsed time: %.6f seconds\n", elapsed);

        print_matrix_to_file_float(T, t_rows, t_cols, "T");
        print_matrix_to_file_float(D, m, n, "D");

        // Free allocated memory
        for (int i = 0; i < m; i++) {
            free(A[i]);
            free(D[i]);
        }
        for (int i = 0; i < k; i++) {
            free(B[i]);
        }
        for (int i = 0; i < l; i++) {
            free(C[i]);
        }
        for (int i = 0; i < t_rows; i++){
            free(T[i]);
        }
        free(A);
        free(B);
        free(C);
        free(D);
        free(T);

    } else if (array_type == DOUBLE){
        printf("type is double\n");
        printf("num_threads is %d\n", num_threads);
        double **A = setup_matrix_double(A, m, k, array_type, RANDOM, "A");
        double **B = setup_matrix_double(B, k, l, array_type, RANDOM, "B");
        double **C = setup_matrix_double(C, l, n, array_type, RANDOM, "C");
        double **D = setup_matrix_double(D, m, n, array_type, ZERO, "D");
        double **T = setup_matrix_double(T, t_rows, t_cols, array_type, ZERO, "T");

        struct timeval start, end;
        gettimeofday(&start, NULL);
        if (t_type == A_TIMES_B){
            printf("Multiplying matrix A * B = T, where A = (%d x %d), B = (%d x %d), T = (%d x %d)\n",
            m, k, k, l, t_rows, t_cols);
            multiply_matrix_double(A, B, T, m, k, l, "T", num_threads, roundedBlockSize);
            printf("Multiplying matrix T * C = D, where T = (%d x %d), C = (%d x %d), D = (%d x %d)\n",
                t_rows, t_cols, l, n, m, n);
            multiply_matrix_double(T, C, D, t_rows, t_cols, n, "D", num_threads, roundedBlockSize);
        } else if (t_type == B_TIMES_C){
            printf("Multiplying matrix B * C = T, where B = (%d x %d), C = (%d x %d), T = (%d x %d)\n",
                k, l, l, n, t_rows, t_cols);
            multiply_matrix_double(B, C, T, k, l, n, "T", num_threads, roundedBlockSize);
            printf("Multiplying matrix A * T = D, where A = (%d x %d), T = (%d x %d), D = (%d x %d)\n",
                m, k, t_rows, t_cols, m, n);
            multiply_matrix_double(A, T, D, m, t_rows, t_cols, "D", num_threads, roundedBlockSize);
        }
        gettimeofday(&end, NULL);
        long seconds = end.tv_sec - start.tv_sec;
        long microseconds = end.tv_usec - start.tv_usec;
        double elapsed = seconds + microseconds * 1e-6;
        printf("Elapsed time: %.6f seconds\n", elapsed);

        print_matrix_to_file_double(T, t_rows, t_cols, "T");
        print_matrix_to_file_double(D, m, n, "D");

        // Free allocated memory
        for (int i = 0; i < m; i++) {
            free(A[i]);
            free(D[i]);
        }
        for (int i = 0; i < k; i++) {
            free(B[i]);
        }
        for (int i = 0; i < l; i++) {
            free(C[i]);
        }
        for (int i = 0; i < t_rows; i++){
            free(T[i]);
        }
        free(A);
        free(B);
        free(C);
        free(D);
        free(T);
    }
}
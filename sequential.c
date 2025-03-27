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

#define FLOAT 1
#define DOUBLE 2

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
    if (array_type == FLOAT){
        printf("type is float\n");
    } else if (array_type == DOUBLE){
        printf("type is double\n");
    }
    printf("num_threads is %d\n", num_threads);

    void **A;
    void **B;
    void **C;
    void **D;
    void **T;
    if (array_type == FLOAT){
    } else if (array_type == DOUBLE){

    }

}

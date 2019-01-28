#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "parasail/parasail.h"
#include "parasail/parasail/io.h"
#include "parasail/parasail/matrices/blosum62.h"

static char* get_alphabet(const parasail_matrix_t *matrix) {
    int i = 0;
    char *alphabet = NULL;

    alphabet = (char*)malloc(sizeof(char)*(matrix->size+1));
    for (i=0; i<matrix->size; ++i) {
        alphabet[i] = '*';
    }
    alphabet[matrix->size+1] = '\0';

    for (i=65; i<91; ++i) {
        if (matrix->mapper[i] < matrix->size) {
            alphabet[matrix->mapper[i]] = i;
        }
    }

    return alphabet;
}

static void print_matrix(const parasail_matrix_t *matrix) {
    int i = 0;
    int j = 0;
    char *alphabet = NULL;

    alphabet = get_alphabet(matrix);

    printf("matrix '%s'\n", matrix->name);
    printf("  ");
    for (i=0; i<matrix->size; ++i) {
        printf("%4c", alphabet[i]);
    }
    printf("\n");
    for (j=0; j<matrix->size; ++j) {
        printf("%c ", alphabet[j]);
        for (i=0; i<matrix->size; ++i) {
            printf("%4d", matrix->matrix[j*matrix->size + i]);
        }
        printf("\n");
    }
    printf("max = %d\n", matrix->max);
    printf("min = %d\n", matrix->min);

    free(alphabet);
}

int main(int argc, char **argv)
{
    parasail_matrix_t *matrix = NULL;
    const parasail_matrix_t *internal_matrix = NULL;
    parasail_matrix_t *user_matrix = NULL;

    if (argc == 2) {
        matrix = parasail_matrix_from_file(argv[1]);
    }
    else {
        printf("missing matrix file argument\n");
        return -1;
    }

    print_matrix(matrix);
    parasail_matrix_free(matrix);

    print_matrix(&parasail_blosum62);

    internal_matrix = parasail_matrix_lookup("blosum62");
    if (NULL == internal_matrix) {
        fprintf(stderr, "matrix lookup failed");
        exit(EXIT_FAILURE);
    }

    user_matrix = parasail_matrix_create("ACGT", 2, -1);
    if (NULL == user_matrix) {
        fprintf(stderr, "matrix create failed");
        exit(EXIT_FAILURE);
    }

    parasail_matrix_free(user_matrix);

    user_matrix = parasail_matrix_copy(internal_matrix);
    if (NULL == user_matrix) {
        fprintf(stderr, "matrix copy failed");
        exit(EXIT_FAILURE);
    }

    parasail_matrix_set_value(user_matrix, 10, 10, 100);
    if (100 != user_matrix->max) {
        fprintf(stderr, "matrix set value failed");
        exit(EXIT_FAILURE);
    }

    parasail_matrix_free(user_matrix);

    return 0;
}


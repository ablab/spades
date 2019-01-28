/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>
#include <string.h>

#include "parasail/parasail.h"
#include "parasail/parasail/matrix_lookup.h"

const parasail_matrix_t* parasail_matrix_lookup(const char *matrixname)
{
    const parasail_matrix_t *matrix = NULL;

    if (matrixname) {
        int index = 0;
        const parasail_matrix_t *current = parasail_matrices[index++];
        while (current) {
            if (0 == strcmp(matrixname, current->name)) {
                matrix = current;
                break;
            }
            current = parasail_matrices[index++];
        }
    }

    return matrix;
}


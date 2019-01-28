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
#include "parasail/parasail/function_lookup.h"

parasail_function_t * parasail_lookup_function(const char *funcname)
{
    const parasail_function_info_t * info = NULL;

    info = parasail_lookup_function_info(funcname);

    if (info && info->pointer) {
        return info->pointer;
    }

    return NULL;
}

const parasail_function_info_t * parasail_lookup_function_info(const char *funcname)
{
    const parasail_function_info_t * f = NULL;

    if (funcname) {
        int index = 0;
        f = &functions[index++];
        while (f->pointer) {
            if (0 == strcmp(funcname, f->name)) {
                break;
            }
            f = &functions[index++];
        }
        if (!f->pointer) {
            /* perhaps caller forgot "parasail/parasail_" prefix? */
            const char *prefix = "parasail/parasail_";
            char *newname = (char*)malloc(strlen(prefix)+strlen(funcname)+1);
            strcpy(newname, prefix);
            strcat(newname, funcname);
            index = 0;
            f = &functions[index++];
            while (f->pointer) {
                if (0 == strcmp(newname, f->name)) {
                    break;
                }
                f = &functions[index++];
            }
            free(newname);
        }
    }

    if (!f->pointer) {
        f = NULL;
    }

    return f;
}

parasail_pfunction_t * parasail_lookup_pfunction(const char *funcname)
{
    const parasail_pfunction_info_t * info = NULL;

    info = parasail_lookup_pfunction_info(funcname);

    if (info && info->pointer) {
        return info->pointer;
    }

    return NULL;
}

parasail_pcreator_t * parasail_lookup_pcreator(const char *funcname)
{
    const parasail_pfunction_info_t * info = NULL;

    info = parasail_lookup_pfunction_info(funcname);

    if (info && info->creator) {
        return info->creator;
    }

    return NULL;
}

const parasail_pfunction_info_t * parasail_lookup_pfunction_info(const char *funcname)
{
    const parasail_pfunction_info_t * f = NULL;

    if (funcname) {
        int index = 0;
        f = &pfunctions[index++];
        while (f->pointer) {
            if (0 == strcmp(funcname, f->name)) {
                break;
            }
            f = &pfunctions[index++];
        }
        if (!f->pointer) {
            /* perhaps caller forgot "parasail/parasail_" prefix? */
            const char *prefix = "parasail/parasail_";
            char *newname = (char*)malloc(strlen(prefix)+strlen(funcname)+1);
            strcpy(newname, prefix);
            strcat(newname, funcname);
            index = 0;
            f = &pfunctions[index++];
            while (f->pointer) {
                if (0 == strcmp(newname, f->name)) {
                    break;
                }
                f = &pfunctions[index++];
            }
            free(newname);
        }
    }

    if (!f->pointer) {
        f = NULL;
    }

    return f;
}


/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifdef PARASAIL_TABLE
#define ENAME parasail_neon_dummy_table
#else
#ifdef PARASAIL_ROWCOL
#define ENAME parasail_neon_dummy_rowcol
#else
#ifdef PARASAIL_TRACE
#define ENAME parasail_neon_dummy_trace
#else
#define ENAME parasail_neon_dummy
#endif
#endif
#endif

extern int ENAME(void);

int ENAME()
{
    return 0;
}


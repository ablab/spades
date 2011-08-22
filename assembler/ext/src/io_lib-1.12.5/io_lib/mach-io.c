/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

/*
 * Machine independant io:
 * For reading and writing to big-endian and little-endian files
 *
 * Routines available:
 *     be_write_int_1()
 *     be_write_int_2()
 *     be_write_int_4()
 *     be_write_int_8()
 *     be_read_int_1()
 *     be_read_int_2()
 *     be_read_int_4()
 *     be_read_int_8()
 *     le_write_int_1()
 *     le_write_int_2()
 *     le_write_int_4()
 *     le_write_int_8()
 *     le_read_int_1()
 *     le_read_int_2()
 *     le_read_int_4()
 *     le_read_int_8()
 *
 * All routine return:
 *    0 - an error has occurred during io operation
 *    1 - value suggessfully read or written
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include "io_lib/stdio_hack.h"
#include "io_lib/mach-io.h"




/**********************************************************************/
/* IO for big-endian files                                            */
/**********************************************************************/

/*
 * Write a big-endian int1
 */
int be_write_int_1(FILE *fp, uint1 *i1)
{
    if (fwrite(i1, sizeof(uint1), 1, fp) != 1) return (0);
    return (1);
}


/*
* Write a big-endian int2
*/
int be_write_int_2(FILE *fp, uint2 *i2)
{
    uint2 i = be_int2(*i2);
    
    if (fwrite(&i, 2, 1, fp) != 1) return (0);
    return (1);
}

/*
 * Write a big-endian int4
 */
int be_write_int_4(FILE *fp, uint4 *i4)
{
    uint4 i = be_int4(*i4);
    
    if (fwrite(&i, 4, 1, fp) != 1) return (0);

    return (1);
}

/*
 * Write a big-endian int8
 */
int be_write_int_8(FILE *fp, uint8 *i8)
{
    uint8 i = be_int8(*i8);
    
    if (fwrite(&i, 8, 1, fp) != 1) return (0);

    return (1);
}


/*
 * Read a big-endian int1
 */
int be_read_int_1(FILE *fp, uint1 *i1)
{
    if (fread(i1, sizeof(uint1), 1, fp) != 1) return (0);
    return (1);
}


/*
 * Read a big-endian int2
 */
int be_read_int_2(FILE *fp, uint2 *i2)
{
    uint2 i;
    
    if (fread(&i, 2, 1, fp) != 1) return (0);
    *i2 = be_int2(i);

    return (1);
}


/*
 * Read a big-endian int4
 */
int be_read_int_4(FILE *fp, uint4 *i4)
{
    uint4 i;
    
    if (fread(&i, 4, 1, fp) != 1) return (0);
    *i4 = be_int4(i);

    return (1);
}


/*
 * Read a big-endian int8
 */
int be_read_int_8(FILE *fp, uint8 *i8)
{
    uint8 i;
    
    if (fread(&i, 8, 1, fp) != 1) return (0);
    *i8 = be_int8(i);

    return (1);
}





/**********************************************************************/
/* IO for little-endian files                                         */
/**********************************************************************/

/*
 * Write a little-endian int1
 */
int le_write_int_1(FILE *fp, uint1 *i1)
{
    if (fwrite(i1, sizeof(uint1), 1, fp) != 1) return (0);
    return (1);
}


/*
 * Write a little-endian int2
 */
int le_write_int_2(FILE *fp, uint2 *i2)
{
    uint2 i = le_int2(*i2);
    
    if (fwrite(&i, 2, 1, fp) != 1) return (0);

    return (1);
}


/*
 * Write a little-endian int4
 */
int le_write_int_4(FILE *fp, uint4 *i4)
{
    uint4 i = le_int4(*i4);
    
    if (fwrite(&i, 4, 1, fp) != 1) return (0);

    return (1);
}


/*
 * Write a little-endian int8
 */
int le_write_int_8(FILE *fp, uint8 *i8)
{
    uint8 i = le_int8(*i8);
    
    if (fwrite(&i, 8, 1, fp) != 1) return (0);

    return (1);
}


/*
 * Read a little-endian int1
 */
int le_read_int_1(FILE *fp, uint1 *i1)
{
    if (fread(i1, sizeof(uint1), 1, fp) != 1) return (0);
    return (1);
}


/*
 * Read a little-endian int2
 */
int le_read_int_2(FILE *fp, uint2 *i2)
{
    uint2 i;
    
    if (fread(&i, 2, 1, fp) != 1) return (0);
    *i2 = le_int2(i);

    return (1);
}

/*
 * Read a little-endian int4
 */
int le_read_int_4(FILE *fp, uint4 *i4)
{
    uint4 i;
    
    if (fread(&i, 4, 1, fp) != 1) return (0);
    *i4 = le_int4(i);

    return (1);
}


/*
 * Read a little-endian int8
 */
int le_read_int_8(FILE *fp, uint8 *i8)
{
    uint8 i;
    
    if (fread(&i, 8, 1, fp) != 1) return (0);
    *i8 = le_int8(i);

    return (1);
}

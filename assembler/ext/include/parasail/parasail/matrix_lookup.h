/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_MATRIX_LOOKUP_H_
#define _PARASAIL_MATRIX_LOOKUP_H_

#include "parasail/parasail/matrices/blosum100.h"
#include "parasail/parasail/matrices/blosum30.h"
#include "parasail/parasail/matrices/blosum35.h"
#include "parasail/parasail/matrices/blosum40.h"
#include "parasail/parasail/matrices/blosum45.h"
#include "parasail/parasail/matrices/blosum50.h"
#include "parasail/parasail/matrices/blosum55.h"
#include "parasail/parasail/matrices/blosum60.h"
#include "parasail/parasail/matrices/blosum62.h"
#include "parasail/parasail/matrices/blosum65.h"
#include "parasail/parasail/matrices/blosum70.h"
#include "parasail/parasail/matrices/blosum75.h"
#include "parasail/parasail/matrices/blosum80.h"
#include "parasail/parasail/matrices/blosum85.h"
#include "parasail/parasail/matrices/blosum90.h"
#include "parasail/parasail/matrices/blosumn.h"
#include "parasail/parasail/matrices/pam10.h"
#include "parasail/parasail/matrices/pam100.h"
#include "parasail/parasail/matrices/pam110.h"
#include "parasail/parasail/matrices/pam120.h"
#include "parasail/parasail/matrices/pam130.h"
#include "parasail/parasail/matrices/pam140.h"
#include "parasail/parasail/matrices/pam150.h"
#include "parasail/parasail/matrices/pam160.h"
#include "parasail/parasail/matrices/pam170.h"
#include "parasail/parasail/matrices/pam180.h"
#include "parasail/parasail/matrices/pam190.h"
#include "parasail/parasail/matrices/pam20.h"
#include "parasail/parasail/matrices/pam200.h"
#include "parasail/parasail/matrices/pam210.h"
#include "parasail/parasail/matrices/pam220.h"
#include "parasail/parasail/matrices/pam230.h"
#include "parasail/parasail/matrices/pam240.h"
#include "parasail/parasail/matrices/pam250.h"
#include "parasail/parasail/matrices/pam260.h"
#include "parasail/parasail/matrices/pam270.h"
#include "parasail/parasail/matrices/pam280.h"
#include "parasail/parasail/matrices/pam290.h"
#include "parasail/parasail/matrices/pam30.h"
#include "parasail/parasail/matrices/pam300.h"
#include "parasail/parasail/matrices/pam310.h"
#include "parasail/parasail/matrices/pam320.h"
#include "parasail/parasail/matrices/pam330.h"
#include "parasail/parasail/matrices/pam340.h"
#include "parasail/parasail/matrices/pam350.h"
#include "parasail/parasail/matrices/pam360.h"
#include "parasail/parasail/matrices/pam370.h"
#include "parasail/parasail/matrices/pam380.h"
#include "parasail/parasail/matrices/pam390.h"
#include "parasail/parasail/matrices/pam40.h"
#include "parasail/parasail/matrices/pam400.h"
#include "parasail/parasail/matrices/pam410.h"
#include "parasail/parasail/matrices/pam420.h"
#include "parasail/parasail/matrices/pam430.h"
#include "parasail/parasail/matrices/pam440.h"
#include "parasail/parasail/matrices/pam450.h"
#include "parasail/parasail/matrices/pam460.h"
#include "parasail/parasail/matrices/pam470.h"
#include "parasail/parasail/matrices/pam480.h"
#include "parasail/parasail/matrices/pam490.h"
#include "parasail/parasail/matrices/pam50.h"
#include "parasail/parasail/matrices/pam500.h"
#include "parasail/parasail/matrices/pam60.h"
#include "parasail/parasail/matrices/pam70.h"
#include "parasail/parasail/matrices/pam80.h"
#include "parasail/parasail/matrices/pam90.h"
#include "parasail/parasail/matrices/nuc44.h"
#include "parasail/parasail/matrices/dnafull.h"
#include "parasail/parasail/matrices/blosum_map.h"
#include "parasail/parasail/matrices/pam_map.h"

#ifdef __cplusplus
extern "C" {
#endif

static const parasail_matrix_t * parasail_matrices[] = {
    &parasail_blosum100,
    &parasail_blosum30,
    &parasail_blosum35,
    &parasail_blosum40,
    &parasail_blosum45,
    &parasail_blosum50,
    &parasail_blosum55,
    &parasail_blosum60,
    &parasail_blosum62,
    &parasail_blosum65,
    &parasail_blosum70,
    &parasail_blosum75,
    &parasail_blosum80,
    &parasail_blosum85,
    &parasail_blosum90,
    &parasail_blosumn,
    &parasail_pam10,
    &parasail_pam100,
    &parasail_pam110,
    &parasail_pam120,
    &parasail_pam130,
    &parasail_pam140,
    &parasail_pam150,
    &parasail_pam160,
    &parasail_pam170,
    &parasail_pam180,
    &parasail_pam190,
    &parasail_pam20,
    &parasail_pam200,
    &parasail_pam210,
    &parasail_pam220,
    &parasail_pam230,
    &parasail_pam240,
    &parasail_pam250,
    &parasail_pam260,
    &parasail_pam270,
    &parasail_pam280,
    &parasail_pam290,
    &parasail_pam30,
    &parasail_pam300,
    &parasail_pam310,
    &parasail_pam320,
    &parasail_pam330,
    &parasail_pam340,
    &parasail_pam350,
    &parasail_pam360,
    &parasail_pam370,
    &parasail_pam380,
    &parasail_pam390,
    &parasail_pam40,
    &parasail_pam400,
    &parasail_pam410,
    &parasail_pam420,
    &parasail_pam430,
    &parasail_pam440,
    &parasail_pam450,
    &parasail_pam460,
    &parasail_pam470,
    &parasail_pam480,
    &parasail_pam490,
    &parasail_pam50,
    &parasail_pam500,
    &parasail_pam60,
    &parasail_pam70,
    &parasail_pam80,
    &parasail_pam90,
    &parasail_nuc44,
    &parasail_dnafull,
    NULL
};

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_MATRIX_LOOKUP_H_ */


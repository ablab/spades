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
 * File:	translate.c
 * Purpose:	Translates between different reading structures.
 * Last update:	01/09/94
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#ifndef NDEBUG
#    define NDEBUG /* disable assertions */
#endif
#include <assert.h>

#include "io_lib/stdio_hack.h"
#include "io_lib/misc.h"
#include "io_lib/scf.h"
#include "io_lib/Read.h"
#include "io_lib/expFileIO.h"
#include "io_lib/traceType.h"
#include "io_lib/translate.h"
#include "io_lib/open_trace_file.h"
#ifdef USE_BIOLIMS
#include "spBiolims.h"
#endif



#ifdef IOLIB_SCF
/*
 * Translates an Scf structure into a Read structure.
 * The Scf structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Read structure upon success.
 *    NULLRead upon failure.
 */
Read *scf2read(Scf *scf) {
    Read *read;
    register int i, i_end;
    TRACE max_val = 0;
    int sections = read_sections(0);
    int nsamples = 0;
    int nbases = 0;

    /* allocate */
    if (sections & READ_SAMPLES)
	nsamples = scf->header.samples;
    if (sections & READ_BASES)
	nbases = scf->header.bases;
    read = read_allocate(nsamples, nbases);

    if (NULLRead == read)
	return NULLRead;

    if (sections & READ_SAMPLES) {
	/* copy the samples */
	i_end = scf->header.samples;
	read->NPoints = i_end;
	
	if (scf->header.sample_size == 1) {
	    for (i = 0; i < i_end; i++) {
		read->traceA[i] = scf->samples.samples1[i].sample_A;
		read->traceC[i] = scf->samples.samples1[i].sample_C;
		read->traceG[i] = scf->samples.samples1[i].sample_G;
		read->traceT[i] = scf->samples.samples1[i].sample_T;
		
		if (read->traceA[i] > max_val) max_val = read->traceA[i];
		if (read->traceC[i] > max_val) max_val = read->traceC[i];
		if (read->traceG[i] > max_val) max_val = read->traceG[i];
		if (read->traceT[i] > max_val) max_val = read->traceT[i];
	    }
	} else { /* sample_size == 2 */
	    for (i = 0; i < i_end; i++) {
		read->traceA[i] = scf->samples.samples2[i].sample_A;
		read->traceC[i] = scf->samples.samples2[i].sample_C;
		read->traceG[i] = scf->samples.samples2[i].sample_G;
		read->traceT[i] = scf->samples.samples2[i].sample_T;
		
		if (read->traceA[i] > max_val) max_val = read->traceA[i];
		if (read->traceC[i] > max_val) max_val = read->traceC[i];
		if (read->traceG[i] > max_val) max_val = read->traceG[i];
		if (read->traceT[i] > max_val) max_val = read->traceT[i];
	    }
	}
	
	read->maxTraceVal = max_val;
    }
    
    if (sections & READ_BASES) {
	/* copy the bases */
	i_end = scf->header.bases;
	read->NBases = i_end;

	for (i = 0; i < i_end; i++) {
	    read->basePos[i] = scf->bases[i].peak_index;
	    read->prob_A[i]  = scf->bases[i].prob_A;
	    read->prob_C[i]  = scf->bases[i].prob_C;
	    read->prob_G[i]  = scf->bases[i].prob_G;
	    read->prob_T[i]  = scf->bases[i].prob_T;
	    read->base[i]    = scf->bases[i].base;
	}
	read->base[i] = 0;
    }
    
    if (sections & READ_COMMENTS) {
	/* allocate and copy the comments */
	if (scf->header.comments_size > 0 && scf->comments) {
	    read->info = (char *)xmalloc(scf->header.comments_size+1);
	    if (NULL == read->info) {
		read_deallocate(read);
		return NULLRead;
	    }

	    memcpy(read->info, scf->comments, scf->header.comments_size);
	    read->info[scf->header.comments_size] = '\0';
	}
    }

    /* other bits and pieces */
    read->leftCutoff = scf->header.bases_left_clip;
    read->rightCutoff = read->NBases - scf->header.bases_right_clip + 1;
    read->format = TT_SCF;

    if (scf->private_data) {
	read->private_data = xmalloc(scf->header.private_size);
	memcpy(read->private_data,scf->private_data, scf->header.private_size);
    }

    return read;
}

/*
 * Translates a Read structure into a Scf structure.
 * The Read structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Scf structure upon success.
 *    NULL upon failure.
 */
Scf *read2scf(Read *read) {
    Scf *scf;
    register int i, i_end;
    int sample_size;

    /* allocate */
    sample_size = read->maxTraceVal >= 0x100 ? 2 : 1;
    scf = scf_allocate(read->NPoints, sample_size, read->NBases, 0, 0);
    if (NULL == scf)
	return NULL;

    /* copy the samples */
    i_end = read->NPoints;
    scf->header.samples = i_end;

    if (sample_size == 1) {
	scf->header.sample_size = 1;
	for (i = 0; i < i_end; i++) {
	    scf->samples.samples1[i].sample_A = (uint_1)read->traceA[i];
	    scf->samples.samples1[i].sample_C = (uint_1)read->traceC[i];
	    scf->samples.samples1[i].sample_G = (uint_1)read->traceG[i];
	    scf->samples.samples1[i].sample_T = (uint_1)read->traceT[i];
	}
    } else {
	scf->header.sample_size = 2;
	for (i = 0; i < i_end; i++) {
	    scf->samples.samples2[i].sample_A = read->traceA[i];
	    scf->samples.samples2[i].sample_C = read->traceC[i];
	    scf->samples.samples2[i].sample_G = read->traceG[i];
	    scf->samples.samples2[i].sample_T = read->traceT[i];
	}
    }

    /* copy the bases */    
    i_end = read->NBases;
    scf->header.bases = i_end;
    
    for (i = 0; i < i_end; i++) {
	scf->bases[i].peak_index = read->basePos ? read->basePos[i] : i;
	scf->bases[i].prob_A     = read->prob_A  ? read->prob_A[i] : 0;
	scf->bases[i].prob_C     = read->prob_A  ? read->prob_C[i] : 0;
	scf->bases[i].prob_G     = read->prob_A  ? read->prob_G[i] : 0;
	scf->bases[i].prob_T     = read->prob_A  ? read->prob_T[i] : 0;
	scf->bases[i].base       = read->base    ? read->base[i] : '-';
    }

    /* allocate and copy the comments */
    if (read->info) {
	scf->header.comments_size = strlen(read->info) + 1;
	scf->comments = (char *)xmalloc(scf->header.comments_size);
	if (NULL == scf->comments) {
	    scf_deallocate(scf);
	    return NULL;
	}

	memcpy(scf->comments, read->info, scf->header.comments_size - 1);

	/* just to make sure */
	scf->comments[scf->header.comments_size-1] = '\0';
    }

    /* other bits and pieces */
    scf->header.bases_left_clip = read->leftCutoff;
    scf->header.bases_right_clip = read->NBases - read->rightCutoff + 1;
    scf->header.code_set = CSET_DEFAULT;
    memcpy(scf->header.version, scf_version_float2str(SCF_VERSION), 4);

    return scf;
}
#endif /* IOLIB_SCF */


#ifdef IOLIB_EXP

#define extend(e, entry, len) \
do { \
    (void)ArrayRef(e->entries[entry],e->Nentries[entry]++); \
    if (NULL == (exp_get_entry(e, entry) = (char *)xmalloc(len))) \
	return NULL; \
} while (0)

/*
 * Translates a Read structure and an Experiment file.
 * The Read structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Exp_info structure upon success.
 *    NULL upon failure.
 */
Exp_info *read2exp(Read *read, char *EN) {
    Exp_info *e;
    char *t = trace_type_int2str(read->format), *p;
    int l = strlen(EN)+1;
    char *sq;
    int i;
    static char valid_bases[256];
    static int valid_setup = 0;

    if (!valid_setup) {
	for (i = 0; i < 256; i++)
	    valid_bases[i] = '-';
	/* IUBC codes */
	for (sq = "acgturymkswbdhvnACGTURYMKSWBDHVN"; *sq; sq++)
	    valid_bases[(unsigned)*sq] = *sq;
	valid_setup = 1;
    }

    if (NULL == (e = exp_create_info()))
	return NULL;

    /* Copy original exp file if present */
    if (read->orig_trace && read->orig_trace_format == TT_EXP) {
	int i, j, k;
	Exp_info *re = (Exp_info *)read->orig_trace;

	for (i = 0; i < MAXIMUM_EFLTS; i++) {
	    if (EFLT_SQ == i ||
		EFLT_QL == i ||
		EFLT_QR == i)
		continue;

	    if (0 == (k = exp_Nentries(re, i)))
		continue;

	    e->Nentries[i] = k;	    
	    ArrayRef(e->entries[i], e->Nentries[i]);
	    for (j = 0; j < k; j++) {
		arr(char *, e->entries[i], j) =
		    strdup(arr(char *, re->entries[i], j));
	    }
	}

    /* Otherwise create our EN, ID, LN and LT lines */
    } else {
	/* Entry name and ID lines */
	if (p = strrchr(EN, '/'))
	    EN = p+1;
	extend(e, EFLT_EN, l);
	sprintf(exp_get_entry(e, EFLT_EN), "%s", EN);
	extend(e, EFLT_ID, l);
	sprintf(exp_get_entry(e, EFLT_ID), "%s", EN);

	/* Trace file & type */
	if (read->trace_name) {
	    char *cp;
	    if (cp = strrchr(read->trace_name, '/'))
		cp++;
	    else
		cp = read->trace_name;
	    extend(e, EFLT_LN, strlen(cp)+1);
	    strcpy(exp_get_entry(e, EFLT_LN), cp);
	}

	if (read->format != TT_ANY && read->format != TT_ANYTR) {
	    extend(e, EFLT_LT, strlen(t)+1);
	    strcpy(exp_get_entry(e, EFLT_LT), t);
	}
    }

    /* Output SQ, QL and QR lines */

    /* Cutoffs */
    if (read->leftCutoff) {
	extend(e, EFLT_QL, 15);
	sprintf(exp_get_entry(e, EFLT_QL), "%d", read->leftCutoff);
    }

    if (read->rightCutoff && read->rightCutoff != read->NBases+1) {
	extend(e, EFLT_QR, 15);
	sprintf(exp_get_entry(e, EFLT_QR), "%d", read->rightCutoff);
    }

    /* Bases */
    extend(e, EFLT_SQ, read->NBases+1);
    sq = exp_get_entry(e, EFLT_SQ);
    for (i = 0; i < read->NBases; i++) {
	sq[i] = valid_bases[(unsigned)read->base[i]];
    }
    sq[i] = 0;

#ifdef USE_BIOLIMS
    /*
     * Johnt:
     * - Added tags below to allow for biolims update
     * - This is all a very big bodge to allow BioLIMS
     *   attributes to be passed through the Read structure
     *   to the Experiment file
     * - Any changes to this should also be mirrored in ../biolims/Exp.cpp
     */
    {
	int1 *qa; /* quality array */
	int i;
	char tmp[1024];
	char *line; /* current line from info */
	char *end;  /* end of current line */
	
	/* AV only supports a single prob value for now */
	qa = (int1 *)malloc((read->NBases+1)*sizeof(int1));

	/* need max 4 bytes per value - see conf2str */
	extend(e,EFLT_AV,(read->NBases+1)*5);

	/* merge into single quality array */
	for(i=0;i<read->NBases;i++){
	    switch(read->base[i]){
	    case 'a':
	    case 'A':
		qa[i] = read->prob_A[i];
		break;
	    case 'c':
	    case 'C':
		qa[i] = read->prob_C[i];
		break;
	    case 'g':
	    case 'G':
		qa[i] = read->prob_G[i];
		break;
	    case 't':
	    case 'T':
		qa[i] = read->prob_T[i];
		break;
	    default:
		qa[i] = 0;
	    }
	}
	conf2str(qa,read->NBases,exp_get_entry(e, EFLT_AV));
	free(qa);
	
	/* Parse the read notes for everything else */
	if( read->info) {
	    for(line=read->info;end=strchr(line,'\n');line=end+1){
		*end='\0'; /* * put back */

		/* look for known tags */
		if(!strncmp(line,EXP_CHEM,EXP_TAGLEN)){
		    /* CH */
		    int chem=0;
		    extend(e, EFLT_CH, 15);
		    if( !strcmp(CH_types[1],line+EXP_TAGLEN))
			chem=1;
		    sprintf(exp_get_entry(e, EFLT_CH), "%d",chem );

		} else if(!strncmp(line,EXP_PRMR,EXP_TAGLEN)){
		    /* PR */
		    int primer=0;
		    extend(e,EFLT_PR,15);
		    for(primer=1;primer<nPR_types;primer++)
			if(!strcmp(PR_types[i],line+EXP_TAGLEN))
			    break;
		    if(primer>=nPR_types)
			primer=1;
		    sprintf(exp_get_entry(e, EFLT_PR), "%d",primer);

		} else if(!strncmp(line,EXP_VECT,EXP_TAGLEN)){
		    /* SV */
		    extend(e,EFLT_SV,strlen(line)-EXP_TAGLEN+1);
		    strcpy(exp_get_entry(e, EFLT_SV),line+EXP_TAGLEN);

		} else if(!strncmp(line,EXP_CLOV,EXP_TAGLEN)){
		    /* CV */
		    extend(e,EFLT_CV,strlen(line)-EXP_TAGLEN+1);
		    strcpy(exp_get_entry(e, EFLT_CV),line+EXP_TAGLEN);

		} else if(!strncmp(line,EXP_CLON,EXP_TAGLEN)){
		    /* CN */
		    extend(e,EFLT_CN,strlen(line)-EXP_TAGLEN+1);
		    strcpy(exp_get_entry(e, EFLT_CN),line+EXP_TAGLEN);

		} else if(!strncmp(line,EXP_FEAT,EXP_TAGLEN)){
		    /* FEAT=start stop key\r\rcomment */
		    /* key and comment have \n encoded as \r */
		    int start,stop,i;
		    char *key; /* biolims feature key */
		    char *comment; /* biolims feature comment */
		    line+=EXP_TAGLEN;
		    start=atoi(line);
		    line=strchr(line,' ')+1;
		    stop=atoi(line);
		    key=strchr(line,' ')+1;
		    comment=strstr(key,"\r\r");
		    *comment='\0'; /* * put back */
		    comment+=2;
		    /* replace \r with \n in key and comment */
		    for(i=0;key[i];i++)
			if(key[i]=='\r') key[i]='\n';
		    for(i=0;comment[i];i++)
			if(comment[i]=='\r') comment[i]='\n';
		    /* could possibly be one of a number of EXP tags excoded
		       as a BioLIMS feature */
		    if(!strncmp(key,featCLON,STADEN_FKEY_LEN)){
			/* CS */
			extend(e, EFLT_CS, 32);
			exp_create_range(exp_get_entry(e, EFLT_CS),start,stop);
		    } else if(!strncmp(key,featVECI,STADEN_FKEY_LEN)){
			/* SI */
			extend(e, EFLT_SI, 32);
			exp_create_range(exp_get_entry(e, EFLT_SI),start,stop);
		    } else if(!strncmp(key,featTEMP,STADEN_FKEY_LEN)){
			/* TN */
			extend(e, EFLT_TN, strlen(comment)+1);
			strcpy(exp_get_entry(e, EFLT_TN),comment);
		    } else if(!strncmp(key,featSTRD,STADEN_FKEY_LEN)){
			/* ST */
			extend(e, EFLT_ST, strlen(comment)+1);
			strcpy(exp_get_entry(e, EFLT_ST),comment);
		    } else if( !strncmp(key,featVECT,STADEN_FKEY_LEN)){
			/* SL and SR */
			extend(e,EFLT_SL,15);
			extend(e,EFLT_SR,15);
			sprintf(exp_get_entry(e, EFLT_SL),"%d",start);
			sprintf(exp_get_entry(e, EFLT_SR),"%d",stop);
		    } else if( !strncmp(key,featGELR,STADEN_FKEY_LEN)){
			/* TG */
			char tag[5]; /* staden note tag (always 4 chars) */
			char strand=*comment; /* first char of comment */
			/* key has format STADEN_GELR:XXXX -
			   where XXXX is staden note tag */
			strncpy(tag,key+STADEN_FKEY_LEN+1,4);
			tag[4]='\0';
			comment+=2; /* skip over strand */
			sprintf(tmp,"%s %c %d..%d\n%s",
				key,strand,start,stop,comment);
			extend(e,EFLT_TG,strlen(tmp)+1);
			strcpy(exp_get_entry(e, EFLT_TG),tmp);
			
		    } else if( !strncmp(key,featCONS,STADEN_FKEY_LEN)){
			/* TC */
			char tag[5]; /* staden note tag (always 4 chars) */
			/* comment has the format Srstart-rend\ncomment
			   S is a single character strand indicator
			   rstart and rend are the real start and end of the
			   tag
			*/
			char strand=*comment; /* first char of comment */
			char *rangestart = comment+2; /*skip over strand*/
			char *rangeend;
			char *emptystring="";
			
			/* key has format STADEN_CONS:XXXX -
			   where XXXX is staden note tag */
			strncpy(tag,key+STADEN_FKEY_LEN+1,4);
			tag[4]='\0';
			/* the REAL range might actually be outside the bases,
			   so this is recorded in the first line of
			   the comment. This is merged with the
			   feature range to allow for complimenting */ 
			comment=strchr(rangestart,'\n');
			if( comment )
			    *(comment++)='\0'; /* *** put back */
			else
			    comment=emptystring;
			
			/* now special processing of bounds */
			rangeend=strchr(rangestart,'-');
			if( rangeend ){
			    long rstart,rend;
			    *(rangeend++)='\0'; /* *** put back */
			    rstart=atol(rangestart);
			    rend=atol(rangeend);
			    *(rangeend-1)='-';
			    
			    /* if start is the same as rstart, just
			       need to extend the stop bounds to the
			       real end If not there has been some
			       complimenting happening so need to
			       extend the start bounds, as its now
			       backwards.  */
			    if( start==rstart)
				stop=start+rend-rstart;
			    else 
				start=stop-rend+rstart;
			} else {
			    /* don't have a range so just use the
			       feature size */
			}
			
			sprintf(tmp,"%s %c %d..%d\n%s",
				key,strand,start,stop,comment);
			extend(e,EFLT_TC,strlen(tmp)+1);
			strcpy(exp_get_entry(e, EFLT_TC),tmp);
			if( comment!=emptystring)
			    *(comment-1)='\n';
			
		    } else {
			/* TG */
			/* a biolims feature, encode this with tag
			   BIOL, and the Biolims Feature key in
			   the first line of the comment*/
			sprintf(tmp,"%s = %d..%d\n%s%s\n%s",/* use strand = */
				tagBIOL,start,stop,featKey,key,comment); 
			extend(e,EFLT_TG,strlen(tmp)+1);
			strcpy(exp_get_entry(e, EFLT_TG),tmp);
		    }
		    *(comment-2)='\r';
		}
		/* else unused value */
		
		*end='\n';
	    }
	}
    }
#else /* USE_BIOLIMS */
    /* Confidence values */
    if (read->prob_A && read->prob_C && read->prob_G && read->prob_T &&
	read->NBases > 0) {
	/* We have some, but are they non zero values? */
	for (i = 0; i < read->NBases; i++) {
	    if (read->prob_A[i] || read->prob_C[i] ||
		read->prob_G[i] || read->prob_T[i])
		break;
	}
	if (i != read->NBases) {
	    int1 *conf = (int1 *)xmalloc(read->NBases);
	    char *cstr = (char *)xmalloc(read->NBases * 5 + 2);

	    for (i = 0; i < read->NBases; i++) {
		switch (read->base[i]) {
		case 'a':
		case 'A':
		    conf[i] = read->prob_A[i];
		    break;
		case 'c':
		case 'C':
		    conf[i] = read->prob_C[i];
		    break;
		case 'g':
		case 'G':
		    conf[i] = read->prob_G[i];
		    break;
		case 't':
		case 'T':
		    conf[i] = read->prob_T[i];
		    break;
		default:
		    conf[i] = (read->prob_A[i] + read->prob_C[i] +
			       read->prob_G[i] + read->prob_T[i]) / 4;
		}
	    }

	    conf2str(conf, read->NBases, cstr);
	    extend(e, EFLT_AV, strlen(cstr)+1);
	    sprintf(exp_get_entry(e, EFLT_AV), "%s", cstr);
	    xfree(conf);
	    xfree(cstr);
	}
    }
#endif /* if USE_BIOLIMS else ... */

    return e;
}


/*
 * Controls the use of the SQ and ON lines when loading an experiment file.
 * The default (value&1 == 1) is to load these into the Read structure.
 * With value&1 == 0 we load the sequence directly from the trace file
 * (LT line).
 * value&2 controls whether to use the SL/SR fields when setting the cutoff.
 * value&2 == 0 implies to do so, and value&2 == 2 implies to not.
 *
 * The default use is to use the SQ and ON lines. Returns the old value.
 */
static int use_experiment_sequence = 1;
int read_experiment_redirect(int value) {
    int old = use_experiment_sequence;

    use_experiment_sequence = value;
    return old;
}


/*
 * Translates an experiment file to a Read structure.
 * The Exp_info structure is left unchanged.
 *
 * Returns:
 *    A pointer to an allocated Read structure upon success.
 *    NULLRead upon failure.
 */
Read *exp2read(Exp_info *e, char *fn) {
    Read *r;
    int q, s, ttype, err = 0;
    char *str;
    int use_exp = use_experiment_sequence;
    FILE *fp;

    if (!exp_Nentries(e, EFLT_LN)) {
	err = 1;
    } else {
	/* Read the trace component of the experiment file */
	ttype = exp_Nentries(e,EFLT_LT)
	    ? trace_type_str2int(exp_get_entry(e, EFLT_LT)) : TT_ANYTR;

	if (fp = open_trace_file(exp_get_entry(e, EFLT_LN), fn)) {
	    if (NULLRead == (r = fread_reading(fp, NULL, ttype)))
		err = 1;
	} else {
	    err = 1;
	}
    }

    if (err) {
	use_exp = 1;
	r = read_allocate(0, 1);
    }
    
    /* Set the left cutoff (QL / SL) */
    q=-1;
    if (exp_Nentries(e, EFLT_QL))
	q = atoi(exp_get_entry(e, EFLT_QL));
    if ((use_exp&2) != 2) {
	s=-1;
	if (exp_Nentries(e, EFLT_SL))
	    s = atoi(exp_get_entry(e, EFLT_SL));
	if (q != -1 || s != -1)
	    r->leftCutoff = MAX(q, s);
    } else {
	r->leftCutoff = q != -1 ? q : 0;
    }

    /* Set the right cutoff (QR / SR) */
    q = INT_MAX;
    if (exp_Nentries(e, EFLT_QR))
	q = atoi(exp_get_entry(e, EFLT_QR));
    if ((use_exp&2) != 2) {
	s = INT_MAX;
	if (exp_Nentries(e, EFLT_SR))
	    s = atoi(exp_get_entry(e, EFLT_SR));
	if (q != INT_MAX || s != INT_MAX)
	    r->rightCutoff = MIN(q, s);
    } else {
	r->rightCutoff = q != INT_MAX ? q : 0;
    }

    if (r->rightCutoff && r->rightCutoff <= r->leftCutoff)
      r->rightCutoff = r->leftCutoff+1;

    /* Bases and base positions, if desired */
    if (use_exp&1) {
	if (exp_Nentries(e, EFLT_SQ) && (str = exp_get_entry(e, EFLT_SQ))) {
	    int slen = strlen(str);
	    
	    if (NULL == (r->base = (char *)xrealloc(r->base, slen+1)))
		return NULLRead;
	    if (NULL == (r->prob_A = (char *)xrealloc(r->prob_A, slen+1)))
		return NULLRead;
	    if (NULL == (r->prob_C = (char *)xrealloc(r->prob_C, slen+1)))
		return NULLRead;
	    if (NULL == (r->prob_G = (char *)xrealloc(r->prob_G, slen+1)))
		return NULLRead;
	    if (NULL == (r->prob_T = (char *)xrealloc(r->prob_T, slen+1)))
		return NULLRead;
	    if (r->basePos) {
		xfree(r->basePos);
		r->basePos = NULL;
	    }

	    /* Clear them */
	    memset(r->prob_A, 0, slen);
	    memset(r->prob_C, 0, slen);
	    memset(r->prob_G, 0, slen);
	    memset(r->prob_T, 0, slen);
		
	    strcpy(r->base, str);
	    r->NBases = slen;

	    /* Copy AV values into prob_* arrays */
	    if (exp_Nentries(e, EFLT_AV) &&
		(str = exp_get_entry(e, EFLT_AV))) {
		int1 *conf = (int1 *)xmalloc((slen+1)*sizeof(*conf));
		int i;

		str2conf(conf, slen, str);
		for (i = 0; i < slen ; i++) {
		    switch(r->base[i]) {
		    case 'a':
		    case 'A':
			r->prob_A[i] = conf[i];
			r->prob_C[i] = 0;
			r->prob_G[i] = 0;
			r->prob_T[i] = 0;
			break;
		    case 'c':
		    case 'C':
			r->prob_A[i] = 0;
			r->prob_C[i] = conf[i];
			r->prob_G[i] = 0;
			r->prob_T[i] = 0;
			break;
		    case 'g':
		    case 'G':
			r->prob_A[i] = 0;
			r->prob_C[i] = 0;
			r->prob_G[i] = conf[i];
			r->prob_T[i] = 0;
			break;
		    case 't':
		    case 'T':
			r->prob_A[i] = 0;
			r->prob_C[i] = 0;
			r->prob_G[i] = 0;
			r->prob_T[i] = conf[i];
			break;
		    default:
			r->prob_A[i] = conf[i];
			r->prob_C[i] = conf[i];
			r->prob_G[i] = conf[i];
			r->prob_T[i] = conf[i];
			break;
		    }
		}

		xfree(conf);
	    } else {
		memset(r->prob_A, 0, slen * sizeof(r->prob_A[0]));
		memset(r->prob_C, 0, slen * sizeof(r->prob_C[0]));
		memset(r->prob_G, 0, slen * sizeof(r->prob_G[0]));
		memset(r->prob_T, 0, slen * sizeof(r->prob_T[0]));
	    }
	}

	r->format = TT_EXP;
    }

    r->orig_trace = e;
    r->orig_trace_format = TT_EXP;
    r->orig_trace_free = (void (*)(void *))exp_destroy_info;

    return r;
}
#endif /* IOLIB_EXP */



/*
 * Takes an original read structure and a set of edit change arrays and
 * produces a new base position array incorporating all the edits. For
 * insertions, interpolation is used to derive a suitable sample position.
 *
 * INPUTS:
 *
 * Read   *r       = The original unedited read structure
 * int     Comp    = 0=Normal sequence, 1=Complemented sequence
 * int     Ned     = Length of edited arrays to follow
 * char   *edBases = Sequence of base characters incorporating ins/del edits
 * uint_2 *edPos   = Corresponding original base numbers, 0 indicates an
 *		     insertion. Base numbers start at 1.
 *
 * OUTPUTS:
 *
 * This array is assumed to be empty with an allocated length of Ned elements.
 *
 * uint_2* basePos = Base positions in samples
 */

void read_update_base_positions( Read *r, int Comp, int Ned, char *edBases,
				 int_2 *edPos, uint_2 *basePos )
{
    int     i, j;
    int     gap;
    int     delta;
    int     o_N;
    int     o_NPoints;
    uint_2* o_basePos;
    int     start;
    int     end;



    /* Check input */
    assert(r);
    assert(edBases);
    assert(edPos);
    assert(basePos);
    assert(Ned>0);
    if( (Ned<=0) || !r || !edBases || !edPos || !basePos )
	return;



    /* Original sequence data */
    o_N       = r->NBases;
    o_NPoints = r->NPoints;
    o_basePos = r->basePos;



    /* Copy original base positions */
    for( i=0; i<Ned; i++ )
    {
	/* If inserted base, set position to zero */
	if( edPos[i] == 0 )
	{
	    basePos[i] = 0;
	}
	else
	{
	    /* Get original position taking into account complementation */
	    basePos[i] = o_basePos[ Comp ? o_N-edPos[i]+1-1 : edPos[i]-1 ];
	}
    }



    /* Do linear interpolation to create positions for inserted bases */
    for( i=0; i<Ned; i++ )
    {
	/* Search for start */
	while( (basePos[i]!=0) && (i<Ned) )
	    i++;
	start = (i==0) ? 0 : basePos[i-1];



	/* Search for end */
	gap = 0;
	while( (basePos[i]==0) && (i<Ned) )
	{
	    gap++;
	    i++;
	}
	if( i == Ned )
	{
	    if( gap == 0 )
		return;
	    end = o_NPoints;
	}
	else
	{
	    end = basePos[i];
	}



	/* Do the interpolation */
	delta = (end - start) / (gap + 1);
	for( j=i-gap; j<i; j++ )
	{
	    /* Watch out for insertions at start (j==0) */
	    if( j==0 )
		basePos[j] = delta;
	    else
		basePos[j] = basePos[j-1] + delta;
	}
    }
}



/*
 * Takes a set of edit change arrays and produces a new set of confidence
 * arrays incorporating all the edits.
 *
 * INPUTS:
 *
 * int    Ned     = Length of edited arrays to follow
 * char*  edBases = Sequence of base characters incorporating ins/del edits
 * int1*  edConf  = Corresponding confidence values, 100 for insertions
 *
 *
 * OUTPUTS:
 *
 * These output arrays are assumed to be empty with an allocated length
 * of Ned elements each. The names and types are identical to the same
 * elements in the Read structure.
 *
 * char*  prob_A  = Base confidence A
 * char*  prob_C  = Base confidence C
 * char*  prob_G  = Base confidence G
 * char*  prob_T  = Base confidence T
 *
 */

void read_update_confidence_values( int Ned, char* edBases, int1* edConf,
                                    char* prob_A, char* prob_C, char* prob_G, char* prob_T )
{
    int  i;
    char c;



    /* Check input */
    assert(Ned>0);
    assert(edBases);
    assert(edConf);
    assert(prob_A);
    assert(prob_C);
    assert(prob_G);
    assert(prob_T);
    if( (Ned<=0) || !edBases || !edConf || !prob_A || !prob_C || !prob_G || !prob_T )
	return;



    /* Copy over confidence values */
    for( i=0; i<Ned; i++ )
    {
	c = (char) edConf[i];
	switch( edBases[i] )
	{
	    case 'A':
	    case 'a':
		prob_A[i] = c;
		prob_C[i] = 0;
		prob_G[i] = 0;
		prob_T[i] = 0;
		break;


	    case 'C':
	    case 'c':
		prob_A[i] = 0;
		prob_C[i] = c;
		prob_G[i] = 0;
		prob_T[i] = 0;
		break;

	    case 'G':
	    case 'g':
		prob_A[i] = 0;
		prob_C[i] = 0;
		prob_G[i] = c;
		prob_T[i] = 0;
		break;

	    case 'T':
	    case 't':
		prob_A[i] = 0;
		prob_C[i] = 0;
		prob_G[i] = 0;
		prob_T[i] = c;
		break;

	    default:
		prob_A[i] = c;
		prob_C[i] = c;
		prob_G[i] = c;
		prob_T[i] = c;
		break;
	}
    }
}



int read_sections(int new_sec) {
    static int sections = READ_ALL;

    if (new_sec)
	sections = new_sec;

    return sections;
}

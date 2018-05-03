/* Function for conversion of binary data to MSA.
 *
 * Contents:
 *    1. The <esl_msa_hmmpgmd2msa> function
 *    2. Test driver
 *    3. Copyright and license information
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <hmmer.h>
#include <hmmpgmd.h>
#include <esl_sqio.h>


/******************************************************************************
 *# 1. The <hmmpgmd2msa> function
 *****************************************************************************/



/* Function:  hmmpgmd2msa()
 * Synopsis:  Convert an HMMPGMD-derived data stream to an MSA, based
 *            on the corresponding hmm
 *
 * Purpose:   Given a data stream from HMMPGMD of the form shown
 *            here, produce an MSA:
 *                 HMMD_SEARCH_STATS
 *                 P7_HITS array of size (nhits) from above?
 *                 then repeats of P7_DOMAIN and P7_ALIDISPLAY data
 *                 for the hits, where each hit with d domains
 *                 produces
 *                   d P7_DOMAINs
 *                   then
 *                   d P7_ALIDISPLAYs
 *            ... optionally adding a sequence with length matching
 *            that of the hmm, which will be included in the alignment.
 *
 *			  A further extension has been the ability to include or exclude
 *            sequences form the list of hits.
 *
 *            This function's expected use is as a helper function for
 *            the hmmer website, which gets the above data stream from
 *            hmmpgmd.
 *
 * Args :     data: a pointer to binary data in the format given above
 *            hmm:  the HMM against which the alidisplay traces and
 *                  additional sequences/traces are threaded to reach
 *                  the returned msa.
 *            qsq : optional sequence to be included in the output msa;
 *                  must have the same number of residues as the hmm
 *                  has states, as each residue i will be aligned to
 *                  state i.
 *            incl: optional array of sequence names, in the case of
 *            		hmmpgmd a list of ints, which are are excluded due
 *            		to the sequence threshold, but have been selected
 *            		to be included in the alignment.  This ties in
 *            		with the way jackhmmer is implemented on the
 *            		HMMER website.
 *       incl_size: required size of the incl array. zero if incl is null.
 *       	  excl: optional array of sequence names, in the case of
 *            		hmmpgmd a list of ints, which are are included as they
 *            		score above threshold, but have been selected
 *            		to be excluded from the alignment.
 *       excl_size: required size of the excl array. zero if excl is null.
 *
 *
 * Returns:   Pointer to completed MSA object. NULL on error
 *
 */
int
hmmpgmd2msa(void *data, P7_HMM *hmm, ESL_SQ *qsq, int *incl, int incl_size, int *excl, int excl_size, ESL_MSA **ret_msa) {
  int i, j;
  int c;
  int status;
  int set_included;

  /* trace of the query sequence with N residues onto model with N match states */
  P7_TRACE          *qtr         = NULL;
  int                extra_sqcnt = 0;

  /* vars used to read from the binary data */
  HMMD_SEARCH_STATS *stats   = NULL;              /* pointer to a single stats object, at the beginning of data */
  P7_HIT            *hits    = NULL;              /* an array of hits, at the appropriate offset in data */

  /* vars used in msa construction */
  P7_TOPHITS         th;
  P7_ALIDISPLAY     *ad, *ad2;
  ESL_MSA           *msa   = NULL;
  P7_DOMAIN         *dom   = NULL;

  char              *p     = (char*)data;        /*pointer used to walk along data, must be char* to allow pointer arithmetic */

  th.N = 0;
  th.unsrt = NULL;
  th.hit   = NULL;

  /* optionally build a faux trace for the query sequence: relative to core model (B->M_1..M_L->E) */
  if (qsq != NULL) {
    if (qsq->n != hmm->M) {
      status = eslFAIL;
      goto ERROR;
    }

    if ((qtr = p7_trace_Create())                      == NULL)  {status = eslFAIL;  goto ERROR; }
    if ((status = p7_trace_Append(qtr, p7T_B, 0, 0))   != eslOK) goto ERROR;
    for (i = 1; i <= qsq->n; i++)
      if ((status = p7_trace_Append(qtr, p7T_M, i, i)) != eslOK) goto ERROR;
    if ((status = p7_trace_Append(qtr, p7T_E, 0, 0))   != eslOK) goto ERROR;
    qtr->M = qsq->n;
    qtr->L = qsq->n;
    extra_sqcnt = 1;
  }

  /* get search stats + hit info */
  stats = (HMMD_SEARCH_STATS*)p;

  /* sanity check */
  if (   ( stats->Z_setby != p7_ZSETBY_NTARGETS    && stats->Z_setby != p7_ZSETBY_OPTION    && stats->Z_setby != p7_ZSETBY_FILEINFO )
      || ( stats->domZ_setby != p7_ZSETBY_NTARGETS && stats->domZ_setby != p7_ZSETBY_OPTION && stats->domZ_setby != p7_ZSETBY_FILEINFO )
      ||   stats->nhits > 10000000
      ||   stats->elapsed > 1000000
  ) {
    status = eslFAIL;
    goto ERROR;
  }

  /* ok, it looks legitimate */
  p    += sizeof(HMMD_SEARCH_STATS);
  hits  = (P7_HIT*)p;
  p    += sizeof(P7_HIT) * stats->nhits;

  /* create a tophits object, to be passed to p7_tophits_Alignment() */
  ESL_ALLOC( th.unsrt, sizeof(P7_HIT) * stats->nhits);
  memcpy( th.unsrt, hits, sizeof(P7_HIT) * stats->nhits);
  ESL_ALLOC( th.hit, sizeof(P7_HIT*) * stats->nhits);
  for (i=0; i<stats->nhits; i++) {
    th.hit[i] = &(th.unsrt[i]);
    if (   th.hit[i]->ndom > 10000
        || th.hit[i]->flags >  p7_IS_INCLUDED + p7_IS_REPORTED + p7_IS_NEW + p7_IS_DROPPED + p7_IS_DUPLICATE
    ) {
      status = eslFAIL;
      goto ERROR;
    }
  }

//  th.unsrt     = NULL;
  th.N         = stats->nhits;
  th.nreported = 0;
  th.nincluded = 0;
  th.is_sorted_by_sortkey = 0;
  th.is_sorted_by_seqidx  = 0;

  for (i = 0; i < th.N; i++) {
    ESL_ALLOC( th.hit[i]->dcl, sizeof(P7_DOMAIN) *  th.hit[i]->ndom);
    /* Go through the hits and set to be excluded or included as necessary */
    set_included = 0;
    if(th.hit[i]->flags & p7_IS_INCLUDED){
      if(excl_size > 0){
        for( c = 0; c < excl_size; c++){
          if(excl[c] == (long)(th.hit[i]->name) ){
            th.hit[i]->flags = p7_IS_DROPPED;
            th.hit[i]->nincluded = 0;
            break;
          }
        }
      }
    }else{
      if(incl_size > 0){
    	for( c = 0; c < incl_size; c++){
          if(incl[c] == (long)th.hit[i]->name ){
            th.hit[i]->flags = p7_IS_INCLUDED;
            set_included = 1;
          }
        }
      }
    }
    /* first grab all the P7_DOMAINs for the hit */
    for (j=0; j < th.hit[i]->ndom; j++) {
      dom = th.hit[i]->dcl + j;
      memcpy(dom , (P7_DOMAIN*)p, sizeof(P7_DOMAIN));
      /* Possibly set domains to be include if being
       * externally set via incl list*/
      if(set_included) th.hit[i]->dcl[j].is_included = 1;
      p += sizeof(P7_DOMAIN);
    }
    /* then grab the P7_ALIDISPLAYs for the hit */
    for (j=0; j < th.hit[i]->ndom; j++) {
      ad = (P7_ALIDISPLAY*)p;
      ESL_ALLOC(th.hit[i]->dcl[j].ad, sizeof(P7_ALIDISPLAY));
      ad2 = th.hit[i]->dcl[j].ad;

      ad2->memsize = ad->memsize;
      ad2->rfline = ad->rfline;
      ad2->mmline = ad->mmline;
      ad2->csline = ad->csline ;
      ad2->model  = ad->model ;
      ad2->mline  = ad->mline ;
      ad2->aseq   = ad->aseq ;
      ad2->ppline = ad->ppline;
      ad2->N      = ad->N;

      ad2->hmmname = ad->hmmname;
      ad2->hmmacc  = ad->hmmacc ;
      ad2->hmmdesc = ad->hmmdesc;
      ad2->hmmfrom = ad->hmmfrom;
      ad2->hmmto   = ad->hmmto;
      ad2->M       = ad->M;

      ad2->sqname  = ad->sqname;
      ad2->sqacc   = ad->sqacc ;
      ad2->sqdesc  = ad->sqdesc;
      ad2->sqfrom  = ad->sqfrom;
      ad2->sqto    = ad->sqto;
      ad2->L       = ad->L;

      p += sizeof(P7_ALIDISPLAY);

      ESL_ALLOC(ad2->mem, ad2->memsize);
      
      memcpy(ad2->mem, p, ad->memsize);
      
      p += ad2->memsize;
      p7_alidisplay_Deserialize(ad2);
    }
  }


  /* use the tophits and trace info above to produce an alignment */
  if ( (status = p7_tophits_Alignment(&th, hmm->abc, &qsq, &qtr, extra_sqcnt, p7_ALL_CONSENSUS_COLS, &msa)) != eslOK) goto ERROR;


  /* free memory */
  if (qtr != NULL) free(qtr);
  for (i = 0; i < th.N; i++) {
    for (j=0; j < th.hit[i]->ndom; j++)
      p7_alidisplay_Destroy(th.hit[i]->dcl[j].ad);

    if (th.hit[i]->dcl != NULL) free (th.hit[i]->dcl);
  }
  if (th.unsrt != NULL) free (th.unsrt);
  if (th.hit != NULL) free (th.hit);

  *ret_msa = msa;
  return eslOK;

ERROR:
  /* free memory */
  if (qtr != NULL) free(qtr);

  for (i = 0; i < th.N; i++) {
    for (j=0; j < th.hit[i]->ndom; j++)
      p7_alidisplay_Destroy(th.hit[i]->dcl[j].ad);

    if (th.hit[i]->dcl != NULL) free (th.hit[i]->dcl);
  }
  if (th.unsrt != NULL) free (th.unsrt);
  if (th.hit != NULL) free (th.hit);

  return status;
}

/******************************************************************************
 *# 2. The <hmmpgmd2stats> function
 *****************************************************************************/



/* Function:  hmmpgmd2stats()
 * Synopsis:  Use a HMMPGMD-derived data stream to extract some simple
 *            statistics regarding its alignment.
 * Purpose:   Given a data stream from HMMPGMD of the form shown
 *            here, produce a vector of floats:
 *                 positions 0 to hmm->M-1 are:
 *                   Fraction of alignments which cover model at that position
 *                 positions hmm->M to hmm->M*2-1
 *                   Fraction of alignments which cover model at that position (mod hmm->M)     
 *                   with a similar residue
 *                 positions hmm->M*2 to hmm->M*3-1
 *                   Fraction of alignments which cover model at that position (mod hmm->M)
 *                   with the consensus residue
 *
 * Args :     data: a pointer to binary data in the format given above
 *            hmm:  the HMM against which the alidisplay traces and
 *                  additional sequences/traces are threaded to reach
 *                  the returned msa.
 *            
 * Returns:   the location where the output vector will be placed.
 *                      caller is responsible for freeing it later   
 *
 */
int hmmpgmd2stats(void *data, P7_HMM *hmm, float** statsOut) 
{
  int i, j, k;
  int status;

  /* trace of the query sequence with N residues onto model with N match states */
  P7_TRACE          *qtr         = NULL;

  /* vars used to read from the binary data */
  HMMD_SEARCH_STATS *stats   = NULL;              /* pointer to a single stats object, at the beginning of data */
  P7_HIT            *hits    = NULL;              /* an array of hits, at the appropriate offset in data */

  P7_TOPHITS         th;
  P7_DOMAIN         *dom;
  P7_ALIDISPLAY     *ad, *ad2;

  int *cover, *id, *similar; //store statistics result per hit
  int readPos, writePos;     //for converting alignment contents into model indexing

  char              *p     = (char*)data;        /*pointer used to walk along data, must be char* to allow pointer arithmetic */

  th.N = 0;
  th.unsrt = NULL;
  th.hit   = NULL;

  //storage for output
  ESL_ALLOC( *statsOut,   sizeof(float) * hmm->M * 3);

  //storage for accumulation per hit
  ESL_ALLOC( cover,   sizeof(int) * hmm->M);
  ESL_ALLOC( id,      sizeof(int) * hmm->M);
  ESL_ALLOC( similar, sizeof(int) * hmm->M);
  for(k = 0; k < hmm->M; k++)
  {
    cover[k] = 0;
    id[k] = 0;
    similar[k] = 0;
    
    (*statsOut)[k         ] = 0;
    (*statsOut)[k+hmm->M  ] = 0;
    (*statsOut)[k+hmm->M*2] = 0;
  }

  /* get search stats + hit info */
  stats = (HMMD_SEARCH_STATS*)p;

  /* sanity check */
  if (   ( stats->Z_setby != p7_ZSETBY_NTARGETS    && stats->Z_setby != p7_ZSETBY_OPTION    && stats->Z_setby != p7_ZSETBY_FILEINFO )
      || ( stats->domZ_setby != p7_ZSETBY_NTARGETS && stats->domZ_setby != p7_ZSETBY_OPTION && stats->domZ_setby != p7_ZSETBY_FILEINFO )
      ||   stats->nhits > 10000000
      ||   stats->elapsed > 1000000
  ) {
    status = eslFAIL;
    goto ERROR;
  }

  /* ok, it looks legitimate */
  p    += sizeof(HMMD_SEARCH_STATS);
  hits  = (P7_HIT*)p;
  p    += sizeof(P7_HIT) * stats->nhits;

  /* create a tophits object, use it to step through the alignments */
  ESL_ALLOC( th.unsrt, sizeof(P7_HIT) * stats->nhits);
  memcpy( th.unsrt, hits, sizeof(P7_HIT) * stats->nhits);
  ESL_ALLOC( th.hit, sizeof(P7_HIT*) * stats->nhits);
  for (i=0; i<stats->nhits; i++) {
    th.hit[i] = &(th.unsrt[i]);
    if (   th.hit[i]->ndom > 10000
        || th.hit[i]->flags >  p7_IS_INCLUDED + p7_IS_REPORTED + p7_IS_NEW + p7_IS_DROPPED + p7_IS_DUPLICATE
    ) {
      status = eslFAIL;
      goto ERROR;
    }
  }

  th.N         = stats->nhits;
  th.nreported = 0;
  th.nincluded = 0;
  th.is_sorted_by_sortkey = 0;
  th.is_sorted_by_seqidx  = 0;

  for (i = 0; i < th.N; i++) 
  {
    ESL_ALLOC( th.hit[i]->dcl, sizeof(P7_DOMAIN) *  th.hit[i]->ndom);
   
    if(th.hit[i]->flags & p7_IS_INCLUDED) th.nincluded++;

    /* first grab all the P7_DOMAINs for the hit */
    for (j=0; j < th.hit[i]->ndom; j++) 
    {
      dom = (P7_DOMAIN*)p;
      th.hit[i]->dcl[j].is_included = dom->is_included;
      p += sizeof(P7_DOMAIN);
    }
    
    /* then grab the P7_ALIDISPLAYs for the hit */
    for (j=0; j < th.hit[i]->ndom; j++) 
    {
      ad = (P7_ALIDISPLAY*)p;
      ESL_ALLOC(th.hit[i]->dcl[j].ad, sizeof(P7_ALIDISPLAY));
      ad2 = th.hit[i]->dcl[j].ad;

      ad2->memsize = ad->memsize;
      ad2->rfline = ad->rfline;
      ad2->mmline = ad->mmline;
      
      ad2->csline = ad->csline ;
      ad2->model  = ad->model ;
      ad2->mline  = ad->mline ;
      ad2->aseq   = ad->aseq ;
      ad2->ppline = ad->ppline;
      ad2->N      = ad->N;

      ad2->hmmname = ad->hmmname;
      ad2->hmmacc  = ad->hmmacc ;
      ad2->hmmdesc = ad->hmmdesc;
      ad2->hmmfrom = ad->hmmfrom;
      ad2->hmmto   = ad->hmmto;
      ad2->M       = ad->M;

      ad2->sqname  = ad->sqname;
      ad2->sqacc   = ad->sqacc ;
      ad2->sqdesc  = ad->sqdesc;
      ad2->sqfrom  = ad->sqfrom;
      ad2->sqto    = ad->sqto;
      ad2->L       = ad->L;    
      
      p += sizeof(P7_ALIDISPLAY);

      ESL_ALLOC(ad2->mem, ad2->memsize);
      
      memcpy(ad2->mem, p, ad->memsize);
      
      p += ad2->memsize;
      
      p7_alidisplay_Deserialize(ad2);
      
      if(th.hit[i]->flags & p7_IS_INCLUDED && th.hit[i]->dcl[j].is_included)
      {
        writePos = ad2->hmmfrom-1;  
        readPos = 0;
        while(readPos < ad2->N)
        {
          //check if model covers residue
          if(isupper(ad2->aseq[readPos]) || ad2->aseq[readPos] == '-')
          {          
            cover[writePos]++;
          
            //check mline for id
            if(isalpha(ad2->mline[readPos]))
            {
              id[writePos]++;
              similar[writePos]++;
            }
            //check mline for not-a-space
            else if(ad2->mline[readPos] == '+')
            {
              similar[writePos]++;
            }
            writePos++;
          }       
          readPos++;
        }
      }
    }
      
    //increment output, adjusting for overlaps
    for(k = 0; k < hmm->M; k++)
    {
      if(cover[k]) (*statsOut)[k]+=1.0;
      
      if(id[k]) (*statsOut)[k+hmm->M]+=(id[k]/cover[k]);
      id[k] = 0;
      
      if(similar[k]) (*statsOut)[k+hmm->M*2]+=(similar[k]/cover[k]);
      similar[k] = 0;
      
      cover[k] = 0;
    }
    
  }

  for(i = 0; i < hmm->M*3; i++)
  {
    (*statsOut)[i] = (*statsOut)[i]/(th.nincluded);
  }

  for(i = hmm->M; i < hmm->M*3; i++)
  {
    if((*statsOut)[i%hmm->M])
    {  
      (*statsOut)[i] = (*statsOut)[i]/(*statsOut)[i%hmm->M];
    }
    else
    {
      (*statsOut)[i] = 0.0;
    }
  }

  
  /* free memory */
  if (qtr != NULL) free(qtr);
  qtr = NULL;
  
  for (i = 0; i < th.N; i++) {
    for (j=0; j < th.hit[i]->ndom; j++)
      if(th.hit[i]->dcl[j].ad)
      {
        p7_alidisplay_Destroy(th.hit[i]->dcl[j].ad);
        th.hit[i]->dcl[j].ad = NULL;
      }
      
    if (th.hit[i]->dcl != NULL) free (th.hit[i]->dcl);
    th.hit[i]->dcl = NULL;
  }
  if (th.unsrt != NULL) free (th.unsrt);
  th.unsrt = NULL;
  if (th.hit != NULL) free (th.hit);
  th.hit = NULL;

  return eslOK;

ERROR:
  /* free memory */
  if (qtr != NULL) free(qtr);
  qtr = NULL;
  
  for (i = 0; i < th.N; i++) {
    for (j=0; j < th.hit[i]->ndom; j++)
      if(th.hit[i]->dcl[j].ad)
      {
        p7_alidisplay_Destroy(th.hit[i]->dcl[j].ad);
        th.hit[i]->dcl[j].ad = NULL;
      }
      
    if (th.hit[i]->dcl != NULL) free (th.hit[i]->dcl);
    th.hit[i]->dcl = NULL;
  }
  if (th.unsrt != NULL) free (th.unsrt);
  th.unsrt = NULL;
  if (th.hit != NULL) free (th.hit);
  th.hit = NULL;
  

  return status;
}




/*****************************************************************
 * 3. Test driver
 *****************************************************************/

//#define hmmpgmd2msa_TESTDRIVE
#ifdef hmmpgmd2msa_TESTDRIVE

//gcc -o alistat_test -msse2 -std=gnu99 -g -O2 -I. -L. -I../easel -L../easel -D hmmpgmd2msa_TESTDRIVE hmmpgmd2msa.c -lhmmer -leasel -lm

/* Test driver. As written, requires files that won't be released with
 * the distribution. So it should be replaced with a tighter test.
 */
int
main(int argc, char **argv) {
  ESL_MSA           *msa   = NULL;
  ESL_SQ            *qsq   = NULL;
  ESL_SQFILE        *qfp   = NULL;              /* open qfile                                      */
  P7_HMMFILE        *hfp   = NULL;              /* open input HMM file                             */
  P7_HMM            *hmm   = NULL;              /* one HMM query                                   */
  ESL_ALPHABET      *abc   = NULL;              /* digital alphabet                                */
  FILE              *fp    = NULL;              /* open file containing the HMMPGMD data */
  void              *data  = NULL;              /* pointer to the full data stream built as if from hmmdmstr */
  long size = 0;

  char   errbuf[eslERRBUFSIZE];
  int status;
  
  float* statsOut;
  int x;
  //char *badstring = "asdlfuhasdfuhasdfhasdfhaslidhflaishdfliasuhdfliasuhdfliasudfh";

  char *hmm_file = "esl_align.hmm";
  char *fa_file  = "esl_align.fa";
  char *dat_file = "esl_align.big.bin";

  if (argc > 1 ) {
    ESL_ALLOC(hmm_file, sizeof(char) * (strlen(argv[1])+1) );
    strcpy(hmm_file, argv[1]);
  }
  if (argc > 2 ) {
    ESL_ALLOC(dat_file, sizeof(char) * (strlen(argv[2])+1) );
    strcpy(dat_file, argv[2]);
  }
  if (argc > 3 ) {
    ESL_ALLOC(fa_file, sizeof(char) * (strlen(argv[3])+1) );
    strcpy(fa_file, argv[3]);
  }

  printf("hmmpgmd2msa:\nhmm: %s\nfa:  %s\ndat: %s\n", hmm_file, fa_file, dat_file);


  /* read the hmm */
  if ( (status = p7_hmmfile_OpenE(hmm_file, NULL, &hfp, errbuf)) != 0 ) goto ERROR;
  if ( (status = p7_hmmfile_Read(hfp, &abc, &hmm)) != 0 ) goto ERROR;

  /* read the query sequence */
  //if ( (status = esl_sqfile_OpenDigital(abc, fa_file, eslSQFILE_UNKNOWN, NULL, &qfp)) != 0) goto ERROR;
  //qsq = esl_sq_CreateDigital(abc);
  //if ( (status = esl_sqio_Read(qfp, qsq)) != eslOK)  goto ERROR;

  //printf("sequence length %d\n", qsq->n);

  /* get stats for the hmmd data */

  if ( (fp = fopen(dat_file, "rb")) == NULL ) goto ERROR;

  fseek (fp , 0 , SEEK_END);
  size = ftell (fp);
  rewind (fp);
  ESL_ALLOC(data, size);
  fread(data, size, 1, fp);

  status = hmmpgmd2stats(data, hmm, &statsOut);
  //status = hmmpgmd2msa(data, hmm, qsq, NULL,0, NULL, 0, &msa);

  for(x = 0; x < hmm->M*3; x++)
  {
    if(statsOut[x] > 1)
    {
      printf("problem x: %d %f\n", x, statsOut[x]);
    }
  }

  for(x = 0; x < hmm->M; x++)
  {
    printf("%d", ((int)(10*statsOut[x])>=10)?9:(int)(10*statsOut[x]));
  }
  printf("\n");
  for(x = hmm->M; x < hmm->M*2; x++)
  {
    printf("%d", ((int)(10*statsOut[x])>=10)?9:(int)(10*statsOut[x]));
  }
  printf("\n");
  for(x = hmm->M*2; x < hmm->M*3; x++)
  {
    printf("%d", ((int)(10*statsOut[x])>=10)?9:(int)(10*statsOut[x]));
  }
  printf("\n");

  if (status != eslOK) goto ERROR;

  //esl_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM);

  exit(0);

ERROR:
  printf ("fail!\n");
  exit(1);

}
#endif

/************************************************************
 * @LICENSE@
 ************************************************************/



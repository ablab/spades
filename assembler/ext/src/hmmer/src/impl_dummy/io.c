/* Saving optimized profiles in two pieces: MSV part and the rest.
 * 
 * To accelerate hmmscan, which is limited by speed of HMM input,
 * hmmpress saves an optimized profile in two pieces. One file gets
 * a bare minimum of information needed to run the MSV filter.
 * The other file gets the rest of the profile. Both files are binary,
 * stored exactly as the <P7_OPROFILE> has the information internally.
 * 
 * By convention, hmmpress calls the two files <hmmfile>.h3f and
 * <hmmfile>.h3p, which nominally stand for "H3 filter" and "H3
 * profile".
 * 
 * Contents:
 *    1. Writing optimized profiles to two files.
 *    2. Reading optimized profiles in two stages.
 *    3. Utility routines.
 *    4. Benchmark driver.
 *    5. Unit tests.
 *    6. Test driver.
 *    7. Example.
 *    8. Copyright and license information.
 *    
 * TODO:
 *    - crossplatform binary compatibility (endedness and off_t)
 *    - Write() could save a tag (model #) instead of name for verifying
 *      that MSV and Rest parts match, saving a malloc for var-lengthed name
 *      in ReadRest().
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"

#include "hmmer.h"
#include "impl_dummy.h"

static uint32_t  v3f_fmagic = 0xb3e6e6e4; /* 3/f binary MSV file, dummy:     "3ffd" = 0x 33 66 66 64  + 0x80808080 */
static uint32_t  v3f_pmagic = 0xb3e6f0e4; /* 3/f binary profile file, dummy: "3fpd" = 0x 33 66 70 64  + 0x80808080 */

static uint32_t  v3e_fmagic = 0xb3e5e6e4; /* 3/e binary MSV file, dummy:     "3efd" = 0x 33 65 66 64  + 0x80808080 */
static uint32_t  v3e_pmagic = 0xb3e5f0e4; /* 3/e binary profile file, dummy: "3epd" = 0x 33 65 70 64  + 0x80808080 */

static uint32_t  v3d_fmagic = 0xb3e4e6e4; /* 3/d binary MSV file, dummy:     "3dfd" = 0x 33 64 66 64  + 0x80808080 */
static uint32_t  v3d_pmagic = 0xb3e4f0e4; /* 3/d binary profile file, dummy: "3dpd" = 0x 33 64 70 64  + 0x80808080 */

static uint32_t  v3c_fmagic = 0xb3e3e6e4; /* 3/c binary MSV file, dummy:     "3cfd" = 0x 33 63 66 64  + 0x80808080 */
static uint32_t  v3c_pmagic = 0xb3e3f0e4; /* 3/c binary profile file, dummy: "3cpd" = 0x 33 63 70 64  + 0x80808080 */

static uint32_t  v3b_fmagic = 0xb3e2e6e4; /* 3/b binary MSV file, dummy:     "3bfd" = 0x 33 62 66 64  + 0x80808080 */
static uint32_t  v3b_pmagic = 0xb3e2f0e4; /* 3/b binary profile file, dummy: "3bpd" = 0x 33 62 70 64  + 0x80808080 */

static uint32_t  v3a_fmagic = 0xe8b3e6e4; /* 3/a binary MSV file, dummy:     "h3fd" = 0x 68 33 66 64  + 0x80808080 */
static uint32_t  v3a_pmagic = 0xe8b3f0e4; /* 3/a binary profile file, dummy: "h3pd" = 0x 68 33 70 64  + 0x80808080 */

/*****************************************************************
 *# 1. Writing optimized profiles to two files.
 *****************************************************************/

/* Function:  p7_oprofile_Write()
 * Synopsis:  Write an optimized profile in two files.
 *
 * Purpose:   Write the MSV filter part of <om> to open binary stream
 *            <ffp>, and the rest of the model to <pfp>. These two
 *            streams will typically be <.h3f> and <.h3p> files 
 *            being created by hmmpress.
 *
 * Args:      ffp  - open binary stream for saving MSV filter part
 *            pfp  - open binary stream for saving rest of profile
 *            om   - optimized profile to save
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEWRITE> on any write failure; for example,
 *            if disk is full. 
 */
int
p7_oprofile_Write(FILE *ffp, FILE *pfp, P7_OPROFILE *om)
{
  /* <ffp> is the part of the oprofile that MSVFilter() needs */
  if (fwrite((char *) &(v3f_fmagic),    sizeof(uint32_t), 1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->M),         sizeof(int),      1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->abc->type), sizeof(int),      1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) om->offs,         sizeof(off_t),    p7_NOFFSETS, ffp) != p7_NOFFSETS) ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(v3f_fmagic),    sizeof(uint32_t), 1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed"); /* sentinel */

  /* <pfp> gets the rest of the oprofile */
  if (fwrite((char *) &(v3f_pmagic),    sizeof(uint32_t), 1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->M),         sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->abc->type), sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(v3f_pmagic),    sizeof(uint32_t), 1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed"); /* sentinel */
  
  return eslOK;
}
/*---------------- end, writing oprofile ------------------------*/


/*****************************************************************
 * 2. Reading optimized profiles in two stages.
 *****************************************************************/

/* Function:  p7_oprofile_ReadMSV()
 * Synopsis:  Read MSV filter part of an optimized profile.
 *
 * Purpose:   Read the MSV filter part of a profile from the
 *            <.h3f> file associated with an open HMM file <hfp>.
 *            Allocate a new model, populate it with this minimal
 *            MSV filter information, and return a pointer to it
 *            in <*ret_om>. 
 *            
 *            Our alphabet may get set by the first HMM we read.  If
 *            <*byp_abc> is <NULL> at start, create a new alphabet and
 *            return a pointer to it in <*byp_abc>. If <*byp_abc> is
 *            non-<NULL>, it is assumed to be a pointer to an existing
 *            alphabet; we verify that the HMM's alphabet matches it
 *            and <*ret_abc> isn't changed.  This is the same
 *            convention used by <p7_hmmfile_Read()>.
 *            
 *            The <.h3f> file was opened automatically, if it existed,
 *            when the HMM file was opened with <p7_hmmfile_OpenE()>.
 *            
 *            When no more HMMs remain in the file, return <eslEOF>.
 *
 * Args:      hfp     - open HMM file, with associated .h3p file
 *            byp_abc - BYPASS: <*byp_abc == ESL_ALPHABET *> if known; 
 *                              <*byp_abc == NULL> if desired; 
 *                              <NULL> if unwanted.
 *            ret_om  - RETURN: newly allocated <om> with MSV filter
 *                      data filled in.
 *            
 * Returns:   <eslOK> on success. <*ret_om> is allocated here;
 *            caller free's with <p7_oprofile_Destroy()>.
 *            <*byp_abc> is allocated here if it was requested;
 *            caller free's with <esl_alphabet_Destroy()>.
 *            
 *            Returns <eslEFORMAT> if <hfp> has no <.h3f> file open,
 *            or on any parsing error.
 *            
 *            Returns <eslEINCOMPAT> if the HMM we read is incompatible
 *            with the existing alphabet <*byp_abc> led us to expect.
 *            
 *            On any returned error, <hfp->errbuf> contains an
 *            informative error message.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_oprofile_ReadMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om)
{
  uint32_t       magic;
  P7_BG         *bg      = NULL;
  P7_PROFILE    *gm      = NULL;
  P7_OPROFILE   *om      = NULL;
  ESL_ALPHABET  *abc     = NULL;
  P7_HMM        *hmm     = NULL;
  int            alphatype;
  off_t          offs[p7_NOFFSETS];
  off_t          roff;
  int            M;
  int            i;
  int            status;

  if (hfp->errbuf != NULL) hfp->errbuf[0] = '\0';
  if (hfp->ffp == NULL) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no MSV profile file; hmmpress probably wasn't run");
  if (feof(hfp->ffp))   { status = eslEOF; goto ERROR; }	/* normal EOF: no more profiles */
  
  /* keep track of the starting offset of the MSV model */
  roff = ftello(hfp->ffp);

  if (! fread( (char *) &magic,     sizeof(uint32_t), 1, hfp->ffp)) { status = eslEOF; goto ERROR; }
  if (magic == v3a_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/a); please hmmpress your HMM file again");
  if (magic == v3b_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/b); please hmmpress your HMM file again");
  if (magic == v3c_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/c); please hmmpress your HMM file again");
  if (magic == v3d_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/d); please hmmpress your HMM file again");
  if (magic == v3e_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/e); please hmmpress your HMM file again");
  if (magic != v3f_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic; not an HMM database?");

  if (! fread( (char *) &M,         sizeof(int),      1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model size M");
  if (! fread( (char *) &alphatype, sizeof(int),      1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read alpha type");
  if (! fread( (char *) offs,       sizeof(off_t),    p7_NOFFSETS, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read hmmpfam offsets");
  if (! fread( (char *) &magic,     sizeof(uint32_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no sentinel magic: .h3f file corrupted?");
  if (magic != v3f_fmagic)                                                    ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad sentinel magic; .h3f file corrupted?");

  /* Set or verify alphabet. */
  if (byp_abc == NULL || *byp_abc == NULL)	{	/* alphabet unknown: whether wanted or unwanted, make a new one */
    if ((abc = esl_alphabet_Create(alphatype)) == NULL)  ESL_XFAIL(eslEMEM, hfp->errbuf, "allocation failed: alphabet");
  } else {			/* alphabet already known: verify it against what we see in the HMM */
    abc = *byp_abc;
    if (abc->type != alphatype) 
      ESL_XFAIL(eslEINCOMPAT, hfp->errbuf, "Alphabet type mismatch: was %s, but current profile says %s", 
		esl_abc_DecodeType(abc->type), esl_abc_DecodeType(alphatype));
  }

  if ((status = p7_hmmfile_Position(hfp, offs[p7_MOFFSET])) != eslOK) goto ERROR;
  if ((status = p7_hmmfile_Read(hfp, &abc, &hmm))           != eslOK) goto ERROR;

  bg  = p7_bg_Create(hmm->abc);
  gm = p7_profile_Create(hmm->M, hmm->abc);
  p7_ProfileConfig(hmm, bg, gm, hmm->M, p7_LOCAL);
  om = p7_oprofile_Create(hmm->M, abc);

  if ((status = p7_oprofile_Convert(gm, om)) != eslOK) goto ERROR;

  om->M     = M;
  om->roff  = roff;
  for (i = 0; i < p7_NOFFSETS; ++i)  om->offs[i] = offs[i];

  /* keep track of the ending offset of the MSV model */
  om->eoff = ftello(hfp->ffp) - 1;;

  if (byp_abc != NULL)  *byp_abc = abc;

  if (gm != NULL) p7_oprofile_Destroy(gm);

  *ret_om = om;
  return eslOK;

 ERROR:
  if (abc != NULL && (byp_abc == NULL || *byp_abc == NULL)) esl_alphabet_Destroy(abc); /* destroy alphabet if we created it here */
  if (gm  != NULL) p7_oprofile_Destroy(gm);
  *ret_om = NULL;
  return status;
}


/* Function:  p7_oprofile_ReadInfoMSV()
 * Synopsis:  Read MSV filter info, but not the scores.
 *
 * Purpose:   Read just enough of the MSV filter header from the
 *            <.h3f> file associated with an open HMM file <hfp>
 *            to skip ahead to the next MSV filter. Allocate a new
 *            model, populate it with just the file offsets of this
 *            model and return a pointer to it in <*ret_om>. 
 *            
 *            The <.h3f> file was opened automatically, if it existed,
 *            when the HMM file was opened with <p7_hmmfile_OpenE()>.
 *            
 *            When no more HMMs remain in the file, return <eslEOF>.
 *
 * Args:      hfp     - open HMM file, with associated .h3p file
 *            byp_abc - BYPASS: <*byp_abc == ESL_ALPHABET *> if known; 
 *                              <*byp_abc == NULL> if desired; 
 *                              <NULL> if unwanted.
 *            ret_om  - RETURN: newly allocated <om> with partial MSV
 *                      filter data filled in.
 *            
 * Returns:   <eslOK> on success. <*ret_om> is allocated here;
 *            caller free's with <p7_oprofile_Destroy()>.
 *            <*byp_abc> is allocated here if it was requested;
 *            caller free's with <esl_alphabet_Destroy()>.
 *            
 *            Returns <eslEFORMAT> if <hfp> has no <.h3f> file open,
 *            or on any parsing error.
 *            
 *            Returns <eslEINCOMPAT> if the HMM we read is incompatible
 *            with the existing alphabet <*byp_abc> led us to expect.
 *            
 *            On any returned error, <hfp->errbuf> contains an
 *            informative error message.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_oprofile_ReadInfoMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om)
{
  uint32_t       magic;
  P7_OPROFILE   *om      = NULL;
  ESL_ALPHABET  *abc     = NULL;
  int            alphatype;
  off_t          offs[p7_NOFFSETS];
  off_t          roff;
  int            M;
  int            i;
  int            status;

  if (hfp->errbuf != NULL) hfp->errbuf[0] = '\0';
  if (hfp->ffp == NULL) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no MSV profile file; hmmpress probably wasn't run");
  if (feof(hfp->ffp))   { status = eslEOF; goto ERROR; }	/* normal EOF: no more profiles */
  
  /* keep track of the starting offset of the MSV model */
  roff = ftello(hfp->ffp);

  if (! fread( (char *) &magic,     sizeof(uint32_t), 1, hfp->ffp)) { status = eslEOF; goto ERROR; }
  if (magic == v3a_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/a); please hmmpress your HMM file again");
  if (magic == v3b_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/b); please hmmpress your HMM file again");
  if (magic == v3c_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/c); please hmmpress your HMM file again");
  if (magic == v3d_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/d); please hmmpress your HMM file again");
  if (magic == v3e_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/e); please hmmpress your HMM file again");
  if (magic != v3f_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic; not an HMM database?");

  if (! fread( (char *) &M,         sizeof(int),      1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model size M");
  if (! fread( (char *) &alphatype, sizeof(int),      1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read alpha type");
  if (! fread( (char *) offs,       sizeof(off_t),    p7_NOFFSETS, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read hmmpfam offsets");
  if (! fread( (char *) &magic,     sizeof(uint32_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no sentinel magic: .h3f file corrupted?");
  if (magic != v3f_fmagic)                                                    ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad sentinel magic; .h3f file corrupted?");

  /* Set or verify alphabet. */
  if (byp_abc == NULL || *byp_abc == NULL)	{	/* alphabet unknown: whether wanted or unwanted, make a new one */
    if ((abc = esl_alphabet_Create(alphatype)) == NULL)  ESL_XFAIL(eslEMEM, hfp->errbuf, "allocation failed: alphabet");
  } else {			/* alphabet already known: verify it against what we see in the HMM */
    abc = *byp_abc;
    if (abc->type != alphatype) 
      ESL_XFAIL(eslEINCOMPAT, hfp->errbuf, "Alphabet type mismatch: was %s, but current profile says %s", 
		esl_abc_DecodeType(abc->type), esl_abc_DecodeType(alphatype));
  }

  om = p7_oprofile_Create(M, abc);

  om->M     = M;
  om->roff  = roff;
  for (i = 0; i < p7_NOFFSETS; ++i)  om->offs[i] = offs[i];

  /* keep track of the ending offset of the MSV model */
  om->eoff = ftello(hfp->ffp) - 1;;

  if (byp_abc != NULL)  *byp_abc = abc;

  *ret_om = om;
  return eslOK;

 ERROR:
  if (abc != NULL && (byp_abc == NULL || *byp_abc == NULL)) esl_alphabet_Destroy(abc); /* destroy alphabet if we created it here */
  *ret_om = NULL;
  return status;
}


/* Function:  p7_oprofile_ReadBlockMSV()
 * Synopsis:  Read the next block of optimized profiles from a hmm file.
 *
 * Purpose:   Reads a block of optimized profiles from open hmm file <hfp> into 
 *            <hmmBlock>.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <sqBlock>.
 * 
 *            Returns <eslEOF> when there is no profiles left in the
 *            file (including first attempt to read an empty file).
 * 
 *            Otherwise return the status of the p7_oprofile_ReadMSV function.
 */
int
p7_oprofile_ReadBlockMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OM_BLOCK *hmmBlock)
{
  int     i;
  int     size = 0;
  int     status = eslOK;

  hmmBlock->count = 0;
  for (i = 0; i < hmmBlock->listSize; ++i)
    {
      status = p7_oprofile_ReadMSV(hfp, byp_abc, &hmmBlock->list[i]);
      if (status != eslOK) break;
      size += hmmBlock->list[i]->M;
      ++hmmBlock->count;
    }

  /* EOF will be returned only in the case were no profiles were read */
  if (status == eslEOF && i > 0) status = eslOK;

  return status;
}

/* Function:  p7_oprofile_ReadRest()
 * Synopsis:  Read the rest of an optimized profile.
 *
 * Purpose:   Read the rest of an optimized profile <om> from
 *            the <.h3p> file associated with an open HMM
 *            file <hfp>. 
 *            
 *            This is the second part of a two-part calling sequence.
 *            The <om> here must be the result of a previous
 *            successful <p7_oprofile_ReadMSV()> call on the same
 *            open <hfp>.
 *
 * Args:      hfp - open HMM file, from which we've previously
 *                  called <p7_oprofile_ReadMSV()>.
 *            om  - optimized profile that was successfully
 *                  returned by  <p7_oprofile_ReadMSV()>.
 *
 * Returns:   <eslOK> on success, and <om> is now a complete
 *            optimized profile.
 *            
 *            Returns <eslEFORMAT> if <hfp> has no <.h3p> file open,
 *            or on any parsing error, and set <hfp->errbuf> to
 *            an informative error message.
 *
 * Throws:    <eslESYS> if an <fseek()> fails to reposition the
 *            binary <.h3p> file.
 *            
 *            <eslEMEM> on allocation error.
 */
int
p7_oprofile_ReadRest(P7_HMMFILE *hfp, P7_OPROFILE *om)
{
  uint32_t       magic;
  int            M;
  int            alphatype;
  int            status = eslOK;

#ifdef HMMER_THREADS
  /* lock the mutex to prevent other threads from reading from the optimized
   * profile at the same time.
   */
  if (hfp->syncRead)
    {
      if (pthread_mutex_lock (&hfp->readMutex) != 0) ESL_EXCEPTION(eslESYS, "mutex lock failed");
    }
#endif

  if (! fread( (char *) &magic,          sizeof(uint32_t), 1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read magic");
  if (magic == v3a_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/a); please hmmpress your HMM file again");
  if (magic == v3b_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/b); please hmmpress your HMM file again");
  if (magic == v3c_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/c); please hmmpress your HMM file again");
  if (magic == v3d_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/d); please hmmpress your HMM file again");
  if (magic == v3e_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/e); please hmmpress your HMM file again");
  if (magic != v3f_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic; not an HMM database file?");

  if (! fread( (char *) &M,              sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model size M");
  if (! fread( (char *) &alphatype,      sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read alphabet type");  
  if (! fread( (char *) &magic,          sizeof(uint32_t), 1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no sentinel magic: .h3p file corrupted?");

  if (magic     != v3f_pmagic)     ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad sentinel magic: .h3p file corrupted?");
  if (M         != om->M)          ESL_XFAIL(eslEFORMAT, hfp->errbuf, "p/f model length mismatch");
  if (alphatype != om->abc->type)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "p/f alphabet type mismatch");

#ifdef HMMER_THREADS
  if (hfp->syncRead)
    {
      if (pthread_mutex_unlock (&hfp->readMutex) != 0) ESL_EXCEPTION(eslESYS, "mutex unlock failed");
    }
#endif

  return eslOK;

 ERROR:
#ifdef HMMER_THREADS
  if (hfp->syncRead)
    {
      if (pthread_mutex_unlock (&hfp->readMutex) != 0) ESL_EXCEPTION(eslESYS, "mutex unlock failed");
    }
#endif

  return status;
}
/*----------- end, reading optimized profiles -------------------*/


/*****************************************************************
 * 3. Utility routines
 *****************************************************************/
/* Function:  p7_oprofile_CreateBlock()
 * Synopsis:  Create a new block of empty <P7_OM_BLOCK>.
 *
 * Purpose:   Creates a block of empty <P7_OM_BLOCK> profile objects.
 *            
 * Returns:   a pointer to the new <P7_OM_BLOCK>. Caller frees this
 *            with <p7_oprofile_DestroyBlock()>.
 *
 * Throws:    <NULL> if allocation fails.
 */
P7_OM_BLOCK *
p7_oprofile_CreateBlock(int count)
{
  int i = 0;

  P7_OM_BLOCK *block = NULL;
  int status = eslOK;

  ESL_ALLOC(block, sizeof(*block));

  block->count = 0;
  block->listSize = 0;
  block->list  = NULL;

  ESL_ALLOC(block->list, sizeof(P7_OPROFILE *) * count);
  block->listSize = count;

  for (i = 0; i < count; ++i)
    {
      block->list[i] = NULL;
    }

  return block;

 ERROR:
  if (block != NULL)
    {
      if (block->list != NULL)  free(block->list);
      free(block);
    }
  
  return NULL;
}

/* Function:  p7_oprofile_DestroyBlock()
 * Synopsis:  Frees an <P7_OM_BLOCK>.
 *
 * Purpose:   Free a Create()'d block of profiles.
 */
void
p7_oprofile_DestroyBlock(P7_OM_BLOCK *block)
{
  int i;

  if (block == NULL) return;

  if (block->list != NULL)
    {
      for (i = 0; i < block->listSize; ++i)
	{
	  if (block->list[i] != NULL) p7_oprofile_Destroy(block->list[i]);
	}
      free(block->list);
    }

  free(block);
  return;
}

/* Function:  p7_oprofile_Position()
 * Synopsis:  Reposition an open hmm file to an offset.
 *
 * Purpose:   Reposition an open <hfp> to offset <offset>.
 *            <offset> would usually be the first byte of a
 *            desired hmm record.
 *            
 * Returns:   <eslOK>     on success;
 *            <eslEOF>    if no data can be read from this position.
 *
 * Throws:    <eslEINVAL>  if the <sqfp> is not positionable.
 *            <eslEFORMAT> if no msv profile opened.
 *            <eslESYS>    if the fseeko() call fails.
 */
int
p7_oprofile_Position(P7_HMMFILE *hfp, off_t offset)
{
  if (hfp->ffp == NULL)  ESL_EXCEPTION(eslEFORMAT, hfp->errbuf, "no MSV profile file; hmmpress probably wasn't run");
  if (hfp->do_stdin)     ESL_EXCEPTION(eslEINVAL, "can't Position() in standard input");
  if (hfp->do_gzip)      ESL_EXCEPTION(eslEINVAL, "can't Position() in a gzipped file");
  if (offset < 0)        ESL_EXCEPTION(eslEINVAL, "bad offset");

  if (fseeko(hfp->ffp, offset, SEEK_SET) != 0) ESL_EXCEPTION(eslESYS, "fseeko() failed");

  return eslOK;
}

/*-------------------- end, utility routines ---------------------*/


/*****************************************************************
 * 4. Benchmark driver.
 *****************************************************************/
#ifdef p7IO_BENCHMARK
/*
  gcc  -g -Wall    -o benchmark-io -I.. -L.. -I../../easel -L../../easel -Dp7IO_BENCHMARK io.c -lhmmer -leasel -lm 
  icc  -O3 -static -o benchmark-io -I.. -L.. -I../../easel -L../../easel -Dp7IO_BENCHMARK io.c -lhmmer -leasel -lm 

  ./benchmark-io Pfam.msv
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_dummy.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <HMM MSV profile file>";
static char banner[] = "benchmark driver for profile input";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_STOPWATCH *w       = esl_stopwatch_Create();
  ESL_ALPHABET  *abc     = NULL;
  char          *msvfile = esl_opt_GetArg(go, 1);
  FILE          *msvfp   = NULL;
  P7_OPROFILE   *om      = NULL;
  int            nmodel  = 0;
  uint64_t       totM    = 0;
  int            status;

  esl_stopwatch_Start(w);

  if ((msvfp = fopen(msvfile, "r")) == NULL) p7_Fail("Failed to open MSV file %s for reading.\n", msvfile);

  while ((status = p7_oprofile_ReadMSV(msvfp, &abc, NULL, &om)) == eslOK)
    {
      nmodel++;
      totM += om->M;

      p7_oprofile_Destroy(om);
    }
  if      (status == eslEFORMAT)   p7_Fail("bad file format in profile file %s",           msvfile);
  else if (status == eslEINCOMPAT) p7_Fail("profile file %s contains different alphabets", msvfile);
  else if (status != eslEOF)       p7_Fail("Unexpected error in reading profiles from %s", msvfile);

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# number of models: %d\n", nmodel);
  printf("# total M:          %" PRId64 "\n", totM);
  
  fclose(msvfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IO_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/


/*****************************************************************
 * 5. Unit tests.
 *****************************************************************/
#ifdef p7IO_TESTDRIVE

static void
utest_ReadWrite(P7_HMM *hmm, P7_OPROFILE *om)
{
  char        *msg         = "oprofile read/write unit test failure";
  ESL_ALPHABET *abc        = NULL;
  P7_OPROFILE *om2         = NULL;
  char         tmpfile[16] = "esltmpXXXXXX";
  char        *mfile       = NULL;
  char        *ffile       = NULL;
  char        *pfile       = NULL;
  char        *ssifile     = NULL;
  FILE        *fp          = NULL;
  FILE        *mfp         = NULL;
  FILE        *ffp         = NULL;
  FILE        *pfp         = NULL;
  ESL_NEWSSI  *nssi        = NULL;
  P7_HMMFILE  *hfp         = NULL;
  uint16_t     fh          = 0;
  float        tolerance   = 0.001;
  char         errbuf[eslERRBUFSIZE];


  /* 1. A mini version of hmmpress: save the test HMM to a file along with its associated .h3{mfpi} files
   */
  if ( esl_tmpfile_named(tmpfile, &fp)          != eslOK) esl_fatal(msg);
  if ( esl_sprintf(&mfile,   "%s.h3m", tmpfile) != eslOK) esl_fatal(msg);
  if ( esl_sprintf(&ffile,   "%s.h3f", tmpfile) != eslOK) esl_fatal(msg);
  if ( esl_sprintf(&pfile,   "%s.h3p", tmpfile) != eslOK) esl_fatal(msg);
  if ( esl_sprintf(&ssifile, "%s.h3i", tmpfile) != eslOK) esl_fatal(msg);

  if ( esl_newssi_Open(ssifile, TRUE, &nssi)    != eslOK) esl_fatal(msg);
  if (( mfp = fopen(mfile, "wb"))               == NULL)  esl_fatal(msg);
  if (( ffp = fopen(ffile, "wb"))               == NULL)  esl_fatal(msg);
  if (( pfp = fopen(pfile, "wb"))               == NULL)  esl_fatal(msg);

  /* the disk offsets are all 0 by construction, if there's only one
   * HMM in the file - but don't want to forget them, if we change the
   * unit test in the future to be multi HMM
   */
  if ((om->offs[p7_MOFFSET] = ftello(mfp))      == -1)    esl_fatal(msg);
  if ((om->offs[p7_FOFFSET] = ftello(ffp))      == -1)    esl_fatal(msg);
  if ((om->offs[p7_POFFSET] = ftello(pfp))      == -1)    esl_fatal(msg);

  if ( p7_hmmfile_WriteASCII(fp,   -1, hmm)     != eslOK) esl_fatal(msg);
  if ( p7_hmmfile_WriteBinary(mfp, -1, hmm)     != eslOK) esl_fatal(msg);
  if ( p7_oprofile_Write(ffp, pfp, om)          != eslOK) esl_fatal(msg);

  if ( esl_newssi_AddFile(nssi, tmpfile, 0, &fh)                           != eslOK) esl_fatal(msg);
  if ( esl_newssi_AddKey (nssi, hmm->name, fh, om->offs[p7_MOFFSET], 0, 0) != eslOK) esl_fatal(msg);
  if ( esl_newssi_Write(nssi)                                              != eslOK) esl_fatal(msg);

  fclose(fp);
  fclose(mfp);
  fclose(ffp); 
  fclose(pfp);
  esl_newssi_Close(nssi);

  /* 2. read the optimized profile back in */
  if ( p7_hmmfile_OpenE(tmpfile, NULL, &hfp, NULL)  != eslOK) esl_fatal(msg);
  if ( p7_oprofile_ReadMSV(hfp, &abc, &om2)         != eslOK) esl_fatal(msg);
  if ( p7_oprofile_ReadRest(hfp, om2)               != eslOK) esl_fatal(msg);

  /* 3. it should be identical to the original  */
  if ( p7_oprofile_Compare(om, om2, tolerance, errbuf) != eslOK) esl_fatal("%s\n%s", msg, errbuf);
       
  p7_oprofile_Destroy(om2);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  remove(ssifile);
  remove(ffile);
  remove(pfile);
  remove(mfile);
  remove(tmpfile);

  free(ssifile);
  free(mfile);
  free(ffile);
  free(pfile);
}
  
#endif /*p7IO_TESTDRIVE*/
/*------------------ end, unit tests ----------------------------*/


/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef p7IO_TESTDRIVE
/* 
   gcc -g -Wall -std=gnu99 -o io_utest -I.. -L.. -I../../easel -L../../easel -Dp7IO_TESTDRIVE io.c -lhmmer -leasel -lm
   ./io_utest
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_dummy.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",        eslARG_INT,     "45", NULL, NULL,  NULL,  NULL, NULL, "size of random model to sample",                 0 },
  { "-L",        eslARG_INT,     "45", NULL, NULL,  NULL,  NULL, NULL, "configure model for length <n>",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for non-optimized Forward, Backward implementations";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  P7_HMM         *hmm  = NULL;
  P7_OPROFILE    *om   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  
  /* Sample a random HMM and optimized profile, in amino acid alphabet.  */
  if ((abc = esl_alphabet_Create(eslAMINO))                    == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))                                 == NULL)  esl_fatal("failed to create null model");
  if (( p7_oprofile_Sample(r, abc, bg, M, L, &hmm, NULL, &om)) != eslOK) esl_fatal("failed to sample HMM and profile");

  /* unit test(s) */
  utest_ReadWrite(hmm, om);

  p7_oprofile_Destroy(om);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}

#endif /*p7IO_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/



/*****************************************************************
 * 7. Example.
 *****************************************************************/
#ifdef p7IO_EXAMPLE
/* gcc -g -Wall -Dp7IO_EXAMPLE -I.. -I../../easel -L.. -L../../easel -o io_example io.c -lhmmer -leasel -lm
 * ./io_example <hmmfile>
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_dummy.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",      0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "verbose: print model info as they're read", 0 }, 
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <HMM file>";
static char banner[] = "example of writing MSV profile part";

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_ALPHABET  *abc     = NULL;
  char          *hmmfile = esl_opt_GetArg(go, 1);
  P7_HMMFILE    *hfp     = NULL;
  P7_HMM        *hmm     = NULL;
  P7_BG         *bg      = NULL;
  P7_PROFILE    *gm      = NULL;
  P7_OPROFILE   *om      = NULL;
  char          *fname   = NULL;
  char          *pname   = NULL;
  FILE          *ffp     = NULL;
  FILE          *pfp     = NULL;
  int            nmodel  = 0;
  uint64_t       totM    = 0;
  int            status;
  char           errbuf[eslERRBUFSIZE];

  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);  

  esl_sprintf(&fname, "%s.h3f", hmmfile);  
  esl_sprintf(&pname, "%s.h3f", hmmfile);  
  if ((ffp = fopen(fname, "wb")) == NULL) p7_Fail("failed to open %s\n", fname);
  if ((pfp = fopen(pname, "wb")) == NULL) p7_Fail("failed to open %s\n", pname);
  free(fname);
  free(pname);

  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK)
    {
      if (nmodel == 0) { 	/* first time initialization, now that alphabet known */
	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, 400);
      }

      if (esl_opt_GetBoolean(go, "-v")) printf("%s\n", hmm->name);
      nmodel++;
      totM += hmm->M;

      gm = p7_profile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, 400, p7_LOCAL);
      om = p7_oprofile_Create(gm->M, abc);
      p7_oprofile_Convert(gm, om);
      
      p7_oprofile_Write(ffp, pfp, om);

      p7_profile_Destroy(gm);
      p7_oprofile_Destroy(om);
      p7_hmm_Destroy(hmm);
    }
  if      (status == eslEFORMAT)   p7_Fail("bad file format in HMM file %s",             hmmfile);
  else if (status == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets",   hmmfile);
  else if (status != eslEOF)       p7_Fail("Unexpected error in reading HMMs from %s",   hmmfile);

  fclose(ffp);
  fclose(pfp);
  p7_bg_Destroy(bg);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IO_EXAMPLE*/

/*****************************************************************
 * @LICENSE@
 *    
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/

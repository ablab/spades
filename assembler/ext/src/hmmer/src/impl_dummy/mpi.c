/* Optional support for MPI parallelization.
 * 
 * Contents:
 *    1. Communicating P7_OPROFILE, a score profile.
 *    2. Benchmark driver.
 *    3. Unit tests.
 *    4. Test driver.
 *    5. Copyright and license information.
 *    
 * SRE, Thu Jun 14 09:59:20 2007 [Janelia] [Tom Waits, Orphans]
 * SVN $Id$
 */
#include "p7_config.h"		

#ifdef HAVE_MPI
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "easel.h"
#include "esl_mpi.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_dummy.h"

/*****************************************************************
 * 1. Communicating P7_OPROFILE, an optimized model.
 *****************************************************************/

/* Function:  p7_oprofile_MPISend()
 * Synopsis:  Send an OPROFILE as an MPI work unit.
 * Incept:    MSF, Wed Oct 21, 2009 [Janelia]
 *
 * Purpose:   Sends an OPROFILE <om> as a work unit to MPI process
 *            <dest> (where <dest> ranges from 0..<nproc-1>), tagged
 *            with MPI tag <tag>, for MPI communicator <comm>, as 
 *            the sole workunit or result. 
 *            
 *            Work units are prefixed by a status code. If <hmm> is
 *            <non-NULL>, the work unit is an <eslOK> code followed by
 *            the packed HMM. If <hmm> is NULL, the work unit is an
 *            <eslEOD> code, which <p7_hmm_MPIRecv()> knows how to
 *            interpret; this is typically used for an end-of-data
 *            signal to cleanly shut down worker processes.
 *            
 *            In order to minimize alloc/free cycles in this routine,
 *            caller passes a pointer to a working buffer <*buf> of
 *            size <*nalloc> characters. If necessary (i.e. if <hmm> is
 *            too big to fit), <*buf> will be reallocated and <*nalloc>
 *            increased to the new size. As a special case, if <*buf>
 *            is <NULL> and <*nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *            
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 * 
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 * 
 * Note:      Compare to p7_hmmfile_WriteBinary(). The two operations (sending
 *            an HMM via MPI, or saving it as a binary file to disk) are
 *            similar.
 */
int
p7_oprofile_MPISend(P7_OPROFILE *om, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  return p7_profile_MPISend(om, dest, tag, comm, buf, nalloc);
}

/* Function:  p7_oprofile_MPIRecv()
 * Synopsis:  Receives an OPROFILE as a work unit from an MPI sender.
 * Incept:    MSF, Wed Oct 21, 2009 [Janelia]
 *
 * Purpose:   Receive a work unit that consists of a single OPROFILE
 *            sent by MPI <source> (<0..nproc-1>, or
 *            <MPI_ANY_SOURCE>) tagged as <tag> for MPI communicator <comm>.
 *            
 *            Work units are prefixed by a status code. If the unit's
 *            code is <eslOK> and no errors are encountered, this
 *            routine will return <eslOK> and a non-<NULL> <*ret_om>.
 *            If the unit's code is <eslEOD> (a shutdown signal), 
 *            this routine returns <eslEOD> and <*ret_om> is <NULL>.
 *   
 *            Caller provides a working buffer <*buf> of size
 *            <*nalloc> characters. These are passed by reference, so
 *            that <*buf> can be reallocated and <*nalloc> increased
 *            if necessary. As a special case, if <*buf> is <NULL> and
 *            <*nalloc> is 0, the buffer will be allocated
 *            appropriately, but the caller is still responsible for
 *            free'ing it.
 *            
 *            Caller may or may not already know what alphabet the OPROFILE
 *            is expected to be in.  A reference to the current
 *            alphabet is passed in <abc>. If the alphabet is unknown,
 *            pass <*abc = NULL>, and when the OPROFILE is received, an
 *            appropriate new alphabet object is allocated and passed
 *            back to the caller via <*abc>.  If the alphabet is
 *            already known, <*ret_abc> is that alphabet, and the new
 *            OPROFILE's alphabet type is verified to agree with it. This
 *            mechanism allows an application to let the first OPROFILE
 *            determine the alphabet type for the application, while
 *            still keeping the alphabet under the application's scope
 *            of control.
 *
 * Returns:   <eslOK> on success. <*ret_om> contains the received OPROFILE;
 *            it is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.  If
 *            <*abc> was passed as <NULL>, it now points to an
 *            <ESL_ALPHABET> object that was allocated here; caller is
 *            responsible for free'ing this.
 *            
 *            Returns <eslEOD> if an end-of-data signal was received.
 *            In this case, <*buf>, <*nalloc>, and <*abc> are left unchanged,
 *            and <*ret_om> is <NULL>.
 *            
 *            Returns <eslEINCOMPAT> if the OPROFILE is in a different alphabet
 *            than <*abc> said to expect. In this case, <*abc> is unchanged,
 *            <*buf> and <*nalloc> may have been changed, and <*ret_om> is
 *            <NULL>.
 *            
 * Throws:    <eslEMEM> on allocation error, in which case <*ret_om> is 
 *            <NULL>.           
 */
int
p7_oprofile_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, P7_OPROFILE **ret_om)
{
  return p7_profile_MPIRecv(source, tag, comm, abc, NULL, buf, nalloc, ret_om);
}

/*----------------- end, P7_OPROFILE communication -------------------*/


/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/

#ifdef p7MPI_BENCHMARK
/* 
  mpicc -g -Wall -L. -I. -L ../easel -I ../easel -D p7MPI_BENCHMARK -o benchmark-mpi mpi.c -lhmmer -leasel -lm
  qsub -N benchmark-mpi -j y -R y -b y -cwd -V -pe lam-mpi-tight 2 'mpirun C ./benchmark-mpi  ~/notebook/1227-msp-statistics/Pfam.hmm > bench.out'
  qsub -N benchmark-mpi -j y -R y -b y -cwd -V -pe lam-mpi-tight 2 'mpirun C ./benchmark-mpi -b ~/notebook/1227-msp-statistics/Pfam.hmm > bench.out'
 */
#include "p7_config.h"

#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "baseline timing: don't send any HMMs",             0 },
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "arrest after start: for debugging MPI under gdb",  0 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for MPI communication";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_ALPHABET   *abc     = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg      = p7_bg_Create(abc);
  int             my_rank;
  int             nproc;
  char           *buf    = NULL;
  int             nbuf   = 0;
  int             subtotalM = 0;
  int             allM   = 0;
  int             stalling = esl_opt_GetBoolean(go, "--stall");

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  while (stalling); 

  /* Master MPI process: */
  if (my_rank == 0) 
    {
      ESL_STOPWATCH  *w       = esl_stopwatch_Create();
      P7_HMMFILE     *hfp     = NULL;
      P7_OPROFILE    *om      = NULL;
      P7_HMM         *hmm     = NULL;

      /* Read HMMs from a file. */
      if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);

      esl_stopwatch_Start(w);
      while (p7_oprofile_ReadMSV(hfp, &abc, &om)  == eslOK &&
	     p7_oprofile_ReadRest(hfp, om)       == eslOK)
	{
	  if (!esl_opt_GetBoolean(go, "-b"))
	    p7_oprofile_MPISend(om, 1, 0, MPI_COMM_WORLD, &buf, &nbuf); /* 1 = dest; 0 = tag */

	  p7_hmm_Destroy(hmm);
	  p7_oprofile_Destroy(om);
	}
      p7_oprofile_MPISend(NULL, 1, 0, MPI_COMM_WORLD, &buf, &nbuf); /* send the "no more HMMs" sign */
      MPI_Reduce(&subtotalM, &allM, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

      printf("total: %d\n", allM);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, "CPU Time: ");
      esl_stopwatch_Destroy(w);
    }
  /* Worker MPI process: */
  else 
    {
      P7_OPROFILE     *om_recd = NULL;      

      while (p7_oprofile_MPIRecv(0, 0, MPI_COMM_WORLD, &buf, &nbuf, &abc, &om_recd) == eslOK) 
	{
	  subtotalM += om_recd->M;
	  p7_oprofile_Destroy(om_recd);  
	}
      MPI_Reduce(&subtotalM, &allM, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

  free(buf);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  MPI_Finalize();
  exit(0);
}

#endif /*p7MPI_BENCHMARK*/
/*---------------------- end, benchmark -------------------------*/


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7MPI_TESTDRIVE

static void
utest_oprofileSendRecv(int my_rank, int nproc)
{
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(42);
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  P7_BG          *bg   = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_OPROFILE    *om   = NULL;
  P7_OPROFILE    *om2  = NULL;
  int             M    = 200;
  int             L    = 400;
  char           *wbuf = NULL;
  int             wn   = 0;
  int             i;
  char            errbuf[eslERRBUFSIZE];

  p7_hmm_Sample(r, M, abc, &hmm); /* master and worker's sampled profiles are identical */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  p7_oprofile_Convert(gm, om);
  p7_bg_SetLength  (bg, L);

  if (my_rank == 0)
    {
      for (i = 1; i < nproc; i++)
	{
	  ESL_DPRINTF1(("Master: receiving test profile\n"));
	  p7_oprofile_MPIRecv(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &wbuf, &wn, &abc, &om2);
	  ESL_DPRINTF1(("Master: test profile received\n"));

	  if (p7_oprofile_Compare(om, om2, 0.001, errbuf) != eslOK) 
	    p7_Die("Received profile not identical to what was sent\n%s", errbuf);

	  p7_oprofile_Destroy(om2);
	}
    }
  else 
    {
      ESL_DPRINTF1(("Worker %d: sending test profile\n", my_rank));
      p7_oprofile_MPISend(om, 0, 0, MPI_COMM_WORLD, &wbuf, &wn);
      ESL_DPRINTF1(("Worker %d: test profile sent\n", my_rank));
    }

  free(wbuf);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return;
}
#endif /*p7MPI_TESTDRIVE*/
/*---------------------- end, unit tests ------------------------*/


/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
#ifdef p7MPI_TESTDRIVE

/* mpicc -o mpi_utest -g -Wall -I../easel -L../easel -I. -L. -Dp7MPI_TESTDRIVE mpi.c -lhmmer -leasel -lm
 * In an MPI environment: (qlogin -pe lam-mpi-tight 2; setenv JOB_ID <jobid>; setenv TMPDIR /tmp/<jobid>....
 *    mpirun C ./mpi_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",              0 },
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "arrest after start: for debugging MPI under gdb",   0 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for mpi.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  int          my_rank;
  int          nproc;

  /* For debugging: stall until GDB can be attached */
  if (esl_opt_GetBoolean(go, "--stall"))  pause();

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  utest_oprofileSendRecv(my_rank, nproc);

  MPI_Finalize();
  return 0;
}

#endif /*p7MPI_TESTDRIVE*/
/*---------------------- end, test driver -----------------------*/


#else /*!HAVE_MPI*/
/* If we don't have MPI compiled in, provide some nothingness to:
 *   a. prevent Mac OS/X ranlib from bitching about .o file that "has no symbols" 
 *   b. prevent compiler from bitching about "empty compilation unit"
 *   c. automatically pass the automated tests.
 */
void p7_mpi_DoAbsolutelyNothing(void) { return; }

#if defined p7MPI_TESTDRIVE || p7MPI_BENCHMARK || p7MPI_EXAMPLE
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

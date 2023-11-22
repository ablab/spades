/* Simple master/worker data parallelization using POSIX threads.
 *
 * The master assigns independent work chunks to workers.  The master
 * blocks until all the workers have started, then tells the workers
 * to go. Meanwhile, workers block as they start, waiting for the go
 * signal. The workers do their work, and exit their threads. The
 * master blocks while it waits for them all to finish and exit.
 *
 * [xref H14/104 for detailed analysis of Farrar's thread synch
 * strategy]
 * 
 * Contents:
 *    1. The <ESL_THREADS> object: a gang of workers.
 *    2. Determining thread number to use.
 *    3. Examples.
 */
#include <esl_config.h>

#ifdef HAVE_PTHREAD

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>   // includes sysconf() and its constants
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#include "easel.h"
#include "esl_threads.h"


/*****************************************************************
 *# 1. The <ESL_THREADS> object: a gang of workers.
 *****************************************************************/ 

/* Function:  esl_threads_Create()
 * Synopsis:  Create a threads object.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Creates an <ESL_THREADS> object, for organizing
 8            a bunch of worker threads that will all run 
 *            the worker function <fnptr>. This object is a shell
 *            for now; the worker threads themselves are 
 *            created individually with <esl_threads_AddThread()>.
 *
 * Returns:   ptr to the new <ESL_THREADS> object.
 * 
 * Throws:    <NULL> on allocation or initialization failure.
 */
ESL_THREADS *
esl_threads_Create(void (*fnptr)(void *))
{
  ESL_THREADS *obj = NULL;
  int          status;

  ESL_ALLOC(obj, sizeof(ESL_THREADS));

  obj->threadCount     = 0;
  obj->threadId        = NULL;
  obj->data            = NULL;
  obj->startThread     = 0;
  obj->func            = fnptr;

  if (pthread_mutex_init(&obj->startMutex, NULL) != 0) ESL_XEXCEPTION(eslESYS, "mutex init failed");
  if (pthread_cond_init (&obj->startCond,  NULL) != 0) ESL_XEXCEPTION(eslESYS, "cond init failed");
  return obj;

 ERROR:
  return NULL;
}

/* Function:  esl_threads_Destroy()
 * Synopsis:  Destroys an <ESL_THREADS> object.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Frees an <ESL_THREADS> object.  
 *
 *            The caller routine must first free the
 *            contents of each <obj->data[]>.
 *
 * Returns:   void
 */
void 
esl_threads_Destroy(ESL_THREADS *obj)
{
  if (obj == NULL) return;

  if (obj->threadId != NULL) free(obj->threadId);
  if (obj->data     != NULL) free(obj->data);
  pthread_mutex_destroy(&obj->startMutex);
  pthread_cond_destroy (&obj->startCond);
  free(obj);
  return;
}

/* Function:  esl_threads_AddThread()
 * Synopsis:  Add a worker thread to the <ESL_THREADS> object.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Create a new worker thread for the <ESL_THREADS> object,
 *            assigning it the work unit pointed to by <data>. 
 *
 *            The caller remains responsible for any memory allocated
 *            to <data>; the <ESL_THREADS> object will only manage
 *            a copy of a pointer to <data>.
 *
 *            Only the master calls this, while it's initializing the
 *            workers; it calls _WaitForStart() after it's done adding
 *            worker threads.
 *
 *            Each worker thread starts immediately, and is supposed
 *            to call _Started() to block soon after it starts, and
 *            tell the master it's ready and waiting.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on allocation failure. 
 *            <eslESYS> if thread creation fails.
 *            <eslEINVAL> if something's wrong with the <obj>.
 */
int 
esl_threads_AddThread(ESL_THREADS *obj, void *data)
{
  int    status;
  void  *p;

  if (obj == NULL) ESL_EXCEPTION(eslEINVAL, "Invalid thread object");

  /* allocate inside the ESL_THREADS object to hold another worker */
  ESL_RALLOC(obj->threadId, p, sizeof(pthread_t) * (obj->threadCount+1));
  ESL_RALLOC(obj->data,     p, sizeof(void *)    * (obj->threadCount+1));
  
  obj->data[obj->threadCount] = data;
  if (pthread_create(&(obj->threadId[obj->threadCount]), NULL, (void *(*)(void *)) obj->func, obj) != 0) ESL_EXCEPTION(eslESYS, "thread creation failed");
  obj->threadCount++;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_threads_GetWorkerCount()
 * Synopsis:  Return the total number of worker threads.
 * Incept:    SRE, Fri Aug 21 13:22:52 2009 [Janelia]
 *
 * Purpose:   Returns the total number of worker threads.
 *
 *            This should only be called in the master after all the
 *            _AddThread() calls are done, or in a worker after the
 *            master's go signal has been received; then this is the
 *            number of parallel workers. 
 */
int
esl_threads_GetWorkerCount(ESL_THREADS *obj)
{
  return obj->threadCount;
}


/* Function:  esl_threads_WaitForStart()
 * Synopsis:  Blocks master until all workers have started.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Make the master thread wait until all the worker threads have
 *            started. When all the worker threads have started and
 *            are blocking at the start mutex, release them.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslESYS> if thread synchronization fails somewhere.
 *            <eslEINVAL> if something is awry with <obj>.
 */
int
esl_threads_WaitForStart(ESL_THREADS *obj)
{
  if (obj == NULL) ESL_EXCEPTION(eslEINVAL, "Invalid thread object");

  if (pthread_mutex_lock (&obj->startMutex) != 0) ESL_EXCEPTION(eslESYS, "mutex lock failed");

  /* wait for all worker threads to start */
  while (obj->startThread < obj->threadCount) {
    if (pthread_cond_wait(&obj->startCond, &obj->startMutex) != 0) ESL_EXCEPTION(eslESYS, "wait cond failed");
  }

  /* release all the worker threads */
  obj->startThread = 0;
  if (pthread_cond_broadcast(&obj->startCond)  != 0) ESL_EXCEPTION(eslESYS, "cond broadcast failed");
  if (pthread_mutex_unlock  (&obj->startMutex) != 0) ESL_EXCEPTION(eslESYS, "mutex unlock failed");
  return eslOK;
}

/* Function:  esl_threads_WaitForFinish()
 * Synopsis:  Blocks master until all workers have completed.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Block the master thread until all the worker threads have
 *            completed. As each worker completes, remove it from the 
 *            <obj>. 
 *            
 *            Upon exit, the <obj> is returned to the same (empty)
 *            state it was in after it was created. It may be reused
 *            for a new problem by adding new workers.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslESYS> if thread synchronization fails somewhere.
 *            <eslEINVAL> if something is awry with <obj>.
 */
int 
esl_threads_WaitForFinish(ESL_THREADS *obj)
{
  int  w;			

  if (obj == NULL) ESL_EXCEPTION(eslEINVAL, "Invalid thread object");

  /* wait for all worker threads to complete */
  for (w = obj->threadCount-1; w >= 0; w--)
    if (pthread_join(obj->threadId[w], NULL) != 0) ESL_EXCEPTION(eslESYS, "pthread join failed");
  obj->threadCount = 0;   // can't touch threadCount until workers are all finished, because they may be using it.

  return eslOK;
}

/* Function:  esl_threads_Started()
 * Synopsis:  Blocks worker until master gives the start signal.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Block a worker thread until master sees that all workers
 *            have started and gives the start signal. Assign the worker
 *            a unique number (0..nworkers-1), and return it in
 *            <*ret_workeridx>. The worker uses this index to 
 *            retrieve its work units.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslESYS> if thread synchronization fails somewhere.
 *            <eslEINVAL> if something is awry with <obj>.
 */
int
esl_threads_Started(ESL_THREADS *obj, int *ret_workeridx)
{
  int           w;
  pthread_t     threadId;
  int           status;

  if (obj == NULL)                                ESL_XEXCEPTION(eslEINVAL, "Invalid thread object");
  if (pthread_mutex_lock (&obj->startMutex) != 0) ESL_XEXCEPTION(eslESYS,   "mutex lock failed");

  /* signal that we're started */
  obj->startThread++;
  if (pthread_cond_broadcast (&obj->startCond) != 0) ESL_XEXCEPTION(eslESYS, "cond broadcast failed");

  /* the logic here is a little inelegant, but it's because we're
   * using the same start mutex/condition for both the worker->master
   * started signal and the master->worker go signal, so maybe it's
   * elegant from that perspective. What's happening is that each time
   * a worker starts, it *broadcasts* the start condition. It's
   * intended for the master; if the master sees that startThread ==
   * threadCount it knows that all workers have started, and it can
   * tell them to go. But since it's broadcast, and workers are also
   * waiting for the "go" start condition (below), workers will wake
   * up below; but they see that startThread > 0, so they go back to
   * sleep.
   *
   * Once we finish starting threads, the master sets startThread back
   * to 0, then broadcasts the start condition again. Now all the
   * workers wake up below, and see startThread=0, break their while
   * loop and proceed.
   *
   * So the cost is a lot of unnecessary worker-waking while we're
   * really just sending a 1:1 signal from a newly started worker to
   * the master. We could instead use two different "started" and "go"
   * mutexes/conditions.
   *
   * [SRE note added; xref H14/104]
   */

  /* wait for the master's go signal to start the calculations */
  while (obj->startThread) {
    if (pthread_cond_wait(&obj->startCond, &obj->startMutex) != 0) ESL_XEXCEPTION(eslESYS, "cond wait failed");
  }

  if (pthread_mutex_unlock (&obj->startMutex) != 0)  ESL_XEXCEPTION(eslESYS, "mutex unlock failed");

  /* Figure out the worker's index */
  threadId = pthread_self();
  for (w = 0; w < obj->threadCount; w++)   // threadCount is stable until thread terminates.
    if (pthread_equal(threadId, obj->threadId[w])) break;
  if (w == obj->threadCount) ESL_XEXCEPTION(eslESYS, "thread not registered");

  *ret_workeridx = w;
  return eslOK;

 ERROR:
  *ret_workeridx = 0;
  return status;
}


/* Function:  esl_threads_GetData()
 * Synopsis:  Return the data associated with this thread.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Return the data pointer associated with the worker thread
 *            <workeridx>. The data pointer was set by the 
 *            <esl_threads_AddThread()> function.
 *
 * Returns:   void *
 */
void *
esl_threads_GetData(ESL_THREADS *obj, int workeridx)
{
  return obj->data[workeridx];
}


/* Function:  esl_threads_Finished()
 * Synopsis:  Terminate the thread.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Terminate a worker thread.
 *            This is currently a no-op, serving as
 *            a placeholder in case we eventually need
 *            any cleanup.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_threads_Finished(ESL_THREADS *obj, int workeridx)
{
  return eslOK;
}


/*****************************************************************
 * 2. Determining thread number to use
 *****************************************************************/

/* Function:  esl_threads_CPUCount()
 * Synopsis:  Figure out how many cores the machine has available.
 * Incept:    SRE, Wed Aug 19 11:31:24 2009 [Janelia]
 *
 * Purpose:   Determine the number of logical processors on this
 *            machine; return that number in <*ret_ncpu>.
 *            
 *            The number of available processors is found by
 *            <sysconf(_SC_NPROCESSORS_ONLN)>, a POSIX standard.  If
 *            sysconf() is somehow unavailable, <*ret_ncpu> is quietly
 *            just set to 1.
 *
 *            This is the number of logical processors, not physical
 *            processors. On systems with hyperthreading, the number
 *            of logical processors is more than the number of
 *            physical cpus. It may or may not be a good thing to
 *            spawn more threads than physical processors.
 *            
 * Args:      ret_ncpu  - RETURN: number of logical cores
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      J5/68; H14/5.
 */
int
esl_threads_CPUCount(int *ret_ncpu)
{
  int   ncpu = 1;

#if defined     (HAVE_SYSCONF) && defined (_SC_NPROCESSORS_ONLN)     // POSIX standard 
  ncpu = sysconf(_SC_NPROCESSORS_ONLN);
#elif defined   (HAVE_SYSCONF) && defined (_SC_NPROC_ONLN)	     // Silicon Graphics IRIX used to be this way 
  ncpu = sysconf(_SC_NPROC_ONLN);                                    // Don't use sysctl() as a fallback. Linux/gcc deprecated it. [H14/5]
#endif
  if (ncpu < 1) ncpu = 1;  // silently ignore any sysconf() problem, including nonexistence; fall back to 1 core.

  *ret_ncpu = ncpu;
  return eslOK;
}


/* Function:  esl_threads_GetCPUCount()
 * Synopsis:  Returns the number of CPU cores on machine.
 * Incept:    SRE, Mon Aug 21 08:52:29 2017
 *
 * Purpose:   Identical to <esl_threads_CPUCount()>, except
 *            it directly returns the result.
 */
int
esl_threads_GetCPUCount(void)
{
  static int ncpu = -1;                         // so we only make system calls once.
  if (ncpu == -1) esl_threads_CPUCount(&ncpu);
  return ncpu;
}


/*****************************************************************
 * 3. Example
 *****************************************************************/

#ifdef eslTHREADS_EXAMPLE
#include "easel.h"
#include "esl_threads.h"

static void 		
worker_thread(void *data)
{
  ESL_THREADS *thr = (ESL_THREADS *) data;
  char  *s         = NULL;
  int    w;

  esl_threads_Started(thr, &w);
  
  s = (char *) esl_threads_GetData(thr, w);
  printf("worker thread %d receives: %s\n", w, s);

  esl_threads_Finished(thr, w);
  return;
}

int
main(void)
{
  ESL_THREADS  *thr  = NULL;
  int           ncpu = 8;
  int           i;
  char        **work = NULL;

  work = malloc(sizeof(char *) * ncpu);
  for (i = 0; i < ncpu; i++) 
    esl_sprintf(&(work[i]), "work packet %d", i);

  thr = esl_threads_Create(&worker_thread);

  for (i = 0; i < ncpu; i++)
    esl_threads_AddThread(thr, (void *) work[i]);

  esl_threads_WaitForStart (thr);
  /* The worker threads now run their work. */
  esl_threads_WaitForFinish(thr);
  esl_threads_Destroy(thr);
  for (i = 0; i < ncpu; i++) free(work[i]);
  free(work);
  return eslOK;
}
#endif /*eslTHREADS_EXAMPLE*/


#ifdef eslTHREADS_EXAMPLE2
#include "easel.h"
#include "esl_threads.h"

/* gcc --std=gnu99 -g -Wall -pthread -o esl_threads_example2 -I. -DeslTHREADS_EXAMPLE2 esl_threads.c easel.c */
int 
main(void)
{
  int ncpu;

  esl_threads_CPUCount(&ncpu);
  printf("Processors: %d\n", ncpu);

  return eslOK;
}
#endif /*eslTHREADS_EXAMPLE2*/
#endif /*HAVE_PTHREAD*/



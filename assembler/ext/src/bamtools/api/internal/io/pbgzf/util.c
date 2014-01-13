#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <sys/stat.h>
#include <unistd.h>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(_SC_NPROCESSORS_ONLN)
#include <unistd.h>
#endif

void
safe_mutex_lock(pthread_mutex_t *mutex)
{
  int ret = pthread_mutex_lock(mutex);
  if(ret != 0) {
      fprintf(stderr, "pthread_mutex_lock error [%d]! Aborting immediately!\n", ret);
      exit(1);
  }
}

void
safe_mutex_unlock(pthread_mutex_t *mutex)
{
  int ret = pthread_mutex_unlock(mutex);
  if(ret != 0) {
      fprintf(stderr, "pthread_mutex_unlock error [%d]! Aborting immediately!\n", ret);
      exit(1);
  }
}

// from pbzip2 version 1.1.6 
/*
   This program, "pbzip2" is copyright (C) 2003-2011 Jeff Gilchrist.
   All rights reserved.

   The library "libbzip2" which pbzip2 uses, is copyright
   (C) 1996-2008 Julian R Seward.  All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   2. The origin of this software must not be misrepresented; you must
   not claim that you wrote the original software.  If you use this
   software in a product, an acknowledgment in the product
   documentation would be appreciated but is not required.

   3. Altered source versions must be plainly marked as such, and must
   not be misrepresented as being the original software.

   4. The name of the author may not be used to endorse or promote
   products derived from this software without specific prior written
   permission.

   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
   OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   Jeff Gilchrist, Ottawa, Canada.
   pbzip2@compression.ca
   pbzip2 version 1.1.6 of Oct 30, 2011
   */
int32_t
detect_cpus()
{
  int32_t ncpu;

  // Set default to 1 in case there is no auto-detect
  ncpu = 1;

  // Autodetect the number of CPUs on a box, if available
#if defined(__APPLE__)
  size_t len = sizeof(ncpu);
  int32_t mib[2];
  mib[0] = CTL_HW;
  mib[1] = HW_NCPU;
  if (sysctl(mib, 2, &ncpu, &len, 0, 0) < 0 || len != sizeof(ncpu))
    ncpu = 1;
#elif defined(_SC_NPROCESSORS_ONLN)
  ncpu = sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(WIN32)
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  ncpu = si.dwNumberOfProcessors;
#else
#warning "CPU autodection is disabled"
#endif

  // Ensure we have at least one processor to use
  if (ncpu < 1)
    ncpu = 1;

  return ncpu;
}

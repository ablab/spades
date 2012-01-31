#pragma once
#include <time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

struct elapsed_timer
{
   elapsed_timer()
   {
      restart();
   }

   double elapsed() const
   {
      timespec ts;
      my_clock_gettime(ts);
      return ts.tv_sec - ts_.tv_sec + (ts.tv_nsec - ts_.tv_nsec) * 1e-9;
   }

   void restart()
   {
      my_clock_gettime(ts_);
   }

private:
   timespec ts_;

   void my_clock_gettime(timespec &ts) const {
      #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
         clock_serv_t cclock;
         mach_timespec_t mts;
         host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
         clock_get_time(cclock, &mts);
         mach_port_deallocate(mach_task_self(), cclock);
         ts.tv_sec = mts.tv_sec;
         ts.tv_nsec = mts.tv_nsec;
      #else
         clock_gettime(CLOCK_REALTIME, &ts);
      #endif      
   }
};

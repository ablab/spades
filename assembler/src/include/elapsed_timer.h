#pragma once
#include <time.h>

struct elapsed_timer
{
   elapsed_timer()
   {
      restart();
   }

   double elapsed() const
   {
      timespec ts;
      clock_gettime(CLOCK_REALTIME, &ts);

      return ts.tv_sec - ts_.tv_sec + (ts.tv_nsec - ts_.tv_nsec) * 1e-9;
   }

   void restart()
   {
      clock_gettime(CLOCK_REALTIME, &ts_);
   }

private:
   timespec ts_;
};

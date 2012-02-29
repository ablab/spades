///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// TaskTimer gathers information on the wall clock time, cpu time, and
/// memory used by a process during some task.  Simply create a
/// TaskTimer object and surround the code in question with a call to
/// that TaskTimer object's Start() and Stop() methods. 
///

/// \class TaskTimer 
///

///   TaskTimer t;
///   t.Start();
///   theTaskToBeTimed();
///   t.Stop();
///

/// When a TaskTimer object has stopped, it can be printed using the <<
/// operator.
///

///   cout << t << endl;
///

/// Repeatedly starting and stopping a timer will accumulate time and
/// resource usage information.  If you want to use the same timer for
/// a different task, you should call the Reset() method between each
/// use of that timer.
///

/// For example, if you want to time two consecutive tasks using the same
/// timer, do something like this:
///

///   TaskTimer t;
///

///   cout << "task 1:" << endl;
///   t.Start();   
///   task1();   
///   t.Stop();   
///   cout << t << endl;
///

///   cout << "task 2:" << endl;
///   t.Reset();
///   t.Start();
///   task2();
///   t.Stop();
///   cout << t << endl;
///

/// If you want to time two tasks in a loop, do something like this:
///

///   TaskTimer t1, t2;
///   for ( int i = 0; i < a_big_number; ++i )
///   {
///     t1.Start();
///     task1(i);
///     t1.Stop();
///

///     t2.Start();
///     task2(i);
///     t2.Stop();
///   }
///

///   cout << "Total time for task1(): " << endl;
///   cout << t1 << endl;
///   cout << "Total time for task2(): " << endl;
///   cout << t2 << endl;
///  
///   By default, a TaskTimer will not track the increase in memory
///   usage during a test.  Use SetTrackMemory( true ) to turn on this
///   feature.


#ifndef INCLUDED_TASKTIMER
#define INCLUDED_TASKTIMER

#include "CoreTools.h"
#include "system/SysConf.h"

#include <sys/times.h>
#include <sys/resource.h>

#include <iostream>

class TaskTimer {
private:
  bool is_running_;
  bool has_run_;
  bool track_memory_;

  struct rusage my_rusage_;

  long start_mem_usage_;
  long diff_mem_usage_;

  struct tms start_tms_;
  struct tms stop_tms_;

  int tasks_timed_;

  clock_t running_time_;
  clock_t u_time_;
  clock_t s_time_;

  clock_t timeout_;

public:
  TaskTimer()
    : is_running_(false), 
      track_memory_(false),
      diff_mem_usage_(0),
      tasks_timed_(0),
      running_time_(0), 
      u_time_(0), 
      s_time_(0),
      timeout_(0)
  { }

  bool IsTrackingMemory() const { return track_memory_; }
  void SetTrackMemory( const bool track_memory ) { track_memory_ = track_memory; }

  bool IsRunning() const { return is_running_; }

  int GetTasksTimed() const { return tasks_timed_; }
  
  // I avoid use of local variables in Start() and Stop() to
  // increase the probability that they will be inlined.

  void Start()
  {
    if ( is_running_ ) {
      cerr << "This timer is already running." << endl;
      return;
    }
    is_running_ = true;
    ++tasks_timed_;

    if ( track_memory_ ) {
      getrusage( RUSAGE_SELF, &my_rusage_ );
      start_mem_usage_ = my_rusage_.ru_maxrss;
    }

    running_time_ -= times( &start_tms_ );
  }

  void Stop()
  {
    running_time_ += times( &stop_tms_ );

    u_time_ += stop_tms_.tms_utime - start_tms_.tms_utime;
    s_time_ += stop_tms_.tms_stime - start_tms_.tms_stime;

    if ( ! is_running_ ) {
      cerr << "This timer is not running." << endl;
      return;
    }

    is_running_ = false;

    if ( track_memory_ ) {
      getrusage( RUSAGE_SELF, &my_rusage_ );
      diff_mem_usage_ += my_rusage_.ru_maxrss - start_mem_usage_;
    }
  }

  void Reset()
  { 
    if ( is_running_ ) {
      cerr << "This timer is running and cannot be reset." << endl;
      return;
    }

    running_time_ = 0;
    tasks_timed_ = 0;
    u_time_ = 0;
    s_time_ = 0;
    timeout_ = 0;
  }

  // ---- added 2009-06-30 ribeiro
  void SetTimeOut(const float timeout_sec)
  {
    const float ticks_per_second = clockTicksPerSecond();
    
    timeout_ = clock_t(timeout_sec * ticks_per_second);
  }

  // ---- added 2009-06-30 ribeiro
  void StartWithTimeOut(const float timeout_sec)
  {
    Start();
    SetTimeOut(timeout_sec);
  }

  // ---- added 2009-06-30 ribeiro
  bool TimedOut()
  {
    if (timeout_ <= 0) return false;

    if (is_running_) {
      struct tms cur_tms;
      if (running_time_ + times(&cur_tms) > timeout_) {
        Stop();
        return true;
      }
      else 
        return false;
    }
    else 
      return (running_time_ > timeout_);
        
  }


  // Print only elapsed time.

  void PrintElapsed( ostream &out )
  {
    if ( is_running_ ) {
      out << "This timer is still running.";
      return;
    }
      
    if ( tasks_timed_ == 0 ) {
      out << "This timer has not been used.";
      return;
    }

    float ticks_per_second = clockTicksPerSecond();
    out << ToString( running_time_ / ticks_per_second, 2 ) << "s";
  }

  // Return elapsed time in seconds as float. Returns -1 on error.

  float Elapsed( )
  {
    if ( is_running_ || tasks_timed_ == 0 )
      return -1;

    float ticks_per_second = clockTicksPerSecond();
    return ( running_time_ / ticks_per_second );
  }

  // Return user CPU time in seconds as float. Returns -1 on error.

  float UserSecs( )
  {
    if ( is_running_ || tasks_timed_ == 0 )
      return -1;

    float ticks_per_second = clockTicksPerSecond();
    return ( u_time_ / ticks_per_second );
  }

  // Return system CPU time in seconds as float. Returns -1 on error.

  float SysSecs( )
  {
    if ( is_running_ || tasks_timed_ == 0 )
      return -1;

    float ticks_per_second = clockTicksPerSecond();
    return ( s_time_ / ticks_per_second );
  }

  friend ostream& operator<<( ostream& out, const TaskTimer& the_timer )
  {
    if ( the_timer.is_running_ ) {
      out << "This timer is still running.";
      return out;
    }

    if ( the_timer.tasks_timed_ == 0 ) {
      out << "This timer has not been used.";
      return out;
    }

    out << "Tasks: " << setw(5) << the_timer.tasks_timed_;
    out << "  ";

    float ticks_per_second = clockTicksPerSecond();

    int old_precision = out.precision( 2 );
    ios::fmtflags old_fmtflags = out.setf( ios::fixed | ios::left );
      
    float elapsed_time = the_timer.running_time_;
    out << "Elapsed: "
        << elapsed_time / ticks_per_second;
    out << setw(4) << "s ";

    float user_time = the_timer.u_time_;
    float sys_time = the_timer.s_time_;

    user_time /= ticks_per_second;
    sys_time /= ticks_per_second;

    out << "CPU (user/sys): "
        << user_time + sys_time << "s (" << user_time << "/" << sys_time << ")";
    out << setw(4) << "s ";

    if ( the_timer.track_memory_ ) {
      out << "Mem: +"
          << the_timer.diff_mem_usage_ / 1024 << "M";
    }

    out.precision( old_precision );
    out.flags( old_fmtflags );

    return out;
  }

  void Write(String filename, bool update = true) {
    float newelapsed = Elapsed();
    float oldelapsed = 0.0;

    ifstream infile(filename.c_str());
    if (infile) 
      infile >> oldelapsed; 
    infile.close();

    ofstream outfile(filename.c_str());
    if (outfile)  
      outfile << (oldelapsed + newelapsed) << endl;
    outfile.close();
  }

};

#endif

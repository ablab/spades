///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Worklist.h
 * \author tsharpe
 * \date Jan 16, 2009
 *
 * \brief A class that uses multiple threads of execution to process a worklist.
 */
#ifndef WORKLIST_H_
#define WORKLIST_H_

#include "system/Assert.h"
#include "system/SysConf.h"
#include "system/WorklistUtils.h"
#include <algorithm>
#include <cstddef>
#include <list>
#include <time.h>

/// A Worklist manages multiple threads to process the Workitems handed to it.
/// It is a templated class with two template parameters, a Workitem type, and
/// a Processor type.
/// <p>
/// See the application WorklistTest for an example of how to use it.
/// <p>
/// Your Workitem must be a class that defines a no-arg constructor, a copy
/// constructor, and allows copying by assignment.  I.e., it must define
/// Workitem(), Workitem(Workitem const&), and
/// Workitem& operator=(Workitem const&).<br>
/// (Primitive types work, too.)<br>
/// It will describe to your Processor what work needs to be done.  (It might
/// be a seed number for a local assembly, for example.)  Since it will be
/// copied onto the worklist when you add it, and copied from the worklist when
/// it comes time to process it, and since there might be oodles of Workitems on
/// the worklist at any one time, you probably want to keep the amount of stuff
/// in the Workitem minimal.  Sharable, thread-safe state, for example, belongs
/// in the Processor, not in the Workitem.
/// <p>
/// Your Processor will be called as follows:<br>
/// p( Workitem& workitem );<br>
/// I.e., it's a functor that processes a Workitem.
/// Therefore, it must either be a function pointer like this:
/// void(*processorFunc)( Workitem& );
/// or a class with this member function defined: void operator()( Workitem& );
/// <p>
/// A separate instance of your Processor object is created (by copying the one
/// you supply to the Worklist constructor) for each thread of execution.
/// Therefore it must have a copy constructor. Your copy constructor should make
/// a separate, private copy of all state that is not thread-safe. And, of
/// course, all global state accessible to your processor needs to be managed in
/// a thread-safe manner.
/// <p>
/// Your processor will be handed multiple Workitems to process, one at a time,
/// and its destructor will be executed when all Workitems have been processed.
/// (Note:  We do not create a Processor for each Workitem.  We create a
/// Processor for each thread.)
/// <p>
/// The Worklist itself is built on an STL list, so you may want to manage the
/// number of pending items in some way if you're concerned about memory.  The
/// add functions return the current size of the list to make this easier.  You
/// can also add a batch of Workitems, and then wait() until the list is
/// drained, and then add another batch.
/// <p>
/// The Worklist destructor waits until all Workitems have been processed, and
/// all threads terminated. The main thread of execution (the one that is
/// executing the Worklist's destructor) will pitch in at this point to help
/// process any remaining Workitems.  Why not, eh?  It's not doing anything
/// else.  And besides, that lets you test your code in a single threaded
/// environment -- just call the Worklist constructor with nThreads = 0.
/// <p><pre>
/// Example:
/// #include "system/Worklist.h"
///
/// class MyProcessor
/// {
/// public:
///   void operator()( char* arg )
///   {
///     // do something with arg
///   };
/// };
///
/// int
/// main( int argc, char** argv )
/// {
///    MyProcessor p;
///    Worklist<char*,MyProcessor> worklist( p, 2 );
///    worklist.add( argv, argv+argc );
///    //Note: we're just letting the worklist evaporate from the stack.  We're
///    // not calling wait() or anything.  When the worklist falls out of scope,
///    // it's destructor is called.  The destructor will wait until all the
///    //work is done.
/// }
/// </pre>
///
/// A special note about exceptions:  Your Processor's operator() function ought
/// to catch every exception from which it can recover.  (As is always true of
/// every function.)  Sadly, you cannot, however, do a catch(...) and just try
/// to ignore all exceptions. The problem is that the pthread library is trying
/// to be helpful by making pthread_exit throw an exception so that the stack
/// gets unwound.  Unfortunately, the exception it throws is anonymous and
/// cannot be ignored. (The exception's destructor terminates the program.)  So
/// you're hosed.  You can't ignore it, because it'll kill the program, but you
/// can't not ignore it if you catch(...) because you don't know its name, and
/// can't rethrow it. So uncaught exceptions will kill the entire process,
/// rather than just the thread that caused the exception, and the author has
/// been unable to figure out a way around this problem.
template<class Workitem, class Processor>
class Worklist
{
public:
    struct ProgressMonitor
    {
        ProgressMonitor() : mTime(0) {}
        Workitem mWorkitem;
        time_t mTime;
    };

    /// Create a Worklist that will be processed by nThreads of execution.
    /// If nThreads is 0, no threads will be created, and the entire
    /// Worklist will be processed by the Worklist destructor.  The default
    /// value for nThreads is the number of processors less one (we assume that
    /// the main thread will continue to do work).
    /// Passing a stack size of 0 (it's default value) will cause threads to be
    /// created with a stack as large as main thread's stack (unless the main
    /// thread's stack size is unlimited, in which case a default value of 8Mb
    /// will be used for each thread's stack size).
    Worklist( Processor const& p,
              size_t nThreads = processorsOnline()-1,
              size_t threadStackSize = 0 )
    : mProcessor(p), mProgressMonitors(new ProgressMonitor[nThreads+1]),
      mTP(nThreads,threadStackSize,threadFunc,this)
    {}

    /// This will wait until all Workitems have been processed, and all threads
    /// have been terminated.  That might be a long time!
    ~Worklist()
    { waitForDone(); delete [] mProgressMonitors; }

    /// Start up a separate thread to monitor thread death.
    /// This thread wakes up every howOftenToCheckInSecs seconds, and replaces
    /// any threads that have died with a new thread.
    /// If you call this more than once, subsequent calls merely adjust the
    /// check interval.
    /// If you call this with howOftenToCheckInSecs == 0, the monitor exits.
    /// If reportRestarts is true, an error message is written to cerr to inform
    /// you of a thread restart.
    void threadDeathMonitor( long howOftenToCheckInSecs,
                             bool reportRestarts = true )
    { mTP.threadDeathMonitor(howOftenToCheckInSecs,reportRestarts); }

    /// Add a single Workitem to the Worklist.
    /// The number of unprocessed Workitems is returned.
    size_t add( Workitem const& workitem )
    { QueueStateManipulator qsm(mQS);
      Assert(!qsm.isDone());
      mList.push_back(workitem);
      return qsm.incSize(); }

    /// Add a mess of Workitems to the Worklist.
    /// The type of (*first) must be a Workitem or a reference to one.
    /// The number of unprocessed Workitems is returned.
    template<class InputIterator>
    size_t add( InputIterator first, InputIterator last )
    { QueueStateManipulator qsm(mQS);
      Assert(!qsm.isDone());
      while ( first != last )
      { mList.push_back(*first); ++first; }
      size_t sz = mList.size(); qsm.setSize(sz); return sz; }

    /// Return the number of unprocessed Workitems.
    size_t size()
    { QueueStateManipulator qsm(mQS);
      return qsm.getSize(); }

    /// Wait until the Worklist is empty.  An empty queue does NOT mean that all
    /// the work is complete.  It's just a good time to add more work if there
    /// is any.  It might be convenient if you can't afford the memory needed to
    /// add all the Workitems at once.
    void waitForEmpty()
    { QueueStateManipulator(mQS).waitForEmpty(); }

    /// Note:  you don't have to call this -- the destructor will automatically
    /// wait until all the work is complete.  You may call this for
    /// entertainment value, if you'd like.
    void waitForDone()
    { mTP.threadDeathMonitor(0); // shut down the monitor.
      QueueStateManipulator(mQS).setDone(); // tell threads to quit
      doWork(true);
      mTP.shutdown(); }

    /// Clear the worklist of all unprocessed items.
    /// If you need to quit early, this is the only way to do it.
    /// The destructor is still going to wait for in-progress work to complete.
    void clear()
    { QueueStateManipulator qsm(mQS); mList.clear(); qsm.setSize(0); }

    /// How many threads are there to process the work.
    size_t getNThreads()
    { return mTP.getNThreads(); }

    /// What happening with a thread.
    /// ThreadIdx runs from 0 to getNThreads()-1.
    ProgressMonitor getProgress( size_t threadIdx )
    { QueueStateManipulator qsm(mQS);
      return mProgressMonitors[threadIdx]; }

private:
    // unimplemented -- no copying
    Worklist( Worklist const& );
    Worklist& operator=( Worklist const& );

    void doWork( bool isMainThread = false );
    static void* threadFunc( void* );


    std::list<Workitem> mList;
    QueueState mQS;
    Processor mProcessor;
    ProgressMonitor* mProgressMonitors;
    ThreadPool mTP; // this needs to be the last member because we don't want
                    // to start threads until everything else is ready
};

template<class Workitem, class Processor>
void Worklist<Workitem,Processor>::doWork( bool isMainThread )
{
    size_t idx = mTP.getNThreads();
    if ( !isMainThread )
        idx = mTP.findThreadIndex(pthread_self());
    ProgressMonitor& progressMonitor = mProgressMonitors[idx];
    Processor processor(mProcessor);

    while ( true )
    {
        if ( true ) // empty block just to control lifetime of qsm
        {
            QueueStateManipulator qsm(mQS);
            qsm.waitForWork();
            if ( !qsm.getSize() ) // wait for work only returns when there's
                                  // work, or when the done flag is set
                break;            // so break if there's no work
                                  // (because that means we're all done)

            time(&progressMonitor.mTime);
            using std::swap;
            swap(progressMonitor.mWorkitem,mList.front());
            mList.pop_front();
            qsm.decSize();
        }
        processor(progressMonitor.mWorkitem);
    }
}

template<class Workitem, class Processor>
void* Worklist<Workitem,Processor>::threadFunc( void* ptr )
{
    reinterpret_cast<Worklist*>(ptr)->doWork();
    return 0;
}

#endif /* WORKLIST_H_ */

///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file WorklistMP.h
 * \author tsharpe
 * \date May 10, 2011
 *
 * \brief
 */
#ifndef WORKLISTMP_H_
#define WORKLISTMP_H_

#include "VecString.h"
#include "feudal/BinaryStream.h"
#include "system/WorklistUtils.h"
#include <cerrno>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <list>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

class SubProcess
{
public:
    SubProcess( char*const* args );

    // compiler-supplied copying and destructor are OK

    pid_t getPID() const { return mPID; }
    int getFD() const { return mFD; }

    bool isIdle() const { return mIsIdle; }
    void setIdle( bool val ) { mIsIdle = val; }

    static void fail( char const* msg )  __attribute__((__noreturn__));

private:
    pid_t mPID;
    int mFD;
    bool mIsIdle;
};

// This is a multi-process worklist.
//
// The parent process creates a bunch of child processes.  Parent and child must
// agree on the nature of the work, i.e., the template args to this class:
//
// ItemMsg is what you send to the child process, explaining what work there is
// to do.  It must be default-constructable, copyable, and
// BinaryStream-serializable.
//
// ResultMsg is what the client sends back, explaining what was done.  It, too,
// must be default-constructable, copyable, and BinaryStream-serializable.
//
// Children will block until the parent dispatches work.  The parent monitors
// each child for death and for results.  After a child returns results, the
// parent may dispatch more work, and repeat the cycle, until no more work
// remains.
//
// This could be made fully distributed, if we ever needed that.
//
template <class ItemMsg, class ResultMsg>
class WorklistMP
{
public:
    class Client
    {
    public:
        Client() : mRdr(3,"<worklist socket>"), mWrtr(3,"<worklist socket>") {}

        // Call these three methods in this order, and repeat until it's time
        // to quit.

        bool isQuittingTime() { return mRdr.atEOF(); }

        ItemMsg getWork() { ItemMsg result; mRdr.read(&result); return result; }

        void reportResults( ResultMsg const& msg )
        { mWrtr.write(msg); mWrtr.flush(); }

    private:
        Client( Client const& ); // unimplemented -- no copying
        Client& operator=( Client const& ); // unimplemented -- no copying

        BinaryReader mRdr;
        BinaryWriter mWrtr;
    };

    WorklistMP( size_t nProcs, vecString const& cmd );
    ~WorklistMP();

    // this method waits until all work is complete, which is, perhaps, not what
    // we want.
    // returns true if, as far as we can tell, everything went OK.
    // the Itr must deref to a ItemMsg const&
    // the ResultsProc is a functor that takes a ResultMsg const&
    template <class Itr, class ResultsProc>
    bool doWork( Itr itr, Itr const& end, ResultsProc proc );

private:
    WorklistMP( WorklistMP const& ); // unimplemented -- no copying
    WorklistMP& operator=( WorklistMP const& ); // unimplemented -- no copying

    void checkStatus( SubProcess const& subProc )
    { int status;
      while ( waitpid(subProc.getPID(),&status,0) == -1 )
        if ( errno != EINTR ) SubProcess::fail("Call to waitpid failed");
      if ( !WIFEXITED(status) || WEXITSTATUS(status) )
      { if ( WIFEXITED(status) )
          std::cout << "Sub-process #" << subProc.getPID()
                      << " exited with exit status: "
                      << WEXITSTATUS(status) << std::endl;
        else if ( WIFSIGNALED(status) )
          std::cout << "Sub-process #" << subProc.getPID()
                      << " exited due to caught signal: "
                      << WTERMSIG(status) << std::endl;
        else
          std::cout << "Sub-process #" << subProc.getPID()
                      << " exited for unknown reasons."
                      << std::endl; } }

    std::list<SubProcess> mSubProcs;
};

template <class ItemMsg, class ResultMsg>
WorklistMP<ItemMsg,ResultMsg>::WorklistMP( size_t nProcs, vecString const& cmd )
{
    size_t nArgs = cmd.size();
    char** args = new char*[nArgs+1];
    for ( size_t idx = 0; idx < nArgs; ++idx )
        args[idx] = const_cast<char*>(cmd[idx].c_str());
    args[nArgs] = 0;
    while ( nProcs-- )
        mSubProcs.push_back(SubProcess(args));
    delete [] args;
}

template <class ItemMsg, class ResultMsg>
WorklistMP<ItemMsg,ResultMsg>::~WorklistMP()
{
    typedef std::list<SubProcess>::iterator Itr;
    for ( Itr itr(mSubProcs.begin()), end(mSubProcs.end()); itr != end; ++itr )
        if ( close(itr->getFD()) )
            SubProcess::fail("Unable to close main process socket");
    for ( Itr itr(mSubProcs.begin()), end(mSubProcs.end()); itr != end; ++itr )
        checkStatus(*itr);
}

template <class ItemMsg, class ResultMsg>
template <class Itr, class ResultsProc>
bool WorklistMP<ItemMsg,ResultMsg>::doWork( Itr itr, Itr const& end,
                                            ResultsProc proc )
{
    bool result = true;
    Dotter dotter(std::distance(itr,end));
    std::cout << "Farming out " << dotter.getNBatches()
              << " jobs to sub-processes." << std::endl;
    typedef std::list<SubProcess>::iterator PItr;
    PItr pend(mSubProcs.end());
    while ( true )
    {
        fd_set fds;
        FD_ZERO(&fds);
        int nfds = 0;
        for ( PItr pitr(mSubProcs.begin()); pitr != pend; ++pitr )
        {
            if ( pitr->isIdle() && itr != end )
            {
                BinaryWriter bw(pitr->getFD(),"<worklist socket>");
                bw.write(*itr);
                ++itr;
                pitr->setIdle(false);
            }
            if ( !pitr->isIdle() )
            {
                int fd = pitr->getFD();
                FD_SET(fd,&fds);
                if ( fd >= nfds )
                    nfds = fd + 1;
            }
        }

        if ( !nfds )
            break; // everyone is idle, and there's no more work

        if ( select(nfds,&fds,0,0,0) == -1 )
            SubProcess::fail("Select failed");

        for ( PItr pitr(mSubProcs.begin()); pitr != pend; ++pitr )
        {
            if ( FD_ISSET(pitr->getFD(),&fds) )
            {
                BinaryReader br(pitr->getFD(),"<worklist socket>");
                if ( br.atEOF() )
                {
                    std::cout << "Sub-process #" << pitr->getPID()
                            << " did not return results.  It probably croaked."
                            << std::endl;
                    result = false;
                    if ( close(pitr->getFD()) )
                        SubProcess::fail("Unable to close main process socket");
                    checkStatus(*pitr);
                    pitr = mSubProcs.erase(pitr);
                }
                else
                {
                    ResultMsg msg;
                    br.read(&msg);
                    proc(msg);
                    pitr->setIdle(true);
                }
                dotter.batchDone();
            }
        }
    }
    return result;
}

#endif /* WORKLISTMP_H_ */

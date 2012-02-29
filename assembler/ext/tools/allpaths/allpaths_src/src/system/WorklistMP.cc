///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file WorklistMP.cc
 * \author tsharpe
 * \date May 13, 2011
 *
 * \brief
 */
#include "system/WorklistMP.h"
#include "system/ErrNo.h"
#include "system/System.h"
#include <cstring>
#include <iostream>
#include <sys/types.h>
#include <sys/socket.h>

SubProcess::SubProcess( char*const* args )
: mIsIdle(true)
{
    int fds[2];
    if ( socketpair(AF_UNIX,SOCK_STREAM,0,fds) )
        fail("Unable to create socketpair for sub-process");
    mFD = fds[0];
    pid_t pid = fork();
    if ( pid )
    {
        mPID = pid;
        if ( close(fds[1]) )
            fail("Unable to close sub-process socket in main");
    }
    else
    {
        if ( close(fds[0]) )
            fail("Unable to close main-process socket in sub-process");
        if ( dup2(fds[1],3) != 3 )
            fail("Unable to dup socket descriptor in sub-process");
        int keepAlive = 1;
        if ( setsockopt(3,SOL_SOCKET,SO_KEEPALIVE,
                            &keepAlive,sizeof(keepAlive)) == -1 )
        {
            ErrNo err;
            std::cout <<
                    "Warning: Unable to set keep alive for sub-process socket"
                    << err << std::endl;
        }
        if ( close(fds[1]) )
            fail("Unable to close original sub-process socket in sub-process");
        execvp(args[0],args);
        fail("Unable to exec sub-process: ");
    }
}

void SubProcess::fail( char const* msg )
{
    ErrNo err;
    FatalErr(msg << err);
}

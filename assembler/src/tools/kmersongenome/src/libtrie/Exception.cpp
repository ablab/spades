/*							-*- C++ -*-
** Exception.cpp
** Login : Julien Lemoine <speedblue@happycoders.org>
** Started on  Thu Jun  1 20:32:08 2006 Julien Lemoine
** $Id$
** 
** Copyright (C) 2006 Julien Lemoine
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
** 
** You should have received a copy of the GNU Lesser General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include <iostream>
#include <assert.h>
#include "Exception.h"

namespace ToolBox
{
  StackTrace::StackTrace(const std::string &error, const std::string &file, unsigned line) :
    _error(error), _file(file), _line(line)
  {
  }

  StackTrace::~StackTrace()
  {
  }

  const std::string& StackTrace::getError() const
  {
    return _error;
  }

  const std::string& StackTrace::getFile() const
  {
    return _file;
  }

  const int StackTrace::getLine() const
  {
    return _line;
  }

  Exception::Exception(const Exception& e, 
		       const std::string &file, int line) throw()
  {
    _trace = e._trace;
    _trace.push_back(StackTrace("", file, line));
  }

  Exception::Exception(const std::string& err, 
		       const std::string &file, int line) throw()
  {
    _trace.push_back(StackTrace(err, file, line));
  }

  Exception::~Exception() throw()
  {
  }

  const char* Exception::what() const throw()
  {
    assert(!_trace.empty());
    return _trace.front().getError().c_str();
  }

  const std::string& Exception::getError() const
  {
    assert(!_trace.empty());
    return _trace.front().getError();
  }

  const std::string& Exception::getFile() const
  {
    assert(!_trace.empty());
    return _trace.front().getFile();
  }

  const int Exception::getLine() const
  {
    assert(!_trace.empty());
    return _trace.front().getLine();
  }

  void Exception::print(std::ostream& stream) const
  {
    bool first = true;
    stream << "Exception : " << std::endl;
    for (std::list<StackTrace>::const_iterator it = _trace.begin();
	 it != _trace.end(); ++it)
      {
	if (it->getError().size())
	  stream << "Error[" << it->getError() << "] ";
	stream << (first ? "throwed at [" : "catched at [") 
	       << it->getFile() << ":" << it->getLine() << "]" << std::endl;
	first = false;
      }
  }

  std::ostream& operator<<(std::ostream& stream, const Exception& e)
  {
    e.print(stream);
    return (stream);
  }

  // all subclass implementation of superclass Exception
  ExceptionCore(EmptyException, Exception)
  ExceptionCore(ExecException, Exception)
  ExceptionCore(MemoryError, Exception)
  ExceptionCore(CancelAlgorithm, Exception)
  ExceptionCore(StrtolError, Exception)
  ExceptionCore(ThreadError, Exception)
  ExceptionCore(ParseError, Exception)
  ExceptionCore(Sentinel, Exception)
  ExceptionCore(NullPointer, Exception)
  ExceptionCore(InputException, Exception)
  ExceptionCore(FileException, InputException)
  ExceptionCore(EOFException, FileException)
  ExceptionCore(MD5Exception, InputException)
  ExceptionCore(XmlException, InputException)
  ExceptionCore(DbException, InputException)
  ExceptionCore(LibraryException, Exception)
  ExceptionCore(ColumnException, Exception)
  ExceptionCore(InvalidIndex, Exception)
  ExceptionCore(InvalidMatrix, Exception)
  ExceptionCore(NetException, Exception)
  ExceptionCore(HostnameError, NetException)
  ExceptionCore(Ipv6SupportError, NetException)
  ExceptionCore(InetptonError, NetException)
  ExceptionCore(InetntopError, NetException)
  ExceptionCore(ConnectionClosed, NetException)
  ExceptionCore(NoConnection, NetException)
  ExceptionCore(Timeout, NetException)
  ExceptionCore(BindError, NetException)
  ExceptionCore(SocketError, NetException)
  ExceptionCore(ListenError, NetException)
  ExceptionCore(SetsockoptError, NetException)
  ExceptionCore(CloseError, NetException)
  ExceptionCore(SelectError, NetException)
  ExceptionCore(ConnectError, NetException)
  ExceptionCore(AcceptError, NetException)
  ExceptionCore(GetpeernameError, NetException)
  ExceptionCore(UnicodeError, Exception)
  ExceptionCore(TokenizerError, Exception)
  ExceptionCore(CompilerError, Exception)
  ExceptionCore(PoolError, Exception)
}

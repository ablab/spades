/*							-*- C++ -*-
** Exception.h
** Login : Julien Lemoine <speedblue@happycoders.org>
** Started on  Thu Jun  1 20:32:07 2006 Julien Lemoine
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

#ifndef   	EXCEPTION_HH_
# define   	EXCEPTION_HH_

#include <iostream>
#include <exception>
#include <string>
#include <list>

/**
 * \file
 * \brief Exception API: provide common exception class
 *
 * <h2>Common exceptions functionalities, tracing for debug...</h2>
 *
 * <p>This API is used to create exception object in the same way
 * everywhere in projects. It provides tracing support and macro to
 * include file position in trace</p>
 */

namespace ToolBox
{
  /**
   * macro used to add file position (filename and line to the
   * exception), it uses preprocessor functions.
   */
#define HERE	__FILE__, __LINE__


  /**
   * macro used to define header of an Exception class
   * @param Name the name of class
   * @param From the name of superclass
   */
#define ExceptionHeader(Name, From)		\
  class  Name : public From			\
  {						\
  public:					\
    Name(const std::string &what,		\
         const std::string &file, int line)	\
          throw();				\
    Name(const Exception &e,			\
	 const std::string &file, int line)	\
	  throw();				\
    virtual ~Name() throw();			\
						\
  private:					\
    Name& operator=(const Name &e);		\
    Name();					\
  }

  /**
   * macro used to define the implementation of an Exception
   * class. The implementation call the superclass constructor.
   * @param Name the name of class
   * @param From the name of superclass
   */
#define ExceptionCore(Name, From)		\
  Name::Name(const std::string &vwhat,		\
         const std::string &file, int line) 	\
    throw() : From(vwhat, file, line)  {}	\
  Name::Name(const Exception &e,		\
         const std::string &file, int line) 	\
    throw() : From(e, file, line)  {}		\
  Name::~Name() throw() {}

  /// store stack of exception
  class StackTrace
  {
  public:
    /// constructor
    StackTrace(const std::string &error, const std::string &file, unsigned line);
    /// destructor
    ~StackTrace();

  private:
    /// avoid default constructor
    StackTrace();

  public:
    /**
     * @brief get the error message of exception
     */
    const std::string& getError() const;
    /**
     * @brief get the filename where exception has been launched
     */
    const std::string& getFile() const;
    /**
     * @brief get the line of the file where exception has been
     * launched
     */
    const int getLine() const;

  protected:
    /// @brief the error message of exception
    std::string	_error;
    /// @brief the filename where exception has been launched
    std::string	_file;
    /// @brief the line of file where exception has been launched
    int		_line;
  };

  /**
   * @brief Generic definition of an exception. This is the
   * superclass of all kind of exception, it contains all
   * implementation.
   * @author Julien Lemoine <speedblue@happycoders.org>
   */
  class Exception : public std::exception
  {
  public:
    /**
     * @brief constructor of exception and add a exception on the stack
     * @param what the description of the error
     * @param file and line are filled using the {@link HERE} macro.
     */
    Exception(const std::string &what,
	      const std::string &file, int line) throw();
    /**
     * @brief constructor of exception from a previous one. This
     * constructor add a trace of current filename and position in the
     * output of exception. This method is usefull for debug purpose.
     * @param e represent the previous exception
     * @param file and line are filled using the {@link HERE} macro.
     */
    Exception(const Exception& e, const std::string &file, 
	      int line) throw();
    /// destructor
    virtual ~Exception() throw();

  public:
    /**
     * @brief compatibility method to be compatible with
     * std::exception. It gives the same result as {@link getError}
     * with char* instead of std::string.
     */
    const char* what() const throw();
    /**
     * @brief get the error message of exception
     */
    const std::string& getError() const;
    /**
     * @brief get the filename where exception has been launched
     */
    const std::string& getFile() const;
    /**
     * @brief get the line of the file where exception has been
     * launched
     */
    const int getLine() const;
    /**
     * @brief display the complet trace of exception.
     */
    void print(std::ostream& stream) const;

  private:
    /// @brief avoid default constructor
    Exception();
    /// @brief avoid affectation operator
    Exception& operator=(const Exception&);

  protected:
    /// @brief stack trace
    std::list<StackTrace>	_trace;
  };

  /**
   * @brief display the complet trace of exception in a
   * std::ostream. This method calls {@link ToolBox::Exception::display}.
   */
  std::ostream& operator<<(std::ostream& stream, const Exception& e);

  // all subclass headers of superclass Exception
  ExceptionHeader(EmptyException, Exception);
  ExceptionHeader(ExecException, Exception);
  ExceptionHeader(MemoryError, Exception);
  ExceptionHeader(CancelAlgorithm, Exception);
  ExceptionHeader(StrtolError, Exception);
  ExceptionHeader(ThreadError, Exception);
  ExceptionHeader(ParseError, Exception);
  ExceptionHeader(Sentinel, Exception);
  ExceptionHeader(NullPointer, Exception);
  ExceptionHeader(InputException, Exception);
  ExceptionHeader(FileException, InputException);
  ExceptionHeader(EOFException, FileException);
  ExceptionHeader(MD5Exception, InputException);
  ExceptionHeader(XmlException, InputException);
  ExceptionHeader(DbException, InputException);
  ExceptionHeader(LibraryException, Exception);
  ExceptionHeader(ColumnException, Exception);
  ExceptionHeader(InvalidIndex, Exception);
  ExceptionHeader(InvalidMatrix, Exception);
  ExceptionHeader(NetException, Exception);
  ExceptionHeader(HostnameError, NetException);
  ExceptionHeader(Ipv6SupportError, NetException);
  ExceptionHeader(InetptonError, NetException);
  ExceptionHeader(InetntopError, NetException);
  ExceptionHeader(ConnectionClosed, NetException);
  ExceptionHeader(NoConnection, NetException);
  ExceptionHeader(Timeout, NetException);
  ExceptionHeader(BindError, NetException);
  ExceptionHeader(SocketError, NetException);
  ExceptionHeader(ListenError, NetException);
  ExceptionHeader(SetsockoptError, NetException);
  ExceptionHeader(CloseError, NetException);
  ExceptionHeader(SelectError, NetException);
  ExceptionHeader(ConnectError, NetException);
  ExceptionHeader(AcceptError, NetException);
  ExceptionHeader(GetpeernameError, NetException);
  ExceptionHeader(UnicodeError, Exception);
  ExceptionHeader(TokenizerError, Exception);
  ExceptionHeader(CompilerError, Exception);
  ExceptionHeader(PoolError, Exception);
}

#endif	    /* !EXCEPTION_HH_ */

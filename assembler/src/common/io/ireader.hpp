/**
 * @file    ireader.hpp
 * @author  Mariya Fomkina
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * IReader is the interface for all other readers and reader wrappers.
 */

#ifndef COMMON_IO_IREADER_HPP_
#define COMMON_IO_IREADER_HPP_

template<typename ReadType>
class IReader {
 public:
  /* 
   * Default destructor.
   */
  virtual ~IReader() {}

  /* 
   * Check whether the stream is opened.
   */
  virtual bool is_open() = 0;

  /* 
   * Check whether we've reached the end of stream.
   */
  virtual bool eof() = 0;

  /*
   * Read single or paired read from stream.
   *
   * @param singleread The single or paired read that will store read
   * data. 
   *
   * @return Reference to this stream.
   */
  virtual IReader& operator>>(ReadType& read) = 0;

  /*
   * Close the stream.
   */
  virtual void close() = 0;

  /* 
   * Close the stream and open it again.
   */
  virtual void reset() = 0;
};

#endif /* COMMON_IO_IREADER_HPP_ */

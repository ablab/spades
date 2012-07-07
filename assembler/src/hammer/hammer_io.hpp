//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_IO_HPP
#define HAMMER_IO_HPP

#include <boost/iostreams/filtering_stream.hpp>
#include <string>

/// structure for boost istreams
struct FIStream {
  std::string fn;
  boost::iostreams::filtering_istream fs;
  std::ifstream stdstream;
  std::vector<char> buffer;
  bool remove_it;
  FIStream(const std::string & fname);
  FIStream(const std::string & fname, bool input_output);
  FIStream(const std::string & fname, bool input_output, uint64_t bufsize);
  ~FIStream();

  static boost::shared_ptr<FIStream> init(const std::string & fname, bool input_output = false);
  static boost::shared_ptr<FIStream> init_buf(const std::string & fname, uint64_t bufsize);
};

/// structure for boost ostreams
struct FOStream {
  std::string fn;
  boost::iostreams::filtering_ostream fs;
  std::ofstream stdstream;
  std::vector<char> buffer;
  bool remove_it;
  FOStream(const std::string & fname);
  FOStream(const std::string & fname, bool input_output);
  FOStream(const std::string & fname, bool input_output, uint64_t bufsize);
  ~FOStream();

  static boost::shared_ptr<FOStream> init(const std::string & fname, bool input_output = false);
  static boost::shared_ptr<FOStream> init_buf(const std::string & fname, uint64_t bufsize);
};


#endif // HAMMER_IO_HPP

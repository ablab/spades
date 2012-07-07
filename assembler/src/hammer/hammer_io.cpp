//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include <string>
#include <vector>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "config_struct_hammer.hpp"
#include "hammer_stats.hpp"
#include "hammer_io.hpp"

FIStream::FIStream(const std::string & fname) :
  fn(fname), stdstream(fname,
                       cfg::get().general_gzip ?
                       (std::ios::in | std::ios::binary) : std::ios::in),
  buffer(0), remove_it(false) {
  if (cfg::get().general_gzip)
    fs.push(boost::iostreams::gzip_decompressor());
  fs.push(stdstream);
}

FIStream::FIStream(const std::string & fname, bool input_output) :
  fn(fname), stdstream(fname,
                       ((!input_output && cfg::get().general_gzip) ||
                        (input_output && cfg::get().input_gzipped)) ?
                       (std::ios::in | std::ios::binary) : std::ios::in),
  buffer(0), remove_it(false) {
  if (cfg::get().general_gzip)
    fs.push(boost::iostreams::gzip_decompressor());
  fs.push(stdstream);
}

FIStream::FIStream(const std::string & fname, bool input_output, uint64_t bufsize) :
  fn(fname), buffer(bufsize), remove_it(false) {
  stdstream.rdbuf()->pubsetbuf(&buffer[0], bufsize);
  stdstream.open(fname, ((!input_output && cfg::get().general_gzip) ||
                         (input_output && cfg::get().input_gzipped)) ?
                 (std::ios::in | std::ios::binary) : std::ios::in);
  if (cfg::get().general_gzip)
    fs.push(boost::iostreams::gzip_decompressor());
  fs.push(stdstream);
}

FIStream::~FIStream() {
  fs.reset();
  if (remove_it && cfg::get().general_remove_temp_files) {
    if (remove( fn.c_str() ) != 0) {
      TIMEDLN("Error deleting file " + fn);
    }
  }
}

boost::shared_ptr<FIStream> FIStream::init(const std::string & fname, bool input_output) {
  boost::shared_ptr<FIStream> p(new FIStream(fname, input_output));
  return p;
}

boost::shared_ptr<FIStream> FIStream::init_buf(const std::string & fname, uint64_t bufsize) {
  boost::shared_ptr<FIStream> p(new FIStream(fname, false, bufsize));
  return p;
}

FOStream::FOStream(const std::string & fname) :
  fn(fname), stdstream(fname,
                       cfg::get().general_gzip ?
                       (std::ios::out | std::ios::binary) : std::ios::out),
  buffer(0), remove_it(false) {
  if (cfg::get().general_gzip)
    fs.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
  fs.push(stdstream);
}

FOStream::FOStream(const std::string & fname, bool input_output) :
  fn(fname), stdstream(fname,
                       ((!input_output && cfg::get().general_gzip) ||
                        (input_output && cfg::get().input_gzipped)) ?
                       (std::ios::out | std::ios::binary) : std::ios::out),
  buffer(0), remove_it(false) {
  if (cfg::get().general_gzip)
    fs.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
  fs.push(stdstream);
}


FOStream::FOStream(const std::string & fname, bool input_output, uint64_t bufsize) : fn(fname), buffer(bufsize), remove_it(false) {
  stdstream.rdbuf()->pubsetbuf(&buffer[0], bufsize);
  stdstream.open(fname, ((!input_output && cfg::get().general_gzip) ||
                         (input_output && cfg::get().input_gzipped)) ?
                 (std::ios::out | std::ios::binary) : std::ios::out);
  if (cfg::get().general_gzip)
    fs.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
  fs.push(stdstream);
}

FOStream::~FOStream() {
  fs.reset();
  if (remove_it && cfg::get().general_remove_temp_files) {
    if (remove( fn.c_str() ) != 0) {
      TIMEDLN("Error deleting file " + fn);
    }
  }
}

boost::shared_ptr<FOStream> FOStream::init(const std::string & fname, bool input_output) {
  boost::shared_ptr<FOStream> p(new FOStream(fname, input_output));
  return p;
}

boost::shared_ptr<FOStream> FOStream::init_buf(const std::string & fname, uint64_t bufsize) {
  boost::shared_ptr<FOStream> p(new FOStream(fname, false, bufsize));
  return p;
}

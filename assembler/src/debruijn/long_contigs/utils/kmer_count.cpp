//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

///*
// * kmer_count.cpp
// *
// *  Created on: Sep 16, 2011
// *      Author: andrey
// */
//
//
//#include "io/single_read.hpp"
//#include "sequence/sequence.hpp"
//#include "sequence/seq.hpp"
//#include "cuckoo.hpp"
//#include <tr1/unordered_map>
//#include <iostream>
//#include <fstream>
//
//const short max = std::numeric_limits<short>::max();
//const size_t K_ = 19;
//
//typedef Seq<K_> Kmer;
//typedef std::tr1::unordered_map<Kmer, short, typename Kmer::hash, typename Kmer::equal_to> MapType;
//
//typedef std::tr1::unordered_map<short, short> HistType;
//
//
//
//template<class ReadStream>
//void ProcessReads(ReadStream &stream, MapType& map) {
//	map.clear();
//	io::SingleRead r;
//	while (!stream.eof()) {
//		stream >> r;
//		Sequence s = r.sequence();
//
//		if (s.size() < K_)
//			continue;
//
//		Seq<K_> kmer = s.start<K_>();
//
//		if (map[kmer] < max) {
//			map[kmer] += 1;
//		}
//		for (size_t j = K_; j < s.size(); ++j) {
//			kmer = kmer << s[j];
//			if (map[kmer] < max) {
//				map[kmer] += 1;
//			}
//		}
//	}
//}
//
//void MakeHist(MapType& map, HistType& hist) {
//	hist.clear();
//	for (auto iter = map.begin(); iter != map.end(); ++ iter) {
//		hist[iter->second] += 1;
//	}
//}
//
//
//int main(int argc, char ** argv) {
//	if (argc < 3) {
//		std::cout << "Specify input file(s) and output" << std::endl;
//		return -1;
//	}
//
//	MapType map;
//
//	for(int i = 1; i < argc - 1; ++i) {
//		//TODO: add filtered stream to exclude Ns
//		typedef io::Reader<io::SingleRead> ReadStream;
//		ReadStream stream(argv[i]);
//
//		ProcessReads(stream, map);
//	}
//
//	HistType hist;
//	MakeHist(map, hist);
//
//	std::ofstream f(argv[argc - 1]);
//	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
//		f << iter->first << " " << iter->second << std::endl;
//	}
//	return 0;
//}

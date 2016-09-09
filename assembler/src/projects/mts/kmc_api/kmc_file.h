/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 2.2.0
  Date   : 2015-04-15
*/

#ifndef _KMC_FILE_H
#define _KMC_FILE_H

#include "kmer_defs.h"
#include "kmer_api.h"
#include <string>
#include <vector>

class CKMCFile
{
	enum open_mode {closed, opened_for_RA, opened_for_listing};
	open_mode is_opened;

	bool end_of_file;

	FILE *file_pre;
	FILE *file_suf;

	uint64* prefix_file_buf;
	uint64 prefix_file_buf_size;
	uint64 prefix_index;			// The current prefix's index in an array "prefix_file_buf", readed from *.kmc_pre
	uint32 single_LUT_size;			// The size of a single LUT (in no. of elements)

	uint32* signature_map;
	uint32 signature_map_size;
	
	uchar* sufix_file_buf;
	uint64 sufix_number;			// The sufix's number to be listed
	uint64 index_in_partial_buf;	// The current byte's number in an array "sufix_file_buf", for listing mode

	uint32 kmer_length;
	uint32 mode;
	uint32 counter_size;
	uint32 lut_prefix_length;
	uint32 signature_len;
	uint32 min_count;
	uint32 max_count;
	uint64 total_kmers;

	uint32 kmc_version;
	uint32 sufix_size;		// sufix's size in bytes 
	uint32 sufix_rec_size;  // sufix_size + counter_size

	uint32 original_min_count;
	uint32 original_max_count;

	static uint64 part_size; // the size of a block readed to sufix_file_buf, in listing mode 
	
	bool BinarySearch(int64 index_start, int64 index_stop, const CKmerAPI& kmer, uint32& counter, uint32 pattern_offset);

	// Open a file, recognize its size and check its marker. Auxiliary function.
	bool OpenASingleFile(const std::string &file_name, FILE *&file_handler, uint64 &size, char marker[]);	

	// Recognize current parameters. Auxiliary function.
	bool ReadParamsFrom_prefix_file_buf(uint64 &size);	

	// Reload a contents of an array "sufix_file_buf" for listing mode. Auxiliary function. 
	void Reload_sufix_file_buf();

	// Implementation of GetCountersForRead for kmc1 database format
	bool GetCountersForRead_kmc1(const std::string& read, std::vector<uint32>& counters);

	// Implementation of GetCountersForRead for kmc2 database format
	bool GetCountersForRead_kmc2(const std::string& read, std::vector<uint32>& counters);
public:
		
	CKMCFile();
	~CKMCFile();

	// Open files *.kmc_pre & *.kmc_suf, read them to RAM, close files. *.kmc_suf is opened for random access
	bool OpenForRA(const std::string &file_name);

	// Open files *kmc_pre & *.kmc_suf, read *.kmc_pre to RAM, *.kmc_suf is buffered
	bool OpenForListing(const std::string& file_name);

	// Return next kmer in CKmerAPI &kmer. Return its counter in float &count. Return true if not EOF
	bool ReadNextKmer(CKmerAPI &kmer, float &count);

	bool ReadNextKmer(CKmerAPI &kmer, uint32 &count);
	// Release memory and close files in case they were opened 
	bool Close();

	// Set the minimal value for a counter. Kmers with counters below this theshold are ignored
	bool SetMinCount(uint32 x);

	// Return a value of min_count. Kmers with counters below this theshold are ignored 
	uint32 GetMinCount(void);

	// Set the maximal value for a counter. Kmers with counters above this theshold are ignored
	bool SetMaxCount(uint32 x);

	// Return a value of max_count. Kmers with counters above this theshold are ignored 
	uint32 GetMaxCount(void);

	// Return the total number of kmers between min_count and max_count
	uint64 KmerCount(void);

	// Return the length of kmers
	uint32 KmerLength(void);

	// Set initial values to enable listing kmers from the begining. Only in listing mode
	bool RestartListing(void);

	// Return true if all kmers are listed
	bool Eof(void);

	// Return true if kmer exists. In this case return kmer's counter in count
	bool CheckKmer(CKmerAPI &kmer, float &count);

	bool CheckKmer(CKmerAPI &kmer, uint32 &count);

	// Return true if kmer exists
	bool IsKmer(CKmerAPI &kmer);

	// Set original (readed from *.kmer_pre) values for min_count and max_count
	void ResetMinMaxCounts(void);

	// Get current parameters from kmer_database
	bool Info(uint32 &_kmer_length, uint32 &_mode, uint32 &_counter_size, uint32 &_lut_prefix_length, uint32 &_signature_len, uint32 &_min_count, uint32 &_max_count, uint64 &_total_kmers);
	
	// Get counters for all k-mers in read
	bool GetCountersForRead(const std::string& read, std::vector<uint32>& counters);
	bool GetCountersForRead(const std::string& read, std::vector<float>& counters);
	private:
		uint32 count_for_kmer_kmc1(CKmerAPI& kmer);
		uint32 count_for_kmer_kmc2(CKmerAPI& kmer, uint32 bin_start_pos);
};

#endif

// ***** EOF

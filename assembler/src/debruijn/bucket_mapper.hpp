#ifndef BUCKET_MAPPER_HPP
#define BUCKET_MAPPER_HPP

#include <iomanip> 
#include "indices/debruijn_kmer_index.hpp"


namespace debruijn_graph {

	template <class Graph>
	class BucketMapper {

		typedef int bucket_id;

		const DeBruijnEdgeIndex<EdgeId>& kmer_index_;
		const Graph& g_;
		const unsigned K_;
		const unsigned bucketNum_;

		std::vector<std::pair<int, int> > coverage_to_multiplicity;
		std::vector<int> buckets;

		int UpdateCoverageCounters( const Sequence& seq ){

			int kmer_counter = 0;
			runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K_);
			std::map<int, int> used_coverages;
			for (size_t j = K_; j < seq.size(); ++j) {

				kmer <<= seq[j];
				int kmer_coverage = kmer_index_[kmer].count_;
				
			//	std::cout << "kmer: " << kmer.str() << " " << kmer_coverage << std::endl;
				if (used_coverages.find(kmer_coverage) != used_coverages.end()){
				
					used_coverages[kmer_coverage] += 1;
				}
				else {
					used_coverages[kmer_coverage] = 1;
				}
				kmer_counter += 1;

			}

			coverage_to_multiplicity.insert(coverage_to_multiplicity.begin(), used_coverages.begin(), used_coverages.end() );
			std::sort(coverage_to_multiplicity.begin(), coverage_to_multiplicity.end());

			return kmer_counter;
		}


		int GetNumberKmersInBucket( int id ) {
			
			int number_of_kmers = 0;
			int lower_coverage = 0, upper_coverage = buckets[id];
			if ( id > 0 ){
				lower_coverage = buckets[id - 1];
			}
			for ( auto it = coverage_to_multiplicity.begin(); it != coverage_to_multiplicity.end(); ++it ){

				if (it->first > lower_coverage && it->first <= upper_coverage) number_of_kmers += it->second;

			}

			std::cout << "number of kmers in bucket "  << number_of_kmers << std::endl;
			return number_of_kmers;

		}


		void CountBuckets( int bucket_size_bound ) {

			int current_bucket_size = 0;
	
			for ( auto it = coverage_to_multiplicity.begin(); it != coverage_to_multiplicity.end(); ++it ) {
				
				current_bucket_size += it->second;
				if ( current_bucket_size > bucket_size_bound ) {
					buckets.push_back(it->first);	
					current_bucket_size = 0;

				}

			}

			if ( current_bucket_size < bucket_size_bound ) {

				buckets.push_back( (coverage_to_multiplicity.end() - 1)->first);
			}

		}

		public:

		BucketMapper( const Graph& g, const DeBruijnEdgeIndex<EdgeId>& kmer_index, unsigned K, unsigned bucketNum ) : g_(g), kmer_index_(kmer_index), K_(K), bucketNum_(bucketNum) {
		
		}

		void InitBuckets( ) {

	
			int kmer_counter = 0;
			for (auto e = g_.SmartEdgeBegin(); !e.IsEnd(); ++e) {

				if (g_.length(*e) >= cfg::get().rr.max_repeat_length) {

					kmer_counter += UpdateCoverageCounters( g_.EdgeNucls(*e) );
				}

			}
		
			std::cout << "kmer_counter: " << kmer_counter << std::endl;
			CountBuckets( kmer_counter / bucketNum_ );
			
			std::cout << "kmer_counter / bucketNum_: " << kmer_counter / bucketNum_ << std::endl;

			std::cout << "Buckets: " << std::endl;
			for (auto it = buckets.begin(); it != buckets.end(); ++it){
				std::cout << *it << " ";
			}
			std::cout << std::endl;
		/*	int distance = 500;
			std::vector< std::vector<double> > bucketToBucket;
			for ( unsigned id = 0; id < bucketNum_; ++id ) {
			
				std::vector<double> histogram(bucketNum_,0);
				SetBucketsForDistance (  id, distance, histogram );
				bucketToBucket.push_back(histogram);
			}

			std::cout << "Histogram:" << std::endl;
			std::cout.precision(4);
			
			for (auto hist = bucketToBucket.begin(); hist != bucketToBucket.end(); ++hist) {
				for (auto it = hist->begin(); it != hist->end(); ++it) {

					std::cout << *it << " ";
				}
				std::cout << std::endl;
			}

			std::cout << std::endl;
		*/
		}

		int GetKmerBucket( const runtime_k::RtSeq& kmer ) {

			 int kmer_coverage = kmer_index_[kmer].count_;
		//	 std::cout << "kmer coverage " << kmer_coverage << std::endl;
			 bucket_id id = 0;
			 while ( id < buckets.size() - 1 && buckets[id] < kmer_coverage ) {
				++id;
			 }

		//	 std::cout << "bucket: " << id << std::endl;
			return id;
		}


		void SetBucketsForDistance ( bucket_id id, int distance, std::vector<double>& histogram ) {


			int kmers_in_bucket_counter = GetNumberKmersInBucket(id);
			for (auto e = g_.SmartEdgeBegin(); !e.IsEnd(); ++e) {

				if (g_.length(*e) >= cfg::get().rr.max_repeat_length) {
 
					Sequence seq =  g_.EdgeNucls(*e) ;
					runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K_);
				
					//std::cout << "Sequence: " << seq.str() << std::endl ;
					//std::cout << "Sequence size: " << seq.size() << std::endl ;
					//std::cout << "j min: " << K_ << " j max: "<< seq.size()  - K_ - distance  << std::endl;
					//std::cout << "Kmer: " << kmer.str() << std::endl ;
					for (size_t j = K_; j < seq.size() - K_ - distance + 1; kmer <<= seq[j], ++j ) {
	
						//std::cout << "Kmer2: " << kmer.str() << std::endl ;
						//std::cout << "j = " << j << std::endl;
						//std::cout << "kmer : " << kmer.str() << std::endl;
						//int kmer_coverage = kmer_index_[kmer].count_;
						bucket_id kmer_bucket_id = GetKmerBucket(kmer);
						//std::cout << "kmer bucket id " << kmer_bucket_id << std::endl;

						if (kmer_bucket_id != id) continue;

						runtime_k::RtSeq kmer_d = seq.start<runtime_k::RtSeq>(K_);
			
						//std::cout << "i:" << std::endl;
						for ( size_t i = j + distance; i < j + K_ + distance; ++i ) {
						
							//std::cout <<  i << ", ";

							kmer_d <<= seq[i];
						}
						//std::cout << "Kmer_d: " << kmer_d.str() << std::endl ;
					
						int kmer_d_bucket_id = GetKmerBucket(kmer_d);
						//std::cout << "kmer d bucket id" << kmer_d_bucket_id << std::endl;

						histogram[kmer_d_bucket_id] += 1; 
						
						//std::cout << "end" << std::endl;
						//std::cout << "Kmer3: " << kmer.str() << std::endl ;
						//
					}
					//std::cout << "Kmer2: " << kmer.str() << std::endl ;
				}
			}


			for (auto it = histogram.begin(); it != histogram.end(); ++it) {

				*it = (double) *it / kmers_in_bucket_counter; 
				//std::cout << *it << " ";
			}
			
			
		}

		void GetProbeblityFromBucketToBucketForDistance ( bucket_id id_from, bucket_id, id_to, int distance ) {

			double probability = 0;
			int kmers_in_bucket_counter = GetNumberKmersInBucket(id);

			for (auto e = g_.SmartEdgeBegin(); !e.IsEnd(); ++e) {

				if (g_.length(*e) >= cfg::get().rr.max_repeat_length) {
 
					Sequence seq =  g_.EdgeNucls(*e) ;
					runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K_);
				
					//std::cout << "Sequence: " << seq.str() << std::endl ;
					//std::cout << "Sequence size: " << seq.size() << std::endl ;
					//std::cout << "j min: " << K_ << " j max: "<< seq.size()  - K_ - distance  << std::endl;
					//std::cout << "Kmer: " << kmer.str() << std::endl ;
					for (size_t j = K_; j < seq.size() - K_ - distance + 1; kmer <<= seq[j], ++j ) {
	
						//std::cout << "Kmer2: " << kmer.str() << std::endl ;
						//std::cout << "j = " << j << std::endl;
						//std::cout << "kmer : " << kmer.str() << std::endl;
						//int kmer_coverage = kmer_index_[kmer].count_;
						bucket_id kmer_bucket_id = GetKmerBucket(kmer);
						//std::cout << "kmer bucket id " << kmer_bucket_id << std::endl;

						if (kmer_bucket_id != id) continue;

						runtime_k::RtSeq kmer_d = seq.start<runtime_k::RtSeq>(K_);
			
						//std::cout << "i:" << std::endl;
						for ( size_t i = j + distance; i < j + K_ + distance; ++i ) {
						
							//std::cout <<  i << ", ";

							kmer_d <<= seq[i];
							if (kmer_bucket_id != id_to) continue;
						}
						//std::cout << "Kmer_d: " << kmer_d.str() << std::endl ;
					
						int kmer_d_bucket_id = GetKmerBucket(kmer_d);
						//std::cout << "kmer d bucket id" << kmer_d_bucket_id << std::endl;

						probability += 1; 
						
						//std::cout << "end" << std::endl;
						//std::cout << "Kmer3: " << kmer.str() << std::endl ;
						//
					}
					//std::cout << "Kmer2: " << kmer.str() << std::endl ;
				}
			}



			return probability / kmers_in_bucket_counter; 
			
		}


	};
}

#endif

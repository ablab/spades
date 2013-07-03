#ifndef BUCKET_MAPPER_HPP
#define BUCKET_MAPPER_HPP

#include <iomanip> 
#include "indices/debruijn_kmer_index.hpp"


namespace debruijn_graph {

	template <class Graph>
	class BucketMapper {

		typedef int bucket_id;

		const Graph& g_;
		const DeBruijnEdgeIndex<EdgeId>& kmer_index_;
		const unsigned K_;
		const unsigned bucketNum_;

		std::vector<std::pair<int, int> > coverage_to_multiplicity;
		std::vector<int> buckets;
   		std::vector<int> number_of_kmers_; 
    		std::set<int> distance_values_;
		std::map<int, std::vector<std::vector<double>> > cache_;

    		double CacheSearch (int distance, int shift, bucket_id id_from, bucket_id id_to) {

			int diff, min_diff = shift + 1;
			int min_dist;
			for (auto cached_distance = distance_values_.begin(); cached_distance != distance_values_.end(); ++cached_distance) {

				diff = fabs(*cached_distance - distance);
				if ( diff < min_diff ) {

					min_diff = diff;
					min_dist = *cached_distance;

				}

			}

			if ( min_diff <= shift ) {

				return cache_[min_dist][id_from][id_to];  
			}

			return -1.0;

		}

		void UpdateCache( int distance, const std::vector<std::vector<double>>& histogram ) {

			distance_values_.insert(distance);
			cache_[distance] = histogram;

		}

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

		

		int CountNumberKmersInBucket( ) {
			
			for ( int id = 0; id < bucketNum_; ++id) {
				int number_of_kmers_i = 0;
				int lower_coverage = 0, upper_coverage = buckets[id];
				if ( id > 0 ){
					lower_coverage = buckets[id - 1];
				}
			////std::cout << lower_coverage << " " << upper_coverage << std::endl;
				for ( auto it = coverage_to_multiplicity.begin(); it != coverage_to_multiplicity.end(); ++it ){

				//std::cout << it->first << std::endl;
					if (it->first > lower_coverage && it->first <= upper_coverage) number_of_kmers_i += it->second;

				}
	
				number_of_kmers_.push_back(number_of_kmers_i);
			}

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


		unsigned getK() {
			return K_;
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
			
			CountNumberKmersInBucket();
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

		int GetCoverageBucket( int kmer_coverage ) {

			 bucket_id id = 0;
			 while ( id < buckets.size() - 1 && buckets[id] < kmer_coverage ) {
				++id;
			 }

		//	 std::cout << "bucket: " << id << std::endl;
			return id;
		}



		int GetNumberKmersInBucket( int id) {
			
			return number_of_kmers_[id];

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

		double GetProbablityFromBucketToBucketForDistance ( bucket_id id_from, bucket_id id_to, int distance, int shift ) {

			double res = 0;
			res = CacheSearch (distance, shift, id_from, id_to);

			if (res != -1.0) return res;

			std::vector< std::vector<double> > bucket_to_bucket;
			for ( unsigned id = 0; id < bucketNum_; ++id ) {
			
				std::vector<double> histogram(bucketNum_,0);
				SetBucketsForDistance (  id, distance, histogram );
				bucket_to_bucket.push_back(histogram);
			}


			UpdateCache(distance,bucket_to_bucket);

			return bucket_to_bucket[id_from][id_to];
/*			int kmers_in_bucket_counter = GetNumberKmersInBucket(id_from);


			for (auto e = g_.SmartEdgeBegin(); !e.IsEnd(); ++e) {

				//std::cout << "next edge" << std::endl ;
				if (g_.length(*e) >= cfg::get().rr.max_repeat_length) {
 
					//std::cout << "length ok" << std::endl ;
					Sequence seq =  g_.EdgeNucls(*e) ;
					//std::cout << "Sequence ok " << seq.size() << std::endl ;
					runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K_);
					//std::cout << "kmer ok" << kmer.str()<< std::endl ;
				
					//std::cout << "Sequence: " << seq.str() << std::endl ;
					//std::cout << "Sequence size: " << seq.size() << std::endl ;
					//std::cout << "j min: " << K_ << " j max: "<< seq.size()  - K_ - distance  << std::endl;
					//std::cout << "Kmer: " << kmer.str() << std::endl ;
					for (size_t j = K_; j < seq.size() - K_ - distance + 1; kmer <<= seq[j], ++j ) {
	
						//std::cout << "Kmer2: " << kmer.str() << std::endl ;
						//std::cout << "j = " << j << std::endl;
						//std::cout << "kmer : " << kmer.str() << std::endl;
						//int kmer_coverage = kmer_index_[kmer].count_;
					//	std::cout << "computing kmer bucket id..." << std::endl;
						bucket_id kmer_bucket_id = GetKmerBucket(kmer);
					//	std::cout << "kmer_bucket_id: " << kmer_bucket_id << std::endl;
						//std::cout << "kmer bucket id " << kmer_bucket_id << std::endl;

						if (kmer_bucket_id != id_from) continue;

						runtime_k::RtSeq kmer_d = seq.start<runtime_k::RtSeq>(K_);
			
						//std::cout << "i:" << std::endl;
						for ( size_t i = j + distance; i < j + K_ + distance; ++i ) {
							//std::cout <<  i << ", ";

							kmer_d <<= seq[i];
						}
							
					//	std::cout << "computing kmer_d bucket id..." << std::endl;
						bucket_id kmer_d_bucket_id = GetKmerBucket(kmer_d);
						if (kmer_d_bucket_id == id_to) probability += 1;
					//	std::cout << probability << std::endl;
					//	std::cout << "kmer_d_ bucket_id: " << kmer_d_bucket_id << std::endl;
						//std::cout << "Kmer_d: " << kmer_d.str() << std::endl ;
					
						//std::cout << "kmer d bucket id" << kmer_d_bucket_id << std::endl;

						
						//std::cout << "end" << std::endl;
						//std::cout << "Kmer3: " << kmer.str() << std::endl ;
						//
					}
					//std::cout << "Kmer2: " << kmer.str() << std::endl ;
				}
				//std::cout << "out of an edge " << g_.length(*e)  << std::endl;
			}


			std::cout << "out of GetNumberKmersInBucket" << std::endl;
			return probability / kmers_in_bucket_counter; 
*/			
		}


	};
}

#endif

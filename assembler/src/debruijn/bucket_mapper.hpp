#ifndef BUCKET_MAPPER_HPP
#define BUCKET_MAPPER_HPP

#include <iomanip> 
#include "indices/debruijn_kmer_index.hpp"


namespace debruijn_graph {

	template <class Graph>
	class BucketMapper {

		typedef unsigned bucket_id;

		const Graph& g_;
		const DeBruijnEdgeIndex<Graph>& kmer_index_;
		const unsigned K_;
		const unsigned bucketNum_;

		FILE* file_;

		std::map<int,int> position_to_coverage_;
		vector<pair<int, int>> coverage_to_multiplicity;
		std::vector<int> buckets;
   		std::vector<int> number_of_kmers_; 
    		std::set<int> distance_values_;
		std::map<int, std::vector<std::vector<double>> > cache_;

    		double CacheSearch (int distance, int shift, bucket_id id_from, bucket_id id_to) {

			int diff, min_diff = shift + 1;
			int min_dist = 0;
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

		double UpdateCoverageCounters(){

			double kmer_counter = 0;
			map<int, int> used_coverages;
			for (auto e = g_.SmartEdgeBegin(); !e.IsEnd(); ++e) {

				if (g_.length(*e) >= cfg::get().rr.max_repeat_length) {
					Sequence seq =  g_.EdgeNucls(*e) ;
					runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K_);
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

					}
			}
			coverage_to_multiplicity.insert(coverage_to_multiplicity.begin(), used_coverages.begin(), used_coverages.end() );
			std::sort(coverage_to_multiplicity.begin(), coverage_to_multiplicity.end());

			return kmer_counter;
		}

		

		void CountNumberCoveragesInBucket( ) {
			
				int number_of_coverages_i = 0;
				int lower_coverage = 0;
				auto upper_coverage = buckets.begin();
				for ( auto it = coverage_to_multiplicity.begin(); it != coverage_to_multiplicity.end(); ++it ){

					if (it->first >= lower_coverage && it->first <= *upper_coverage) {
						number_of_coverages_i += it->second;
					}
					else if (it->first > *upper_coverage) {
						number_of_kmers_.push_back(number_of_coverages_i);
						number_of_coverages_i = it->second;
						lower_coverage = it->first+1;
						++upper_coverage;

					}

				}
	
				number_of_kmers_.push_back(number_of_coverages_i);

		}

		int UpdateCoverageCountersFromFile( FILE* file ){

			if (file == NULL) return 0;
			cout << "file != null" << endl;
			int position_counter = 0;
			int position, coverage;
			std::map<int, int> used_coverages;
			while ( fscanf(file, "%d%d", &position, &coverage) != EOF ) {
				//cout << position << " " << coverage << endl;
				
				position_to_coverage_[position] = coverage;
				if (used_coverages.find(coverage) != used_coverages.end()){
				
					used_coverages[coverage] += 1;
				}
				else {
					used_coverages[coverage] = 1;
				}
				position_counter += 1;
			}

			coverage_to_multiplicity.insert(coverage_to_multiplicity.begin(), used_coverages.begin(), used_coverages.end() );
			std::sort(coverage_to_multiplicity.begin(), coverage_to_multiplicity.end());
			return position_counter;
		}

		void CountBuckets( double bucket_size_bound ) {

			double current_bucket_size = 0;
			cout << "Counting buckets..." << endl;
			int i = 0;
			for ( auto it = coverage_to_multiplicity.begin(); it != coverage_to_multiplicity.end(); ++it ) {
				
				current_bucket_size += it->second;
				if ( current_bucket_size > bucket_size_bound ) {
					buckets.push_back(it->first);	
					current_bucket_size = 0;
					++i;

				}

			}

			if ( current_bucket_size <= bucket_size_bound ) {

				cout << i << ": " << current_bucket_size << "coverages" << endl;
				buckets.push_back( current_bucket_size );
			}

		}

		public:

		BucketMapper( const Graph& g, const DeBruijnEdgeIndex<Graph>& kmer_index, unsigned K, unsigned bucketNum ) : g_(g), kmer_index_(kmer_index), K_(K), bucketNum_(bucketNum) {
		
		}

		~BucketMapper() {
			
			if (file_ != NULL)
				fclose(file_);
		}

		unsigned K() {
			return K_;
		}

		void LoadProbabilities() {

			INFO("Loading probabilities...");
			if (file_ != NULL) {
	
					int distance, dimension;	
					while ( fscanf(file_, "%d %d", &distance, &dimension) != EOF ) {
						vector<vector<double>> histogram;
						for (auto i = 0; i < dimension; ++i){
							histogram.push_back(vector<double>(dimension));
							for (auto j = 0; j < dimension; ++j) {
								auto res = fscanf(file_, "%lf", &histogram[i][j] );
								VERIFY(res != EOF);
							}
							UpdateCache( distance, histogram ); 
						}
					}
				

			}

		
		}

		void SetBucketsForDistanceFromFile ( bucket_id id, int distance, int genome_size, std::vector<int>& histogram ) {

			for (unsigned i = 0; i < genome_size - distance; ++i ) {
				auto position_i_bucket_id = GetCoverageBucket(position_to_coverage_[i]);
				if ( position_i_bucket_id != id ) continue;
				unsigned j =  i + distance;
				auto position_j_bucket_id = GetCoverageBucket(position_to_coverage_[j]);
				histogram[position_j_bucket_id] += 1;
			}
		}


		double CheckChiSquareForGenomePositions( unsigned id, int distance, int genome_size ) {
		
			vector<double> counter_ij(bucketNum_, 0);
			vector<double> counter_jk(bucketNum_, 0);
			vector<vector<double>> counter_ijk;
			for (unsigned i = 0; i < bucketNum_; ++i) {
				counter_ijk.push_back(vector<double>(bucketNum_,0));
			}
			int n_j(0);
			for (unsigned i = 0; i < genome_size - 2 * distance; ++i ) {
				
				auto kmer_i_bucket_id = GetCoverageBucket(position_to_coverage_[i]);
				unsigned j =  i + distance;
				auto kmer_j_bucket_id = GetCoverageBucket(position_to_coverage_[j]);
				if (kmer_j_bucket_id != id) continue;
				n_j++;
				unsigned k =  i + 2 * distance;
				auto kmer_k_bucket_id = GetCoverageBucket(position_to_coverage_[k]);
				
				counter_ijk[kmer_i_bucket_id][kmer_k_bucket_id] += 1; 
				counter_ij[kmer_i_bucket_id] += 1; 
				counter_jk[kmer_k_bucket_id] += 1; 
		
			}

		/*	for (unsigned k = 0; k < bucketNum_; ++k) {
				printf("%5.4f ", counter_jk[k] / n_j );
			}
			printf("\n"); */
			double val = 0.0;
		
			for (unsigned i = 0; i < bucketNum_; ++i) {
				if (counter_ij[i] == 0) continue;
				for (unsigned k = 0; k < bucketNum_; ++k) {
					if ( counter_ijk[i][k] == 0 || counter_ij[i] == 0 || counter_jk[k] == 0 ) continue;
					val += counter_ij[i] * pow(counter_ijk[i][k] / counter_ij[i] - counter_jk[k] / n_j, 2) / ( counter_jk[k] / n_j );
				}
			}

			return val;

		}

		double CheckChiSquare( int id, int distance ) {
		
			vector<double> counter_ij(bucketNum_, 0);
			vector<double> counter_jk(bucketNum_, 0);
			vector<vector<double>> counter_ijk;
			for (unsigned i = 0; i < bucketNum_; ++i) {
				counter_ijk.push_back(vector<double>(bucketNum_,0));
			}
			double n_j(0);
			for (auto e = g_.SmartEdgeBegin(); !e.IsEnd(); ++e) {
				if (g_.length(*e) >= cfg::get().rr.max_repeat_length) { 
					Sequence seq =  g_.EdgeNucls(*e) ;
					runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K_);
					for (size_t i = K_; i < seq.size() - K_ - 2 * distance + 1; kmer <<= seq[i], ++i ) {
						bucket_id kmer_bucket_id = GetKmerBucket(kmer);
						runtime_k::RtSeq kmer_d_1 = seq.start<runtime_k::RtSeq>(K_);
						for ( size_t j = i + distance; j < i + K_ + distance; ++j ) {
							kmer_d_1 <<= seq[j];
						}
						int kmer_d_1_bucket_id = GetKmerBucket(kmer_d_1);
						if (kmer_d_1_bucket_id != id) continue;
						++n_j;
						runtime_k::RtSeq kmer_d_2 = seq.start<runtime_k::RtSeq>(K_);
						for ( size_t k = i + 2 * distance; k < i + K_ + 2 * distance; ++k ) {
							kmer_d_2 <<= seq[k];
						}

						int kmer_d_2_bucket_id = GetKmerBucket(kmer_d_2);
						counter_ijk[kmer_bucket_id][kmer_d_2_bucket_id] += 1; 
						counter_ij[kmer_bucket_id] += 1; 
						counter_jk[kmer_d_2_bucket_id] += 1; 
					}
				}
			}

			for (unsigned k = 0; k < bucketNum_; ++k) {
				printf("%5.4f ", counter_jk[k] / n_j );
			}
			printf("\n");
			double val = 0.0;
			for (unsigned i = 0; i < bucketNum_; ++i) {
				for (unsigned k = 0; k < bucketNum_; ++k) {
					val += counter_ij[i] * pow(counter_ijk[i][k] / counter_ij[i] - counter_jk[k] / n_j, 2) / ( counter_jk[k] / n_j );
				}
			}

			return val;

		}


		void InitBucketsFromFile() {

			FILE* file = fopen("/johnny/ksenia/ECOLI_LANE1_raw.cov", "r");
			int position_counter = UpdateCoverageCountersFromFile(file);

			CountBuckets(position_counter/bucketNum_);	

			CountNumberCoveragesInBucket();
			std::cout << "kmer_counter / bucketNum_: " << position_counter / bucketNum_ << std::endl;

			std::cout << "Buckets: " << std::endl;
			int i = 0;
			for (auto it = buckets.begin(); it != buckets.end(); ++it, ++i){
				std::cout << i << ": max coverage: " << *it << " number of kmers in bucket: " << GetNumberKmersInBucket(i) << endl;
			}
			std::cout << std::endl;

			int distance = 500;
			//std::cout << "Check Chi Square. \n  Distance: " << distance << std::endl;
			for (unsigned i = 0; i < bucketNum_; ++i) {
				//printf("%d %4.10f\n", i, CheckChiSquareForGenomePositions( i, distance, position_counter ));
				//CheckChiSquareForGenomePositions( i, distance, position_counter );
					
				std::vector<int> histogram(bucketNum_,0);
				SetBucketsForDistanceFromFile (  i, position_counter, distance, histogram );
				for (auto it = histogram.begin(); it != histogram.end(); ++it) {

					printf("d ",*it);
				}
				std::cout << std::endl;

			}

			fclose(file);

		}


		void InitBuckets( ) {

			file_ = fopen("/home/ksenia/probabilities.store", "aw");
			LoadProbabilities();

			double kmer_counter = UpdateCoverageCounters( );
					

			std::cout << "kmer_counter: " << kmer_counter << std::endl;
			CountBuckets( kmer_counter / bucketNum_ );
			
			CountNumberCoveragesInBucket();
			std::cout << "kmer_counter / bucketNum_: " << kmer_counter / bucketNum_ << std::endl;

			std::cout << "Buckets: " << std::endl;
			int i = 0;
			for (auto it = buckets.begin(); it != buckets.end(); ++it, ++i){
				std::cout << i << ": max coverage: " << *it << " number of coverages in bucket: " << GetNumberKmersInBucket(i) << endl;
			}
			std::cout << std::endl;
		
			//int distance = 500;
			/*std::cout << "Check Chi Square. \n  Distance: " << distance << std::endl;
			for (unsigned i = 0; i < bucketNum_; ++i) {
				printf("%d %4.10f\n", i, CheckChiSquare( i, distance ));
			}*/

			for ( int distance = 500; distance <= 500; distance += 1) {
			cout << distance << endl;
			std::vector< std::vector<double> > bucketToBucket;
			for ( unsigned id = 0; id < bucketNum_; ++id ) {
			
				std::vector<double> histogram(bucketNum_,0);
				SetBucketsForDistance (  id, distance, histogram );
				bucketToBucket.push_back(histogram);
			}

			std::cout << "Histogram:" << std::endl;
			
			for (auto hist = bucketToBucket.begin(); hist != bucketToBucket.end(); ++hist) {
				for (auto it = hist->begin(); it != hist->end(); ++it) {
						printf("%5.2f ",*it);
					}
					std::cout << std::endl;
				}

				std::cout << std::endl;
			}
		
		}


		int GetKmerBucket( const runtime_k::RtSeq& kmer ) {

			 int kmer_coverage = kmer_index_[kmer].count_;
			 bucket_id id = 0;
			 while ( id < buckets.size() - 1 && buckets[id] <= kmer_coverage ) {
				++id;
			 }
			return id;
		}

		int GetCoverageBucket( int kmer_coverage ) {

			 bucket_id id = 0;
			 while ( id < buckets.size() - 1 && buckets[id] <= kmer_coverage ) {
				++id;
			 }
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
					for (size_t j = K_; j < seq.size() - K_ - distance + 1; kmer <<= seq[j], ++j ) {
						bucket_id kmer_bucket_id = GetKmerBucket(kmer);
						if (kmer_bucket_id != id) continue;
						runtime_k::RtSeq kmer_d = seq.start<runtime_k::RtSeq>(K_);
						for ( size_t i = j + distance; i < j + K_ + distance; ++i ) {
							kmer_d <<= seq[i];
						}
						int kmer_d_bucket_id = GetKmerBucket(kmer_d);
						histogram[kmer_d_bucket_id] += 1; 
					}
				}
			}
			/*for (auto it = histogram.begin(); it != histogram.end(); ++it) {
				*it = (double) *it / kmers_in_bucket_counter; 
			}*/
		}


		void SaveProbability( const vector<vector<double>>& histogram, int distance ) {

			fprintf(file_, "%d %lu\n", distance,  histogram.size());
			for ( unsigned i = 0; i < histogram.size(); ++i) {
				for ( unsigned j = 0; j < histogram[i].size(); ++j) {
					fprintf(file_, "%4.5f\n", histogram[i][j] );
				}
			}
			
		}
		double GetProbabilityFromBucketToBucketForDistance ( bucket_id id_from, bucket_id id_to, int distance, int shift ) {

			double res = 0;

			//std::cout << "search in cache" << std::endl;
			res = CacheSearch (distance, shift, id_from, id_to);

			if (res != -1.0) return res;
			//std::cout << "not found" << std::endl;

			std::vector< std::vector<double> > bucket_to_bucket;
			std::cout << "buckets are updated..." << std::endl;
			for ( unsigned id = 0; id < bucketNum_; ++id ) {
				
				//std::cout << "id: " << id << std::endl;
				std::vector<double> histogram(bucketNum_,0);
				SetBucketsForDistance (  id, distance, histogram );
				bucket_to_bucket.push_back(histogram);
			}

			//std::cout << "buckets updated" << std::endl;

			UpdateCache(distance,bucket_to_bucket);

			//SaveProbability(bucket_to_bucket, distance);
			

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

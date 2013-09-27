#ifndef BUCKET_MAPPER_HPP
#define BUCKET_MAPPER_HPP

#include <iomanip> 
#include "indices/debruijn_kmer_index.hpp"


namespace debruijn_graph {


	class MarkovChain {
			
		vector<vector<int>> counter_ij_;
		vector<vector<int>> counter_jk_;
		vector<vector<vector<int>>> counter_ijk_;
		const int dim_;
		vector<int> n_j_;
		vector<vector<double>> transition_probabilities_;
		public:

			explicit MarkovChain( int dim ) : dim_(dim), n_j_(dim_,0) {
			
				vector<int> tmp(dim_,0);
				for (int i = 0; i < dim_; ++i) {
					counter_ij_.push_back(tmp);
					counter_jk_.push_back(tmp);
				}
				vector<vector<int>> tmp2;
				for (int i = 0; i < dim_; ++i) {
					tmp2.push_back(tmp);
				}
				for (int i = 0; i < dim_; ++i) {
					counter_ijk_.push_back(tmp2);
				}

			
			}

			void GenerateRandomTransitionMatrix() {
				
				for ( int i = 0; i < dim_; ++i ) {
					vector<double>random_vector;
					for ( int j = 0; j < dim_-1; ++j ){
						random_vector.push_back((double)rand() / RAND_MAX);
					}
					sort(random_vector.begin(), random_vector.end());
					transition_probabilities_.push_back(vector<double>());
					double current, prev(0.0);
					for ( int j = 0; j < dim_; ++j ){
						if ( j < dim_ - 1)
							current = random_vector[j]; 
						else 
							current = 1;
						transition_probabilities_[i].push_back(current - prev);
						if ( j < dim_ - 1)
							prev = random_vector[j];
					}
				}
				for (int i = 0; i < dim_; ++i) {
					for (int j = 0; j < dim_; ++j) {
						DEBUG(transition_probabilities_[i][j] << " ");
					}
					DEBUG("\n");
				}
			}

			void Sample(int number_of_observations){

				FILE *file = fopen("/johnny/ksenia/markov.sample","w");
				srand(time(NULL));
				for (int i = 0; i < number_of_observations; ++i) {
					int state = rand() % dim_;
					double cube_value = (double)rand() / RAND_MAX;
					double sum = 0;
					unsigned j = 0;
					while ( cube_value > sum + transition_probabilities_[state][j] ) {
						sum += transition_probabilities_[state][j];
						++j;
					}
					cube_value = (double)rand() / RAND_MAX;
					sum = 0;
					unsigned k = 0;
					while ( cube_value > sum + transition_probabilities_[j][k] ) {
						sum += transition_probabilities_[j][k];
						++k;
					}
					n_j_[j] += 1;
					counter_ij_[state][j] += 1;
					counter_jk_[j][k] += 1;
					counter_ijk_[state][j][k] += 1;
					fprintf(file, "%d %d %d\n", state, j, k);
				}
				fclose(file);
			}

			double CheckChiSquare (unsigned j) {

				double val = 0.0;
				for (int i = 0; i < dim_; ++i) {
					for (int k = 0; k < dim_; ++k) {
						if ( counter_ijk_[i][j][k] == 0 || counter_ij_[i][j] == 0 || counter_jk_[j][k] == 0 ) continue;
						val += counter_ij_[i][j] * pow( (double) counter_ijk_[i][j][k] / counter_ij_[i][j] - (double) counter_jk_[j][k] / n_j_[j], 2) / 
							( (double) counter_jk_[j][k] / n_j_[j] );
					}
				}

				return val;

			}



	};

	template <class Graph, class KmerIndex>
	class BucketMapper {

		typedef unsigned bucket_id;

		const Graph& g_;
		const KmerIndex & kmer_index_;
		const unsigned K_;
		const unsigned number_of_buckets_;

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
						int kmer_coverage = kmer_index_[kmer].count;
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
			int position_counter = 0;
			int position, coverage;
			std::map<int, int> used_coverages;
			while ( fscanf(file, "%d%d", &position, &coverage) != EOF ) {
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
			sort(coverage_to_multiplicity.begin(), coverage_to_multiplicity.end());
			return position_counter;
		}

		void CountBuckets( double bucket_size_bound ) {

			double current_bucket_size = 0;
			DEBUG("Counting buckets...\n");
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
				buckets.push_back( current_bucket_size );
			}

		}

		public:

		BucketMapper(const Graph& g,
				const KmerIndex& kmer_index,
				unsigned K, unsigned bucketNum ) : g_(g), kmer_index_(kmer_index), K_(K), number_of_buckets_(bucketNum) {
			file_ = NULL;
		}


		~BucketMapper() {
			if (file_ != NULL)
				fclose(file_);
		}

		unsigned K() const {
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
								VERIFY(fscanf(file_, "%lf", &histogram[i][j] ) != EOF);
							}
							UpdateCache( distance, histogram ); 
						}
					}
				

			}
		}
	

		void SetBucketsForDistanceFromFile3D ( int distance, int genome_size, vector<vector<vector<int>>>& histogram ) {

			for (int i = 0; i < genome_size - 2* distance; ++i ) {
				auto position_i_bucket_id = GetCoverageBucket(position_to_coverage_[i]);
				//if ( position_i_bucket_id != id ) continue;
				int j =  i + distance;
				auto position_j_bucket_id = GetCoverageBucket(position_to_coverage_[j]);
				int k =  j + distance;
				auto position_k_bucket_id = GetCoverageBucket(position_to_coverage_[k]);
				histogram[position_i_bucket_id][position_j_bucket_id][position_k_bucket_id] += 1;
			}
		}


		void SetBucketsForDistanceFromFile ( bucket_id id, int distance, int genome_size, std::vector<int>& histogram ) {

			for (int i = 0; i < genome_size - distance; ++i ) {
				auto position_i_bucket_id = GetCoverageBucket(position_to_coverage_[i]);
				if ( position_i_bucket_id != id ) continue;
				int j =  i + distance;
				auto position_j_bucket_id = GetCoverageBucket(position_to_coverage_[j]);
				histogram[position_j_bucket_id] += 1;
			}
		}


		double CheckChiSquareForGenomePositions( unsigned id, int distance, int genome_size ) {
		
			vector<double> counter_ij(number_of_buckets_, 0);
			vector<double> counter_jk(number_of_buckets_, 0);
			vector<vector<double>> counter_ijk;
			for (unsigned i = 0; i < number_of_buckets_; ++i) {
				counter_ijk.push_back(vector<double>(number_of_buckets_,0));
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

		/*	for (unsigned k = 0; k < number_of_buckets_; ++k) {
				printf("%5.4f ", counter_jk[k] / n_j );
			}
			printf("\n"); */
			double val = 0.0;
			for (unsigned i = 0; i < number_of_buckets_; ++i) {
				if (counter_ij[i] == 0) continue;
				for (unsigned k = 0; k < number_of_buckets_; ++k) {
					if ( counter_ijk[i][k] == 0 || counter_ij[i] == 0 || counter_jk[k] == 0 ) continue;
					val += counter_ij[i] * pow((double)counter_ijk[i][k] / counter_ij[i] - (double)counter_jk[k] / n_j, 2) / ( (double)counter_jk[k] / n_j );
				}
			}

			return val;

		}
		double CheckChiSquare( int id, int distance ) const {
		
			vector<double> counter_ij(number_of_buckets_, 0);
			vector<double> counter_jk(number_of_buckets_, 0);
			vector<vector<double>> counter_ijk;
			for (unsigned i = 0; i < number_of_buckets_; ++i) {
				counter_ijk.push_back(vector<double>(number_of_buckets_,0));
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
			for (unsigned k = 0; k < number_of_buckets_; ++k) {
				printf("%5.4f ", counter_jk[k] / n_j );
			}
			printf("\n");
			double val = 0.0;
			for (unsigned i = 0; i < number_of_buckets_; ++i) {
				for (unsigned k = 0; k < number_of_buckets_; ++k) {
					val += (double)counter_ij[i] * pow((double)counter_ijk[i][k] / counter_ij[i] - (double)counter_jk[k] / n_j, 2) / ( (double)counter_jk[k] / n_j );
				}
			}

			return val;

		}


		void InitBucketsFromFile() {

			FILE* file = fopen("/johnny/ksenia/ECOLI_LANE1_raw.cov", "r");
			//FILE* file = fopen("/johnny/ksenia/fake.cov", "r");
			int position_counter = UpdateCoverageCountersFromFile(file);
			CountBuckets(position_counter/number_of_buckets_);	
			CountNumberCoveragesInBucket();
			DEBUG("kmer_counter / number_of_buckets_: " << position_counter / number_of_buckets_ << "\n" );
			DEBUG("Buckets:\n");
			int i = 0;
			for (auto it = buckets.begin(); it != buckets.end(); ++it, ++i){
				DEBUG(": max coverage: " << *it << " number of kmers in bucket: " << GetNumberKmersInBucket(i) << "\n");
			}
			std::cout << std::endl;
			int distance = 500;
			vector<vector<vector<int>>> histogram;
			for (unsigned i = 0; i < number_of_buckets_; ++i) {
				histogram.push_back(vector<vector<int>>());
				for (unsigned j = 0; j < number_of_buckets_; ++j) {
					histogram[i].push_back(vector<int>(number_of_buckets_,0));
				}
			}
			SetBucketsForDistanceFromFile3D(distance, position_counter, histogram);
			for (unsigned i = 0; i < number_of_buckets_; ++i) {
				cout << i << endl;
				for (unsigned j = 0; j < number_of_buckets_; ++j) {
					for (unsigned k = 0; k < number_of_buckets_; ++k) {
			
						cout << histogram[i][j][k] << " ";
					}
					cout << endl;
				}
				cout << endl;
			}
			for (unsigned i = 0; i < number_of_buckets_; ++i) {
				printf("%d %4.10f\n", i, CheckChiSquareForGenomePositions( i, distance, position_counter ));
			}
			fclose(file);
		}

		
		void InitBuckets( ) {

			file_ = fopen("/home/ksenia/probabilities.store", "aw");
		//	LoadProbabilities();
		/*	int dim = 10;
			INFO("Initializing markov chain..");
			auto mc = MarkovChain(dim);
			INFO("Generating random matrix..");
			mc.GenerateRandomTransitionMatrix();
			INFO("Sampling...");
			mc.Sample(40000);
			for (int i = 0; i < dim; ++i) {
				INFO("Checking chi-square");
				printf("%d %4.10f\n", i, mc.CheckChiSquare(i));
			}
		*/
			double kmer_counter = UpdateCoverageCounters( );
			DEBUG("kmer_counter: " << kmer_counter << "\n");
			CountBuckets( kmer_counter / number_of_buckets_ );
			CountNumberCoveragesInBucket();
			DEBUG( "kmer_counter / number_of_buckets_: " << kmer_counter / number_of_buckets_ << "\nBuckets:\n");
			int i = 0;
			for (auto it = buckets.begin(); it != buckets.end(); ++it, ++i){
				DEBUG( i << ": max coverage: " << *it << " number of coverages in bucket: " << GetNumberKmersInBucket(i) << "\n");
			}
			/*std::cout << "Check Chi Square. \n  Distance: " << distance << std::endl;
			for (unsigned i = 0; i < number_of_buckets_; ++i) {
				printf("%d %4.10f\n", i, CheckChiSquare( i, distance ));
			}*/
		/*	std::vector< std::vector<double> > bucket_to_bucket;
			for ( unsigned id = 0; id < number_of_buckets_; ++id ) {

			
				std::vector<double> histogram(number_of_buckets_,0);
				SetBucketsForDistance (  id, distance, histogram );
				bucket_to_bucket.push_back(histogram);
			}
			std::cout << "Histogram:" << std::endl;
			for (auto hist = bucket_to_bucket.begin(); hist != bucket_to_bucket.end(); ++hist) {
				for (auto it = hist->begin(); it != hist->end(); ++it) {
						printf("%5.2f ",*it);
					}
					std::cout << std::endl;
				}

			std::cout << std::endl;
			*/
		}
		
		int GetKmerBucket( const runtime_k::RtSeq& kmer ) const {
			 int kmer_coverage = kmer_index_[kmer].count;
			 bucket_id id = 0;
			 while ( id < buckets.size() - 1 && buckets[id] <= kmer_coverage ) {
				++id;
			 }
			return id;
		}

		int GetCoverageBucket( int kmer_coverage ) const {
			 bucket_id id = 0;
			 while ( id < buckets.size() - 1 && buckets[id] <= kmer_coverage ) {
				++id;
			 }
			return id;
		}

		int GetNumberKmersInBucket( int id) const {
			return number_of_kmers_[id];

		}
		void SetBucketsForDistance ( int id, int distance, vector<double>& histogram ) const {

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
			int kmers_in_bucket_counter = GetNumberKmersInBucket(id);
			for (auto it = histogram.begin(); it != histogram.end(); ++it) {
				*it = (double) *it / kmers_in_bucket_counter; 
			}
		}


		void SetBucketsForDistance ( int distance, vector<vector<double>>& histogram ) const {
			for (auto e = g_.SmartEdgeBegin(); !e.IsEnd(); ++e) {
				if (g_.length(*e) >= cfg::get().rr.max_repeat_length) {
					Sequence seq =  g_.EdgeNucls(*e) ;
					runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K_);
					for (size_t j = K_; j < seq.size() - K_ - distance + 1; kmer <<= seq[j], ++j ) {
						bucket_id kmer_bucket_id = GetKmerBucket(kmer);
						//if (kmer_bucket_id != id) continue;
						runtime_k::RtSeq kmer_d = seq.start<runtime_k::RtSeq>(K_);
						for ( size_t i = j + distance; i < j + K_ + distance; ++i ) {
							kmer_d <<= seq[i];
						}
						int kmer_d_bucket_id = GetKmerBucket(kmer_d);
						histogram[kmer_bucket_id][kmer_d_bucket_id] += 1; 
					}
				}
			}
			for (unsigned id = 0; id < number_of_buckets_; ++id) {
				int kmers_in_bucket_counter = GetNumberKmersInBucket(id);
				for (auto it = histogram[id].begin(); it != histogram[id].end(); ++it) {
					*it = (double) *it / kmers_in_bucket_counter; 
				}
			}
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
			res = CacheSearch (distance, shift, id_from, id_to);
			if (res != -1.0) return res;
			std::vector< std::vector<double> > bucket_to_bucket;
			for ( unsigned id = 0; id < number_of_buckets_; ++id ) {
				std::vector<double> histogram(number_of_buckets_,0);
				bucket_to_bucket.push_back(histogram);
			}
			SetBucketsForDistance ( distance, bucket_to_bucket );
			UpdateCache(distance,bucket_to_bucket);
			//SaveProbability(bucket_to_bucket, distance);
			return bucket_to_bucket[id_from][id_to];
		}
	};
}

#endif

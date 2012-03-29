#pragma once 

template<class T, class Value, class Hash>
struct parallel_unordered_map
{
private:
	struct bucket_hasher
	{
		bucket_hasher(size_t bucket_size, Hash const& hash = Hash())
			: bucket_size_	(bucket_size)
			, hash_    		(hash)
		{
		}

		size_t operator()(T const& t) const
		{
			return hash_(t) % bucket_size_;
		}

	private:
		size_t 	bucket_size_;
		Hash	hash_;
	};
	//struct bucket_hasher;

	typedef	std::tr1::unordered_map<T, Value, Hash>				origin_map_t;
	typedef	std::tr1::unordered_map<T, Value, bucket_hasher>	bucket_map_t;
	typedef	std::vector<bucket_map_t>							map_arr_t;
	typedef	std::vector<Hash>									hash_arr_t;

    typedef typename origin_map_t::value_type      value_type;
    typedef typename origin_map_t::size_type       size_type;
    typedef typename origin_map_t::hasher          hasher;
    typedef typename origin_map_t::key_equal       key_equal;

	parallel_unordered_map(size_t nthreads, Hash const& hash = Hash())
		: nthreads_		(nthreads)
		, bucket_size_	(bucket_size(nthreads))
		, buckets_		(nthreads, origin_map_t(10, bucket_hasher(bucket_size(nthreads), hash)))
    	, bucket_hashs_	(nthreads, hash)
	{
	}

	bool insert(const value_type& value, const size_t thread_num)
	{
		size_t bucket_num = thread_num;
		size_t h = bucket_hashs_[bucket_num](value.first);
		return buckets_[bucket_num].insert(value);
	}

	void merge_buckets(origin_map_t& om)
	{

	}


private:
	static size_t bucket_size(size_t nthreads)
	{
		return size_t(
				std::ceil(
					(1 << sizeof(size_t)) /
					double(nthreads)));
	}

private:
	size_t 		nthreads_;
	size_t 		bucket_size_;
	map_arr_t 	buckets_;
	hash_arr_t	bucket_hashs_;
};

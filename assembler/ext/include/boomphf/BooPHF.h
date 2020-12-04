// BooPHF library
// intended to be a minimal perfect hash function with fast and low memory
// construction, at the cost of (slightly) higher bits/elem than other state of
// the art libraries once built.  should work with arbitray large number of
// elements, based on a cascade of "collision-free" bit arrays

#pragma once
#include <cstdio>
#include <climits>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>

#include <array>
#include <unordered_map>
#include <vector>
#include <memory> // for make_shared
#include <iosfwd>
#include <unistd.h>

namespace boomphf {

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark utils
////////////////////////////////////////////////////////////////

#define L8 0x0101010101010101ULL // Every lowest 8th bit set: 00000001...
#define G2 0xAAAAAAAAAAAAAAAAULL // Every highest 2nd bit: 101010...
#define G4 0x3333333333333333ULL // 00110011 ... used to group the sum of 4 bits.
#define G8 0x0F0F0F0F0F0F0F0FULL

static inline unsigned popcount_64(uint64_t x) {
    // Step 1:  00 - 00 = 0;  01 - 00 = 01; 10 - 01 = 01; 11 - 01 = 10;
    x = x - ((x & G2) >> 1);
    // step 2:  add 2 groups of 2.
    x = (x & G4) + ((x >> 2) & G4);
    // 2 groups of 4.
    x = (x + (x >> 4)) & G8;
    // Using a multiply to collect the 8 groups of 8 together.
    x = x * L8 >> 56;
    return x;
}

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark hasher
////////////////////////////////////////////////////////////////

typedef std::array<uint64_t,2> hash_pair_t;
typedef hash_pair_t internal_hash_t; // ou hash_pair_t directement ?  __uint128_t
typedef std::vector<internal_hash_t>::iterator vectorit_hash128_t;

struct InternalHasher {
    uint64_t operator()(const internal_hash_t& key) const {
        uint64_t s0 = key[0];
        uint64_t s1 = key[1];
        s1 ^= s1 << 23;
        return  (s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26)) + s0;
    }
};

template<class InnerHasher> class XorshiftHashFunctors {
    /*  Xorshift128*
        Written in 2014 by Sebastiano Vigna (vigna@acm.org)

        To the extent possible under law, the author has dedicated all copyright
        and related and neighboring rights to this software to the public domain
        worldwide. This software is distributed without any warranty.

        See <http://creativecommons.org/publicdomain/zero/1.0/>.

        This is the fastest generator passing BigCrush without
        systematic failures, but due to the relatively short period it is
        acceptable only for applications with a mild amount of parallelism;
        otherwise, use a xorshift1024* generator.

        The state must be seeded so that it is not everywhere zero. If you have
        a nonzero 64-bit seed, we suggest to pass it twice through
        MurmurHash3's avalanching function. */
  public:
    template<class Item>
    hash_pair_t hashpair128(const Item& key) const {
        auto h = inner_hasher_(key);
        return { h.first, h.second };
    }

    hash_pair_t hashpair128(const internal_hash_t &key) const {
        return key;
    }

    //return next hash an update state s
    uint64_t next(hash_pair_t &s) const {
        uint64_t s1 = s[0];
        const uint64_t s0 = s[1];
        s[0] = s0;
        s1 ^= s1 << 23; // a
        return (s[1] = (s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26))) + s0; // b, c
    }

  private:
    InnerHasher inner_hasher_;
};


////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark iterators
////////////////////////////////////////////////////////////////

template <typename Iterator>
struct iter_range {
    iter_range(Iterator b, Iterator e)
            : m_begin(std::move(b)), m_end(std::move(e)) {}

    Iterator begin() const { return m_begin; }
    Iterator end() const { return m_end; }

    Iterator m_begin, m_end;
};

template <typename Iterator>
iter_range<Iterator> range(Iterator begin, Iterator end) {
    return iter_range<Iterator>(std::move(begin), std::move(end));
}

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark BitVector
////////////////////////////////////////////////////////////////

class bitVector {

  public:

    bitVector()
            : _size(0) {
        _bitArray = nullptr;
    }

    bitVector(uint64_t n)
            : _size(n) {
        _nchar  = (1ULL+n/64ULL);
        _bitArray =  (uint64_t *) calloc(_nchar,sizeof(uint64_t));
    }

    ~bitVector() {
        if (_bitArray != nullptr)
            free(_bitArray);
    }

    //copy constructor
    bitVector(bitVector const &r) {
        _size =  r._size;
        _nchar = r._nchar;
        _ranks = r._ranks;
        _bitArray = nullptr;
        if (r._bitArray) {
            _bitArray = (uint64_t *) calloc(_nchar,sizeof(uint64_t));
            memcpy(_bitArray, r._bitArray, _nchar*sizeof(uint64_t) );
        }
    }

    // Copy assignment operator
    bitVector &operator=(bitVector const &r) {
        if (&r != this) {
            _size =  r._size;
            _nchar = r._nchar;
            _ranks = r._ranks;
            if (_bitArray != nullptr)
                free(_bitArray);
            _bitArray = nullptr;
            if (r._bitArray) {
                _bitArray = (uint64_t *) calloc(_nchar, sizeof(uint64_t));
                memcpy(_bitArray, r._bitArray, _nchar*sizeof(uint64_t) );
            }
        }
        return *this;
    }

    // Move assignment operator
    bitVector &operator=(bitVector &&r) noexcept {
        if (&r != this) {
            if (_bitArray != nullptr)
                free(_bitArray);

            _size =  r._size;
            _nchar = r._nchar;
            _ranks = std::move(r._ranks);
            _bitArray = r._bitArray;
            r._bitArray = nullptr;
        }
        return *this;
    }

    // Move constructor
    bitVector(bitVector &&r) noexcept
            : _bitArray( nullptr) ,_size(0)  {
        *this = std::move(r);
    }


    void resize(uint64_t newsize) {
        _nchar  = (1ULL+newsize/64ULL);
        _bitArray = (uint64_t *) realloc(_bitArray,_nchar*sizeof(uint64_t));
        _size = newsize;
    }

    size_t size() const { return _size; }
    uint64_t bitSize() const {return (_nchar*64ULL + _ranks.capacity()*64ULL );}

    //clear whole array
    void clear() {
        memset(_bitArray,0,_nchar*sizeof(uint64_t));
    }

    //clear collisions in interval, only works with start and size multiple of 64
    void clearCollisions(uint64_t start, size_t size, bitVector * cc) {
        assert( (start & 63) ==0);
        assert( (size & 63) ==0);
        uint64_t ids = (start/64ULL);
        for (uint64_t ii =0;  ii< (size/64ULL); ii++) {
            _bitArray[ids+ii] =  _bitArray[ids+ii] & (~ (cc->get64(ii)) );
        }

        cc->clear();
    }


    //clear interval, only works with start and size multiple of 64
    void clear(uint64_t start, size_t size) {
        assert( (start & 63) ==0);
        assert( (size & 63) ==0);
        memset(_bitArray + (start/64ULL),0,(size/64ULL)*sizeof(uint64_t));
    }

    //for debug purposes
    void print() const {
        printf("bit array of size %lli: \n", _size);
        for (uint64_t ii = 0; ii< _size; ii++) {
            if (ii%10==0)
                printf(" (%llu) ",ii);
            int val = (_bitArray[ii >> 6] >> (ii & 63 ) ) & 1;
            printf("%i",val);
        }
        printf("\n");

        printf("rank array : size %lu \n",_ranks.size());
        for (uint64_t ii = 0; ii< _ranks.size(); ii++) {
            printf("%llu :  %lli,  ",ii,_ranks[ii]);
        }
        printf("\n");
    }

    // return value at pos
    uint64_t operator[](uint64_t pos) const {
        return (_bitArray[pos >> 6ULL] >> (pos & 63)) & 1;
    }

    //atomically   return old val and set to 1
    uint64_t atomic_test_and_set(uint64_t pos) {
        uint64_t oldval = __sync_fetch_and_or(_bitArray + (pos >> 6), (uint64_t) (1ULL << (pos & 63)) );
        return (oldval >> (pos & 63)) & 1;
    }


    uint64_t get(uint64_t pos) const {
        return (*this)[pos];
    }

    uint64_t get64(uint64_t cell64) const {
        return _bitArray[cell64];
    }

    //set bit pos to 1
    void set(uint64_t pos) {
        assert(pos<_size);
        __sync_fetch_and_or(_bitArray + (pos >> 6ULL), (1ULL << (pos & 63)) );
    }

    //set bit pos to 0
    void reset(uint64_t pos) {
        __sync_fetch_and_and(_bitArray + (pos >> 6ULL), ~(1ULL << (pos & 63) ));
    }

    // return value of last rank
    // add offset to all ranks computed
    uint64_t build_ranks(uint64_t offset = 0) {
        _ranks.reserve(2 + _size/_nb_bits_per_rank_sample);

        uint64_t curent_rank = offset;
        for (size_t ii = 0; ii < _nchar; ii++) {
            if (((ii*64) % _nb_bits_per_rank_sample) == 0) {
                _ranks.push_back(curent_rank);
            }
            curent_rank +=  popcount_64(_bitArray[ii]);
        }

        return curent_rank;
    }

    uint64_t rank(uint64_t pos) const {
        uint64_t word_idx = pos / 64ULL;
        uint64_t word_offset = pos % 64;
        uint64_t block = pos / _nb_bits_per_rank_sample;
        uint64_t r = _ranks[block];
        for (uint64_t w = block * _nb_bits_per_rank_sample / 64; w < word_idx; ++w)
            r += popcount_64(_bitArray[w]);
        uint64_t mask = (uint64_t(1) << word_offset ) - 1;
        r += popcount_64( _bitArray[word_idx] & mask);

        return r;
    }

    void save(std::ostream& os) const {
        os.write(reinterpret_cast<char const*>(&_size), sizeof(_size));
        os.write(reinterpret_cast<char const*>(&_nchar), sizeof(_nchar));
        os.write(reinterpret_cast<char const*>(_bitArray), (std::streamsize)(sizeof(uint64_t) * _nchar));
        size_t sizer = _ranks.size();
        os.write(reinterpret_cast<char const*>(&sizer),  sizeof(size_t));
        os.write(reinterpret_cast<char const*>(_ranks.data()), (std::streamsize)(sizeof(_ranks[0]) * _ranks.size()));
    }

    void load(std::istream& is) {
        is.read(reinterpret_cast<char*>(&_size), sizeof(_size));
        is.read(reinterpret_cast<char*>(&_nchar), sizeof(_nchar));
        this->resize(_size);
        is.read(reinterpret_cast<char *>(_bitArray), (std::streamsize)(sizeof(uint64_t) * _nchar));

        size_t sizer;
        is.read(reinterpret_cast<char *>(&sizer),  sizeof(size_t));
        _ranks.resize(sizer);
        is.read(reinterpret_cast<char*>(_ranks.data()), (std::streamsize)(sizeof(_ranks[0]) * _ranks.size()));
    }


  protected:
    uint64_t*  _bitArray;
    uint64_t _size;
    uint64_t _nchar;

    // epsilon =  64 / _nb_bits_per_rank_sample   bits
    // additional size for rank is epsilon * _size
    static constexpr uint64_t _nb_bits_per_rank_sample = 512; //512 seems ok
    std::vector<uint64_t> _ranks;
};

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark level
////////////////////////////////////////////////////////////////

static inline uint64_t fastrange64(uint64_t word, uint64_t p) {
    return (uint64_t)(((__uint128_t)word * (__uint128_t)p) >> 64);
}

class level{
  public:
    level() {}

    ~level() {}

    uint64_t get(uint64_t hash_raw) const {
        uint64_t hashi = fastrange64(hash_raw, hash_domain);
        return bitset.get(hashi);
    }

    uint64_t hash_domain;
    bitVector bitset;
};


////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark mphf
////////////////////////////////////////////////////////////////

/* Hasher_t returns a single hash when operator()(elem_t key) is called.
   if used with XorshiftHashFunctors, it must have the following operator: operator()(elem_t key, uint64_t seed) */
template<typename Hasher_t>
class mphf {
    /* this mechanisms gets P hashes out of Hasher_t */
    typedef XorshiftHashFunctors<Hasher_t> MultiHasher_t ;

  public:
    static constexpr uint64_t NOT_FOUND = -1ULL;

    enum ConflictPolicy {
        Error,
        Warning,
        Ignore
    };

    mphf()
            : _built(false) {}

    ~mphf() {}

    // allow perc_elem_loaded  elements to be loaded in ram for faster construction (default 3%), set to 0 to desactivate

    mphf(size_t n,
         ConflictPolicy policy = ConflictPolicy::Warning,
         double gamma = 2.0, float perc_elem_loaded = 0.03f, unsigned nb_levels = 25)
            : _built(false) {
        init(n, policy, gamma, perc_elem_loaded, nb_levels);
    }

    void init(size_t n,
              ConflictPolicy policy = ConflictPolicy::Warning,
              double gamma = 2.0, float perc_elem_loaded = 0.03f, unsigned nb_levels = 25) {
        _nb_levels = nb_levels;
        _gamma = gamma;
        _hash_domain = size_t(ceil(double(n) * gamma));
        _nelem = n;
        _policy = policy;
        _percent_elem_loaded_for_fastMode = perc_elem_loaded;
        _fastmode = _percent_elem_loaded_for_fastMode > 0.0;

        setup();
    }

    template<typename Range>
    void build(const Range &input_range) {
        if (_nelem == 0)
            return;

        uint64_t offset = 0;
        for (unsigned i_level = 0; i_level < _nb_levels; ++i_level) {
            auto &level = _levels[i_level];
            bitVector collisions(level.hash_domain); // temp collision bitarray for this level
            processLevel(input_range, i_level, collisions);
            level.bitset.clearCollisions(0, level.hash_domain, &collisions);
            offset = level.bitset.build_ranks(offset);
        }

        _lastbitsetrank = offset;
        std::vector<internal_hash_t>().swap(setLevelFastmode);   // clear setLevelFastmode reallocating

        _built = true;
    }

    template<typename Range>
    void build(const std::vector<Range> &ranges,
               unsigned nthreads = 1) {
        if (_nelem == 0)
            return;

        uint64_t offset = 0;
        for (unsigned i_level = 0; i_level < _nb_levels; ++i_level) {
            auto &level = _levels[i_level];
            bitVector collisions(level.hash_domain); // temp collision bitarray for this level
            processLevel(ranges, i_level, collisions, nthreads);
            level.bitset.clearCollisions(0, level.hash_domain, &collisions);
            offset = level.bitset.build_ranks(offset);
        }

        _lastbitsetrank = offset;
        std::vector<internal_hash_t>().swap(setLevelFastmode);   // clear setLevelFastmode reallocating

        _built = true;
    }


    template<class elem_t>
    uint64_t lookup(const elem_t &elem) const {
        if (!_built) return NOT_FOUND;

        uint64_t non_minimal_hp;
        unsigned level;

        hash_pair_t bbhash = _hasher.hashpair128(elem);
        uint64_t level_hash = getLevel(bbhash, &level, _nb_levels);

        if (level == (_nb_levels-1)) {
            auto in_final_map = _final_hash.find(bbhash);
            if (in_final_map == _final_hash.end())
                return NOT_FOUND;

            return (in_final_map->second != NOT_FOUND ?
                    in_final_map->second + _lastbitsetrank : NOT_FOUND);
        } else {
            non_minimal_hp = fastrange64(level_hash, _levels[level].hash_domain);
        }

        return _levels[level].bitset.rank(non_minimal_hp); // minimal_hp
    }

    uint64_t size() const {
        return _nelem;
    }

    uint64_t mem_size() const {
        if (!_built)
            return 0;

        uint64_t totalsizeBitset = 0;
        for (unsigned ii = 0; ii < _nb_levels; ii++)
            totalsizeBitset += _levels[ii].bitset.bitSize();

        uint64_t totalsize = totalsizeBitset +  _final_hash.size()*42*8 ;  // unordered map takes approx 42B per elem [personal test] (42B with uint64_t key, would be larger for other type of elem)

        return totalsize / 8;
    }

    double prob_collision() const {
        return _proba_collision;
    }

    uint64_t last_level_size() const {
        return _final_hash.size();
    }

    void save(std::ostream& os) const {
        os.write(reinterpret_cast<char const*>(&_gamma), sizeof(_gamma));
        os.write(reinterpret_cast<char const*>(&_nb_levels), sizeof(_nb_levels));
        os.write(reinterpret_cast<char const*>(&_lastbitsetrank), sizeof(_lastbitsetrank));
        os.write(reinterpret_cast<char const*>(&_nelem), sizeof(_nelem));
        for (int ii=0; ii<_nb_levels; ii++) {
            _levels[ii].bitset.save(os);
        }

        //save final hash
        size_t final_hash_size = _final_hash.size();

        os.write(reinterpret_cast<char const*>(&final_hash_size), sizeof(size_t));
        for (auto it = _final_hash.begin(); it != _final_hash.end(); ++it) {
            os.write(reinterpret_cast<char const*>(&(it->first)), sizeof(internal_hash_t));
            os.write(reinterpret_cast<char const*>(&(it->second)), sizeof(uint64_t));
        }

    }

    void load(std::istream& is) {
        is.read(reinterpret_cast<char*>(&_gamma), sizeof(_gamma));
        is.read(reinterpret_cast<char*>(&_nb_levels), sizeof(_nb_levels));
        is.read(reinterpret_cast<char*>(&_lastbitsetrank), sizeof(_lastbitsetrank));
        is.read(reinterpret_cast<char*>(&_nelem), sizeof(_nelem));

        _levels.resize(_nb_levels);
        for (int ii=0; ii<_nb_levels; ii++)
            _levels[ii].bitset.load(is);

        // mini setup, recompute size of each level
        _proba_collision = 1.0 -  pow(((_gamma*(double)_nelem -1 ) / (_gamma*(double)_nelem)),_nelem-1);
        uint64_t previous_idx =0;
        _hash_domain = (size_t)(ceil(double(_nelem) * _gamma)) ;
        for (int ii=0; ii<_nb_levels; ii++) {
            _levels[ii].hash_domain =  ((uint64_t(_hash_domain * pow(_proba_collision,ii)) + 63) / 64) * 64;
            if (_levels[ii].hash_domain == 0)
                _levels[ii].hash_domain = 64;
        }

        //restore final hash

        _final_hash.clear();
        size_t final_hash_size ;

        is.read(reinterpret_cast<char *>(&final_hash_size), sizeof(size_t));

        for (unsigned int ii=0; ii<final_hash_size; ii++) {
            internal_hash_t key;
            uint64_t value;

            is.read(reinterpret_cast<char *>(&key), sizeof(internal_hash_t));
            is.read(reinterpret_cast<char *>(&value), sizeof(uint64_t));

            _final_hash[key] = value;
        }
        _built = true;
    }


  private:
    void setup() {
        if (_fastmode)
            setLevelFastmode.resize(_percent_elem_loaded_for_fastMode * (double)_nelem);

        _proba_collision = 1.0 -  pow(((_gamma*(double)_nelem -1 ) / (_gamma*(double)_nelem)),_nelem-1);
        _levels.resize(_nb_levels);

        // build levels
        for (unsigned ii = 0; ii<_nb_levels; ii++) {
            // round size to nearest superior multiple of 64, makes it easier to clear a level
            _levels[ii].hash_domain = ((uint64_t(_hash_domain * pow(_proba_collision, ii)) + 63) / 64) * 64;
            if (_levels[ii].hash_domain == 0)
                _levels[ii].hash_domain = 64;
        }

        _fastModeLevel = _nb_levels;
        for (unsigned ii = 0; ii < _nb_levels; ii++) {
            if (pow(_proba_collision, ii) < _percent_elem_loaded_for_fastMode) {
                _fastModeLevel = ii;
                break;
            }
        }
    }

    constexpr uint64_t iterate_hash(hash_pair_t &bbhash, unsigned level) const {
        if (level == 0)
            return bbhash[0];
        else if (level == 1)
            return bbhash[1];

        return _hasher.next(bbhash);
    }

    // compute level and returns hash of last level reached
    uint64_t getLevel(internal_hash_t bbhash, unsigned *res_level, unsigned maxlevel) const {
        unsigned level = 0;
        uint64_t hash_raw = 0;

        for (level = 0; level < _nb_levels - 1 && level < maxlevel; ++level) {
            hash_raw = iterate_hash(bbhash, level);
            if (_levels[level].get(hash_raw)) {
                *res_level = level;
                return hash_raw;
            }
        }

        *res_level = level;
        return iterate_hash(bbhash, level);
    }

    // insert into bitarray
    void insertIntoLevel(uint64_t level_hash, int level,
                         bitVector &collisions) {
        uint64_t hashl = fastrange64(level_hash, _levels[level].hash_domain);

        if (_levels[level].bitset.atomic_test_and_set(hashl))
            collisions.atomic_test_and_set(hashl);
    }

    void processHash(internal_hash_t val, unsigned i, bitVector &collisions) {
        unsigned level; uint64_t level_hash;
        level_hash = getLevel(val, &level, i);

        if (level != i)
            return;

        // insert into lvl i
        if (_fastmode && i == _fastModeLevel) {
            uint64_t idxl2 = __sync_fetch_and_add(&_idxLevelsetLevelFastmode,1);
            //si depasse taille attendue pour setLevelFastmode, fall back sur slow mode mais devrait pas arriver si hash ok et proba avec nous
            if (idxl2 >= setLevelFastmode.size())
                _fastmode = false;
            else
                setLevelFastmode[idxl2] = val; // create set for fast mode
        }

        // insert to level i+1 : either next level of the cascade or final hash if last level reached
        if (i == _nb_levels-1) { //stop cascade here, insert into exact hash
            uint64_t hashidx =  __sync_fetch_and_add(&_final_hashidx, 1);

            // calc rank de fin  precedent level qq part, puis init hashidx avec ce rank, direct minimal, pas besoin inser ds bitset et rank
#           pragma omp critical
            {
                if (_final_hash.count(val)) { // key already in final hash
                    if (_policy == ConflictPolicy::Ignore) {
                        _final_hash[val] = NOT_FOUND;
                    } else {
                        fprintf(stderr,"The impossible happened : collision on 128 bit hashes... please switch to safe branch, and play the lottery.");
                        fprintf(stderr,"Another more likely explanation might be that you have duplicate keys in your input.\
                                        If so, you can ignore this message, but be aware that too many duplicate keys will increase ram usage\n");
                        if (_policy == ConflictPolicy::Error)
                            abort();
                    }
                } else {
                    _final_hash[val] = hashidx;
                }
            }
        } else {
            insertIntoLevel(level_hash, i, collisions); //should be safe
        }
    }

    template<typename Range>
    void processLevel(Range const& input_range,
                      unsigned level, bitVector &collisions) {
        _levels[level].bitset = bitVector(_levels[level].hash_domain);

        _final_hashidx = 0;
        _idxLevelsetLevelFastmode = 0;

        if (_fastmode && level > _fastModeLevel) {
            for (const auto &entry : setLevelFastmode)
                processHash(_hasher.hashpair128(entry), level, collisions);
        } else {
            for (const auto &entry : input_range)
                processHash(_hasher.hashpair128(entry), level, collisions);
        }

        if (_fastmode && level == _fastModeLevel) { //shrink to actual number of elements in set
            setLevelFastmode.resize(_idxLevelsetLevelFastmode);
        }
    }

    template<typename Range>
    void processLevel(std::vector<Range> const& ranges,
                      unsigned level, bitVector &collisions,
                      unsigned nthreads) {
        _levels[level].bitset = bitVector(_levels[level].hash_domain);

        _final_hashidx = 0;
        _idxLevelsetLevelFastmode = 0;

        if (_fastmode && level > _fastModeLevel) {
#           pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < setLevelFastmode.size(); ++i) {
                processHash(_hasher.hashpair128(setLevelFastmode[i]), level, collisions);
            }
        } else {
#           pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < ranges.size(); ++i) {
                for (const auto &entry : ranges[i])
                    processHash(_hasher.hashpair128(entry), level, collisions);
            }
        }

        if (_fastmode && level == _fastModeLevel) { //shrink to actual number of elements in set
            setLevelFastmode.resize(_idxLevelsetLevelFastmode);
        }
    }

  private:
    std::vector<level> _levels;
    int _nb_levels;

    MultiHasher_t _hasher;

    double _gamma;
    uint64_t _hash_domain;
    uint64_t _nelem;
    double _proba_collision;
    uint64_t _lastbitsetrank;
    ConflictPolicy _policy;
    std::unordered_map<internal_hash_t,uint64_t, InternalHasher> _final_hash; // InternalHasher   Hasher_t
    uint64_t _final_hashidx;

    // fast build mode , requires  that _percent_elem_loaded_for_fastMode %   elems are loaded in ram
    float _percent_elem_loaded_for_fastMode;
    bool _fastmode;
    uint64_t _idxLevelsetLevelFastmode;
    std::vector<internal_hash_t> setLevelFastmode;
    int _fastModeLevel;

    bool _built;
};

}

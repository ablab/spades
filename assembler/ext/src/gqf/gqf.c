#include <stdlib.h>
#if 0
#undef NDEBUG
# include <assert.h>
#else
# define assert(x)
#endif
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <gqf/gqf.h>

/* Can be 
	 0 (choose size at run-time), 
	 8, 16, 32, or 64 (for optimized versions),
	 or other integer <= 56 (for compile-time-optimized bit-shifting-based versions)
	 */
#define BITS_PER_SLOT 0

/******************************************************************
 * Code for managing the metadata bits and slots w/o interpreting *
 * the content of the slots.
 ******************************************************************/

/* Must be >= 6.  6 seems fastest. */
#define BLOCK_OFFSET_BITS (6)

#define SLOTS_PER_BLOCK (1ULL << BLOCK_OFFSET_BITS)
#define METADATA_WORDS_PER_BLOCK ((SLOTS_PER_BLOCK + 63) / 64)

#define NUM_SLOTS_TO_LOCK (1ULL<<12)
#define CLUSTER_SIZE (1ULL<<8)

#define METADATA_WORD(qf,field,slot_index) (get_block((qf), (slot_index) / SLOTS_PER_BLOCK)->field[((slot_index) % SLOTS_PER_BLOCK) / 64])
#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) - 1ULL)
#define MAX_VALUE(nbits) ((1ULL << (nbits)) - 1)
#define BILLION 1000000000L

typedef struct __attribute__ ((__packed__)) qfblock {
	uint8_t offset; /* Code works with uint16_t, uint32_t, etc, but uint8_t seems just as fast as anything else */
	uint64_t occupieds[METADATA_WORDS_PER_BLOCK];
	uint64_t runends[METADATA_WORDS_PER_BLOCK];

#if BITS_PER_SLOT == 8
	uint8_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT == 16
	uint16_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT == 32
	uint32_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT == 64
	uint64_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT != 0
	uint8_t   slots[SLOTS_PER_BLOCK * BITS_PER_SLOT / 8];
#else
	uint8_t   slots[];
#endif
} qfblock;

static __inline__ unsigned long long rdtsc(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#ifdef LOG_WAIT_TIME
static inline bool qf_spin_lock(QF *qf, volatile int *lock, uint64_t idx, bool
																flag_spin)
{
	struct timespec start, end;
	bool ret;

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
	if (!flag_spin) {
		ret = !__sync_lock_test_and_set(lock, 1);
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
		qf->wait_times[idx].locks_acquired_single_attempt++;
		qf->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
																												start.tv_sec) +
			end.tv_nsec - start.tv_nsec;
	} else {
		if (!__sync_lock_test_and_set(lock, 1)) {
			clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
			qf->wait_times[idx].locks_acquired_single_attempt++;
			qf->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
																													start.tv_sec) +
			end.tv_nsec - start.tv_nsec;
		} else {
			while (__sync_lock_test_and_set(lock, 1))
				while (*lock);
			clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
			qf->wait_times[idx].total_time_spinning += BILLION * (end.tv_sec -
																														start.tv_sec) +
				end.tv_nsec - start.tv_nsec;
		}
		ret = true;
	}
	qf->wait_times[idx].locks_taken++;

	return ret;

	/*start = rdtsc();*/
	/*if (!__sync_lock_test_and_set(lock, 1)) {*/
		/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);*/
		/*qf->wait_times[idx].locks_acquired_single_attempt++;*/
		/*qf->wait_times[idx].total_time_single += BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;*/
	/*} else {*/
		/*while (__sync_lock_test_and_set(lock, 1))*/
			/*while (*lock);*/
		/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);*/
		/*qf->wait_times[idx].total_time_spinning += BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;*/
	/*}*/

	/*end = rdtsc();*/
	/*qf->wait_times[idx].locks_taken++;*/
	/*return;*/
}
#else
/**
 * Try to acquire a lock once and return even if the lock is busy.
 * If spin flag is set, then spin until the lock is available.
 */
static inline bool qf_spin_lock(volatile int *lock, bool flag_spin)
{
	if (!flag_spin) {
		return !__sync_lock_test_and_set(lock, 1);
	} else {
		while (__sync_lock_test_and_set(lock, 1))
			while (*lock);
		return true;
	}

	return false;
}
#endif

/*static inline void qf_spin_lock(volatile int *lock)*/
/*{*/
	/*while (__sync_lock_test_and_set(lock, 1))*/
		/*while (*lock);*/
	/*return;*/
/*}*/

static inline void qf_spin_unlock(volatile int *lock)
{
	__sync_lock_release(lock);
	return;
}

/*
static void modify_metadata(QF *qf, uint64_t *metadata, int cnt)
{
#ifdef LOG_WAIT_TIME
	qf_spin_lock(qf, &qf->mem->metadata_lock,qf->num_locks, true);
#else
	qf_spin_lock(&qf->mem->metadata_lock, true);
#endif
	*metadata = *metadata + cnt;
	qf_spin_unlock(&qf->mem->metadata_lock);
	return;
}
*/

static void modify_metadata(QF *qf, uint64_t *metadata, int cnt) {
    __sync_fetch_and_add(metadata, cnt);
}


static inline int popcnt(uint64_t val) {
       asm("popcnt %[val], %[val]"
                       : [val] "+r" (val)
                       :
                       : "cc");
       return val;
}

static inline int popcntv(const uint64_t val, int ignore)
{
	if (ignore % 64)
		return popcnt(val & ~BITMASK(ignore % 64));
	else
		return popcnt(val);
}

// Returns the number of 1s up to (and including) the pos'th bit
// Bits are numbered from 0
static inline int bitrank(uint64_t val, int pos) {
	val = val & ((2ULL << pos) - 1);
    return popcnt(val);
}

static const uint8_t kSelectInByte[2048] = {
  8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0,
  1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0,
  2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0,
  1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0,
  3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0,
  1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0,
  2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0,
  1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
  4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0,
  1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8, 8, 8, 1,
  8, 2, 2, 1, 8, 3, 3, 1, 3, 2, 2, 1, 8, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2,
  2, 1, 8, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1,
  4, 3, 3, 1, 3, 2, 2, 1, 8, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4,
  4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1,
  3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 7, 7, 1, 7, 2,
  2, 1, 7, 3, 3, 1, 3, 2, 2, 1, 7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
  7, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3,
  3, 1, 3, 2, 2, 1, 7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1,
  4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2,
  2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 8, 8, 8, 8, 8, 8, 2,
  8, 8, 8, 3, 8, 3, 3, 2, 8, 8, 8, 4, 8, 4, 4, 2, 8, 4, 4, 3, 4, 3, 3, 2, 8, 8,
  8, 5, 8, 5, 5, 2, 8, 5, 5, 3, 5, 3, 3, 2, 8, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3,
  4, 3, 3, 2, 8, 8, 8, 6, 8, 6, 6, 2, 8, 6, 6, 3, 6, 3, 3, 2, 8, 6, 6, 4, 6, 4,
  4, 2, 6, 4, 4, 3, 4, 3, 3, 2, 8, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2,
  6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 7, 8, 7, 7, 2, 8, 7,
  7, 3, 7, 3, 3, 2, 8, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2, 8, 7, 7, 5,
  7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3,
  3, 2, 8, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4, 6, 4, 4, 2,
  6, 4, 4, 3, 4, 3, 3, 2, 7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6, 5,
  5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 4, 8, 4, 4, 3, 8, 8, 8, 8, 8, 8,
  8, 5, 8, 8, 8, 5, 8, 5, 5, 3, 8, 8, 8, 5, 8, 5, 5, 4, 8, 5, 5, 4, 5, 4, 4, 3,
  8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 3, 8, 8, 8, 6, 8, 6, 6, 4, 8, 6,
  6, 4, 6, 4, 4, 3, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6, 5, 6, 5, 5, 3, 8, 6, 6, 5,
  6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7,
  7, 3, 8, 8, 8, 7, 8, 7, 7, 4, 8, 7, 7, 4, 7, 4, 4, 3, 8, 8, 8, 7, 8, 7, 7, 5,
  8, 7, 7, 5, 7, 5, 5, 3, 8, 7, 7, 5, 7, 5, 5, 4, 7, 5, 5, 4, 5, 4, 4, 3, 8, 8,
  8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 3, 8, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4,
  6, 4, 4, 3, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 3, 7, 6, 6, 5, 6, 5,
  5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 4, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6,
  6, 4, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 6, 8, 6, 6, 5,
  8, 6, 6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8,
  8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 4, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7,
  8, 7, 7, 5, 8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 4, 8, 8, 8, 8, 8, 8,
  8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 4,
  8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6,
  6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6,
  8, 6, 6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7,
  8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 8,
  8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6,
  6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7
};

static inline uint64_t bitselect(uint64_t x, uint64_t k) {
  if (k >= popcnt(x)) { return 64; }

  const uint64_t kOnesStep4  = 0x1111111111111111ULL;
  const uint64_t kOnesStep8  = 0x0101010101010101ULL;
  const uint64_t kMSBsStep8  = 0x80ULL * kOnesStep8;

  uint64_t s = x;
  s = s - ((s & 0xA * kOnesStep4) >> 1);
  s = (s & 0x3 * kOnesStep4) + ((s >> 2) & 0x3 * kOnesStep4);
  s = (s + (s >> 4)) & 0xF * kOnesStep8;
  uint64_t byteSums = s * kOnesStep8;

  uint64_t kStep8 = k * kOnesStep8;
  uint64_t geqKStep8 = (((kStep8 | kMSBsStep8) - byteSums) & kMSBsStep8);
  uint64_t place = popcnt(geqKStep8) * 8;
  uint64_t byteRank = k - (((byteSums << 8) >> place) & (uint64_t)(0xFF));
  return place + kSelectInByte[((x >> place) & 0xFF) | (byteRank << 8)];
}

static inline uint64_t bitselectv(const uint64_t val, int ignore, int rank)
{
	return bitselect(val & ~BITMASK(ignore % 64), rank);
}

#if BITS_PER_SLOT > 0
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
	return &qf->blocks[block_index];
}
#else
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
	return (qfblock *)(((char *)qf->blocks) + block_index * (sizeof(qfblock) +
																													 SLOTS_PER_BLOCK *
																													 qf->metadata->bits_per_slot /
																													 8));
}
#endif

static inline int is_runend(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, runends, index) >> ((index % SLOTS_PER_BLOCK) % 64)) & 1ULL;
}

static inline int is_occupied(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, occupieds, index) >> ((index % SLOTS_PER_BLOCK) % 64)) & 1ULL;
}

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->metadata->xnslots);
	return get_block(qf, index / SLOTS_PER_BLOCK)->slots[index % SLOTS_PER_BLOCK];
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	get_block(qf, index / SLOTS_PER_BLOCK)->slots[index % SLOTS_PER_BLOCK] =
		value & BITMASK(qf->metadata->bits_per_slot);
}

#elif BITS_PER_SLOT > 0

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems to generate buggy code.  :/  */
	assert(index < qf->metadata->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * BITS_PER_SLOT / 8];
	return (uint64_t)(((*p) >> (((index % SLOTS_PER_BLOCK) * BITS_PER_SLOT) %
															8)) & BITMASK(BITS_PER_SLOT));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems to generate buggy code.  :/  */
	assert(index < qf->metadata->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * BITS_PER_SLOT / 8];
	uint64_t t = *p;
	uint64_t mask = BITMASK(BITS_PER_SLOT);
	uint64_t v = value;
	int shift = ((index % SLOTS_PER_BLOCK) * BITS_PER_SLOT) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

#else

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	return (uint64_t)(((*p) >> (((index % SLOTS_PER_BLOCK) * qf->metadata->bits_per_slot)
															% 8)) & BITMASK(qf->metadata->bits_per_slot));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	uint64_t t = *p;
	uint64_t mask = BITMASK(qf->metadata->bits_per_slot);
	uint64_t v = value;
	int shift = ((index % SLOTS_PER_BLOCK) * qf->metadata->bits_per_slot) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

#endif

static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index);

static inline uint64_t block_offset(const QF *qf, uint64_t blockidx)
{
	/* If we have extended counters and a 16-bit (or larger) offset
		 field, then we can safely ignore the possibility of overflowing
		 that field. */
	if (sizeof(qf->blocks[0].offset > 1) || 
			get_block(qf, blockidx)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
		return get_block(qf, blockidx)->offset;

	return run_end(qf, SLOTS_PER_BLOCK * blockidx - 1) - SLOTS_PER_BLOCK *
		blockidx + 1;
}

static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index)
{
	uint64_t bucket_block_index       = hash_bucket_index / SLOTS_PER_BLOCK;
	uint64_t bucket_intrablock_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t bucket_blocks_offset = block_offset(qf, bucket_block_index);

	// uint64_t bucket_intrablock_rank   = bitrankv(get_block(qf, bucket_block_index)->occupieds, METADATA_WORDS_PER_BLOCK, bucket_intrablock_offset);
	uint64_t bucket_intrablock_rank   = bitrank(get_block(qf,
																												bucket_block_index)->occupieds[0],
																							bucket_intrablock_offset);

	if (bucket_intrablock_rank == 0) {
		if (bucket_blocks_offset <= bucket_intrablock_offset)
			return hash_bucket_index;
		else
			return SLOTS_PER_BLOCK * bucket_block_index + bucket_blocks_offset - 1;
	}

	/* FIXME: Should handle offset overflow  --- Ah screw it, just use 16-bit offset.  That'll never overflow.  */
	uint64_t runend_block_index  = bucket_block_index + bucket_blocks_offset /
		SLOTS_PER_BLOCK;
	uint64_t runend_ignore_bits  = bucket_blocks_offset % SLOTS_PER_BLOCK;
	uint64_t runend_rank         = bucket_intrablock_rank - 1;
	// uint64_t runend_block_offset = bitselectv(get_block(qf, runend_block_index)->runends, METADATA_WORDS_PER_BLOCK, runend_ignore_bits, runend_rank);
	uint64_t runend_block_offset = bitselectv(get_block(qf,
																											runend_block_index)->runends[0],
																						runend_ignore_bits, runend_rank);
	if (runend_block_offset == SLOTS_PER_BLOCK) {
		if (bucket_blocks_offset == 0 && bucket_intrablock_rank == 0) {
			/* The block begins in empty space, and this bucket is in that region of empty space */
			return hash_bucket_index;
		} else {
			do {
				// runend_rank        -= popcntv(get_block(qf, runend_block_index)->runends, METADATA_WORDS_PER_BLOCK, runend_ignore_bits);
				runend_rank        -= popcntv(get_block(qf,
																								runend_block_index)->runends[0],
																			runend_ignore_bits);
				runend_block_index++;
				runend_ignore_bits  = 0;
				// runend_block_offset = bitselectv(get_block(qf, runend_block_index)->runends, METADATA_WORDS_PER_BLOCK, runend_ignore_bits, runend_rank);
				runend_block_offset = bitselectv(get_block(qf,
																									 runend_block_index)->runends[0],
																				 runend_ignore_bits, runend_rank);
			} while (runend_block_offset == SLOTS_PER_BLOCK);
		}
	}

	uint64_t runend_index = SLOTS_PER_BLOCK * runend_block_index + runend_block_offset;
	if (runend_index < hash_bucket_index)
		return hash_bucket_index;
	else
		return runend_index;
}

static inline int offset_lower_bound(const QF *qf, uint64_t slot_index)
{
	const qfblock * b = get_block(qf, slot_index / SLOTS_PER_BLOCK);
	const uint64_t slot_offset = slot_index % SLOTS_PER_BLOCK;
	const uint64_t boffset = b->offset;
	const uint64_t occupieds = b->occupieds[0] & BITMASK(slot_offset+1);
	assert(SLOTS_PER_BLOCK == 64);
	if (boffset <= slot_offset) {
		const uint64_t runends = (b->runends[0] & BITMASK(slot_offset)) >> boffset;
		return popcnt(occupieds) - popcnt(runends);
	}
	return boffset - slot_offset + popcnt(occupieds);
}

static inline int is_empty(const QF *qf, uint64_t slot_index)
{
	return offset_lower_bound(qf, slot_index) == 0;
}

/* static inline int is_empty(const QF *qf, uint64_t slot_index) */
/* { */
/* 	return get_slot(qf, slot_index) == 0 */
/* 		&& !is_occupied(qf, slot_index)  */
/* 		&& !is_runend(qf, slot_index) */
/* 		&& run_end(qf, slot_index) == slot_index; */
/* } */

static inline int might_be_empty(const QF *qf, uint64_t slot_index)
{
	return !is_occupied(qf, slot_index)
		&& !is_runend(qf, slot_index);
}

static inline int probably_is_empty(const QF *qf, uint64_t slot_index)
{
	return get_slot(qf, slot_index) == 0
		&& !is_occupied(qf, slot_index)
		&& !is_runend(qf, slot_index);
}

static inline uint64_t find_first_empty_slot(QF *qf, uint64_t from)
{
	do {
		/* while (!probably_is_empty(qf, from)) */
		/*   from++; */
		int t = offset_lower_bound(qf, from);
		assert(t>=0);
		if (t == 0)
			break;
		from = from + t;
	} while(1);
	return from;
}

static inline uint64_t shift_into_b(const uint64_t a, const uint64_t b,
																		const int bstart, const int bend,
																		const int amount)
{
	const uint64_t a_component = bstart == 0 ? (a >> (64 - amount)) : 0;
	const uint64_t b_shifted_mask = BITMASK(bend - bstart) << bstart;
	const uint64_t b_shifted = ((b_shifted_mask & b) << amount) & b_shifted_mask;
	const uint64_t b_mask = ~b_shifted_mask;
	return a_component | b_shifted | (b & b_mask);
}

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64

static inline void shift_remainders(QF *qf, uint64_t start_index, uint64_t
																		empty_index)
{
	uint64_t start_block  = start_index / SLOTS_PER_BLOCK;
	uint64_t start_offset = start_index % SLOTS_PER_BLOCK;
	uint64_t empty_block  = empty_index / SLOTS_PER_BLOCK;
	uint64_t empty_offset = empty_index % SLOTS_PER_BLOCK;

	assert (start_index <= empty_index && empty_index < qf->metadata->xnslots);

	while (start_block < empty_block) {
		memmove(&get_block(qf, empty_block)->slots[1], 
						&get_block(qf, empty_block)->slots[0],
						empty_offset * sizeof(qf->blocks[0].slots[0]));
		get_block(qf, empty_block)->slots[0] = get_block(qf,
																										 empty_block-1)->slots[SLOTS_PER_BLOCK-1];
		empty_block--;
		empty_offset = SLOTS_PER_BLOCK-1;
	}

	memmove(&get_block(qf, empty_block)->slots[start_offset+1], 
					&get_block(qf, empty_block)->slots[start_offset],
					(empty_offset - start_offset) * sizeof(qf->blocks[0].slots[0]));
}

#else

#define REMAINDER_WORD(qf, i) ((uint64_t *)&(get_block(qf, (i)/qf->metadata->bits_per_slot)->slots[8 * ((i) % qf->metadata->bits_per_slot)]))

static inline void shift_remainders(QF *qf, const uint64_t start_index, const
																		uint64_t empty_index)
{
	uint64_t last_word = (empty_index + 1) * qf->metadata->bits_per_slot / 64;
	const uint64_t first_word = start_index * qf->metadata->bits_per_slot / 64;
	int bend = ((empty_index + 1) * qf->metadata->bits_per_slot) % 64;
	const int bstart = (start_index * qf->metadata->bits_per_slot) % 64;

        assert(start_index <= empty_index && empty_index < qf->metadata->xnslots);
        assert(first_word <= last_word);
	while (last_word != first_word) {
		*REMAINDER_WORD(qf, last_word) = shift_into_b(*REMAINDER_WORD(qf, last_word-1),
																									*REMAINDER_WORD(qf, last_word),
																									0, bend, qf->metadata->bits_per_slot);
		last_word--;
		bend = 64;
	}
	*REMAINDER_WORD(qf, last_word) = shift_into_b(0, *REMAINDER_WORD(qf, last_word),
																								bstart, bend, qf->metadata->bits_per_slot);
}

#endif

static inline void qf_dump_block(const QF *qf, uint64_t i)
{
	uint64_t j;

	printf("%-192d", get_block(qf, i)->offset);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02llx ", j);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->occupieds[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->runends[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#elif BITS_PER_SLOT == 64
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02lx ", get_block(qf, i)->slots[j]);
#else
	for (j = 0; j < SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#endif

	printf("\n");

	printf("\n");
}

void qf_dump(const QF *qf)
{
	uint64_t i;

	printf("%llu %llu %llu\n",
				 qf->metadata->nblocks,
				 qf->metadata->ndistinct_elts,
				 qf->metadata->nelts);

	for (i = 0; i < qf->metadata->nblocks; i++) {
		qf_dump_block(qf, i);
	}

}

static inline void find_next_n_empty_slots(QF *qf, uint64_t from, uint64_t n,
																					 uint64_t *indices)
{
	while (n) {
		indices[--n] = find_first_empty_slot(qf, from);
		from = indices[n] + 1;
	}
}

static inline void shift_slots(QF *qf, int64_t first, uint64_t last, uint64_t
															 distance)
{
	int64_t i;
	if (distance == 1)
		shift_remainders(qf, first, last+1);
	else
		for (i = last; i >= first; i--)
			set_slot(qf, i + distance, get_slot(qf, i));
}

static inline void shift_runends(QF *qf, int64_t first, uint64_t last,
																 uint64_t distance)
{
	assert(last < qf->metadata->xnslots && distance < 64);
	uint64_t first_word = first / 64;
	uint64_t bstart = first % 64;
	uint64_t last_word = (last + distance + 1) / 64;
	uint64_t bend = (last + distance + 1) % 64;

	if (last_word != first_word) {
		METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)), 
																														METADATA_WORD(qf, runends, 64*last_word),
																														0, bend, distance);
		bend = 64;
		last_word--;
		while (last_word != first_word) {
			METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)), 
																															METADATA_WORD(qf, runends, 64*last_word),
																															0, bend, distance);
			last_word--;
		}
	}
	METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(0, METADATA_WORD(qf,
																																					 runends,
																																					 64*last_word),
																													bstart, bend, distance);

}

static inline void insert_replace_slots_and_shift_remainders_and_runends_and_offsets(QF		*qf, 
																																										 int		 operation, 
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders, 
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 noverwrites)
{
	uint64_t empties[67];
	uint64_t i, j;
	uint64_t ninserts = total_remainders - noverwrites;
	uint64_t insert_index = overwrite_index + noverwrites;

	if (ninserts > 0) {
		/* First, shift things to create n empty spaces where we need them. */
		find_next_n_empty_slots(qf, insert_index, ninserts, empties);

		for (i = 0; i < ninserts - 1; i++)
			shift_slots(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		if (ninserts > 0)
			shift_slots(qf, insert_index, empties[ninserts - 1] - 1, ninserts);

		for (i = 0; i < ninserts - 1; i++)
			shift_runends(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		if (ninserts > 0)
			shift_runends(qf, insert_index, empties[ninserts - 1] - 1, ninserts);


		for (i = noverwrites; i < total_remainders - 1; i++)
			METADATA_WORD(qf, runends, overwrite_index + i) &= ~(1ULL <<
																													 (((overwrite_index
																															+ i) %
																														 SLOTS_PER_BLOCK)
																														% 64));

		switch (operation) {
			case 0: /* insert into empty bucket */
				assert (noverwrites == 0);
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										SLOTS_PER_BLOCK) % 64);
				break;
			case 1: /* append to bucket */
				METADATA_WORD(qf, runends, overwrite_index + noverwrites - 1)      &=
					~(1ULL << (((overwrite_index + noverwrites - 1) % SLOTS_PER_BLOCK) %
										 64));
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										SLOTS_PER_BLOCK) % 64);
				break;
			case 2: /* insert into bucket */
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) &=
					~(1ULL << (((overwrite_index + total_remainders - 1) %
											SLOTS_PER_BLOCK) % 64));
				break;
			default:
				fprintf(stderr, "Invalid operation %d\n", operation);
				abort();
		}

		if (ninserts > 0) {
			for (i = bucket_index / SLOTS_PER_BLOCK + 1; i <= empties[ninserts -
					 1]/SLOTS_PER_BLOCK; i++) {
				if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
					get_block(qf, i)->offset += ninserts;
			}
			for (j = 0; j < ninserts -1 ; j++)
				for (i = empties[ninserts - j - 1] / SLOTS_PER_BLOCK + 1; i <=
						 empties[ninserts - j - 2]/SLOTS_PER_BLOCK; i++) {
					if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
						get_block(qf, i)->offset += ninserts - j - 1;
				}
		}
	}
	for (i = 0; i < total_remainders; i++)
		set_slot(qf, overwrite_index + i, remainders[i]);

	modify_metadata(qf, &qf->metadata->noccupied_slots, ninserts);
}

static inline void remove_replace_slots_and_shift_remainders_and_runends_and_offsets(QF	*qf, 
                                                                                     int		 operation, 
                                                                                     uint64_t		 bucket_index,
                                                                                     uint64_t		 overwrite_index,
                                                                                     const uint64_t	*remainders, 
                                                                                     uint64_t		 total_remainders,
                                                                                     uint64_t		 nremovals)
{
	uint64_t i;

	for (i = 0; i < total_remainders; i++)
		set_slot(qf, overwrite_index + i, remainders[i]);
	if (is_runend(qf, overwrite_index + total_remainders - 1) != 
			is_runend(qf, overwrite_index + total_remainders + nremovals - 1))
		METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) ^= 1ULL
			<< ((overwrite_index + total_remainders - 1) % 64);

	uint64_t current_bucket = bucket_index;
	uint64_t current_slot = overwrite_index + total_remainders;
	uint64_t current_distance = nremovals;

	/* FIXME: update block offsets */

	while (current_distance > 0) {
		if (is_runend(qf, current_slot + current_distance - 1))
			do {
				current_bucket++;
			} while (current_bucket < current_slot + current_distance &&
							 !is_occupied(qf, current_bucket));

		if (current_bucket <= current_slot) {
			set_slot(qf, current_slot, get_slot(qf, current_slot + current_distance));
			if (is_runend(qf, current_slot) != 
					is_runend(qf, current_slot + current_distance))
				METADATA_WORD(qf, runends, current_slot) ^= 1ULL << (current_slot % 64);
			current_slot++;

		} else if (current_bucket <= current_slot + current_distance) {
			uint64_t i;
			for (i = current_slot; i < current_bucket; i++) {
				set_slot(qf, i, 0);
				METADATA_WORD(qf, runends, current_slot) &= ~(1ULL << (current_slot % 64));
			}

			current_distance = current_slot + current_distance - current_bucket;
			current_slot = current_bucket;

		} else {
			current_distance = 0;
		}

	}

	modify_metadata(qf, &qf->metadata->noccupied_slots, -nremovals);
}


/*****************************************************************************
 * Code that uses the above to implement a QF with keys and inline counters. *
 *****************************************************************************/

/* 
	 Counter format:
	 0 xs:    <empty string>
	 1 x:     x
	 2 xs:    xx
	 3 0s:    000
	 >2 xs:   xbc...cx  for x != 0, b < x, c != 0, x
	 >3 0s:   0c...c00  for c != 0
	 */
static inline uint64_t *encode_counter(QF *qf, uint64_t remainder, uint64_t
																			 counter, uint64_t *slots)
{
	uint64_t digit = remainder;
	uint64_t base = (1ULL << qf->metadata->bits_per_slot) - 1;
	uint64_t *p = slots;

	if (counter == 0)
		return p;

	*--p = remainder;

	if (counter == 1)
		return p;

	if (counter == 2) {
		*--p = remainder;
		return p;
	}

	if (counter == 3 && remainder == 0) {
		*--p = remainder;
		*--p = remainder;
		return p;    
	}

	if (counter == 3 && remainder > 0) {
		*--p = 0;
		*--p = remainder;
		return p;
	}

	if (remainder == 0)
		*--p = remainder;
	else 
		base--;

	if (remainder)
		counter -= 3;
	else
		counter -= 4;
	do {
		digit = counter % base;
		digit++; /* Zero not allowed */
		if (remainder && digit >= remainder)
			digit++; /* Cannot overflow since digit is mod 2^r-2 */
		*--p = digit;
		counter /= base;
	} while (counter);

	if (remainder && digit >= remainder)
		*--p = 0;

	*--p = remainder;

	return p;
}

/* Returns the length of the encoding. 
REQUIRES: index points to first slot of a counter. */
static inline uint64_t decode_counter(const QF *qf, uint64_t index, uint64_t
																			*remainder, uint64_t *count)
{
	uint64_t base;
	uint64_t rem;
	uint64_t cnt;
	uint64_t digit;
	uint64_t end;

	*remainder = rem = get_slot(qf, index);

	if (is_runend(qf, index)) { /* Entire run is "0" */
		*count = 1; 
		return index;
	}

	digit = get_slot(qf, index + 1);

	if (is_runend(qf, index + 1)) {
		*count = digit == rem ? 2 : 1;
		return index + (digit == rem ? 1 : 0);
	}

	if (rem > 0 && digit >= rem) {
		*count = digit == rem ? 2 : 1;
		return index + (digit == rem ? 1 : 0);
	}

	if (rem > 0 && digit == 0 && get_slot(qf, index + 2) == rem) {
		*count = 3;
		return index + 2;
	}

	if (rem == 0 && digit == 0) {
		if (get_slot(qf, index + 2) == 0) {
			*count = 3;
			return index + 2;
		} else {
			*count = 2;
			return index + 1;
		}
	}

	cnt = 0;
	base = (1ULL << qf->metadata->bits_per_slot) - (rem ? 2 : 1);

	end = index + 1;
	while (digit != rem && !is_runend(qf, end)) {
		if (digit > rem)
			digit--;
		if (digit && rem)
			digit--;
		cnt = cnt * base + digit;

		end++;
		digit = get_slot(qf, end);
	}

	if (rem) {
		*count = cnt + 3;
		return end;
	}

	if (is_runend(qf, end) || get_slot(qf, end + 1) != 0) {
		*count = 1;
		return index;
	}

	*count = cnt + 4;
	return end + 1;
}

/* return the next slot which corresponds to a 
 * different element 
 * */
static inline uint64_t next_slot(QF *qf, uint64_t current) 
{
	uint64_t rem = get_slot(qf, current);
	current++;

	while (get_slot(qf, current) == rem && current <= qf->metadata->nslots) {
		current++;
	}
	return current;
}

static inline bool insert1(QF *qf, uint64_t hash, bool lock, bool spin)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;

	if (lock) {
#ifdef LOG_WAIT_TIME
		if (!qf_spin_lock(qf, &qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK],
											hash_bucket_index/NUM_SLOTS_TO_LOCK, spin))
			return false;
		if (hash_bucket_lock_offset > (NUM_SLOTS_TO_LOCK - CLUSTER_SIZE)) {
			if (!qf_spin_lock(qf, &qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												hash_bucket_index/NUM_SLOTS_TO_LOCK+1, spin)) {
				qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#else
		if (!qf_spin_lock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], spin))
			return false;
		if (hash_bucket_lock_offset > (NUM_SLOTS_TO_LOCK - CLUSTER_SIZE)) {
			if (!qf_spin_lock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												spin)) {
				qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#endif
	}
	if (is_empty(qf, hash_bucket_index) /* might_be_empty(qf, hash_bucket_index) && runend_index == hash_bucket_index */) {
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);

		modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
		modify_metadata(qf, &qf->metadata->noccupied_slots, 1);
		modify_metadata(qf, &qf->metadata->nelts, 1);
	} else {
		uint64_t runend_index              = run_end(qf, hash_bucket_index);
		int operation = 0; /* Insert into empty bucket */
		uint64_t insert_index = runend_index + 1;
		uint64_t new_value = hash_remainder;

		/* printf("RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

		uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																	 hash_bucket_index
																																	 - 1) + 1;

		if (is_occupied(qf, hash_bucket_index)) {

			/* Find the counter for this remainder if it exists. */
			uint64_t current_remainder = get_slot(qf, runstart_index);
			uint64_t zero_terminator = runstart_index;

			/* The counter for 0 is special. */
			if (current_remainder == 0) {
				uint64_t t = runstart_index + 1;
				while (t < runend_index && get_slot(qf, t) != 0)
					t++;
				if (t < runend_index && get_slot(qf, t+1) == 0)
					zero_terminator = t+1; /* Three or more 0s */
				else if (runstart_index < runend_index && get_slot(qf, runstart_index
																													 + 1) == 0)
					zero_terminator = runstart_index + 1; /* Exactly two 0s */
				/* Otherwise, exactly one 0 (i.e. zero_terminator == runstart_index) */

				/* May read past end of run, but that's OK because loop below
					 can handle that */
				if (hash_remainder != 0) {
					runstart_index = zero_terminator + 1;
					current_remainder = get_slot(qf, runstart_index);
				}
			}

			/* Skip over counters for other remainders. */
			while (current_remainder < hash_remainder && runstart_index <=
						 runend_index) {
				/* If this remainder has an extended counter, skip over it. */
				if (runstart_index < runend_index && 
						get_slot(qf, runstart_index + 1) < current_remainder) {
					runstart_index = runstart_index + 2;
					while (runstart_index < runend_index &&
								 get_slot(qf, runstart_index) != current_remainder)
						runstart_index++;
					runstart_index++;

					/* This remainder has a simple counter. */
				} else {
					runstart_index++;
				}

				/* This may read past the end of the run, but the while loop
					 condition will prevent us from using the invalid result in
					 that case. */
				current_remainder = get_slot(qf, runstart_index);
			}

			/* If this is the first time we've inserted the new remainder,
				 and it is larger than any remainder in the run. */
			if (runstart_index > runend_index) {
				operation = 1;
				insert_index = runstart_index;
				new_value = hash_remainder;
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);

				/* This is the first time we're inserting this remainder, but
					 there are larger remainders already in the run. */
			} else if (current_remainder != hash_remainder) {
				operation = 2; /* Inserting */
				insert_index = runstart_index;
				new_value = hash_remainder;
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);

				/* Cases below here: we're incrementing the (simple or
					 extended) counter for this remainder. */

				/* If there's exactly one instance of this remainder. */
			} else if (runstart_index == runend_index || 
								 (hash_remainder > 0 && get_slot(qf, runstart_index + 1) >
									hash_remainder) ||
								 (hash_remainder == 0 && zero_terminator == runstart_index)) {
				operation = 2; /* Insert */
				insert_index = runstart_index;
				new_value = hash_remainder;

				/* If there are exactly two instances of this remainder. */
			} else if ((hash_remainder > 0 && get_slot(qf, runstart_index + 1) ==
									hash_remainder) ||
								 (hash_remainder == 0 && zero_terminator == runstart_index + 1)) {
				operation = 2; /* Insert */
				insert_index = runstart_index + 1;
				new_value = 0;

				/* Special case for three 0s */
			} else if (hash_remainder == 0 && zero_terminator == runstart_index + 2) {
				operation = 2; /* Insert */
				insert_index = runstart_index + 1;
				new_value = 1;

				/* There is an extended counter for this remainder. */
			} else {

				/* Move to the LSD of the counter. */
				insert_index = runstart_index + 1;
				while (get_slot(qf, insert_index+1) != hash_remainder)
					insert_index++;

				/* Increment the counter. */
				uint64_t digit, carry;
				do {
					carry = 0;
					digit = get_slot(qf, insert_index);
					// Convert a leading 0 (which is special) to a normal encoded digit
					if (digit == 0) {
						digit++;
						if (digit == current_remainder)
							digit++;
					}

					// Increment the digit
					digit = (digit + 1) & BITMASK(qf->metadata->bits_per_slot);

					// Ensure digit meets our encoding requirements
					if (digit == 0) {
						digit++;
						carry = 1;
					}
					if (digit == current_remainder)
						digit = (digit + 1) & BITMASK(qf->metadata->bits_per_slot);
					if (digit == 0) {
						digit++;
						carry = 1;
					}

					set_slot(qf, insert_index, digit);
					insert_index--;
				} while(insert_index > runstart_index && carry);

				/* If the counter needs to be expanded. */
				if (insert_index == runstart_index && (carry > 0 || (current_remainder
																														 != 0 && digit >=
																														 current_remainder)))
				{
					operation = 2; /* insert */
					insert_index = runstart_index + 1;
					if (!carry)						/* To prepend a 0 before the counter if the MSD is greater than the rem */
						new_value = 0;
					else if (carry) {			/* Increment the new value because we don't use 0 to encode counters */
						new_value = 2;
						/* If the rem is greater than or equal to the new_value then fail*/
						assert(!current_remainder || new_value < current_remainder);
					}
				} else {
					operation = -1;
				}
			}
		}

		if (operation >= 0) {
			uint64_t empty_slot_index = find_first_empty_slot(qf, runend_index+1);

			assert(insert_index <= empty_slot_index);
			shift_remainders(qf, insert_index, empty_slot_index);

			set_slot(qf, insert_index, new_value);

			shift_runends(qf, insert_index, empty_slot_index-1, 1);
			switch (operation) {
				case 0:
					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index
																																	%
																																	SLOTS_PER_BLOCK)
																																 % 64);
					break;
				case 1:
					METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL <<
																													(((insert_index-1) %
																														SLOTS_PER_BLOCK) %
																													 64));
					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index
																																	%
																																	SLOTS_PER_BLOCK)
																																 % 64);
					break;
				case 2:
					METADATA_WORD(qf, runends, insert_index)   &= ~(1ULL <<
																													((insert_index %
																														SLOTS_PER_BLOCK) %
																													 64));
					break;
				default:
					fprintf(stderr, "Invalid operation %d\n", operation);
					abort();
			}
			/* 
			 * Increment the offset for each block between the hash bucket index
			 * and block of the empty slot  
			 * */
			uint64_t i;
			for (i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
					 empty_slot_index/SLOTS_PER_BLOCK; i++) {
				if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
					get_block(qf, i)->offset++;
				//assert(get_block(qf, i)->offset != 0);
			}
			modify_metadata(qf, &qf->metadata->noccupied_slots, 1);
		}
		modify_metadata(qf, &qf->metadata->nelts, 1);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
	}

	if (lock) {
		qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
		if (hash_bucket_lock_offset > (NUM_SLOTS_TO_LOCK - CLUSTER_SIZE))
			qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
	}

	return true;
}

static inline bool insert(QF *qf,uint64_t hash, uint64_t count, bool lock,
                          bool spin)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t runend_index             = run_end(qf, hash_bucket_index);
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;

	if (lock) {
		/* take the lock for two lock-blocks; the lock-block in which the
		 * hash_bucket_index falls and the next lock-block */

#ifdef LOG_WAIT_TIME
		if (!qf_spin_lock(qf, &qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK],
											hash_bucket_index/NUM_SLOTS_TO_LOCK, spin))
			return false;
		if (hash_bucket_lock_offset > (NUM_SLOTS_TO_LOCK - CLUSTER_SIZE)) {
			if (!qf_spin_lock(qf, &qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												hash_bucket_index/NUM_SLOTS_TO_LOCK+1, spin)) {
				qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#else
		if (!qf_spin_lock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], spin))
			return false;
		if (hash_bucket_lock_offset > (NUM_SLOTS_TO_LOCK - CLUSTER_SIZE)) {
			if (!qf_spin_lock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1], spin)) {
				qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#endif
	}
	/* Empty slot */
	if (might_be_empty(qf, hash_bucket_index) && runend_index ==
			hash_bucket_index) {
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		
		modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
		modify_metadata(qf, &qf->metadata->noccupied_slots, 1);
		modify_metadata(qf, &qf->metadata->nelts, 1);
		/* This trick will, I hope, keep the fast case fast. */
		if (count > 1) {
			qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
			if (hash_bucket_lock_offset > (NUM_SLOTS_TO_LOCK - CLUSTER_SIZE))
				qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
			insert(qf, hash, count - 1, lock, spin);
		}
	} else { /* Non-empty slot */
		uint64_t new_values[67];
		int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
                                                                      hash_bucket_index
                                                                      - 1) + 1;

		if (!is_occupied(qf, hash_bucket_index)) { /* Empty bucket, but its slot is occupied. */
			uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
			insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
                                                                              0, 
                                                                              hash_bucket_index, 
                                                                              runstart_index, 
                                                                              p, 
                                                                              &new_values[67] - p, 
                                                                              0);
			modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
		} else { /* Non-empty bucket */

			uint64_t current_remainder, current_count, current_end;

			/* Find the counter for this remainder, if one exists. */
			current_end = decode_counter(qf, runstart_index, &current_remainder,
																	 &current_count);
			while (current_remainder < hash_remainder && !is_runend(qf, current_end)) {
				runstart_index = current_end + 1;
				current_end = decode_counter(qf, runstart_index, &current_remainder,
																		 &current_count);	
			}

			/* If we reached the end of the run w/o finding a counter for this remainder,
				 then append a counter for this remainder to the run. */
			if (current_remainder < hash_remainder) {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																					1, /* Append to bucket */
																																					hash_bucket_index, 
																																					current_end + 1, 
																																					p, 
																																					&new_values[67] - p, 
																																					0);
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
				/* Found a counter for this remainder.  Add in the new count. */
			} else if (current_remainder == hash_remainder) {
				uint64_t *p = encode_counter(qf, hash_remainder, current_count + count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																					is_runend(qf, current_end) ? 1 : 2, 
																																					hash_bucket_index, 
																																					runstart_index, 
																																					p, 
																																					&new_values[67] - p, 
																																					current_end - runstart_index + 1);
				/* No counter for this remainder, but there are larger
					 remainders, so we're not appending to the bucket. */
			} else {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																					2, /* Insert to bucket */
																																					hash_bucket_index, 
																																					runstart_index, 
																																					p, 
																																					&new_values[67] - p, 
																																					0);	
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
			}
		}
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);
		
		modify_metadata(qf, &qf->metadata->nelts, count);
	}

	if (lock) {
		qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
		if (hash_bucket_lock_offset > (NUM_SLOTS_TO_LOCK - CLUSTER_SIZE))
			qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
	}

	return true;
}

/* inline static void _qf_remove(QF *qf, __uint128_t hash, uint64_t count) */
/* { */
/* 	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot); */
/* 	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot; */
/* 	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK; */
/* 	uint64_t current_remainder, current_count, current_end; */
/* 	uint64_t new_values[67]; */

/* 	/\* Empty bucket *\/ */
/* 	if (is_occupied(qf, hash_bucket_index)) */
/* 		return; */

/* 	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index - 1) + 1; */

/* 	/\* Find the counter for this remainder, if one exists. *\/ */
/* 	current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count); */
/* 	while (current_remainder < hash_remainder && !is_runend(qf, current_end)) { */
/* 		runstart_index = current_end + 1; */
/* 		current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);	 */
/* 	} */
/* 	if (current_remainder != hash_remainder) */
/* 		return; */

/* 	uint64_t *p = encode_counter(qf, hash_remainder,  */
/* 															 count > current_count ? 0 : current_count - count,  */
/* 															 &new_values[67]); */
/* } */

/***********************************************************************
 * Code that uses the above to implement key-value-counter operations. *
 ***********************************************************************/

void qf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t value_bits,
             uint32_t seed)
{
	uint64_t num_slots, xnslots, nblocks;
	uint64_t key_remainder_bits, bits_per_slot;
	uint64_t size;

	assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */
	num_slots = nslots;
	assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */
	xnslots = nslots + 10*sqrt((double)nslots);
	nblocks = (xnslots + SLOTS_PER_BLOCK - 1) / SLOTS_PER_BLOCK;
	key_remainder_bits = key_bits;
	while (nslots > 1) {
		assert(key_remainder_bits > 0);
		key_remainder_bits--;
		nslots >>= 1;
	}

	fprintf(stderr, "%llu %llu %llu\n", key_bits, key_remainder_bits, value_bits);
	bits_per_slot = key_remainder_bits + value_bits;
	assert (BITS_PER_SLOT == 0 || BITS_PER_SLOT == bits_per_slot);
	assert(bits_per_slot > 1);
#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
	size = nblocks * sizeof(qfblock);
#else
	size = nblocks * (sizeof(qfblock) + SLOTS_PER_BLOCK * bits_per_slot / 8);
#endif

	qf->mem = (qfmem *)calloc(sizeof(qfmem), 1);
	
	if (true) {
		qf->metadata = (qfmetadata *)calloc(sizeof(qfmetadata), 1);

		qf->metadata->size = size;
		qf->metadata->seed = seed;
		qf->metadata->nslots = num_slots;
		qf->metadata->xnslots = qf->metadata->nslots +
			10*sqrt((double)qf->metadata->nslots);
		qf->metadata->key_bits = key_bits;
		qf->metadata->value_bits = value_bits;
		qf->metadata->key_remainder_bits = key_remainder_bits;
		qf->metadata->bits_per_slot = bits_per_slot;

		qf->metadata->range = qf->metadata->nslots;
		fprintf(stderr, "nslots: %llu\n", qf->metadata->range);
		qf->metadata->range <<= qf->metadata->bits_per_slot;
        fprintf(stderr, "bits per slot: %llu range: %016llx\n", qf->metadata->bits_per_slot, (uint64_t)qf->metadata->range);
		qf->metadata->nblocks = (qf->metadata->xnslots + SLOTS_PER_BLOCK - 1) /
			SLOTS_PER_BLOCK;
		qf->metadata->nelts = 0;
		qf->metadata->ndistinct_elts = 0;
		qf->metadata->noccupied_slots = 0;
		qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;

		qf->blocks = (qfblock *)calloc(size, 1);

	}
	
	/* initialize all the locks to 0 */
	qf->mem->metadata_lock = 0;
	qf->mem->locks = (volatile int *)calloc(qf->metadata->num_locks,
																					sizeof(volatile int));
#ifdef LOG_WAIT_TIME
	qf->wait_times = (wait_time_data* )calloc(qf->metadata->num_locks+1,
																						sizeof(wait_time_data));
#endif
}

/* The caller should call qf_init on the dest QF before calling this function. 
 */
void qf_copy(QF *dest, QF *src)
{
	memcpy(dest->mem, src->mem, sizeof(qfmem));
	memcpy(dest->metadata, src->metadata, sizeof(qfmetadata));
	memcpy(dest->blocks, src->blocks, src->metadata->size);
}

/* free up the memory if the QF is in memory.
 * else unmap the mapped memory from pagecache.
 *
 * It does not delete the file on disk for on-disk QF.
 */
void qf_destroy(QF *qf)
{
	assert(qf->blocks != NULL);
    free(qf->mem);
    free(qf->metadata);
    free(qf->blocks);
}

void qf_close(QF *qf)
{
	assert(qf->blocks != NULL);
	munmap(qf->metadata, qf->metadata->size + sizeof(qfmetadata));
	close(qf->mem->fd);
}

/* 
 * Will read the on-disk QF using mmap.
 * Data won't be copied in memory.
 *
 */
void qf_read(QF *qf, const char *path)
{
	struct stat sb;
	int ret;

	qf->mem = (qfmem *)calloc(sizeof(qfmem), 1);
	qf->mem->fd = open(path, O_RDWR, S_IRWXU);
	if (qf->mem->fd < 0) {
		perror("Couldn't open file:\n");
		exit(EXIT_FAILURE);
	}

	ret = fstat (qf->mem->fd, &sb);
	if ( ret < 0) {
		perror ("fstat");
		exit(EXIT_FAILURE);
	}

	if (!S_ISREG (sb.st_mode)) {
		fprintf (stderr, "%s is not a file\n", path);
		exit(EXIT_FAILURE);
	}

	qf->metadata = (qfmetadata *)mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED,
																qf->mem->fd, 0);

	qf->blocks = (qfblock *)(qf->metadata + sizeof(qfmetadata));
}

void qf_reset(QF *qf)
{
	qf->metadata->nelts = 0;
	qf->metadata->ndistinct_elts = 0;
	qf->metadata->noccupied_slots = 0;

#ifdef LOG_WAIT_TIME
	memset(qf->wait_times, 0, (qf->metadata->num_locks+1)*sizeof(wait_time_data));
#endif
#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
	memset(qf->blocks, 0, qf->metadata->nblocks* sizeof(qfblock));
#else
	memset(qf->blocks, 0, qf->metadata->nblocks*(sizeof(qfblock) + SLOTS_PER_BLOCK *
																		 qf->metadata->bits_per_slot / 8));
#endif
}

void qf_serialize(const QF *qf, const char *filename)
{
	FILE *fout;
	fout = fopen(filename, "wb+");
	if (fout == NULL) {
		perror("Error opening file for serializing\n");
		exit(EXIT_FAILURE);
	}

	fwrite(qf->metadata, sizeof(qfmetadata), 1, fout);

	/* we don't serialize the locks */
	fwrite(qf->blocks, qf->metadata->size, 1, fout);

	fclose(fout);
}

void qf_deserialize(QF *qf, const char *filename)
{
	FILE *fin;
	fin = fopen(filename, "rb");
	if (fin == NULL) {
		perror("Error opening file for deserializing\n");
		exit(EXIT_FAILURE);
	}

	qf->mem = (qfmem *)calloc(sizeof(qfmem), 1);
	qf->metadata = (qfmetadata *)calloc(sizeof(qfmetadata), 1);

	fread(qf->metadata, sizeof(qfmetadata), 1, fin);

	/* initlialize the locks in the QF */
	qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;
	qf->mem->metadata_lock = 0;
	/* initialize all the locks to 0 */
	qf->mem->locks = (volatile int *)calloc(qf->metadata->num_locks, sizeof(volatile int));

	qf->blocks = (qfblock *)calloc(qf->metadata->size, 1);
	fread(qf->blocks, qf->metadata->size, 1, fin);

	fclose(fin);
}

bool qf_insert(QF *qf, uint64_t key, uint64_t value, uint64_t count, bool
							 lock, bool spin)
{
	/*uint64_t hash = (key << qf->metadata->value_bits) | (value & BITMASK(qf->metadata->value_bits));*/
	if (count == 1)
		return insert1(qf, key, lock, spin);
	else
		return insert(qf, key, count, lock, spin);
}

uint64_t qf_count_key_value(const QF *qf, uint64_t key, uint64_t value,
                            bool lock)
{
	uint64_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->bits_per_slot);
	int64_t hash_bucket_index = hash >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
    uint64_t res = 0;

	if (lock) {
		/* take the lock for two lock-blocks; the lock-block in which the
		 * hash_bucket_index falls and the next lock-block */

		if (!qf_spin_lock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], lock))
			return 0;
		if (hash_bucket_lock_offset > (NUM_SLOTS_TO_LOCK - CLUSTER_SIZE)) {
			if (!qf_spin_lock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1], lock)) {
				qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return 0;
			}
		}
	}
    
	if (!is_occupied(qf, hash_bucket_index)) {
        res = 0;
        goto end;
    }

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
                                                                  hash_bucket_index-1)
                             + 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder) {
            res = current_count;
            goto end;
        }
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));

    res = 0;

end:
    if (lock) {
		qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
		if (hash_bucket_lock_offset > (NUM_SLOTS_TO_LOCK - CLUSTER_SIZE))
			qf_spin_unlock(&qf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
	}

    return res;
}

/* initialize the iterator at the run corresponding
 * to the position index
 */
bool qf_iterator(QF *qf, QFi *qfi, uint64_t position)
{
	assert(position < qf->metadata->nslots);
	if (!is_occupied(qf, position)) {
		uint64_t block_index = position;
		uint64_t idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
		if (idx == 64) {
			while(idx == 64 && block_index < qf->metadata->nblocks) {
				block_index++;
				idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
			}
		}
		position = block_index * SLOTS_PER_BLOCK + idx;
	}

	qfi->qf = qf;
	qfi->num_clusters = 0;
	qfi->run = position;
	qfi->current = position == 0 ? 0 : run_end(qfi->qf, position-1) + 1;
	if (qfi->current < position)
		qfi->current = position;

#ifdef LOG_CLUSTER_LENGTH
	qfi->c_info = (cluster_data* )calloc(qf->metadata->nslots/32, sizeof(cluster_data));
	qfi->cur_start_index = position;
	qfi->cur_length = 1;
#endif

	if (qfi->current >= qf->metadata->nslots)
		return false;
	return true;
}

int qfi_get(QFi *qfi, uint64_t *key, uint64_t *value, uint64_t *count)
{
	assert(qfi->current < qfi->qf->metadata->nslots);

	uint64_t current_remainder, current_count;
	decode_counter(qfi->qf, qfi->current, &current_remainder, &current_count);

	*key = (qfi->run << qfi->qf->metadata->bits_per_slot) | current_remainder;
	*value = 0;   // for now we are not using value
	*count = current_count; 
	
	qfi->qf->metadata->ndistinct_elts++;
	qfi->qf->metadata->nelts += current_count;

	/*qfi->current = end_index;*/ 		//get should not change the current index
																		//of the iterator
	return 0;
}

int qfi_next(QFi *qfi)
{
	if (qfi_end(qfi))
		return 1;
	else {
		/* move to the end of the current counter*/
		uint64_t current_remainder, current_count;
		qfi->current = decode_counter(qfi->qf, qfi->current, &current_remainder,
																	&current_count);
		
		if (!is_runend(qfi->qf, qfi->current)) {
			qfi->current++;
#ifdef LOG_CLUSTER_LENGTH
			qfi->cur_length++;
#endif
			if (qfi->current > qfi->qf->metadata->nslots)
				return 1;
			return 0;
		}
		else {
#ifdef LOG_CLUSTER_LENGTH
			/* save to check if the new current is the new cluster. */
			uint64_t old_current = qfi->current;
#endif
			uint64_t block_index = qfi->run / SLOTS_PER_BLOCK;
			uint64_t rank = bitrank(get_block(qfi->qf, block_index)->occupieds[0],
															qfi->run % SLOTS_PER_BLOCK);
			uint64_t next_run = bitselect(get_block(qfi->qf,
																							block_index)->occupieds[0],
																		rank);
			if (next_run == 64) {
				rank = 0;
				while (next_run == 64 && block_index < qfi->qf->metadata->nblocks) {
					block_index++;
					next_run = bitselect(get_block(qfi->qf, block_index)->occupieds[0],
															 rank);
				}
			}
			if (block_index == qfi->qf->metadata->nblocks) {
				/* set the index values to max. */
				qfi->run = qfi->current = qfi->qf->metadata->xnslots;
				return 1;
			}
			qfi->run = block_index * SLOTS_PER_BLOCK + next_run;
			qfi->current++;
			if (qfi->current < qfi->run)
				qfi->current = qfi->run;
#ifdef LOG_CLUSTER_LENGTH
			if (qfi->current > old_current + 1) { /* new cluster. */
				if (qfi->cur_length > 10) {
					qfi->c_info[qfi->num_clusters].start_index = qfi->cur_start_index;
					qfi->c_info[qfi->num_clusters].length = qfi->cur_length;
					qfi->num_clusters++;
				}
				qfi->cur_start_index = qfi->run;
				qfi->cur_length = 1;
			} else {
				qfi->cur_length++;
			}
#endif
			return 0;
		}
	}
}

inline int qfi_end(QFi *qfi)
{
	if (qfi->current >= qfi->qf->metadata->xnslots /*&& is_runend(qfi->qf, qfi->current)*/)
		return 1;
	else
		return 0;
}

/*
 * Merge qfa and qfb into qfc 
 */
/*
 * iterate over both qf (qfa and qfb)
 * simultaneously 
 * for each index i 
 * min(get_value(qfa, ia) < get_value(qfb, ib))  
 * insert(min, ic) 
 * increment either ia or ib, whichever is minimum.  
 */
void qf_merge(QF *qfa, QF *qfb, QF *qfc) 
{
	QFi qfia, qfib;
	qf_iterator(qfa, &qfia, 0);
	qf_iterator(qfb, &qfib, 0);

	uint64_t keya, valuea, counta, keyb, valueb, countb;
	qfi_get(&qfia, &keya, &valuea, &counta);
	qfi_get(&qfib, &keyb, &valueb, &countb);
	do {
		if (keya < keyb) {
			qf_insert(qfc, keya, valuea, counta, true, true);
			qfi_next(&qfia);
			qfi_get(&qfia, &keya, &valuea, &counta);
		}
		else {
			qf_insert(qfc, keyb, valueb, countb, true, true);
			qfi_next(&qfib);
			qfi_get(&qfib, &keyb, &valueb, &countb);
		}
	} while(!qfi_end(&qfia) && !qfi_end(&qfib));

	if (!qfi_end(&qfia)) {
		do {
			qfi_get(&qfia, &keya, &valuea, &counta);
			qf_insert(qfc, keya, valuea, counta, true, true);
		} while(!qfi_next(&qfia));
	}
	if (!qfi_end(&qfib)) {
		do {
			qfi_get(&qfib, &keyb, &valueb, &countb);
			qf_insert(qfc, keyb, valueb, countb, true, true);
		} while(!qfi_next(&qfib));
	}

	return;
}

/*
 * Merge an array of qfs into the resultant QF
 */
void qf_multi_merge(QF *qf_arr[], int nqf, QF *qfr)
{
	int i;
	QFi qfi_arr[nqf];
	int flag = 0;
	int smallest_i = 0;
	uint64_t smallest_key = UINT64_MAX;
	for (i=0; i<nqf; i++) {
		qf_iterator(qf_arr[i], &qfi_arr[i], 0);
	}

	while (!flag) {
		uint64_t keys[nqf];
		uint64_t values[nqf];
		uint64_t counts[nqf];
		for (i=0; i<nqf; i++)
			qfi_get(&qfi_arr[i], &keys[i], &values[i], &counts[i]);
		
		do {
			smallest_key = UINT64_MAX;
			for (i=0; i<nqf; i++) {
				if (keys[i] < smallest_key) {
					smallest_key = keys[i]; smallest_i = i;
				}
			}
			qf_insert(qfr, keys[smallest_i], values[smallest_i], counts[smallest_i],
								true, true);
			qfi_next(&qfi_arr[smallest_i]);
			qfi_get(&qfi_arr[smallest_i], &keys[smallest_i], &values[smallest_i],
							&counts[smallest_i]);
		} while(!qfi_end(&qfi_arr[smallest_i]));

		/* remove the qf that is exhausted from the array */
		if (smallest_i < nqf-1)
			memmove(&qfi_arr[smallest_i], &qfi_arr[smallest_i+1],
							(nqf-smallest_i-1)*sizeof(qfi_arr[0]));
		nqf--;
		if (nqf == 1)
			flag = 1;
	}
	if (!qfi_end(&qfi_arr[0])) {
		do {
			uint64_t key, value, count;
			qfi_get(&qfi_arr[0], &key, &value, &count);
			qf_insert(qfr, key, value, count, true, true);
		} while(!qfi_next(&qfi_arr[0]));
	}

	return;
}

/* find cosine similarity between two QFs. */
uint64_t qf_inner_product(QF *qfa, QF *qfb)
{
	uint64_t acc = 0;
	QFi qfi;
	QF *qf_mem, *qf_disk;

	// create the iterator on the larger QF.
	if (qfa->metadata->size > qfb->metadata->size) {
		qf_mem = qfb;
		qf_disk = qfa;
	} else {
		qf_mem = qfa;
		qf_disk = qfb;
	}

	qf_iterator(qf_disk, &qfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		uint64_t count_mem;
		qfi_get(&qfi, &key, &value, &count);
		if ((count_mem = qf_count_key_value(qf_mem, key, 0, false)) > 0) {
			acc += count*count_mem;
		}
	} while (!qfi_next(&qfi));

	return acc;
}

/* magnitude of a QF. */
uint64_t qf_magnitude(QF *qf)
{
	return sqrt(qf_inner_product(qf, qf));
}


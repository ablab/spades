#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "io_lib/os.h"
#include "io_lib/hash_table.h"
#include "io_lib/jenkins_lookup3.h"

/* =========================================================================
 * TCL's hash function. Basically hash*9 + char.
 * =========================================================================
 */

uint32_t HashTcl(uint8_t *data, int len) {
    uint32_t hash = 0;
    int i;

    for (i = 0; i < len; i++) {
	hash += (hash<<3) + data[i];
    }

    return hash;
}

/* =========================================================================
 * Paul Hsieh's hash function
 * http://www.azillionmonkeys.com/qed/hash.html
 * =========================================================================
 */

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((const uint8_t *)(d))[1] << 8UL)\
                      +((const uint8_t *)(d))[0])
#endif

uint32_t HashHsieh(uint8_t *data, int len) {
    uint32_t hash = 0, tmp;
    int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
    case 3: hash += get16bits (data);
	hash ^= hash << 16;
	hash ^= data[sizeof (uint16_t)] << 18;
	hash += hash >> 11;
	break;
    case 2: hash += get16bits (data);
	hash ^= hash << 11;
	hash += hash >> 17;
	break;
    case 1: hash += *data;
	hash ^= hash << 10;
	hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 2;
    hash += hash >> 15;
    hash ^= hash << 10;

    return hash;
}

/* =========================================================================
 * Bob Jenkins' hash function
 * http://burtleburtle.net/bob/hash/doobs.html
 *
 * See jenkins_lookup3.c for a new version of this that has good hash
 * characteristics for a full 64-bit hash value.
 * =========================================================================
 */

#define hashsize(n) ((uint32_t)1<<(n))
#define hashmask(n) (hashsize(n)-1)

/*
--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a 
  structure that could supported 2x parallelism, like so:
      a -= b; 
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------
*/
#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

/*
--------------------------------------------------------------------
hash() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  len     : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6*len+35 instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (uint8_t **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------
*/

uint32_t HashJenkins(uint8_t *k, int length /*, uint32_t initval */)
{
   register uint32_t a,b,c,len;

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = 0; /* initval; */        /* the previous hash value */

   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a += (k[0] +((uint32_t)k[1]<<8) +((uint32_t)k[2]<<16) +((uint32_t)k[3]<<24));
      b += (k[4] +((uint32_t)k[5]<<8) +((uint32_t)k[6]<<16) +((uint32_t)k[7]<<24));
      c += (k[8] +((uint32_t)k[9]<<8) +((uint32_t)k[10]<<16)+((uint32_t)k[11]<<24));
      mix(a,b,c);
      k += 12; len -= 12;
   }

   /*------------------------------------- handle the last 11 bytes */
   c += length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c+=((uint32_t)k[10]<<24);
   case 10: c+=((uint32_t)k[9]<<16);
   case 9 : c+=((uint32_t)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b+=((uint32_t)k[7]<<24);
   case 7 : b+=((uint32_t)k[6]<<16);
   case 6 : b+=((uint32_t)k[5]<<8);
   case 5 : b+=k[4];
   case 4 : a+=((uint32_t)k[3]<<24);
   case 3 : a+=((uint32_t)k[2]<<16);
   case 2 : a+=((uint32_t)k[1]<<8);
   case 1 : a+=k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}

/*
 * An interface to the above hash functions.
 * Returns:
 *    A 32-bit hash key, suitable for masking down to smaller bit sizes
 */
uint32_t hash(int func, uint8_t *key, int key_len) {
    switch (func) {
    case HASH_FUNC_HSIEH:
	return HashHsieh(key, key_len);

    case HASH_FUNC_TCL:
	return HashTcl(key, key_len);
	
    case HASH_FUNC_JENKINS:
	return HashJenkins(key, key_len);

    case HASH_FUNC_JENKINS3:
      {
	uint32_t pc = 0, pb = 0;
	HashJenkins3(key, key_len, &pc, &pb);
	return pc;
      }
    }
    
    return 0;
}

/*
 * As per hash() above but returns a 64-bit key. For 32-bit hash functions
 * this is simply a duplication of the 32-bit value.
 */
uint64_t hash64(int func, uint8_t *key, int key_len) {
    uint32_t pc = 0, pb = 0;

    switch (func) {
    case HASH_FUNC_HSIEH:
	pb = pc = HashHsieh(key, key_len);
	break;

    case HASH_FUNC_TCL:
	pb = pc = HashTcl(key, key_len);
	break;
	
    case HASH_FUNC_JENKINS:
	pb = pc = HashJenkins(key, key_len);
	break;

    case HASH_FUNC_JENKINS3:
	HashJenkins3(key, key_len, &pc, &pb);
	break;
    }
    
    return pc + (((uint64_t)pb)<<32);
}

/* =========================================================================
 * Hash Table handling code
 * =========================================================================
 */

/* Multiplicative factors indicating when to grow or shrink the hash table */
#define HASH_TABLE_RESIZE 3

/*
 * Creates a HashItem for use with HashTable h.
 *
 * Returns:
 *    A pointer to new HashItem on success
 *    NULL on failure.
 */
static HashItem *HashItemCreate(HashTable *h) {
    HashItem *hi;

    hi = (h->options & HASH_POOL_ITEMS
	  ? pool_alloc(h->hi_pool) : malloc(sizeof(*hi)));
    if (NULL == hi) return NULL;

    hi->data.p    = NULL;
    hi->data.i    = 0;
    hi->next      = NULL;
    hi->key       = NULL;
    hi->key_len   = 0;

    h->nused++;
    
    return hi;
}

/*
 * Deallocates a HashItem created via HashItemCreate.
 *
 * This function will not remove the item from the HashTable so be sure to
 * call HashTableDel() first if appropriate.
 */
static void HashItemDestroy(HashTable *h, HashItem *hi, int deallocate_data) {
    if (!hi) return;

    if (!(h->options & HASH_NONVOLATILE_KEYS) || (h->options & HASH_OWN_KEYS))
	if (hi->key)
	    free(hi->key);

    if (deallocate_data && hi->data.p)
	free(hi->data.p);

    if (h->options & HASH_POOL_ITEMS) {
        pool_free(h->hi_pool, hi);
    } else {
	free(hi);
    }

    h->nused--;
}

/*
 * Creates a new HashTable object. Size will be rounded up to the next
 * power of 2. It is a starting point and hash tables may be grown or shrunk
 * as needed (if HASH_DYNAMIC_SIZE is used).
 *
 * If HASH_POOL_ITEMS is used, HashItems will be allocated in blocks to reduce
 * malloc overhead in the case where a large number of items is required.
 * HashItems allocated this way will be put on a free list when destroyed; the
 * memory will only be reclaimed when the entire hash table is destroyed.
 *
 * Options are as defined in the header file (see HASH_* macros).
 *
 * Returns:
 *    A pointer to a HashTable on success
 *    NULL on failure
 */
HashTable *HashTableCreate(int size, int options) {
    HashTable *h;
    int i, bits;
    uint32_t mask;

    if (!(h = (HashTable *)malloc(sizeof(*h))))
	return NULL;

    if (options & HASH_POOL_ITEMS) {
        h->hi_pool = pool_create(sizeof(HashItem));
	if (NULL == h->hi_pool) {
	    free(h);
	    return NULL;
	}
    } else {
        h->hi_pool = NULL;
    }

    if (size < 4)
	size = 4; /* an inconsequential minimum size */

    /* Round the requested size to the next power of 2 */
    bits = 0;
    size--;
    while (size) {
	size /= 2;
	bits++;
    }
    size = 1<<bits;
    mask = size-1;

    h->nbuckets = size;
    h->mask = mask;
    h->options = options;
    h->nused = 0;
    h->bucket = (HashItem **)malloc(sizeof(*h->bucket) * size);
    if (NULL == h->bucket) {
        HashTableDestroy(h, 0);
        return NULL;
    }

    for (i = 0; i < size; i++) {
	h->bucket[i] = NULL;
    }

    return h;
}

/*
 * Deallocates a HashTable object (created by HashTableCreate).
 *
 * The deallocate_data parameter is a boolean to indicate whether the
 * data attached to the hash table should also be free()d. DO NOT USE
 * this if the HashData attached was not a pointer allocated using
 * malloc().
 */
void HashTableDestroy(HashTable *h, int deallocate_data) {
    int i;

    if (!h)
	return;

    if (h->bucket) {
        for (i = 0; i < h->nbuckets; i++) {
	    HashItem *hi = h->bucket[i], *next = NULL;
	    for (hi = h->bucket[i]; hi; hi = next) {
	        next = hi->next;
		HashItemDestroy(h, hi, deallocate_data);
	    }
	}

	free(h->bucket);
    }

    if (h->hi_pool) pool_destroy(h->hi_pool);

    free(h);
}

/*
 * Resizes a HashTable to have 'newsize' buckets.
 * This is called automatically when adding or removing items so that the
 * hash table keeps at a sensible scale.
 *
 * FIXME: Halving the size of the hash table is simply a matter of coaelescing
 * every other bucket. Instead we currently rehash (which is slower).
 * Doubling the size of the hash table currently requires rehashing, but this
 * too could be optimised by storing the full 32-bit hash of the key along
 * with the key itself. This then means that it's just a matter of seeing what
 * the next significant bit is. It's a memory vs speed tradeoff though and
 * re-hashing is pretty quick.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int HashTableResize(HashTable *h, int newsize) {
    HashTable *h2;
    int i;

    /* fprintf(stderr, "Resizing to %d\n", newsize); */

    /* Create a new hash table and rehash everything into it */
    h2 = HashTableCreate(newsize, h->options);

    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi, *next;
	for (hi = h->bucket[i]; hi; hi = next) {
	    uint64_t hv = hash64(h2->options & HASH_FUNC_MASK,
				 (uint8_t *)hi->key, hi->key_len) & h2->mask;
	    next = hi->next;
	    hi->next = h2->bucket[hv];
	    h2->bucket[hv] = hi;
	}
    }

    /* Swap the links over & free */
    free(h->bucket);
    h->bucket   = h2->bucket;
    h->nbuckets = h2->nbuckets;
    h->mask     = h2->mask;

    if (h2->hi_pool)
	pool_destroy(h2->hi_pool);
    free(h2);

    return 0;
}

/*
 * Adds a HashData item to HashTable h with a specific key. Key can be binary
 * data, but if key_len is passed as zero then strlen() will be used to
 * determine the key length.
 *
 * The "new" pointer may be passed as NULL. When not NULL it is filled out
 * as a boolean to indicate whether the key is already in this hash table.
 *
 * The HASH_ALLOW_DUP_KEYS option (specified when using HashTableCreate)
 * will allow duplicate keys to be stored, and hence *new is also zero.
 * By default duplicate keys are disallowed.
 *
 * Keys are considered to be volatile memory (ie temporary storage) and so the
 * hash table takes separate copies of them. To avoid this use the
 * HASH_NONVOLATILE_KEYS option.
 *
 * If the HASH_OWN_KEYS option was specified when creating the table then
 * keys will be considered to be owned by the hash table. In this case
 * the key will be freed when the table is destroyed regardless of
 * whether the HASH_NONVOLATILE_KEYS option was used to allocate its
 * own private copy.
 *
 * Returns:
 *    The HashItem created (or matching if a duplicate) on success
 *    NULL on failure
 */
HashItem *HashTableAdd(HashTable *h, char *key, int key_len, HashData data,
		       int *new) {
    uint64_t hv;
    HashItem *hi;

    if (!key_len)
	key_len = strlen(key);

    hv = hash64(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len)
	& h->mask;

    /* Already exists? */
    if (!(h->options & HASH_ALLOW_DUP_KEYS)) {
	for (hi = h->bucket[hv]; hi; hi = hi->next) {
	    if (key_len == hi->key_len && key[0] == hi->key[0] &&
		memcmp(key, hi->key, key_len) == 0) {
		if (new) *new = 0;
		return hi;
	    }
	}
    }

    /* No, so create a new one and link it in */
    if (NULL == (hi = HashItemCreate(h)))
	return NULL;

    if (h->options & HASH_NONVOLATILE_KEYS)
	hi->key = key;
    else {
	hi->key = (char *)malloc(key_len+1);
	memcpy(hi->key, key, key_len);
	hi->key[key_len] = 0; /* null terminate incase others print keys */
    }
    hi->key_len = key_len;
    hi->data = data;
    hi->next = h->bucket[hv];
    h->bucket[hv] = hi;

    if ((h->options & HASH_DYNAMIC_SIZE) &&
	h->nused > HASH_TABLE_RESIZE * h->nbuckets)
	HashTableResize(h, h->nbuckets*4);

    if (new) *new = 1;

    return hi;
}


/*
 * Removes a specified HashItem from the HashTable. (To perform this it needs
 * to rehash based on the hash key as hash_item only has a next pointer and
 * not a previous pointer.)
 * 
 * The HashItem itself is also destroyed (by an internal call to
 * HashItemDestroy). The deallocate_data parameter controls whether the data
 * associated with the HashItem should also be free()d.
 *
 * See also the HashTableRemove() function to remove by key instead of
 * HashItem.
 *
 * Returns 0 on success
 *        -1 on failure (eg HashItem not in the HashTable);
 */
int HashTableDel(HashTable *h, HashItem *hi, int deallocate_data) {
    uint64_t hv;
    HashItem *next, *last;

    hv = hash64(h->options & HASH_FUNC_MASK,
		(uint8_t *)hi->key, hi->key_len) & h->mask;

    for (last = NULL, next = h->bucket[hv]; next;
	 last = next, next = next->next) {
	if (next == hi) {
	    /* Link last to next->next */
	    if (last)
		last->next = next->next;
	    else
		h->bucket[hv] = next->next;

	    HashItemDestroy(h, hi, deallocate_data);

	    return 0;
	}
    }

    return -1;
}


/*
 * Searches the HashTable for the data registered with 'key' and removes
 * these items from the HashTable. In essence this is a combination of
 * HashTableSearch and HashTableDel functions.
 *
 * If HASH_ALLOW_DUP_KEYS is used this will remove all items matching 'key',
 * otherwise just a single item will be removed.
 *
 * If 'deallocate_data' is true the data associated with the HashItem will
 * be free()d.
 *
 * Returns
 *    0 on success (at least one item found)
 *   -1 on failure (no items found).
 */
int HashTableRemove(HashTable *h, char *key, int key_len,
		    int deallocate_data) {
    uint64_t hv;
    HashItem *last, *next, *hi;
    int retval = -1;

    if (!key_len)
	key_len = strlen(key);

    hv = hash64(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len)
	& h->mask;

    last = NULL;
    next = h->bucket[hv];

    while (next) {
	hi = next;
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0) {
	    /* An item to remove, adjust links and destroy */
	    if (last)
		last->next = hi->next;
	    else
		h->bucket[hv] = hi->next;

	    next = hi->next;
	    HashItemDestroy(h, hi, deallocate_data);

	    retval = 0;
	    if (!(h->options & HASH_ALLOW_DUP_KEYS))
		break;

	} else {
	    /* We only update last when it's something we haven't destroyed */
	    last = hi;
	    next = hi->next;
	}
    }

    return retval;
}

/*
 * Searches the HashTable for the data registered with 'key'.
 * If HASH_ALLOW_DUP_KEYS is used this will just be the first one found.
 * You will then need to use HashTableNext to iterate through the matches.
 *
 * Returns
 *    HashItem if found
 *    NULL if not found
 */
HashItem *HashTableSearch(HashTable *h, char *key, int key_len) {
    uint64_t hv;
    HashItem *hi;

    if (!key_len)
	key_len = strlen(key);

    hv = hash64(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len)
	& h->mask;
    for (hi = h->bucket[hv]; hi; hi = hi->next) {
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0)
	    return hi;
    }

    return NULL;
}

/*
 * Find the next HashItem (starting from 'hi') to also match this key.
 * This is only valid when the HASH_ALLOW_DUP_KEYS is in use.
 *
 * Returns
 *    HashItem if found
 *    NULL if not found
 */
HashItem *HashTableNext(HashItem *hi, char *key, int key_len) {
    if (!hi)
	return NULL;

    for (hi = hi->next; hi; hi = hi->next) {
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0)
	    return hi;
    }

    return NULL;
}

/*
 * Dumps a textual represenation of the hash table to stdout.
 */
void HashTableDump(HashTable *h, FILE *fp) {
    int i;
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    fprintf(fp, "%.*s\n", hi->key_len, hi->key);
	}
    }
}

/*
 * Produces some simple statistics on the hash table population.
 */
void HashTableStats(HashTable *h, FILE *fp) {
    int i;
    double avg = (double)h->nused / h->nbuckets;
    double var = 0;
    int maxlen = 0;
    int filled = 0;
    int clen[51];

    for (i = 0; i <= 50; i++)
	clen[i] = 0;

    for (i = 0; i < h->nbuckets; i++) {
	int len = 0;
	HashItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    len++;
	}
	if (len > 0) {
	    filled++;
	    if (len > maxlen)
		maxlen = len;
	}
	clen[len <= 50 ? len : 50]++;
	var += (len-avg) * (len-avg);
    }
    var /= h->nbuckets;
    /* sd = sqrt(var); */

    fprintf(fp, "Nbuckets  = %d\n", h->nbuckets);
    fprintf(fp, "Nused     = %d\n", h->nused);
    fprintf(fp, "Avg chain = %f\n", avg);
    fprintf(fp, "Chain var.= %f\n", var);
    fprintf(fp, "%%age full = %f\n", (100.0*filled)/h->nbuckets);
    fprintf(fp, "max len   = %d\n", maxlen);
    for (i = 0; i <= maxlen; i++) {
	fprintf(fp, "Chain %2d   = %d\n", i, clen[i]);
    }
}

/*
 * --------------------------------------------------------------------
 * Below we have a specialisation of the HashTable code where the data
 * attached to the hash table is a position,size pair. This allows for the
 * hash table to encode positions and sizes of items within a file archive.
 * --------------------------------------------------------------------
 */

/*
 * Writes the HashTable structures to 'fp'.
 * This is a specialisation of the HashTable where the HashData is a
 * position,size tuple.
 *
 * This consists of the following format:
 * Header:
 *    ".hsh" (magic numebr)
 *    x4     (1-bytes of version code, eg "1.00")
 *    x1     (HASH_FUNC_? function used)
 *    x1     (number of file headers: FH. These count from 1 to FH inclusive)
 *    x1     (number of file footers: FF. These count from 1 to FF inclusive)
 *    x1     (number of archives indexed: NA)
 *    x4     (4-bytes big-endian; number of hash buckets)
 *    x8     (offset to add to item positions. eg size of this index)
 *    x4     (4-bytes big-endian; number of bytes in hash file, inc. header)
 * Archive name: (NH copies of, or just 1 zero byte if none)
 *    x1     (length, zero => no name, eg when same file as hash index)
 *    ?      (archive filename)
 * File headers (FH copies of):
 *    x1     (archive no.)
 *    x7     (position)
 *    x4     (size)
 * File footers (FH copies of):
 *    x1     (archive no.)
 *    x7     (position)
 *    x4     (size)
 * Buckets (multiples of)
 *    x4     (4-byte offset of linked list pos,  rel. to the start of the hdr)
 * Items (per bucket chain, not written if Bucket[?]==0)
 *    x1     (key length, zero => end of chain)
 *    ?      (key)
 *    x0.5   (File header to use. zero => none) top 4 bits
 *    x0.5   (File footer to use. zero => none) bottom 4 bits
 *    x8     (position)
 *    x4     (size)
 * ... arbitrary gap (but likely none)
 * Index footer:
 *    ".hsh" (magic number)
 *    x8     (offset to Hash Header. +ve = absolute, -ve = relative to end)
 *
 * It is designed such that on-disk querying of the hash table can be done
 * purely by forward seeks. (This is generally faster due to pre-fetching of
 * the subsequent blocks by many disk controllers.)
 *
 * Returns: the number of bytes written on success
 *         -1 for error
 */
uint64_t HashFileSave(HashFile *hf, FILE *fp, int64_t offset) {
    int i;
    HashItem *hi;
    uint32_t *bucket_pos;
    uint64_t hfsize = 0, be_hfsize;
    HashTable *h = hf->h;
    HashFileFooter foot;
   
    /* Compute the coordinates of the hash items */
    hfsize = HHSIZE;				 /* header */
    hfsize += h->nbuckets * 4; 			 /* buckets */
    for (i = 0; i < hf->nheaders; i++)		 /* headers */
	hfsize += 12;
    for (i = 0; i < hf->nfooters; i++)		 /* footers */
	hfsize += 12;
    if (hf->narchives) {
	for (i = 0; i < hf->narchives; i++)
	    hfsize += strlen(hf->archives[i])+1; /* archive filename */
    } else {
	hfsize++;
    }
    bucket_pos = (uint32_t *)calloc(h->nbuckets, sizeof(uint32_t));
    for (i = 0; i < h->nbuckets; i++) {
	bucket_pos[i] = hfsize;

	if (!(hi = h->bucket[i]))
	    continue;
	for (; hi; hi = hi->next) {
	    hfsize += 1 + 1 + hi->key_len + 8 + 4; /* keys, pos, size */
	}
	hfsize++;				/* list-end marker */
    }
    hfsize += sizeof(foot);

    /* Write the header: */
    memcpy(hf->hh.magic, HASHFILE_MAGIC, 4);
    if (hf->narchives > 1)
	memcpy(hf->hh.vers,  HASHFILE_VERSION, 4);
    else
	memcpy(hf->hh.vers,  HASHFILE_VERSION100, 4);
    hf->hh.hfunc     = h->options & HASH_FUNC_MASK;
    hf->hh.nheaders  = hf->nheaders;
    hf->hh.nfooters  = hf->nfooters;
    hf->hh.narchives = hf->narchives == 1 ? 0 : hf->narchives;
    hf->hh.nbuckets  = be_int4(h->nbuckets);
    hf->hh.offset    = offset == HASHFILE_PREPEND
	? be_int8(hfsize) /* archive will be appended to this file */
	: be_int8(offset);
    hf->hh.size     = be_int4(hfsize);
    fwrite(&hf->hh, HHSIZE, 1, fp);

    /* Write the archive filename, if known */
    if (hf->narchives) {
	for (i = 0; i < hf->narchives; i++) {
	    fputc(strlen(hf->archives[i]), fp);
	    fputs(hf->archives[i], fp);
	}
    } else {
	/* Compatibility with v1.00 file format */
	fputc(0, fp);
    }

    /* Write out the headers and footers */
    for (i = 0; i < hf->nheaders; i++) {
	HashFileSection hs;
	hs.pos  = be_int8(hf->headers[i].pos);
	*(char *)&hs.pos = hf->headers[i].archive_no;
	fwrite(&hs.pos, 8, 1, fp);
	hs.size = be_int4(hf->headers[i].size);
	fwrite(&hs.size, 4, 1, fp);
    }

    for (i = 0; i < hf->nfooters; i++) {
	HashFileSection hs;
	hs.pos  = be_int8(hf->footers[i].pos);
	*(char *)&hs.pos = hf->footers[i].archive_no;
	fwrite(&hs.pos, 8, 1, fp);
	hs.size = be_int4(hf->footers[i].size);
	fwrite(&hs.size, 4, 1, fp);
    }

    /* Write out hash buckets */
    for (i = 0; i < h->nbuckets; i++) {
	uint32_t zero = 0;
	uint32_t be32;

	if (!(hi = h->bucket[i])) {
	    fwrite(&zero, 4, 1, fp);
	    continue;
	}

	be32 = be_int4(bucket_pos[i]);
	fwrite(&be32, 4, 1, fp);
    }
    free(bucket_pos);

    /*
     * Write the hash_item linked lists. The first item is the
     * hash key length. We append a zero to the end of the list so we
     * can check this key length to determine the end of this hash
     * item list.
     */
    for (i = 0; i < h->nbuckets; i++) {
	if (!(hi = h->bucket[i]))
	    continue;
	for (; hi; hi = hi->next) {
	    uint64_t be64;
	    uint32_t be32;
	    HashFileItem *hfi = (HashFileItem *)hi->data.p;
	    unsigned char headfoot = 0;

	    fprintf(fp, "%c%.*s", hi->key_len,
		    hi->key_len, hi->key);
	    headfoot = (((hfi->header) & 0xf) << 4) | ((hfi->footer) & 0xf);
	    fwrite(&headfoot, 1, 1, fp);
	    be64 = be_int8(hfi->pos);
	    *(char *)&be64 = hfi->archive;
	    fwrite(&be64, 8, 1, fp);
	    be32 = be_int4(hfi->size);
	    fwrite(&be32, 4, 1, fp);
	}
        fputc(0, fp);
    }

    /* Finally write the footer referencing back to the header start */
    memcpy(foot.magic, HASHFILE_MAGIC, 4);
    be_hfsize = be_int8(-hfsize);
    memcpy(foot.offset, &be_hfsize, 8);
    fwrite(&foot, sizeof(foot), 1, fp);

    return hfsize;
}

#if 0
/*
 * Reads an entire HashTable from fp.
 *
 * Returns:
 *    A filled out HashTable pointer on success
 *    NULL on failure   
 */
HashFile *HashFileLoad_old(FILE *fp) {
    int i;
    HashTable *h;
    HashItem *hi;
    HashFile *hf;
    uint32_t *bucket_pos;
    unsigned char *htable;
    int htable_pos;
    int fnamelen;

    if (NULL == (hf = (HashFile *)calloc(1, sizeof(HashFile))))
	return NULL;
    if (NULL == (htable = (unsigned char *)malloc(HHSIZE)))
	return NULL;

    /* Read and create the hash table header */
    if (HHSIZE != fread(htable, 1, HHSIZE, fp))
	return NULL;
    memcpy(&hf->hh, htable, HHSIZE);
    hf->hh.nbuckets = be_int4(hf->hh.nbuckets);
    hf->hh.offset = be_int8(hf->hh.offset);
    hf->hh.size = be_int4(hf->hh.size);
    hf->h = h = HashTableCreate(hf->hh.nbuckets, hf->hh.hfunc);
    bucket_pos = (uint32_t *)calloc(h->nbuckets, sizeof(uint32_t));

    /* Load the archive filename */
    if (hf->narchives) {
	hf->archives = (char **)malloc(hf->narchives * sizeof(char *));
    } else {
	hf->archives = NULL;
    }

    if (hf->narchives) {
	for (i = 0; i < hf->narchives; i++) {
	    fnamelen = fgetc(fp);
	    hf->archives[i] = malloc(fnamelen+1);
	    fread(hf->archives[i], 1, fnamelen, fp);
	    hf->archives[i][fnamelen] = 0;
	}
    } else {
	/* Consume 0 byte for v1.00 format */
	fgetc(fp);
    }

    /* Load the rest of the hash table to memory */
    htable_pos = HHSIZE + fnamelen + 1;
    if (NULL == (htable = (unsigned char *)realloc(htable, hf->hh.size)))
	return NULL;
    if (hf->hh.size-htable_pos !=
	fread(&htable[htable_pos], 1, hf->hh.size-htable_pos, fp))
	return NULL;

    /* Read the header / footer items */
    for (i = 0; i < hf->hh.nheaders; i++)
	htable_pos += 8; /* skip them for now */
    for (i = 0; i < hf->hh.nfooters; i++)
	htable_pos += 8; /* skip them for now */

    /* Identify the "bucket pos". Detemines which buckets have data */
    for (i = 0; i < h->nbuckets; i++) {
	memcpy(&bucket_pos[i], &htable[htable_pos], 4);
	bucket_pos[i] = be_int4(bucket_pos[i]);
	htable_pos += 4;
    }

    /* Read the hash table items */
    for (i = 0; i < h->nbuckets; i++) {
	if (!bucket_pos[i])
	    continue;
	for (;;) {
	    int c;
	    unsigned char uc;
	    char key[256];
	    uint64_t pos;
	    uint32_t size;
	    HashFileItem *hfi;

	    c = htable[htable_pos++];
	    if (c == EOF || !c)
		break;

	    /* key */
	    memcpy(key, &htable[htable_pos], c);
	    htable_pos += c;

	    /* header/footer */
	    uc = htable[htable_pos++];
	    hfi = (HashFileItem *)malloc(sizeof(*hfi));
	    hfi->header = (uc >> 4) & 0xf;
	    hfi->footer = uc & 0xf;

	    /* pos */
	    memcpy(&pos, &htable[htable_pos], 8);
	    htable_pos += 8;
	    hfi->pos = be_int8(pos) + hf->hh.offset;

	    /* size */
	    memcpy(&size, &htable[htable_pos], 4);
	    htable_pos += 4;
	    hfi->size = be_int4(size);

	    /* Add to the hash table */
	    hi = HashItemCreate(h);
	    hi->next = h->bucket[i];
	    h->bucket[i] = hi;
	    hi->key_len = c;
	    hi->key = (char *)malloc(c+1);
	    memcpy(hi->key, key, c);
	    hi->key[c] = 0; /* For debugging convenience only */
	    hi->data.p = hfi;
	}
    }

    fprintf(stderr, "done\n");
    fflush(stderr);
    free(bucket_pos);

    return hf;
}
#endif

/*
 * Opens a stored hash table file. It also internally keeps an open file to
 * hash and the archive files.
 *
 * Returns the HashFile pointer on success
 *         NULL on failure
 */
HashFile *HashFileFopen(FILE *fp) {
    HashFile *hf = HashFileCreate(0, 0);
    int archive_len;
    int i, fnamelen;

    /* Set the stdio buffer to be small to avoid massive I/O wastage */

    /* Read the header */
    hf->hfp = fp;
    hf->hf_start = ftello(hf->hfp);

    if (HHSIZE != fread(&hf->hh, 1, HHSIZE, hf->hfp)) {
	HashFileDestroy(hf);
	return NULL;
    }
    if (memcmp(HASHFILE_MAGIC, &hf->hh, 4) != 0) {
	HashFileFooter foot;
	int64_t offset;

	/* Invalid magic number, try other end of file! */
	fseeko(hf->hfp, -(off_t)sizeof(HashFileFooter), SEEK_END);
	if (sizeof(foot) != fread(&foot, 1, sizeof(foot), hf->hfp)) {
	    HashFileDestroy(hf);
	    return NULL;
	}
	if (memcmp(HASHFILE_MAGIC, &foot.magic, 4) != 0) {
	    HashFileDestroy(hf);
	    return NULL;
	}
	memcpy(&offset, foot.offset, 8);
	offset = be_int8(offset);
	fseeko(hf->hfp, offset, SEEK_CUR);
	hf->hf_start = ftello(hf->hfp);
	if (HHSIZE != fread(&hf->hh, 1, HHSIZE, hf->hfp)) {
	    HashFileDestroy(hf);
	    return NULL;
	}
    }
    if (memcmp(hf->hh.vers, HASHFILE_VERSION,    4) != 0 &&
	memcmp(hf->hh.vers, HASHFILE_VERSION100, 4) != 0) {
	/* incorrect version */
	HashFileDestroy(hf);
	return NULL;
    }
    
    hf->hh.nbuckets = be_int4(hf->hh.nbuckets);
    hf->hh.offset   = be_int8(hf->hh.offset);
    hf->hh.size     = be_int4(hf->hh.size);

    /* Load the archive filename(s) */
    hf->narchives = hf->hh.narchives;

    /* Old archives had narchives fixed as zero, so check again */
    if (!hf->narchives) {
	int n = fgetc(fp);
	if (!n) {
	    archive_len = 1;
	} else {
	    ungetc(n, fp);
	    hf->narchives = 1;
	}
    }

    if (hf->narchives) {
	hf->archives = (char **)malloc(hf->narchives * sizeof(char *));
	hf->afp = calloc(hf->narchives, sizeof(FILE *));
    } else {
	hf->archives = NULL;
	hf->afp = &hf->hfp;
    }

    if (hf->narchives) {
	archive_len = 0;
	for (i = 0; i < hf->narchives; i++) {
	    fnamelen = fgetc(fp);
	    hf->archives[i] = malloc(fnamelen+1);
	    fread(hf->archives[i], 1, fnamelen, fp);
	    hf->archives[i][fnamelen] = 0;
	    archive_len += fnamelen+1;
	}
    }

    hf->header_size = HHSIZE + archive_len +
	12 * (hf->hh.nheaders + hf->hh.nfooters);
    hf->nheaders = hf->hh.nheaders;
    hf->nfooters = hf->hh.nfooters;

    /* Load the header and footer locations */
    hf->headers = hf->nheaders
	? (HashFileSection *)malloc(hf->nheaders * sizeof(HashFileSection))
	: NULL;
    for (i = 0; i < hf->nheaders; i++) {
	fread(&hf->headers[i].pos,  8, 1, hf->hfp);
	fread(&hf->headers[i].size, 4, 1, hf->hfp);
	hf->headers[i].archive_no = *(char *)&hf->headers[i].pos;
	*(char *)&hf->headers[i].pos = 0;
	hf->headers[i].pos  = be_int8(hf->headers[i].pos) + hf->hh.offset;
	hf->headers[i].size = be_int4(hf->headers[i].size);
	hf->headers[i].cached_data = NULL;
    }

    hf->footers = hf->nfooters
	? (HashFileSection *)malloc(hf->nfooters * sizeof(HashFileSection))
	: NULL;
    for (i = 0; i < hf->nfooters; i++) {
	fread(&hf->footers[i].pos,  8, 1, hf->hfp);
	fread(&hf->footers[i].size, 4, 1, hf->hfp);
	hf->footers[i].archive_no = *(char *)&hf->footers[i].pos;
	*(char *)&hf->footers[i].pos = 0;
	hf->footers[i].pos  = be_int8(hf->footers[i].pos) + hf->hh.offset;
	hf->footers[i].size = be_int4(hf->footers[i].size);
	hf->footers[i].cached_data = NULL;
    }

    return hf;
}

HashFile *HashFileOpen(char *fname) {
    FILE *fp;
    HashFile *hf;

    /* Open the hash and read the header */
    if (NULL == (fp = fopen(fname, "rb")))
	return NULL;

    if (!(hf = HashFileFopen(fp)))
	return NULL;

    /* Open the main archive too? Usually deferred */
    if (hf->narchives) {
	int i;

	hf->afp = malloc(hf->narchives * sizeof(FILE *));
	if (hf->afp == NULL)
	    return NULL;

	/* Delay opening the main archive until required */
	for (i = 0; i < hf->narchives; i++) {
	    hf->afp[i] = NULL;
	}

#if 0
	    if (NULL == (hf->afp[i] = fopen(hf->archives[i], "rb"))) {
		/* Possibly done via a relative pathname (optimal infact) */
		char *cp;
		char aname[1024];
		if (NULL == (cp = strrchr(fname, '/'))) {
		    HashFileDestroy(hf);
		    return NULL;
		}
		sprintf(aname, "%.*s%s", (int)(cp-fname+1), fname,
			hf->archives[i]);
		if (NULL == (hf->afp[i] = fopen(aname, "rb"))) {

		    return NULL;
		}
	    }
#endif
    } else {
	hf->afp = &hf->hfp;
    }

    return hf;
}

HashFile *HashFileLoad(FILE *fp) {
    HashFile *hf;
    char *htable;
    off_t htable_pos;
    int i;
    HashItem *hi;
    HashTable *h;
    uint32_t *bucket_pos;
    uint32_t hsize;

    /* Open the hash table */
    fseeko(fp, 0, SEEK_SET);
    if (NULL == (hf = HashFileFopen(fp)))
	return NULL;

    HashTableDestroy(hf->h, 1);
    h = hf->h = HashTableCreate(hf->hh.nbuckets, hf->hh.hfunc);
    bucket_pos = (uint32_t *)calloc(h->nbuckets, sizeof(uint32_t));

    /* Also load in the entire thing to memory */
    htable = (char *)malloc(hf->hh.size);
    fseeko(fp, hf->hf_start, SEEK_SET);
    hsize = fread(htable, 1, hf->hh.size, fp);
    if (hf->hh.size != hsize) {
	free(htable);
	return NULL;
    }

    /*
     * HashFileOpen has already decoded the headers up to and including the
     * individual file header/footer sections, but not the buckets and item
     * lists, so we start from there.
     */
    htable_pos = hf->header_size;
    
    /* Identify the "bucket pos". Detemines which buckets have data */
    for (i = 0; i < h->nbuckets; i++) {
	memcpy(&bucket_pos[i], &htable[htable_pos], 4);
	bucket_pos[i] = be_int4(bucket_pos[i]);
	htable_pos += 4;
    }

    /* Read the hash table items */
    for (i = 0; i < h->nbuckets; i++) {
	if (!bucket_pos[i])
	    continue;
	for (;;) {
	    int c;
	    unsigned char uc;
	    char key[256];
	    uint64_t pos;
	    uint32_t size;
	    HashFileItem *hfi;

	    c = htable[htable_pos++];
	    if (c == EOF || !c)
		break;

	    /* key */
	    memcpy(key, &htable[htable_pos], c);
	    htable_pos += c;

	    /* header/footer */
	    uc = htable[htable_pos++];
	    hfi = (HashFileItem *)malloc(sizeof(*hfi));
	    hfi->header = (uc >> 4) & 0xf;
	    hfi->footer = uc & 0xf;

	    /* archive no. + pos */
	    memcpy(&pos, &htable[htable_pos], 8);
	    htable_pos += 8;
	    hfi->archive = *(char *)&pos;
	    *(char *)&pos = 0;
	    hfi->pos = be_int8(pos) + hf->hh.offset;

	    /* size */
	    memcpy(&size, &htable[htable_pos], 4);
	    htable_pos += 4;
	    hfi->size = be_int4(size);

	    /* Add to the hash table */
	    hi = HashItemCreate(h);
	    hi->next = h->bucket[i];
	    h->bucket[i] = hi;
	    hi->key_len = c;
	    hi->key = (char *)malloc(c+1);
	    memcpy(hi->key, key, c);
	    hi->key[c] = 0; /* For debugging convenience only */
	    hi->data.p = hfi;
	}
    }

    fflush(stderr);
    free(bucket_pos);
    free(htable);

    return hf;
}

/*
 * Searches the named HashFile for a specific key.
 * When found it returns the position and size of the object in pos and size.
 *
 * Returns
 *    0 on success (pos & size updated)
 *   -1 on failure
 */
int HashFileQuery(HashFile *hf, uint8_t *key, int key_len,
		  HashFileItem *item) {
    uint64_t hval;
    uint32_t pos;
    int klen;
    int cur_offset = 0;

    /* Hash 'key' to compute the bucket number */
    hval = hash64(hf->hh.hfunc, key, key_len) & (hf->hh.nbuckets-1);

    /* Read the bucket to find the first linked list item location */
    if (-1 == fseeko(hf->hfp, hf->hf_start + 4*hval + hf->header_size,SEEK_SET))
	return -1;
    if (4 != fread(&pos, 1, 4, hf->hfp))
	return -1;
    pos = be_int4(pos);
    cur_offset = 4*hval + 4 + hf->header_size;

    if (0 == pos)
	/* No bucket pos => key not present */
	return -1;

    /* Jump to the HashItems list and look through for key */
    if (-1 == fseeko(hf->hfp, pos - cur_offset, SEEK_CUR))
	return -1;

    for (klen = fgetc(hf->hfp); klen; klen = fgetc(hf->hfp)) {
	char k[256];
	unsigned char headfoot;
	uint64_t pos;
	uint32_t size;

	fread(k, klen, 1, hf->hfp);
	fread(&headfoot, 1, 1, hf->hfp);
	item->header = (headfoot >> 4) & 0xf;
	item->footer = headfoot & 0xf;
	fread(&pos, 8, 1, hf->hfp);
	item->archive = *(char *)&pos;
	*(char *)&pos = 0;
	pos = be_int8(pos) + hf->hh.offset;
	fread(&size, 4, 1, hf->hfp);
	size = be_int4(size);
	if (klen == key_len && 0 == memcmp(key, k, key_len)) {
	    item->pos = pos;
	    item->size = size;
	    return 0;
	}
    }

    return -1;
}

HashFile *HashFileCreate(int size, int options) {
    HashFile *hf;

    if (NULL == (hf = (HashFile *)calloc(1, sizeof(*hf))))
	return NULL;
    if (NULL == (hf->h = HashTableCreate(size, options)))
	return NULL;

    return hf;
}

void HashFileDestroy(HashFile *hf) {
    if (!hf)
	return;

    if (hf->h)
	HashTableDestroy(hf->h, 1);

    if (hf->narchives) {
	int i;
	for (i = 0; i < hf->narchives; i++)
	    if (hf->archives[i])
		free(hf->archives[i]);
	free(hf->archives);
    }

    if (hf->headers) {
	int i;
	for (i = 0; i < hf->nheaders; i++) {
	    if (hf->headers[i].cached_data)
		free(hf->headers[i].cached_data);
	}
	free(hf->headers);
    }

    if (hf->footers) {
	int i;
	for (i = 0; i < hf->nfooters; i++) {
	    if (hf->footers[i].cached_data)
		free(hf->footers[i].cached_data);
	}
	free(hf->footers);
    }

    if (hf->afp) {
	int i;

	for (i = 0; i < hf->narchives; i++)
	    if (hf->afp[i] && hf->afp[i] != hf->hfp)
		fclose(hf->afp[i]);

	if (hf->afp != &hf->hfp)
	    free(hf->afp);
    }

    if (hf->hfp)
	fclose(hf->hfp);

    free(hf);
}


/*
 * Opens a specific archive number.
 * Returns 0 on success,
 *        -1 on failure
 */
static int HashFileOpenArchive(HashFile *hf, int archive_no) {
    if (hf->narchives && archive_no > hf->narchives)
	return -1;

    if (hf->afp[archive_no])
	return 0;

    if (NULL == (hf->afp[archive_no] = fopen(hf->archives[archive_no], "rb")))
	return -1;

    return 0;
}


/*
 * Extracts the contents for a file out of the HashFile.
 */
char *HashFileExtract(HashFile *hf, char *fname, size_t *len) {
    HashFileItem hfi;
    size_t sz, pos;
    char *data;
    HashFileSection *head = NULL, *foot = NULL;

    /* Find out if and where the item is in the archive */
    if (-1 == HashFileQuery(hf, (uint8_t *)fname, strlen(fname), &hfi))
	return NULL;

    /* Work out the size including header/footer and allocate */
    sz = hfi.size;
    if (hfi.header) {
	head = &hf->headers[hfi.header-1];
	sz += head->size;
    }
    if (hfi.footer) {
	foot = &hf->footers[hfi.footer-1];
	sz += foot->size;
    }
    *len = sz;

    if (NULL == (data = (char *)malloc(sz+1)))
	return NULL;
    data[sz] = 0;

    /* Header */
    pos = 0;
    if (head) {
	HashFileOpenArchive(hf, head->archive_no);
	if (!hf->afp[head->archive_no])
	    return NULL;

	fseeko(hf->afp[head->archive_no], head->pos, SEEK_SET);
	fread(&data[pos], head->size, 1, hf->afp[head->archive_no]);
	pos += head->size;
    }

    /* Main file */
    HashFileOpenArchive(hf, hfi.archive);
    if (!hf->afp[hfi.archive])
	return NULL;

    fseeko(hf->afp[hfi.archive], hfi.pos, SEEK_SET);
    fread(&data[pos], hfi.size, 1, hf->afp[hfi.archive]);
    pos += hfi.size;

    /* Footer */
    if (foot) {
	HashFileOpenArchive(hf, foot->archive_no);
	if (!hf->afp[foot->archive_no])
	    return NULL;

	fseeko(hf->afp[foot->archive_no], foot->pos, SEEK_SET);
	fread(&data[pos], foot->size, 1, hf->afp[foot->archive_no]);
	pos += foot->size;
    }

    return data;
}

/*
 * Iterates through members of a hash table returning items sequentially.
 *
 * Returns the next HashItem on success
 *         NULL on failure.
 */
HashItem *HashTableIterNext(HashTable *h, HashIter *iter) {
    do {
	if (iter->hi == NULL) {
	    if (++iter->bnum >= h->nbuckets)
		break;
	    iter->hi = h->bucket[iter->bnum];
	} else {
	    iter->hi = iter->hi->next;
	}
    } while (!iter->hi);
    
    return iter->hi;
}

void HashTableIterReset(HashIter *iter) {
    if (iter) {
	iter->bnum = -1;
	iter->hi = NULL;
    }
}

HashIter *HashTableIterCreate(void) {
    HashIter *iter = (HashIter *)malloc(sizeof(*iter));

    HashTableIterReset(iter);
    return iter;
}

void HashTableIterDestroy(HashIter *iter) {
    if (iter)
	free(iter);
}

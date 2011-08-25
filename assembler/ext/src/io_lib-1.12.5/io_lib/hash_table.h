#ifndef _HASH_TABLE_H_
#define _HASH_TABLE_H_

#include <inttypes.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>

#include "io_lib/pooled_alloc.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The data referenced by the hash table */
typedef union {
    uint64_t i;
    void *p;
} HashData;

/* A hash item with "next" pointer to use in a linked list */
typedef struct HashItemStruct {
    HashData data;        /* user defined data attached to this key */
    char    *key;         /* key we hashed on */
    int      key_len;     /* and its length */
    struct HashItemStruct *next;
} HashItem;

/* The main hash table structure itself */
typedef struct {
    int       options;  /* HASH_FUNC & HASH_OPT macros */
    uint32_t  nbuckets; /* Number of hash buckets; power of 2 */
    uint32_t  mask;	/* bit-mask equiv of nbuckets */
    int       nused;    /* How many hash entries we're storing */
    HashItem **bucket;  /* The bucket "list heads" themselves */
    pool_alloc_t *hi_pool; /* Pool of allocated HashItem structs */
} HashTable;

/* An iterator on HashTable items */
typedef struct {
    int bnum;
    HashItem *hi;
} HashIter;

#define HASHFILE_MAGIC ".hsh"
#define HASHFILE_VERSION100 "1.00"
#define HASHFILE_VERSION "1.01"
#define HASHFILE_PREPEND -1

/* File format: the hash table header */
typedef struct {
    char magic[4];
    char vers[4];
    char hfunc;
    unsigned char nheaders;
    unsigned char nfooters;
    unsigned char narchives;
    uint32_t nbuckets;
    int64_t offset;
    uint32_t size;
} HashFileHeader;

/* sizeof(HashFileHeader) minus terminal padding */
#define HHSIZE 28

typedef struct {
    char magic[4];
    char offset[8];
} HashFileFooter;

/* The data block attached to the hash table */
typedef struct {
    uint64_t pos;
    uint32_t size;
    unsigned char archive;
    unsigned char header; /* zero if not set */
    unsigned char footer; /* zero if not set */
} HashFileItem;

/* Common headers or footers to prepend to the archive contents */
typedef struct {
    unsigned char archive_no;
    uint64_t pos;
    uint32_t size;
    unsigned char *cached_data;
} HashFileSection;

/*
 * The main structure for the HashFile functions.
 *
 * We obtain an existing HashFile by opening a stored hash file or by
 * loading the entire thing.
 * New empty ones can be created using HashFileCreate.
 */
typedef struct {
    HashFileHeader hh;		/* on-disk file header */
    HashTable *h;		/* the in-memory hash table */
    int nheaders;		/* number of common file headers */
    HashFileSection *headers;	/* on-disk common file headers struct */
    int nfooters;		/* number of common file footers */
    HashFileSection *footers;	/* on-disk common file footers struct */
    int narchives;		/* number of archive files, 0 if inline file */
    char **archives;		/* archive filenames */
    FILE *hfp;			/* hash FILE */
    FILE **afp;			/* archive FILE(s) */
    int header_size;		/* size of header + filename + N(head/feet) */
    off_t hf_start;		/* location of HashFile header in file */
} HashFile;

/* Functions to to use HashTable.options */
#define HASH_FUNC_HSIEH       0
#define HASH_FUNC_TCL         1
#define HASH_FUNC_JENKINS     2
#define HASH_FUNC_JENKINS3    3
#define HASH_FUNC_MASK        7

/* Other HashTable.options values */
#define HASH_NONVOLATILE_KEYS (1<<3)
#define HASH_ALLOW_DUP_KEYS   (1<<4)
#define HASH_DYNAMIC_SIZE     (1<<5)
#define HASH_OWN_KEYS	      (1<<6)
#define HASH_POOL_ITEMS       (1<<7)

/* Hashing prototypes */
uint32_t hash(int func, uint8_t *key, int key_len);
uint64_t hash64(int func, uint8_t *key, int key_len);
uint32_t HashJenkins(uint8_t *k, int length);
uint32_t HashTcl(uint8_t *data, int len);
uint32_t HashHsieh(uint8_t *k, int length);

/* HashTable management prototypes */
HashTable *HashTableCreate(int size, int options);
void HashTableDestroy(HashTable *h, int deallocate_date);
int HashTableResize(HashTable *h, int newsize);
HashItem *HashTableAdd(HashTable *h, char *key, int key_len,
		       HashData data, int *added);
int HashTableDel(HashTable *h, HashItem *hi, int deallocate_data);
int HashTableRemove(HashTable *h, char *key, int key_len, int deallocate_data);
HashItem *HashTableSearch(HashTable *h, char *key, int key_len);
HashItem *HashTableNext(HashItem *hi, char *key, int key_len);

void HashTableStats(HashTable *h, FILE *fp);
void HashTableDump(HashTable *h, FILE *fp);

/* Iterator prototypes */
HashIter *HashTableIterCreate(void);
void HashTableIterDestroy(HashIter *iter);
HashItem *HashTableIterNext(HashTable *h, HashIter *iter);
void HashTableIterReset(HashIter *iter);

/* HashFile prototypes */
uint64_t HashFileSave(HashFile *hf, FILE *fp, int64_t offset);
HashFile *HashFileLoad(FILE *fp);
int HashFileQuery(HashFile *hf, uint8_t *key, int key_len, HashFileItem *item);
char *HashFileExtract(HashFile *hf, char *fname, size_t *len);


HashFile *HashFileCreate(int size, int options);
void HashFileDestroy(HashFile *hf);
HashFile *HashFileOpen(char *fname);
HashFile *HashFileFopen(FILE *fp);

#ifdef __cplusplus
}
#endif

#endif /* _HASH_TABLE_H_ */

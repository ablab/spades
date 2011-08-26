#ifndef _TAR_FORMAT_H
#define _TAR_FORMAT_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Our own tar block defines - we cannot rely on UNIX to provide these for us
 * as the sun tar.h is minimal and Alliant's does not even exist.
 */
#define TBLOCK  512
#define NAMSIZ  100

/* Values used in typeflag field. */
#define REGTYPE         '0'             /* Regular File */
#define AREGTYPE        '\0'            /* Regular File */
#define LNKTYPE         '1'             /* Hard Link */
#define SYMTYPE         '2'             /* Symbolic Link */
#define CHRTYPE         '3'             /* Character Special File */
#define BLKTYPE         '4'             /* Block Special File */
#define DIRTYPE         '5'             /* Directory */
#define FIFOTYPE        '6'             /* FIFO */
#define CONTTYPE        '7'             /* Reserved */

/*
 * There will usually be more data than this in a tar header - but we don't
 * need to concern ourselves with it.
 */
typedef union hblock {
    char data[TBLOCK];
    struct header {
	char name[NAMSIZ];
	char mode[8];
	char uid[8];
	char gid[8];
	char size[12];
	char mtime[12];
	char chksum[8];
	char typeflag;
	char linkname[NAMSIZ];
	char magic[6];
	char version[2];
	char uname[32];
	char gname[32];
	char devmajor[8];
	char devminor[8];
	char prefix[155];
    } header;
} tar_block;

#ifdef __cplusplus
}
#endif

#endif

#ifndef _VLEN_H_
#define _VLEN_H_

#ifdef __cplusplus
extern "C" {
#endif

extern int vflen(char *fmt, va_list ap);
extern int flen(char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif /* _VLEN_H_ */

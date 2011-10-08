#ifndef TRACE_H__
#define TRACE_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#ifdef HAVEVAM  /* we use variadic macros if possible */

#define wtrace(fmt, ...)     wtrace_x(0, fmt, ## __VA_ARGS__)
#define wtrace_buf(fmt, ...) wtrace_x(1, fmt, ## __VA_ARGS__)
int wtrace_x(int, const char*, ...);

#else /* otherwise default to the buffered version */

#define wtrace wtrace_buf
int wtrace_buf(const char *, ...);

#endif /* HAVEVAM */

#endif


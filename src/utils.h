#ifndef UTILS_H
#define UTILS_H 1   // Pourquoi 1? 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <errno.h> 


void* _safeCalloc(const char* src,int line,size_t num_elements, size_t size);
#define safeCalloc(N,SIZE) _safeCalloc(__FILE__, __LINE__, N, SIZE);

void* _safeRealloc(const char* src, int line, void *ptr, size_t size);
#define safeRealloc(PTR,SIZE) _safeRealloc(__FILE__, __LINE__, PTR, SIZE);

FILE* _safeOpen(const char* src, int line, const char* filename, const char* mode);
#define safeOpen(FILENAME,MODE) _safeOpen(__FILE__, __LINE__, FILENAME, MODE);


#endif

#ifndef UTILS_H
#define UTILS_H 1
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>   // C Error Codes




int* safeReallocInt(int* tLen, unsigned int length);
char* safeReallocChar(char* bases, int unsigned position);    


void* _safeCalloc(const char* source,int line,size_t nmemb, size_t size)
#define safeCalloc(N,SIZE) _safeCalloc(__FILE__,__LINE__,N,SIZE)   ;

void* _safeRealloc(const char* source,int line,void *ptr, size_t size);
#define safeRealloc(PTR,SIZE) _safeRealloc(__FILE__,__LINE__,PTR,SIZE)   ;

FILE* _safeOpen(const char* source,int line,const char* filename,const char* mode)   ;

#define safeOpen(FILENAME,MODE) _safeOpen(__FILE__,__LINE__,FILENAME,MODE)   ;


#endif



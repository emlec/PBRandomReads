#include "utils.h"
#include <zlib.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>


void* _safeCalloc(const char* src, int line,size_t num_elements, size_t size)
    {
    void* ptr = calloc(num_elements,size);      
    if(ptr==NULL) {
        fprintf(stderr,"[%s:%d] OUT OF MEMORY", src,line); 
        exit(EXIT_FAILURE);
        }
    return ptr;
    }

void* _safeRealloc(const char* src, int line,void* ptr, size_t size)
    {
    void* ptr_temp = realloc(ptr,size);
    if(ptr==NULL) {
        fprintf(stderr,"[%s:%d] OUT OF MEMORY", src,line); 
        exit(EXIT_FAILURE);
        }
    return ptr_temp;
    }

gzFile _safeOpen(const char* src, int line, const char* filename, const char* mode)  //gzopen can be used to read a file which is not in gzip format; in this case gzread will directly read from the file without decompression. 
    {
    gzFile file = gzopen(filename, mode);
    if (file == NULL)                               // gzopen returns NULL if the file could not be opened
        {
        fprintf(stderr,"[%s:%d] Error opening file %s : %s.\n", src, line, filename, strerror(errno));
        exit(EXIT_FAILURE);
        }
    return file;
    }





#include "utils.h"


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

FILE* _safeOpen(const char* src, int line, const char* filename, const char* mode)  
    {
    FILE*  fic = fopen(filename, mode);
    if (fic == NULL) 
        {
        fprintf(stderr,"[%s:%d] Error opening file %s : %s.\n", src, line, filename, strerror(errno));
        exit(EXIT_FAILURE);
        }
    return fic;
    }



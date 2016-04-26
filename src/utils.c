#include "utils.h"



    
int* safeReallocInt(int* tLen, unsigned int length)
    {
    int *temp = realloc(tLen, ((length+1) * sizeof(int)));      
    if (temp == NULL) {fprintf(stderr,"%s", strerror(12)); exit(EXIT_FAILURE);};  //Out of memory
    return temp;
    }
    
void* _safeCalloc(const char* source,int line,size_t nmemb, size_t size)
    {
	void* ptr = calloc(nmemb,size);      
	if(ptr==NULL) {
        fprintf(stderr,"[%s:%d] OUT OF MEMORY", src,line); 
        exit(EXIT_FAILURE);
        }
    return ptr;
    }

void* _safeRealloc(const char* source,int line,void *ptr, size_t size) {
	void* ptr = realloc(ptr,size);      
	if(ptr==NULL) {
        fprintf(stderr,"[%s:%d] OUT OF MEMORY", src,line); 
        exit(EXIT_FAILURE);
        }
    return ptr;
	}


void* _safeRealloc(const char* source,int line,void *ptr, size_t size) {
	void* ptr = realloc(ptr,size);      
	if(ptr==NULL) {
        fprintf(stderr,"[%s:%d] OUT OF MEMORY", src,line); 
        exit(EXIT_FAILURE);
        }
    return ptr;
	}

      
FILE* _safeOpen(const char* src,int line,const char* filename,const char* mode)  
    {
    FILE*  fic = fopen(filename,mode);
    if (fic == NULL) 
        {
        fprintf(stderr,"[%s:%d]Error opening file \"%s\" : %s.\n", src,line,filename,strerror(errno));  //No such file or directory
        exit(EXIT_FAILURE);
        }
    }



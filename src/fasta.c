#include "fasta.h"
#include "utils.h"
#include <stdlib.h>
#include <zlib.h>
#include <errno.h>
#include <stdio.h>

SequencePtr _initStructSequence(const char* src, int line){
    
    int buffer_size=1024;
    
    SequencePtr seq = (SequencePtr)safeCalloc(1,sizeof(Sequence));
    if(seq==NULL) {
        fprintf(stderr,"[%s:%d] OUT OF MEMORY",src ,line); 
        exit(EXIT_FAILURE);
        }
    seq->name = (char*)safeCalloc(1, sizeof(char));
    seq->name[0] = 0;
    seq->bases = (char*)safeCalloc(buffer_size, sizeof(char));
    seq->bases[0] = 0;
    seq->length = 0;
    return seq;
}

void readOneSequenceFromFile(gzFile file, SequencePtr seq) {
    
    
    int buffer_size=1024;
    int c =0;
    unsigned int current_bases_added = 0;
    unsigned int current_header_size = 1;
    
    while ((c=gzgetc(file))!=EOF)  // gzgetc does not check to see if file is NULL
        {  
        if (c=='>') { 
            if(seq->name[0]!=0) {
                fprintf(stderr,"BE CAREFUL! The file contains more than one sequence");
                exit(EXIT_FAILURE);
            }           
            while ((c=gzgetc(file))!=EOF && c!='\n') {  // ecrire le nom de la séquence dans name
                
                seq->name = safeRealloc(seq->name, current_header_size*sizeof(char)); 
                seq->name[current_header_size-1]=c;
                current_header_size++;
            }
            fprintf(stderr, "Name of the sequence : %s\n", seq->name);
            continue; 
            }
            
        else {
            if(isspace(c)) continue;
            else {
                if (current_bases_added+1>buffer_size) {
                    buffer_size*=2;
                    seq->bases = safeRealloc(seq->bases, buffer_size*sizeof(char));}
                seq->bases[current_bases_added]=c;
                current_bases_added++;
            } 
        }
    }
    seq->length=current_bases_added;
    gzclose(file);
}



void check_reference (SequencePtr seq){
    
    unsigned int i; 
    fprintf(stderr,"Checking sequence of reference :\n");
    for (i=0; i<seq->length; i++) {
        fprintf(stderr, "Sequence %c\t position %d\n", seq->bases[i],i);
        }
    fprintf(stderr,"Size of the sequence of reference : %d\n", seq->length);
    }


void SequenceFree(SequencePtr seq) {
    if(seq==NULL) return;   // void donc return rien accepté
    free(seq->name);
    free(seq->bases);
    free(seq);
    }

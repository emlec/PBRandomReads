#include "fasta.h"
#include "utils.h"
#include <stdlib.h>
#include <zlib.h>


SequencePtr readOneSequenceFromFile(gzFile file) {
    int c =0;
    unsigned int current_bases_added = 0;
    /*unsigned int current_header_size = 0;*/
    int buffer_size=1024;
    SequencePtr seq = (SequencePtr)safeCalloc(1,sizeof(Sequence));   //void* calloc(size_t num_elements, size_t size)
    seq->name = (char*)safeCalloc(1, sizeof(char));
    seq->name[0] = 0;                           // '\0'
    seq->bases = (char*)safeCalloc(buffer_size, sizeof(char));
    seq->bases[0] = 0;
    seq->length = 0;
    
    
    while ((c=gzgetc(file))!=EOF)  // gzgetc does not check to see if file is NULL
        {  
        if (c=='>')  
            {
            if(seq->name[0]!=0) {
            fprintf(stderr,"BE CAREFUL! The file contains more than one sequence");
            exit(EXIT_FAILURE);
            }/*
            while ((c=gzgetc(file))!=EOF && c!='\n') {  // ecrire le nom de la séquence dans name
                current_header_size++;
                seq->name = safeRealloc(seq->name, current_header_size*size(char));} 
                seq->name[*/
            }
        else {
            if(isspace(c)) continue;
            if (current_bases_added+1>buffer_size) {
                buffer_size*=2;
                seq->bases = safeRealloc(seq->bases, buffer_size*sizeof(char));}
            seq->bases[current_bases_added]=c;
            current_bases_added++;
            }
        }

        seq->length=current_bases_added+1;
        fprintf(stderr,"Size of the sequence of reference : %d\n", seq->length);
        fprintf(stderr, "Part2 successfull\n\n");
        
        
return seq;
}


void SequenceFree(SequencePtr seq) {
    if(seq==NULL) return;   // void donc return rien accepté
    free(seq->name);
    free(seq->bases);
    free(seq);
    }

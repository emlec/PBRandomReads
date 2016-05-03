#ifndef FASTA_H
#define FASTA_H
#include <stdio.h>
#include <zlib.h>

typedef struct sequence_t {
    char* bases;
    char* name;
    unsigned int length;
} Sequence,*SequencePtr;


SequencePtr readOneSequenceFromFile(gzFile file);
void SequenceFree(SequencePtr seq);




#endif


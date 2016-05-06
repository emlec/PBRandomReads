#ifndef FASTA_H
#define FASTA_H
#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h> 

typedef struct sequence_t {
    char* bases;
    char* name;
    unsigned int length;
} Sequence,*SequencePtr;

SequencePtr _initStructSequence(const char* src, int line);
#define initStructSequence() _initStructSequence(__FILE__, __LINE__);
void SequenceFree(SequencePtr seq);

void readOneSequenceFromFile(gzFile file, SequencePtr seq);
void check_reference (SequencePtr seq);


#endif

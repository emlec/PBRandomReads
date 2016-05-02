#ifndef FASTA_H
#define FASTA_H
#include <stdio.h>

typedef struct sequence_t {
	char* bases;
	char* name;
	unsigned int length;
} Sequence,*SequencePtr;

SequencePtr readOneSequenceFromFile(FILE* in);
void SequenceFree(SequencePtr seq);




#endif


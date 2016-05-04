#ifndef PBINSILICO_H
#define PBINSILICO_H
#include <stdio.h>
#include <zlib.h>
#include "fasta.h"
#include "reads_distribution.h"
#include "utils.h"
#include "inSilico_reads.h"

#define CASE_MUT(B, OPT) case B : mute = OPT[rand()%3]; gzprintf(inSilicoDataSet->inSilicoReadsFile,"%c", tolower(mute)); break

typedef struct InSilicoReads_t {
    SequencePtr seq;
    gzFile inSilicoReadsFile;
    Read_DistributionPtr reads;
} InSilicoReads, *InSilicoReadsPtr;

InSilicoReadsPtr _initStructInSilicoReads(const char* src, int line, SequencePtr seq, gzFile inSilicoReadsFile, Read_DistributionPtr reads);
#define initStructInSilicoReads(SEQ_PTR, FIC, READS_PTR) _initStructInSilicoReads(__FILE__, __LINE__, SEQ_PTR, FIC, READS_PTR);

void create_PBRandomReads (InSilicoReadsPtr inSilicoDataSet);
void check_PBDataSet (InSilicoReadsPtr inSilicoDataSet);
void inSilicoDataSet_Free(InSilicoReadsPtr inSilicoDataSet);

#endif


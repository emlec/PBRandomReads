#ifndef PBINSILICO_H
#define PBINSILICO_H
#include <stdio.h>
#include <zlib.h>
#include "fasta.h"
#include "reads_distribution.h"
#include "utils.h"
#include "inSilico_reads.h"

#define CASE_MUT(B, OPT) case B : mute = OPT[rand()%3]; gzprintf(inSilicoDataSet->inSilicoReadsFile,"%c", tolower(mute)); fprintf(stderr, "%s : from %c to %c\n", "Substitution", B, tolower(mute)); break

typedef struct MuteStats_t {
    float cptINS;
    float cptDEL;
    float cptSUBST;
    float cptNOMUTE;
} MuteStats, *MuteStatsPtr;

typedef struct InSilicoReads_t {
    SequencePtr seq;
    gzFile inSilicoReadsFile;
    Read_DistributionPtr reads;
    MuteStatsPtr muteStats;
} InSilicoReads, *InSilicoReadsPtr;


MuteStatsPtr _initStructMuteStats(const char* src, int line); 
#define initStructMuteStats() _initStructMuteStats(__FILE__, __LINE__);
void inSilicoDataSet_Free(InSilicoReadsPtr inSilicoDataSet);

InSilicoReadsPtr _initStructInSilicoReads(const char* src, int line, SequencePtr seq, gzFile inSilicoReadsFile, Read_DistributionPtr reads,  MuteStatsPtr muteStats);
#define initStructInSilicoReads(SEQ_PTR, FIC, READS_PTR, STATS_PTR) _initStructInSilicoReads(__FILE__, __LINE__, SEQ_PTR, FIC, READS_PTR, STATS_PTR);
void MuteStats_Free(MuteStatsPtr muteStats);


void create_PBRandomReads (InSilicoReadsPtr inSilicoDataSet);
void check_PBDataSet(InSilicoReadsPtr inSilicoDataSet, char* filename_InSilicoReads);
void inSilicoDataSet_Free(InSilicoReadsPtr inSilicoDataSet);

void write_OriginalStrand_Read(InSilicoReadsPtr inSilicoDataSet, unsigned int i, int unsigned readNb);
void write_ReverseComplementaryStrand_Read(InSilicoReadsPtr inSilicoDataSet, unsigned int i, int unsigned readNb);
//char complementaryBase(char base);
void assembly_length_cutoff(InSilicoReadsPtr inSilicoDataSet,  Read_DistributionPtr InSilicoDistr, unsigned int coverage);


#endif


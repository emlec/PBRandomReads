#ifndef READS_DISTRIBUTION_H
#define READS_DISTRIBUTION_H
#include <stdio.h>
#include <zlib.h>

typedef struct ReadDistribution_t
    {
    int* num_elements;
    unsigned int max_length;
    } Read_Distribution, *Read_DistributionPtr; 
 
Read_DistributionPtr _initStructReadDistribution(const char* src, int line);
#define initStructReadDistribution() _initStructReadDistribution(__FILE__, __LINE__);
void ReadFree(Read_DistributionPtr reads);

void make_distribution(gzFile file, Read_DistributionPtr reads, int seq_type);
void check_distribution (Read_DistributionPtr reads);
void hist_to_R (Read_DistributionPtr reads);

#endif

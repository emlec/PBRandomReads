#ifndef READS_DISTRIBUTION_H
#define READS_DISTRIBUTION_H
#include <stdio.h>
#include <zlib.h>

typedef struct LengthToOccurence_t
    {
    int* num_elements;
    unsigned int max_length;
    } Read_Distribution, *Read_DistributionPtr; 
    
Read_DistributionPtr make_distribution(gzFile file);
void check_distribution (Read_DistributionPtr read);
void ReadFree(Read_DistributionPtr read);

#endif

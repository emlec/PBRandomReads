#include "reads_distribution.h"
#include "utils.h"
#include <stdlib.h>
#include <zlib.h>


Read_DistributionPtr _initStructReadDistribution(const char* src, int line){
    
    Read_DistributionPtr reads = (Read_DistributionPtr)safeCalloc(1, sizeof(Read_Distribution));
    if(reads==NULL) {
        fprintf(stderr,"[%s:%d] OUT OF MEMORY",src ,line); 
        exit(EXIT_FAILURE);
        }
    reads->num_elements=(int*)safeCalloc(1, sizeof(int));
    reads->num_elements[0]=0;
    reads->max_length=0;
    
    return reads;
}

void make_distribution(gzFile file, Read_DistributionPtr reads, int seq_type){   // fasta seq_type = 2, fastq seq_type = 4
    
    int current_read_base = 0;
    unsigned int current_read_length = 1;
    unsigned int table_size = 1;
    int nLine = 0;                            

    fprintf(stderr,"\nChecking ReadDistribution :\n");
    do 
        {
        current_read_base=gzgetc(file);
        if (current_read_base=='\n')            
            {   
            if (nLine%seq_type==1)          // 1==sequence  
                {
                if (current_read_length>table_size)
                    {
                    reads->num_elements=safeRealloc(reads->num_elements, current_read_length*sizeof(int));
                    while(table_size<current_read_length) {reads->num_elements[table_size]=0; table_size++;}; 
                    reads->max_length=current_read_length;
                    }
                (reads->num_elements[current_read_length-1])++;
                }
            current_read_length=0;
            nLine++;
            } 
        else 
            {
            current_read_length++;
            }
        }
    while (current_read_base!=EOF);
}


void check_distribution (Read_DistributionPtr reads){
    
    unsigned int i;
    unsigned int cptBases = 0;
        
    for (i=0; i<reads->max_length; i++) {
        if (reads->num_elements[i] != 0) {
        fprintf(stderr, "Size of read : %d bases\tnumber of reads : %d\n", i+1, reads->num_elements[i]);
        cptBases+=((i+1)*(reads->num_elements[i]));
        }
    }
    fprintf(stderr,"The max size of read is : %d bases\t Total number of bases %d\n", reads->max_length, cptBases);
}

void hist_to_R (Read_DistributionPtr reads){
    unsigned int i;
    gzFile output_hist_file = NULL;
    output_hist_file = safeOpen("output_hist","w");
    gzprintf(output_hist_file,"%s\t%s\n", "size of read", "number of read");
    for (i=0; i<reads->max_length; i++) {
        gzprintf(output_hist_file,"%d\t%d\n", i+1, reads->num_elements[i]);
    }
    gzclose(output_hist_file);
}


void ReadFree(Read_DistributionPtr reads) {
    if(reads==NULL) return;   // void donc return rien accepté
    free(reads->num_elements);
    free(reads);
    }

#include "reads_distribution.h"
#include "utils.h"
#include <stdlib.h>
#include <zlib.h>


Read_DistributionPtr make_distribution(gzFile file){

    int current_read_base = 0;
    unsigned int current_read_length = 0;
    unsigned int table_size = 1;
    int type_of_line = 0;                            

        
    Read_DistributionPtr read = (Read_DistributionPtr)safeCalloc(1, sizeof(Read_Distribution));
    read->num_elements=(int*)safeCalloc(1, sizeof(int));
    read->num_elements[0]=0;
    read->max_length=0;
            

    while ((current_read_base=gzgetc(file))!=EOF) 
        {
        if (current_read_base=='\n')            
            {   
            if (type_of_line%4==1)          // 1==sequence  
                {
                if (current_read_length>=table_size)
                    {
                    read->num_elements=safeRealloc(read->num_elements, (current_read_length+1)*sizeof(int));
                    while(table_size<=current_read_length) {read->num_elements[table_size]=0; table_size++;}; 
                    table_size = current_read_length + 1;
                    }
                read->num_elements[current_read_length-1]++;
                if (current_read_length>read->max_length) { read->max_length=current_read_length; printf("The max size of read is : %d bases\n", read->max_length);}
                }
            current_read_length=0;
            type_of_line++;
            } 
         else 
            {
            current_read_length++;
            }
        }
    return read;
    }


void check_distribution (Read_DistributionPtr read){
    
    unsigned int i;    
    for (i=0; i<=read->max_length; i++) {
        if (read->num_elements != 0) {
        fprintf(stdout, "size of read %dbases\tnumber of reads%d\n", i+1, read->num_elements[i]);
        }
    fprintf(stdout, "Part1 successfull\n\n");
    }
}


void ReadFree(Read_DistributionPtr read) {
    if(read==NULL) return;   // void donc return rien accepté
    free(read->num_elements);
    free(read);
    }

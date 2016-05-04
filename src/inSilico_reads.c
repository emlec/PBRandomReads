#include <stdlib.h>
#include <zlib.h>
#include "inSilico_reads.h"
#include "reads_distribution.h"
#include "utils.h"
#include "fasta.h"


InSilicoReadsPtr _initStructInSilicoReads(const char* src, int line, SequencePtr seqPtr, gzFile inSilicoReadsFile, Read_DistributionPtr readsPtr){
    InSilicoReadsPtr inSilicoDataSet =(InSilicoReadsPtr)safeCalloc(1, sizeof(InSilicoReads));
    if(inSilicoDataSet==NULL) {
        fprintf(stderr,"[%s:%d] OUT OF MEMORY",src ,line); 
        exit(EXIT_FAILURE);
    }
    inSilicoDataSet->seq= seqPtr;   //je prends la valeur contenue dans seqPtr c'est à dire l'adresse de seq
    inSilicoDataSet->inSilicoReadsFile=inSilicoReadsFile;
    inSilicoDataSet->reads=readsPtr; 
    return inSilicoDataSet;
}

void create_PBRandomReads (InSilicoReadsPtr inSilicoDataSet){

    unsigned int i; 
    int j;
    unsigned int k=0;
    int start_position = 0;
    int prob;
    char ins;               // the inserted base
    char mute;              // the mutated base
    int readNb = 1;         // sequence generated (print in the header of the output)


    for (i=0; i<=inSilicoDataSet->reads->max_length; i++)        // for each size of read
            {                  
            if(inSilicoDataSet->reads->num_elements[i]!=0)                      // if read exist
                {
                for (j=0; j<inSilicoDataSet->reads->num_elements[i]; j++)           // for each read
                    {  
                    start_position = 1 + rand()%inSilicoDataSet->seq->length;      // select the beginning in the reference genome
                    gzprintf(inSilicoDataSet->inSilicoReadsFile, ">m160129_165300_42263_c100880132550000001823194304021670_s1_X0/%d/0_%d\n", readNb, i);  // For FALCON compatibility
                    
                    for (k=1; k<=i;k++)                                  // for each base 
                        {
                        if ((start_position+k)>inSilicoDataSet->seq->length) {start_position=start_position-inSilicoDataSet->seq->length;}    
                        prob = (rand()%101);                               // select the probability of mutation [1;100]
                        
                        // INSERTION : 11% [0-11]
                        if (prob<=11) 
                            {            
                            ins = "ATCG"[rand()%4];
                            gzprintf(inSilicoDataSet->inSilicoReadsFile, "%c%c", toupper(inSilicoDataSet->seq->bases[start_position+k]), tolower(ins));
                            }
                                                  
                        // SUBSTITUTION : 1% [12]
                        else if (prob<=12) 
                            {
                            switch(toupper(inSilicoDataSet->seq->bases[start_position+k])) 
                                {
                                CASE_MUT('A', "TCG");      
                                CASE_MUT('T', "ACG");
                                CASE_MUT('C', "TAG");
                                CASE_MUT('G', "TCA");
                                };
                            }
                            
                        // CORRECT : 84% [13-96]
                        else if (prob<=96) 
                            {
                            gzprintf(inSilicoDataSet->inSilicoReadsFile, "%c", toupper(inSilicoDataSet->seq->bases[start_position+k]));
                            }
                        }
                    gzprintf(inSilicoDataSet->inSilicoReadsFile, "\n");
                    readNb++;                                   // Incrémente le nombre de reads Insilico crées
                    }
                }
                else {continue;}
            }
        
    }

void check_PBDataSet(InSilicoReadsPtr inSilicoDataSet){
    
    Read_DistributionPtr InSilicoDistr = initStructReadDistribution();
  
    make_distribution(inSilicoDataSet->inSilicoReadsFile,InSilicoDistr, 2);
    
    fprintf(stderr, "check length max %d", InSilicoDistr->max_length);
    check_distribution(InSilicoDistr);
    ReadFree(InSilicoDistr);
}



void inSilicoDataSet_Free(InSilicoReadsPtr inSilicoDataSet){
    if(inSilicoDataSet==NULL) return;
    free(inSilicoDataSet);
}

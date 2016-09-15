#include <stdlib.h>
#include <zlib.h>
#include "inSilico_reads.h"
#include "reads_distribution.h"
#include "utils.h"
#include "fasta.h"


static char complementaryBase(char base){     // STATIC : used only in this file 
    //char base_complement;
    switch (toupper(base)) {
        case 'A' :  return 'T'; break;
        case 'T' :  return 'A'; break;
        case 'C' :  return 'G'; break;
        case 'G' :  return 'C'; break;
        case 'N' :  return 'N'; break;
        default : fprintf(stderr, "Debug the complementary base : [%s:%d] \n", __FILE__, __LINE__); exit(EXIT_FAILURE);
    };
    //return base_complement;
}


/*static makeRnd(SequencePtr seqPtr )                       TODO  +  fputc instead of fprintf
    {
    int plus_strand=(rand...);
    int start = ..
    int end = start + length...
    int fragment_length= end-start;
    for(i=0; i < fragment_length ;++i)  {
        char c;
        if(plus_strand)
            {
            c= seqPtr->bases[(start+i)%(seqPtr->length)];
            }
        else
            {
            c = seqPtr->bases[(end-i)%(seqPtr->length)];
            switch(c)
                {
                case 'A' : case 'a': c= 't'; 
                }
            }
        
    
    
        }
    
    
    
    
    }
*/
MuteStatsPtr _initStructMuteStats(const char* src, int line) {
    MuteStatsPtr muteStats = (MuteStatsPtr)safeCalloc(1, sizeof(MuteStats));
    muteStats->cptINS=0;
    muteStats->cptDEL=0;
    muteStats->cptSUBST=0;
    muteStats->cptNOMUTE=0;
    return muteStats;
}

InSilicoReadsPtr _initStructInSilicoReads(const char* src, int line, SequencePtr seqPtr, gzFile inSilicoReadsFile, Read_DistributionPtr readsPtr, MuteStatsPtr muteStatsPtr){
    InSilicoReadsPtr inSilicoDataSet =(InSilicoReadsPtr)safeCalloc(1, sizeof(InSilicoReads));
    if(inSilicoDataSet==NULL) {
        fprintf(stderr,"[%s:%d] OUT OF MEMORY",src ,line); 
        exit(EXIT_FAILURE);
    }
    inSilicoDataSet->seq= seqPtr;   //je prends la valeur contenue dans seqPtr c'est à dire l'adresse de seq
    inSilicoDataSet->muteStats = muteStatsPtr;
    inSilicoDataSet->inSilicoReadsFile=inSilicoReadsFile;
    inSilicoDataSet->reads=readsPtr; 
    return inSilicoDataSet;
}


void create_PBRandomReads (InSilicoReadsPtr inSilicoDataSet){

    unsigned int i; 
    int j;
    int unsigned prob_strand=0;
    unsigned int readNb = 1;         // sequence generated (print in the header of the output)


    fprintf(stderr,"In Silico read writing :\n");
    for (i=0; i<=inSilicoDataSet->reads->max_length; i++)        // for each size of read
            {
            if(inSilicoDataSet->reads->num_elements[i]!=0)                      // if read exist
                {
                fprintf(stderr, "\nSize of input reads : %d bases\n", i+1);
                for (j=0; j<inSilicoDataSet->reads->num_elements[i]; j++)           // for each read
                    {
                    fprintf(stderr, "Read number %d to %d\n", j+1, inSilicoDataSet->reads->num_elements[i]);
                    prob_strand=rand()%2 + 1;                                       // select the strand
                    fprintf(stderr, "prob_strand : %d\n", prob_strand);
                    switch(prob_strand){
                        case 1 : 
                            write_OriginalStrand_Read(inSilicoDataSet,i, readNb); 
                            readNb++;
                            break;
                        case 2 : 
                            write_ReverseComplementaryStrand_Read(inSilicoDataSet, i, readNb); 
                            readNb++;
                            break;
                        default : 
                            fprintf(stderr, "Debug the probability of strand : [%s:%d] \n", __FILE__, __LINE__); 
                            exit(EXIT_FAILURE);
                        }
                    }
                    //readNb++;  //Number of reads
                }
            else {continue;}
            }
    gzclose(inSilicoDataSet->inSilicoReadsFile);
    }


void check_PBDataSet(InSilicoReadsPtr inSilicoDataSet, char* filename_InSilicoReads){
    
    float cptTOT = inSilicoDataSet->muteStats->cptINS + inSilicoDataSet->muteStats->cptDEL + inSilicoDataSet->muteStats->cptSUBST + inSilicoDataSet->muteStats->cptNOMUTE;
    
    Read_DistributionPtr InSilicoDistr = initStructReadDistribution();
    inSilicoDataSet->inSilicoReadsFile = safeOpen(filename_InSilicoReads,"r");
    make_distribution(inSilicoDataSet->inSilicoReadsFile,InSilicoDistr, 4);
    check_distribution(InSilicoDistr);
    fprintf(stderr, "Statistiques of mutations :\nInsertion %.1f%%\t Deletion %.1f%%\t Substitution %.1f%%\t  No mutation %.1f%% \n", (inSilicoDataSet->muteStats->cptINS/cptTOT)*100, (inSilicoDataSet->muteStats->cptDEL/cptTOT)*100, (inSilicoDataSet->muteStats->cptSUBST/cptTOT)*100, (inSilicoDataSet->muteStats->cptNOMUTE/cptTOT)*100);
    // Analysis
    hist_to_R (InSilicoDistr);
    assembly_length_cutoff(inSilicoDataSet, InSilicoDistr,20); 
    ReadFree(InSilicoDistr);
    gzclose(inSilicoDataSet->inSilicoReadsFile);
}

void inSilicoDataSet_Free(InSilicoReadsPtr inSilicoDataSet){
    if(inSilicoDataSet==NULL) return;
    free(inSilicoDataSet);
}

void MuteStats_Free(MuteStatsPtr muteStats){
    if(muteStats==NULL) return;
    free(muteStats);
}



void write_OriginalStrand_Read(InSilicoReadsPtr inSilicoDataSet, unsigned int i, int unsigned readNb){
    int k;
    int start_position = 0;
    int unsigned prob_mutation=0;
    //char* strand;
    char ins;               // the inserted base
    char mute;              // the mutated base
    float cptIns=0;
    float cptDel=0;
    float cptSubst=0;
    float cptNoMute=0;
    int cpt_quality=0;
    int j=0;
    
    
    //strand = "originalStrand"; 
    gzprintf(inSilicoDataSet->inSilicoReadsFile, "@m160129_165300_42263_c100880132550000001823194304021670_s1_X0/%d/0_%d\n", readNb, i);  // For FALCON compatibility      //Delete RCS/OS
   
    start_position = 1 + rand()%inSilicoDataSet->seq->length;      // select the beginning in the reference genome
    for (k=0; k<=i;k++){
        if ((start_position+k)>=inSilicoDataSet->seq->length) {start_position=start_position-inSilicoDataSet->seq->length;};    // circularize the genome of reference  strandOriginal
        prob_mutation = (rand()%100 + 1);                               // select the probability of mutation [1;100]
        fprintf(stderr, "%s\tstart_position %d\tposition %d\tinitial base %c\tProbability of mutation : %d\t", "Ostrand", start_position, start_position+k, inSilicoDataSet->seq->bases[start_position+k],prob_mutation);
        
        // INSERTION : 11% [0-11]
        if (prob_mutation<=11) 
            {            
            ins = "ATCG"[rand()%4];
            gzprintf(inSilicoDataSet->inSilicoReadsFile, "%c%c", toupper(inSilicoDataSet->seq->bases[start_position+k]), tolower(ins));
            fprintf(stderr, "%s : from %c to %c%c\n", "Insertion", toupper(inSilicoDataSet->seq->bases[start_position+k]), toupper(inSilicoDataSet->seq->bases[start_position+k]), tolower(ins));
            cptIns++;
            }
                                  
        // SUBSTITUTION : 1% [12]
        else if (prob_mutation<=12) 
            {
            switch(inSilicoDataSet->seq->bases[start_position+k])
                {
                CASE_MUT('A', "TCG");      
                CASE_MUT('T', "ACG");
                CASE_MUT('C', "TAG");
                CASE_MUT('G', "TCA"); 
                default :  fprintf(stderr, "Debug the substitution case : [%s:%d] \n", __FILE__, __LINE__); exit(EXIT_FAILURE);
                };
            cptSubst++;
            }
            
        // CORRECT : 84% [13-96]
        else if (prob_mutation<=96) 
            {
            gzprintf(inSilicoDataSet->inSilicoReadsFile, "%c", toupper(inSilicoDataSet->seq->bases[start_position+k]));
            fprintf(stderr, "%s : from %c to %c\n", "No mutation", toupper(inSilicoDataSet->seq->bases[start_position+k]),toupper(inSilicoDataSet->seq->bases[start_position+k]));
            cptNoMute++;
            }
            
        // DELETION : 4% [97-100]
        else if (prob_mutation<=100) 
            {
            gzprintf(inSilicoDataSet->inSilicoReadsFile, "", "");
            fprintf(stderr, "%s : from %c to nothing\n", "Deletion", toupper(inSilicoDataSet->seq->bases[start_position+k]));
            cptDel++;
            }
        }
        // For artificial read_quality (to control Insilico Data Set Distribution with fastqc)
        cpt_quality=((1*cptSubst)+(1*cptNoMute)+(2*cptIns));

        //fprintf(stderr, "quality_seq %s", quality_seq);
        
        // Newline
        gzprintf(inSilicoDataSet->inSilicoReadsFile, "\n+\n"); // gzprintf : buffer size 8192bytes
        for (j=0; j<cpt_quality;j++) {
            gzprintf(inSilicoDataSet->inSilicoReadsFile, "%c", '+');
            }
        gzprintf(inSilicoDataSet->inSilicoReadsFile, "\n");

        
        fprintf(stderr, "Insertion %.1f%%\t Deletion %.1f%%\t Substitution %.1f%%\t  No mutation %.1f%% \n", ((cptIns/(i+1))*100), ((cptDel/(i+1))*100), ((cptSubst/(i+1))*100), ((cptNoMute/(i+1)*100)));
        //fprintf(stderr, "\n+\n%s\n", quality_seq);
        inSilicoDataSet->muteStats->cptINS+=cptIns;
        inSilicoDataSet->muteStats->cptSUBST+=cptSubst;
        inSilicoDataSet->muteStats->cptNOMUTE+=cptNoMute;
        inSilicoDataSet->muteStats->cptDEL+=cptDel;
        
}


void write_ReverseComplementaryStrand_Read(InSilicoReadsPtr inSilicoDataSet, unsigned int i, int unsigned readNb){


    int k;
    int start_position = 0;
    int unsigned prob_mutation=0;
    //char* strand;
    char ins;               // the inserted base
    char mute;              // the mutated base
    char base_complement;
    float cptIns=0;
    float cptDel=0;
    float cptSubst=0;
    float cptNoMute=0;
    int cpt_quality=0;
    int j=0;
    
    //strand = "reverseComplementaryStrand"; 
    gzprintf(inSilicoDataSet->inSilicoReadsFile, "@m160129_165300_42263_c100880132550000001823194304021670_s1_X0/%d/0_%d\n", readNb,i);  // For FALCON compatibility  Delete RCS/OS
    start_position = 1 + rand()%inSilicoDataSet->seq->length;      // select the beginning in the reference genome
    
    
    for (k=0; k<=i;k++){
        if ((start_position-k)<0) {
            start_position=start_position+inSilicoDataSet->seq->length;  // circularize the genome of reference  strandOriginal
            }
        base_complement = complementaryBase(inSilicoDataSet->seq->bases[start_position-k]);
        prob_mutation = (rand()%100 + 1);                               // select the probability of mutation [1;100]
        fprintf(stderr, "%s\tstart_position %d\tposition %d\tinitial base %c\tprobability of mutation : %d\t","RCStrand",start_position, start_position-k,inSilicoDataSet->seq->bases[start_position-k],prob_mutation);

        // INSERTION : 11% [1-11]
        if (prob_mutation<=11) 
            {            
            ins = "ATCG"[rand()%4];
            gzprintf(inSilicoDataSet->inSilicoReadsFile, "%c%c", base_complement, tolower(ins));
            fprintf(stderr, "%s : from %c to %c%c\n", "Insertion", inSilicoDataSet->seq->bases[start_position-k], base_complement, tolower(ins));
            cptIns++;
            }
                                  
        // SUBSTITUTION : 1% [12]
        else if (prob_mutation<=12) 
            {
            switch(toupper(base_complement)) 
                {
                CASE_MUT('A', "TCG");      
                CASE_MUT('T', "ACG");
                CASE_MUT('C', "TAG");
                CASE_MUT('G', "TCA"); 
                default :  fprintf(stderr, "Debug the substitution case : [%s:%d] \n", __FILE__, __LINE__); exit(EXIT_FAILURE);
                };
            cptSubst++;
            }
          
        // CORRECT : 84% [13-96]
        else if (prob_mutation<=96) 
            {
            gzprintf(inSilicoDataSet->inSilicoReadsFile, "%c", toupper(base_complement));
            fprintf(stderr, "%s : from %c to %c\n", "No mutation",inSilicoDataSet->seq->bases[start_position-k], base_complement);
            cptNoMute++;
            }
        
        // DELETION : 4% [97-100]
        else if (prob_mutation<=100) 
            {
            gzprintf(inSilicoDataSet->inSilicoReadsFile, "", "");
            fprintf(stderr, "%s : from %c to nothing\n", "Deletion", inSilicoDataSet->seq->bases[start_position-k]);
            cptDel++;
            }

        }
        // For artificial read_quality (to control Insilico Data Set Distribution with fastqc)
        cpt_quality=((1*cptSubst)+(1*cptNoMute)+(2*cptIns));


        //fprintf(stderr, "quality_seq %s", quality_seq);
       
        // Newline
        gzprintf(inSilicoDataSet->inSilicoReadsFile, "\n+\n"); // gzprintf : buffer size 8192bytes
        for (j=0; j<cpt_quality;j++) {
            gzprintf(inSilicoDataSet->inSilicoReadsFile, "%c", '+');
            }
        gzprintf(inSilicoDataSet->inSilicoReadsFile, "\n");
        
        
        fprintf(stderr, "Insertion %.1f%%\t Deletion %.1f%%\t Substitution %.1f%%\t  No mutation %.1f%% \n", ((cptIns/(i+1))*100), ((cptDel/(i+1))*100), ((cptSubst/(i+1))*100), ((cptNoMute/(i+1)*100)));
        //fprintf(stderr, "\n+\n%s\n", quality_seq);
        inSilicoDataSet->muteStats->cptINS+=cptIns;
        inSilicoDataSet->muteStats->cptSUBST+=cptSubst;
        inSilicoDataSet->muteStats->cptNOMUTE+=cptNoMute;
        inSilicoDataSet->muteStats->cptDEL+=cptDel;
        
        
                                       
}


void assembly_length_cutoff(InSilicoReadsPtr inSilicoDataSet,  Read_DistributionPtr InSilicoDistr, unsigned int coverage){
    unsigned int i = InSilicoDistr->max_length; 
    unsigned int cpt_bases = 0;
    fprintf(stderr,"The max inSilico read length is %d\n", InSilicoDistr->max_length);
    while (cpt_bases<((inSilicoDataSet->seq->length)*coverage)){
        if (i==0) {
            fprintf(stderr, "Insufficient reads to cover %d X the genome\n", coverage);
            return;
        }
        else {
            cpt_bases+=((InSilicoDistr->num_elements[i-1])*i);
            i--;
        }
    }
    fprintf(stdout, "\nTHE LENGTH CUTOFF FOR THIS ASSEMBLY USING %d COVERAGE IS %d BASES (TOTAL OF BASES USED %d).\n", coverage, i+1, cpt_bases);
}
    


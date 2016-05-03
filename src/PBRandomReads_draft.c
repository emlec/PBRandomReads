#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include "utils.h"
#include <zlib.h>
#include "fasta.h"
#include "reads_distribution.h"


#define CASE_MUT(B, OPT) case B : mute = OPT[rand()%3]; gzprintf(inSilicoReads,"%c", tolower(mute)); break
 
/*struct PBRandomReads {
    struct sequence_t seq;
    gzfile output;
    struct LengthToOccurence;
    };
*/
 
struct sequence {char* bases; unsigned int seqLength;};

int main (int argc, char** argv) {
    char* filename_InSilicoReads = NULL;
    char* filename_reference = NULL;
    int opt;
    time_t begin, end; 

     
// Part1 : a table (tLen) containing the number of read for each size 
          
    gzFile PBReads_fastqFile = NULL;        // argv[1] = fastq file containing the initial PBReads  
    Read_DistributionPtr reads = (Read_DistributionPtr)safeCalloc(1, sizeof(Read_Distribution));
    reads->num_elements=(int*)safeCalloc(1, sizeof(int));
    reads->num_elements[0]=0;
    reads->max_length=0;
    

// Part 2 :  ref is a sequence structure containing informations (bases and length) of the reference
// Example : a genome of 8 bases AATTCCGG, ref.bases[0]='A', seqLength = 8

    gzFile reference_file = NULL;        // fasta file containing the sequence of reference
    SequencePtr reference = (SequencePtr)safeCalloc(1,sizeof(Sequence));   //void* calloc(size_t num_elements, size_t size)
    //reference->name = (char*)safeCalloc(1, sizeof(char));
    //reference->name[0] = 0;                           // '\0'
    //reference->bases = (char*)safeCalloc(buffer_size, sizeof(char));
    //reference->bases[0] = 0;
    //reference->length = 0;


// Part 3 : Create a fasta file containing the InSilico reads 
// Example : >Read1
//           ATCCc (lowerletter for insertion or substitution)

    gzFile inSilicoReads = NULL;                // Will contain the PB InSilico reads          
    int j;
    unsigned int k=0;
    int start_position = 0;
    int prob;
    char ins;               // the inserted base
    char mute;              // the mutated base
    int readNb = 1;         // sequence generated (print in the header of the output)
    unsigned int i; 


// BEGIN

    time(&begin);

//Parsing arguments using getopt
    while ((opt = getopt(argc, argv, "o:r:s:")) != -1) {
        switch (opt) {
            case 'o':
                filename_InSilicoReads = optarg;
                break;
            case 'r':
                filename_reference = optarg;
                break;
            case 's': 
                srand(time(NULL));
            default: 
                fprintf(stderr, "error %s\n", argv[0]); // TODO 
                exit(EXIT_FAILURE);
           }
        }

    if (filename_reference == NULL) {
        fprintf(stderr, "Error undefined reference\n");
        return EXIT_FAILURE;
        }
    else {
        reference_file = safeOpen(filename_reference,"rb")
        }

    // PROCESS INPUT
    if (optind +1 == argc) {   //UN SEUL FICHIER SPECIFIE
        PBReads_fastqFile = safeOpen(argv[optind], "rb");
        }
    else if (optind==argc) //PAS DE FICHIER SPECIFIE, ex : gunzip -c file.fastq.gz | program -o *.fasta -r *.fasta
        {
        fprintf(stderr,"Reading fastq from STDIN\n");
        PBReads_fastqFile = gzdopen(fileno(stdin), "rb");
        }
    else
        {
        fprintf(stderr,"%s : Too many input files specified\n",argv[0]);
        return EXIT_FAILURE;
        }

    //OUTPUT
    if (filename_InSilicoReads == NULL) {
        fprintf(stderr,"Writing output to STDOUT\n");
        inSilicoReads = gzdopen(fileno(stdout), "wb");
        }
    else {
        inSilicoReads = safeOpen(filename_InSilicoReads,"wb");
        }


// Part1 : a table (tLen) containing the number of read for each size 

    reads=make_distribution(PBReads_fastqFile);
    check_distribution (reads);
    gzclose(PBReads_fastqFile);


// Part 2 :  ref is a sequence structure containing informations (bases and length) of the genome of reference
// Example : a genome of 8 bases AATTCCGG, ref.bases[0]='A', seqLength = 8
    
    reference=readOneSequenceFromFile(reference_file);
    gzclose(reference_file);

// Part 3 : Create a fasta file containing the InSilico reads 
// Example : >Read1
//           ATCCc (lowerletter for insertion or substitution)


    for (i=0; i<=reads->max_length; i++)        // for each size of read
        {                  
        if(reads->num_elements[i]!=0)                      // if read exist
            {
            for (j=0; j<reads->num_elements[i]; j++)           // for each read
                {  
                start_position = 1 + rand()%reference->length;      // select the beginning in the reference genome
                gzprintf(inSilicoReads, ">m160129_165300_42263_c100880132550000001823194304021670_s1_X0/%d/0_%d\n", readNb, i);  // For FALCON compatibility
                
                for (k=1; k<=i;k++)                                  // for each base 
                    {
                    if ((start_position+k)>reference->length) {start_position=start_position-reference->length;}    
                    prob = (rand()%101);                               // select the probability of mutation [1;100]
                    
                    // INSERTION : 11% [0-11]
                    if (prob<=11) 
                        {            
                        ins = "ATCG"[rand()%4];
                        gzprintf(inSilicoReads, "%c%c", toupper(reference->bases[start_position+k]), tolower(ins));
                        }
                                              
                    // SUBSTITUTION : 1% [12]
                    else if (prob<=12) 
                        {
                        switch(toupper(reference->bases[start_position+k])) 
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
                        gzprintf(inSilicoReads, "%c", toupper(reference->bases[start_position+k]));
                        }
                    }
                gzprintf(inSilicoReads, "\n");
                readNb++;                                   // Incrémente le nombre de reads Insilico crées
                }
            }
            else {continue;}
        }
    gzclose(inSilicoReads);
    fprintf(stderr, "Part3 successfull\n\n");
    SequenceFree(reference);
    ReadFree(reads);
    time(&end);

    printf("Elapsed time %f sec.\n", difftime(end, begin));

    return 0;
    }

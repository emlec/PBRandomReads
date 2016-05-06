#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <zlib.h>
#include "utils.h"
#include "fasta.h"
#include "reads_distribution.h"
#include "inSilico_reads.h"

// Part1 : a table (tLen) containing the number of read for each size 
// Part 2 :  ref is a sequence structure containing informations (bases and length) of the reference
// Example : a genome of 8 bases AATTCCGG, ref.bases[0]='A', seqLength = 8
// Part 3 : Create a fasta file containing the InSilico reads 
// Example : >Read1
//           ATCCc (lowerletter for insertion or substitution)



int main (int argc, char** argv) {
    char* filename_InSilicoReads = NULL;
    char* filename_reference = NULL;
    int opt;
    time_t begin, end; 

          
    gzFile PBReads_fastqFile = NULL;        // argv[1] = fastq file containing the initial PBReads  
    Read_DistributionPtr readsDistrPtr = initStructReadDistribution();
    
    gzFile reference_file = NULL;        // fasta file containing the sequence of reference
    SequencePtr referencePtr = initStructSequence();

    gzFile inSilicoReadsFile = NULL;                // Will contain the PB InSilico reads   
    MuteStatsPtr muteStatsPtr = initStructMuteStats();      
    InSilicoReadsPtr inSilicoDataSet_Ptr = NULL;
   
    time(&begin);

//Parsing arguments using getopt
    while ((opt = getopt(argc, argv, "o:r:s")) != -1) {
        switch (opt) {
            case 'o':
                filename_InSilicoReads = optarg;
                break;
            case 'r':
                filename_reference = optarg;
                break;
            case 's': 
                srand(time(NULL));
                break;
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
        reference_file = safeOpen(filename_reference,"rb");
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
        inSilicoReadsFile = gzdopen(fileno(stdout), "wb");
        }
    else {
        inSilicoReadsFile = safeOpen(filename_InSilicoReads,"wb");
        }


// Part1
    fprintf(stderr,"##########\nPART 1\n##########\n");
    make_distribution(PBReads_fastqFile, readsDistrPtr,4);
    check_distribution(readsDistrPtr);
    fprintf(stderr, "\nPart 1 successfull\n\n");

// Part 2
    fprintf(stderr,"##########\nPART 2\n##########\n");
    readOneSequenceFromFile(reference_file, referencePtr);
    check_reference (referencePtr);
    fprintf(stderr, "\nPart 2 successfull\n\n");
    
// Part 3
    fprintf(stderr,"##########\nPART 3\n##########\n");
    inSilicoDataSet_Ptr = initStructInSilicoReads(referencePtr, inSilicoReadsFile, readsDistrPtr, muteStatsPtr);
    create_PBRandomReads(inSilicoDataSet_Ptr);
    check_PBDataSet(inSilicoDataSet_Ptr, filename_InSilicoReads);
    fprintf(stderr, "\nPart 3 successfull\n\n");



// End

    time(&end);
    fprintf(stderr, "Elapsed time %f sec.\n", difftime(end, begin));

    gzclose(PBReads_fastqFile);
    gzclose(reference_file);
    gzclose(inSilicoReadsFile);
    
    SequenceFree(referencePtr);
    ReadFree(readsDistrPtr);
    MuteStats_Free(muteStatsPtr);
    inSilicoDataSet_Free(inSilicoDataSet_Ptr);

    return 0;
    }

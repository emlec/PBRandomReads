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

    make_distribution(PBReads_fastqFile, readsDistrPtr,4);
    check_distribution(readsDistrPtr);
    fprintf(stderr, "Part1 successfull\n\n");


// Part 2
    
    readOneSequenceFromFile(reference_file, referencePtr);
    check_reference (referencePtr);
    fprintf(stderr, "Part2 successfull\n\n");
    
// Part 3

    inSilicoDataSet_Ptr = initStructInSilicoReads(referencePtr, inSilicoReadsFile, readsDistrPtr);
    create_PBRandomReads(inSilicoDataSet_Ptr);
    check_PBDataSet(inSilicoDataSet_Ptr);
    fprintf(stderr, "Part3 successfull\n\n");




// End

    gzclose(PBReads_fastqFile);
    gzclose(inSilicoReadsFile);
    SequenceFree(referencePtr);
    ReadFree(readsDistrPtr);
    inSilicoDataSet_Free(inSilicoDataSet_Ptr);

    time(&end);

    fprintf(stderr, "Elapsed time %f sec.\n", difftime(end, begin));

    return 0;
    }

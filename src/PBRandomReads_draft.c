#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include "utils.h"


#define CASE_MUT(B, OPT) case B : mute = OPT[rand()%3]; fprintf(output_file,"%c", tolower(mute)); break
 
 
typedef struct LengthToOccurence
	{
	int* lengths;
	unsigned int num_elements;
	}; 
 
struct sequence {char* bases; unsigned int seqLength;};

int main (int argc, char** argv) {
    char* filename_out = NULL;
    char* filename_reference = NULL;
    int opt;
    time_t begin, end; 

     
// Part1 : a table (tLen) containing the number of read for each size 
          
    FILE* input_file = NULL;        // argv[1] = fastq file containing the initial PBReads  
    
    int c = 0;                                
    unsigned int length = 0;                  // Read Length
    int nLine = 0;                            // Number of line (1==sequence)  
    int* tLen = NULL;
    unsigned int table_size = 1;              
    unsigned int max_length=0;
    unsigned int i;

// Part 2 :  ref is a sequence structure containing informations (bases and length) of the reference
// Example : a genome of 8 bases AATTCCGG, ref.bases[0]='A', seqLength = 8

    FILE* reference_file = NULL;        // fasta file containing the sequence of reference

    int b = 0;                     
    struct sequence ref;
    unsigned int position = 1;

// Part 3 : Create a fasta file containing the InSilico reads 
// Example : >Read1
//           ATCCc (lowerletter for insertion or substitution)

    FILE* output_file = NULL;                // Will contain the PB InSilico reads          
    int j;
    unsigned int k=0;
    int start_position = 0;
    int prob;
    char ins;               // the inserted base
    char mute;              // the mutated base
    int readNb = 1;         // sequence generated (print in the header of the output)


// BEGIN

    time(&begin);

//Parsing arguments using getopt
    while ((opt = getopt(argc, argv, "o:r:")) != -1) {
          switch (opt) {
           case 'o':
               filename_out = optarg;
               break;
           case 'r':
                filename_reference = optarg;
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
        reference_file = safeOpen(filename_reference,"r")
        }

    // PROCESS INPUT
    if (optind +1 == argc) {   //UN SEUL FICHIER SPECIFIE
        input_file = safeOpen(argv[optind], "r");
        }
    else if (optind==argc) //PAS DE FICHIER SPECIFIE, ex : gunzip -c file.fastq.gz | program -o *.fasta -r *.fasta
        {
        fprintf(stderr,"Reading fastq from STDIN\n");
        input_file = stdin; //et avec gzip ??? :-)
          
        }
    else
        {
        fprintf(stderr,"%s : Too many input files specified\n",argv[0]);
        return EXIT_FAILURE;
        }

    //OUTPUT
    if (filename_out == NULL) {
        fprintf(stderr,"Writing output to STDOUT\n");
        output_file = stdout;
        }
    else {
        output_file = safeOpen(filename_out,"w");
        }


    

// Part1 : a table (tLen) containing the number of read for each size 

    tLen = safeCalloc(table_size, sizeof(int));   
    while ((c=fgetc(input_file))!=EOF) 
        {
        if (c=='\n')
            {   
            if (nLine%4==1) 
                {
                if (length>=table_size)
                    {
                    tLen = safeRealloc(tLen, (length+1)*sizeof(int));
                    while(table_size<=length) {tLen[table_size]=0; table_size++;}; 
                    table_size = length + 1;
                    }
                tLen[length]++;
                if (length>max_length) { max_length=length ; printf("The max size of read is : %d bases\n", max_length);}
                }
            length=0;
            nLine++;                // Nb of characters for each line
            } 
         else 
            {
            length++;
            }
        }
    // To check the script    
    for (i=0; i<=max_length; i++) 
        {
        fprintf(stderr, "tLen  : size of reads : %d\t number of reads %d\n", i, tLen[i]);
        }
    fclose (input_file);
    fprintf(stderr, "Part1 successfull\n\n");


// Part 2 :  ref is a sequence structure containing informations (bases and length) of the genome of reference
// Example : a genome of 8 bases AATTCCGG, ref.bases[0]='A', seqLength = 8
    
    
    ref.bases = safeCalloc(position, sizeof(char));
    while ((b=fgetc(reference_file))!=EOF)
        {  
        if (b=='>')  
            {
            while ((b=fgetc(reference_file))!=EOF && b!='\n') { continue; }
            }
        else {
            if(isspace(b)) continue;
            ref.bases = safeRealloc(ref.bases, (position+1)*sizeof(char));
            ref.bases[position]=b;
            position++;
            }
        }
        ref.bases[position]='\0';
        ref.seqLength=position-1;
        fprintf(stderr,"Size of the sequence of reference : %d\n", ref.seqLength);
        fprintf(stderr, "Part2 successfull\n\n");
    
    fclose(reference_file);

// Part 3 : Create a fasta file containing the InSilico reads 
// Example : >Read1
//           ATCCc (lowerletter for insertion or substitution)


    srand(time(NULL));
    for (i=1; i<=max_length; i++)        // for each size of read
        {                  
        if(tLen[i]!=0)                      // if read exist
            {
            for (j=0; j<tLen[i]; j++)           // for each read
                {  
                start_position = 1 + rand()%ref.seqLength;      // select the beginning in the reference genome
                fprintf(output_file, ">m160129_165300_42263_c100880132550000001823194304021670_s1_X0/%d/0_%d\n", readNb, i);  // For FALCON compatibility
                
                for (k=1; k<=i;k++)                                  // for each base 
                    {
                    if ((start_position+k)>ref.seqLength) {start_position=start_position-ref.seqLength;}    
                    prob = (rand()%101);                               // select the probability of mutation [1;100]
                    
                    // INSERTION : 11% [0-11]
                    if (prob<=11) 
                        {            
                        ins = "ATCG"[rand()%4];
                        fprintf(output_file, "%c%c", toupper(ref.bases[start_position+k]), tolower(ins));
                        }
                                              
                    // SUBSTITUTION : 1% [12]
                    else if (prob<=12) 
                        {
                        switch(toupper(ref.bases[start_position+k])) 
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
                        fprintf(output_file, "%c", toupper(ref.bases[start_position+k]));
                        }
                    }
                fprintf(output_file, "\n");
                readNb++;                                   // Incrémente le nombre de reads Insilico crées
                }
            }
            else {continue;}
        }
    fclose (output_file);
    fprintf(stderr, "Part3 successfull\n\n");
    free (tLen);
    time(&end);

    printf("Elapsed time %f sec.\n", difftime(end, begin));

    return 0;
    }

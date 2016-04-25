#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>   // C Error Codes

#define CASE_MUT(B, OPT) case B : mute = OPT[rand()%3]; fprintf(output,"%c", tolower(mute)); break

struct sequence { char bases[200000]; unsigned int seqLength;};  // je n'arrive pas à le mettre en dynamique 


void safeCalloc (int* tLen)
    {
    if (tLen == NULL) {printf("\ntLen : %s", strerror(12)); exit(EXIT_FAILURE);};  //Out of memory
    }
    
int* safeRealloc(int* tLen, unsigned int length)
    {
    int *temp = realloc(tLen, ((length+1) * sizeof(int)));      
    if (temp == NULL) {printf("\ntLen : %s", strerror(12)); exit(EXIT_FAILURE);};  //Out of memory
    return temp;
    }
        
void checkFile(FILE* fic)   
    {
    if (fic == NULL) 
        {
        printf ("Error opening file : %s\n", strerror(2));  //No such file or directory
        exit(EXIT_FAILURE);;
        }
    }
    
void checkArg(int argc)   
    {
        if (argc != 4)
        {
        printf ("The number of argument is not correct : %s\n", strerror(22));
        exit(EXIT_FAILURE);
        } 
    }



/*struct sequence* safeRealloc (struct bases* ref, int structCount)
    {
    struct bases* updateRef;    
    updateRef = realloc(ref, (structCount+1)* sizeof(struct bases));
    if (updateRef == NULL) {printf("%s", strerror(12)); exit(EXIT_FAILURE);};  //Out of memory
    return updateRef;
    }
*/

int main (int argc, char** argv) {
    
    time_t begin, end; 

     
// Part1 : a table (tLen) containing the number of read for each size 
          
    FILE* input = fopen(argv[1], "r");        // argv[1] = fastq file containing the initial PBReads  
    
    int c = 0;                                
    unsigned int length = 0;                  // Read Length
    int nLine = 0;                            // Number of line (1==sequence)  
    int* tLen = NULL;
    unsigned int table_size = 1;              
    unsigned int max_length=0;
    unsigned int i;

// Part 2 :  ref is a sequence structure containing informations (bases and length) of the reference
// Example : a genome of 8 bases AATTCCGG, ref.bases[0]='A', seqLength = 8

    FILE* reference_in = fopen (argv[2], "r");         // fasta file containing the sequence of reference

    int b = 0;                     
    struct sequence ref;
    int unsigned position = 1;

// Part 3 : Create a fasta file containing the InSilico reads 
// Example : >Read1
//           ATCCc (lowerletter for insertion or substitution)

    FILE* output = fopen (argv[3], "w+");                // Will contain the PB InSilico reads          
    int j;
    unsigned int k=0;
    int start_position = 0;
    int prob;
    char ins;               // the inserted base
    char mute;              // the mutated base
    int readNb = 1;         // sequence generated (print in the header of the output)


// BEGIN

    time(&begin);
    checkArg(argc);

// Part1 : a table (tLen) containing the number of read for each size 

    // Open the fastq file containing the reads from PacBio sequencer
    checkFile(input);
    tLen = calloc(table_size, sizeof(int));   
    safeCalloc(tLen);

    while ((c=fgetc(input))!=EOF) 
        {
        if (c=='\n')
            {   
            if (nLine%4==1) 
                {
                if (length>=table_size)
                    {
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
        printf("tLen  : taille du read  : %d\t nombre de read %d\n", i, tLen[i]);
        }
    fclose (input);
    printf("Part1 successfull\n\n");


// Part 2 :  ref is a sequence structure containing informations (bases and length) of the reference
// Example : a genome of 8 bases AATTCCGG, ref.bases[0]='A', seqLength = 8
    
    checkFile (reference_in);

    while ((b=fgetc(reference_in))!=EOF)
        {  
        if (b=='>')  
            {
            while ((b=fgetc(reference_in))!=EOF && b!='\n') { continue; }
            }
        else {
            if(isspace(b)) continue;
            ref.bases[position]=b;
            printf("char %c, position %d\n", ref.bases[position], position);
            position++;
            }
        }
        ref.bases[position]='\0';
        printf("char %c, position %d\n", ref.bases[position], position);
        ref.seqLength=position-1;
        printf("Size of the sequence of reference : %d\n", ref.seqLength);
        printf("Part2 successfull\n\n");

// Part 3 : Create a fasta file containing the InSilico reads 
// Example : >Read1
//           ATCCc (lowerletter for insertion or substitution)
        
    checkFile (output);
    srand(time(NULL));
    for (i=0; i<=max_length; i++)        // for each size of read
        {                  
        if(tLen[i]!=0)                      // if read exist
            {
            for (j=0; j<tLen[i]; j++)           // for each read
                {
                start_position = 1 + rand()%ref.seqLength;      // select the beginning in the reference genome
                printf("Ref seq length : %d\tstart position %d\t Base : %c\n", ref.seqLength,start_position,ref.bases[start_position]);
                fprintf(output, ">InSilico read %d\n", readNb); 
                
                for (k=0; k<i;k++)                                  // for each base 
                    {
                    if ((start_position + k)>ref.seqLength) {start_position=start_position-ref.seqLength;}    
                    prob = (rand()%101);                               // select the probability of mutation [1;100]
                    printf("Base : %c\t probabilité de mutation %d\n", ref.bases[start_position+k], prob);
                    
                    // INSERTION : 11% [0-11]
                    if (prob<=11) 
                        {            
                        ins = "ATCG"[rand()%4];
                        fprintf(output, "%c%c", toupper(ref.bases[start_position+k]), tolower(ins));
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
                        fprintf(output, "%c", toupper(ref.bases[start_position+k]));
                        }
                    }
                fprintf(output, "\n");
                readNb++;                                   // Incrémente le nombre de reads Insilico crées
                }
            }
        }
    fclose (output);
    printf("Part3 successfull\n\n");
    free (tLen);
    time(&end);

    printf("Elapsed time %f sec.\nPB reads stored in %s \nEND of the program\n", difftime(end, begin),argv[3]);

    return 0;
    }

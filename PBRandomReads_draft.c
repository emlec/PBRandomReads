#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

#define CASE_MUT(B, OPT) case B : mute = OPT[rand()%3]; fprintf(output,"%c", tolower(mute)); break

// TODO : modularité


int main () {
    
    time_t begin, end; 
    
// Part1 : a table (tLen) containing the number of read for each size 
// Example : 2 reads of 100 bases and 1 read of 10 bases : tLen[100]=2 and tLen[10]=1

    char inputFileName[30];        // fastq file containing the initial PBReads     
    FILE* input = NULL; 
    
    char c = 0;                     // lecture fastq file
    unsigned int length = 0;                 // Read Length
    int nLine = 0;                  // Number of line (1==sequence)  
    int* tLen = NULL;
    unsigned int table_size = 1;           // taille pour stocker les séquences de 0 à 100. 
    unsigned int max_length=0;


// Part 2 :  a table (tRef) containing  the sequence of reference
// Example : a genome of 8 bases AATTCCGG, tRef[0]='A', tRef[1]='A', ..., tRef[8]='\0' 

    char referenceFileName[30];    // fasta file containing the sequence of reference
    FILE* reference = NULL; 

    char b = 0;                     // lecture fasta reference file
    char* tRef = NULL;              
    unsigned int position = 0;

// Part 3 : Create a fasta file containing the InSilico reads 
// Example : >Read1
//           ATCCc (lowerletter for insertion or substitution)


    char outputFileName[30];       // Will contain the PB InSilico reads           
    FILE* output = NULL;
    unsigned int i;
    int j;
    unsigned int k=0;
    int start_position = 0;
    int prob;
    char ins;               // the inserted base
    char mute;              // the mutated base
    int readNb = 1;



// BEGIN

    time(&begin);       // TODO décaler une fois modulaire


// Part1 : a table (tLen) containing the number of read for each size 
// Example : 2 reads of 100 bases and 1 read of 10 bases : tLen[100]=2 and tLen[10]=1
    
    
    // Open the fastq file containing the reads from PacBio sequencer
    printf ("Input file name :");
    scanf ("%s", inputFileName);
    input = fopen (inputFileName, "r");
    if (input == NULL) 
        {
        printf ("Cannot open the file %s\n", inputFileName);
        return 0;
        }
    
    // tLen[length] = nb of reads
    tLen = calloc(table_size, sizeof(int));    
    if (tLen == NULL) {printf("\ntLen : Failed to allocate memory"); EXIT_FAILURE;};
    
    while ((c=fgetc(input))!=EOF) 
        {
        // At the end of a line
        if (c=='\n')
            {   
            // If line = sequence
            if (nLine%4==1) 
                {
                printf("table_size %d\n", table_size);
                if (length>=table_size)
                    {
                    // realloc memory
                    tLen = realloc(tLen, ((length-table_size+1) * sizeof(int)));   //TODO passer par un pointeur intermédiaire pour éviter fuite de mémoire
                    if (tLen == NULL) {printf("\ntLen : Failed to allocate memory"); EXIT_FAILURE;};
                    table_size = length + 1;
                    }
                tLen[length]++;
                printf(" length : %d, nb of reads : %d\n", length, tLen[length]);
                if (length>max_length) { max_length=length;}
                }
            
            length=0;
            // Je compte le nombre de caractère de la linge quelle quelle soit
            nLine++;
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
    printf ("max length = %d\n", max_length);   // Pour la boucle de l'étape 3 
    fclose (input);
    printf("Part1 successfull\n");


// Part 2 :  a table (tRef) containing  the sequence of reference
// Example : a genome of 8 bases AATTCCGG, tRef[0]='A', tRef[1]='A', ..., tRef[8]='\0' 

    // Open the reference file 
    printf ("Reference file name :");
    scanf ("%s", referenceFileName);
    reference = fopen(referenceFileName, "r");
    if (reference == NULL) 
        {
        printf ("Cannot open the file %s\n", referenceFileName);
        return 0;
        }

    // tRef[position]=base     
    tRef = malloc(position * sizeof(char));
    if (tRef == NULL) {printf("\ntRef : Failed to allocate memory"); EXIT_FAILURE;};
  
    while ((b=fgetc(reference))!=EOF)
        {  
        if (b=='>')   // On passe le header 
            {
            while ((b=fgetc(reference))!=EOF && b!='\n')   // Pourquoi pas uniquement while (c!='\n') ?
                {
                continue;
                }
            }
        else {
            if(isspace(b)) continue;
            tRef = realloc(tRef, ((position+1) * sizeof(char)));
            tRef[position]=b;
            position++;
            }
        }
         
    tRef[position]='\0';
    
    // To check
    for (i=0; i<=position ;i++)
        {
            printf("position %d, base %c\n", i, tRef[i]);
        }
    fclose (reference);
    printf("Part2 successfull\n");
    
    

// Part 3 : Create a fasta file containing the InSilico reads 
// Example : >Read1
//           ATCCc (lowerletter for insertion or substitution)


    // Create the output file that will contain the in-silico reads
    printf ("Output file name :");
    scanf ("%s", outputFileName);
    output = fopen (outputFileName, "w+");
    
    srand(time(NULL));      

    for (i=0; i<=max_length; i++)        // parcours le tableau des tailles 
        {                  
        if(tLen[i]!=0)                  // s'il y a des reads
            {
            for (j=0; j<tLen[i]; j++)        // pour chaque read de taille n
                {            
                printf("position %d", position);
                start_position = rand()%position;         // on choisit un début sur la référence
                printf("start position %d\n", start_position);
                fprintf(output, ">InSilico read %d\n", readNb); 
                
                for (k=0; k<i;k++)                       // parcours le nb de base de la séquence
                    {                       
                    prob = (rand()%101);                    // probabilité de mutation [0;100]
                    printf("probabilité de mutation %d\n", prob);
                    
                    // INSERTION : 11% [0-11]
                    if (prob<=11) 
                        {            
                        ins = "ATCG"[rand()%4];
                        fprintf(output, "%c%c", toupper(tRef[start_position]), tolower(ins));
                        start_position+=2;
                        }
                        
                    // SUBSTITUTION : 1% [12]
                    else if (prob<=12) 
                        {
                        switch(toupper(tRef[start_position])) 
                            {
                            CASE_MUT('A', "TCG");           // eq to case 'A': printf("%c","TGC"[rand()%3]);break;
                            CASE_MUT('T', "ACG");
                            CASE_MUT('C', "TAG");
                            CASE_MUT('G', "TCA");
                            };
                        start_position+=1;
                        }
                    // CORRECT : 84% [13-96]
                    else if (prob<=96) 
                        {
                        fprintf(output, "%c", toupper(tRef[start_position]));
                        start_position+=1;
                        }
                    }
                fprintf(output, "\n");
                readNb++;                                   // Incrémente le nombre de reads Insilico crées
                }
            }
        }
    fclose (output);
    free (tLen);
    free (tRef);
    
    time(&end);
    printf("Elapsed time %f sec.\nPB reads stored in %s \nEND of the program\n", difftime(end, begin),outputFileName);
    
    return 0;
    }

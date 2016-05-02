#include "fasta.h"
#include "utils.h"

void SequenceFree(SequencePtr seq) {
if(seq==NULL) return;
free(seq->name);
free(seq->bases);
free(seq);
}



SequencePtr readOneSequenceFromFile(FILE* in) {
SequencePtr seq = (SequencePtr)safeCalloc(1,sizeof(Sequence));
seq->name = safeCalloc(1, sizeof(char));
seq->name[0] = 0;
seq->bases =  safeCalloc(1, sizeof(char));
seq->bases[0] = 0;


    while ((b=fgetc(reference_file))!=EOF)
        {  
        if (b=='>')  
            {
            if(seq->name[0]!=0) {
            fprintf(stderr,"ATTENTION VOTRE FASTA CONTIENT DEUX SEQUCNE");
            exit(EXIT_FAILURES);
            }
            while ((b=fgetc(reference_file))!=EOF && b!='\n') { continue; }
            }
        else {
            if(isspace(b)) continue;
            seq->bases = safeRealloc(seq->bases, (position+1)*sizeof(char));
            ref.bases[position]=b;
            position++;
            }
        }
        ref.bases[position]='\0';
        ref.seqLength=position-1;
        fprintf(stderr,"Size of the sequence of reference : %d\n", ref.seqLength);
        fprintf(stderr, "Part2 successfull\n\n");



return seq;
}







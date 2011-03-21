//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#pragma hdrstop

//---------------------------------------------------------------------------

#define READ_LENGTH 100
//#define INSERT_LENGTH 200
#pragma argsused
int main(int argc, char* argv[])
{
char *Genome;
FILE *in;
FILE *out1;
char ch;
long i;
char ch1, ch2;
long j;
int t;
int INSERT_LENGTH;

i=0;
Genome = (char *) malloc(6000000);
if ((in = fopen("MG1655-K12_cut.fasta", "r"))
           == NULL) return 1;
if ((out1 = fopen("MG1655-K12_emulpaired.fasta", "w"))
           == NULL) return 1;

Genome[i] = fgetc(in);
t=0;
while (t<5) {
  if (Genome[i]>20) {t=0;i++;}
  else t++;
  Genome[i] = fgetc(in);
}
//i=1000;
int genome_length =i;
INSERT_LENGTH = 200;
for(j=0; j<i-(2*READ_LENGTH+INSERT_LENGTH); j++) {
 // printf(out1, "%i ",INSERT_LENGTH);
  ch1 = Genome[j+READ_LENGTH];
  ch2 = Genome[j+2*READ_LENGTH+INSERT_LENGTH];
  Genome[j+READ_LENGTH]=0;
  Genome[j+2*READ_LENGTH+INSERT_LENGTH]=0;
  fprintf(out1, "%s ",Genome+j);
  fprintf(out1, "%s ",Genome+j+READ_LENGTH+INSERT_LENGTH);
  Genome[j+READ_LENGTH]=ch1;
  Genome[j+2*READ_LENGTH+INSERT_LENGTH]=ch2;

/*  switch(rand()%11) {
    case 0:
    case 1:
    case 2: INSERT_LENGTH = 100; break;
    case 3:
    case 4: INSERT_LENGTH = 99; break;
    case 5:
    case 6: INSERT_LENGTH = 101; break;
    case 7: INSERT_LENGTH = 102; break;
    case 8: INSERT_LENGTH = 98; break;
    case 9: INSERT_LENGTH = 97; break;
    case 10:INSERT_LENGTH = 103; break;
  }
*/
}

fclose(out1);


        return 0;
}
//---------------------------------------------------------------------------

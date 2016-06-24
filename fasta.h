/* fasta.h
 * Declarations for simple FASTA i/o library
 * SRE, Sun Sep  8 05:37:38 2002 [AA2721, transatlantic]
 * CVS $Id: fasta.h,v 1.1 2003/10/05 18:43:39 eddy Exp $
 */

#include <stdio.h>

#define FASTA_MAXLINE 512
uint32_t counter;

typedef struct fastafile {
  FILE *fp;
  char  buffer[FASTA_MAXLINE];
} fasta_file_t;

extern fasta_file_t *open_fasta(char *seqfile);
extern int        read_fasta(fasta_file_t *ffp, char **ret_seq, char **ret_name, int *ret_L, uint32_t* pos_seq);
extern void       close_fasta(fasta_file_t *ffp);

/*
 * fasta.c
 *
 *  Created on: 25/02/2016
 *      Author: galvez
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>

#include "fasta.h"
/*
 *
 *
 * Args:
 *           seqfile   - name of a FASTA file to open.
 *           seq       - RETURN: one sequence
 *           name      - RETURN: name of the sequence
 *           seqlen    - RETURN: length of the sequence in residues
 *           ffp       - ptr to a FASTAFILE object.
 *
 * Returns:
 *           OpenFASTA() returns a FASTAFILE pointer, or NULL on failure (for
 *           instance, if the file doesn't exist, or isn't readable).
 *
 *           ReadFASTA() returns 1 on success, or a 0 if there are no
 *           more sequences to read in the file.
 *
 *           CloseFASTA() "always succeeds" and returns void.
 */


fasta_file_t * open_fasta(char *seqfile)
{
  fasta_file_t *ffp;

  ffp = malloc(sizeof(fasta_file_t));
  ffp->fp = fopen(seqfile, "r");              /* Assume seqfile exists & readable!   */
  if (ffp->fp == NULL) { free(ffp); return NULL; }
  if ((fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) == NULL)
    { free(ffp); return NULL; }
  return ffp;
}

int read_fasta(fasta_file_t *ffp, char **ret_seq, char **ret_name, int *ret_L, uint32_t *pos_seq)
{
  char *s;
  char *name;
  char *seq;
  int   n;
  int   nalloc;

  /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
   */
  if (ffp->buffer[0] != '>') return 0;

  /* Parse out the name: the first non-whitespace token after the >
   */
  *pos_seq = counter;
  s  = strtok(ffp->buffer+1, " \t\n");
  name = malloc(sizeof(char) * (strlen(s)+1));
  strcpy(name, s);
  counter++;

  /* Everything else 'til the next descline is the sequence.
   * Note the idiom for dynamic reallocation of seq as we
   * read more characters, so we don't have to assume a maximum
   * sequence length.
   */
  seq = malloc(sizeof(char) * 128);     /* allocate seq in blocks of 128 residues */
  nalloc = 128;
  n = 0;

  while (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp))
    {

      if (ffp->buffer[0] == '>') break;	/* a-ha, we've reached the next descline */
      counter++;
      for (s = ffp->buffer; *s != '\0'; s++)
	{
	  if (! isalpha(*s)) continue;  /* accept any alphabetic character */

	  seq[n] = *s;                  /* store the character, bump length n */
	  n++;
	  if (nalloc == n)	        /* are we out of room in seq? if so, expand */
	    {			        /* (remember, need space for the final '\0')*/
	      nalloc += 128;
	      seq = realloc(seq, sizeof(char) * nalloc);
	    }
	}
    }
  seq[n] = '\0';

  *ret_name = name;
  *ret_seq  = seq;
  *ret_L    = n;
  return 1;
}

void
close_fasta(fasta_file_t *ffp)
{
  fclose(ffp->fp);
  free(ffp);
}


#ifdef TEST_FASTA
/* Test the fasta parsing API.
 *  to compile:  gcc -o test -DTEST_FASTA_STUFF -Wall -g fasta.c
 *  to run:      ./test myseqs.fa
 */
int
main(int argc, char **argv)
{
  fasta_file_t *ffp;
  char *seq;
  char *name;
  int   L;
				/* argv[1] is the name of a FASTA file */
  ffp = open_fasta(argv[1]);
  while (read_fasta(ffp, &seq, &name, &L))
    {
      printf(">%s\n", name);
      printf("%s\n",  seq);

      free(seq);
      free(name);
    }
  close_fasta(ffp);
  exit(0);
}
#endif /*TEST_FASTA*/

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>

void merge_scores(int * scores, char ** titles, unsigned long int size);

void sort_scores(int * scores, char ** titles, unsigned long int size);

double getTimestamp();

void merge_sequences(char ** sequences, char ** titles, unsigned short int * sequences_lengths, unsigned long int size);
void sort_sequences (char ** sequences, char ** titles, unsigned short int * sequences_lengths, unsigned long int size);

void sort_score_simple(int* score, int* index, unsigned long int size);
void merge_score_simple(int*score, int*index,unsigned long int size);

#endif

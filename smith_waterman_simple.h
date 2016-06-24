/*
 * smith_waterman_simple.h
 *
 *  Created on: 12/05/2016
 *      Author: galvez
 */

#ifndef SMITH_WATERMAN_SIMPLE_H_
#define SMITH_WATERMAN_SIMPLE_H_

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GAP -3
#define MATCH 3
#define MISMATCH -1

typedef enum {
	true, false
} bool;

typedef struct {
	char *a;
	unsigned int alen;
	char *b;
	unsigned int blen;
} seq_pair;
typedef seq_pair *seq_pair_t;

typedef struct {
	int score;
	unsigned int prev[2];
} entry;
typedef entry *entry_t;

typedef struct {
	unsigned int m;
	unsigned int n;
	entry_t **mat;
} matrix;
typedef matrix *matrix_t;

 char* reverse(char *str);

 int get_score(seq_pair_t problem, matrix_t S);

 matrix_t create_matrix(unsigned int m, unsigned int n);

 seq_pair_t traceback(seq_pair_t problem, matrix_t S, bool local);

void destroy_matrix(matrix_t S);

void destroy_seq_pair(seq_pair_t pair);

 int smith_waterman(seq_pair_t problem);




#endif /* SMITH_WATERMAN_SIMPLE_H_ */

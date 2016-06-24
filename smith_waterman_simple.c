/*
 * smith_waterman_simple.c
 *
 *  Created on: 12/05/2016
 *      Author: Felipe Sulser
 */
#include "smith_waterman_simple.h"

/* reverse a string in place, return str */
char* reverse(char *str) {
	char *left = str;
	char *right = left + strlen(str) - 1;
	char tmp;

	while (left < right) {
		tmp = *left;
		*(left++) = *right;
		*(right--) = tmp;
	}

	return str;
}
 seq_pair_t traceback(seq_pair_t problem, matrix_t S, bool local) {
	seq_pair_t result = malloc(sizeof(seq_pair));
	unsigned int i = S->m - 1;
	unsigned int j = S->n - 1;
	unsigned int k = 0;
	char c[S->m + S->n + 1];
	char d[S->m + S->n + 1];

	memset(c, '\0', sizeof(c));
	memset(d, '\0', sizeof(d));

	if (local == true) {
		unsigned int l, m;
		double max = FLT_MIN;

		for (l = 0; l < S->m; l++) {
			for (m = 0; m < S->m; m++) {
				if (S->mat[l][m]->score > max) {
					i = l;
					j = m;
				}
			}
		}
	}

	if (S->mat[i][j]->prev[0] != 0 && S->mat[i][j]->prev[1] != 0) {
		while (i > 0 || j > 0) {
			unsigned int new_i = S->mat[i][j]->prev[0];
			unsigned int new_j = S->mat[i][j]->prev[1];

			if (new_i < i)
				*(c + k) = *(problem->a + i - 1);
			else
				*(c + k) = '-';

			if (new_j < j)
				*(d + k) = *(problem->b + j - 1);
			else
				*(d + k) = '-';

			k++;

			i = new_i;
			j = new_j;
		}
	}

	result->a = malloc(sizeof(char) * k + 1);
	result->b = malloc(sizeof(char) * k + 1);

	memset(result->a, '\0', sizeof(result->a));
	memset(result->b, '\0', sizeof(result->b));

	reverse(c);
	reverse(d);

	strcpy(result->a, c);
	strcpy(result->b, d);

	result->alen = k;
	result->blen = k;

	return result;
}


// works globally
int get_score(seq_pair_t problem, matrix_t S) {
	entry_t** val = S->mat;
	return val[problem->alen - 1][problem->blen - 1]->score;
}

matrix_t create_matrix(unsigned int m, unsigned int n) {
	matrix_t S = malloc(sizeof(matrix));
	unsigned int i, j;

	S->m = m;
	S->n = n;

	S->mat = malloc(sizeof(entry_t) * m * n);

	for (i = 0; i < m; i++) {
		S->mat[i] = malloc(sizeof(entry_t) * n);
	}

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			S->mat[i][j] = malloc(sizeof(entry));
		}
	}

	return S;
}

void destroy_matrix(matrix_t S) {
	unsigned int i, j;

	for (i = 0; i < S->m; i++) {
		for (j = 0; j < S->n; j++) {
			free(S->mat[i][j]);
		}
		free(S->mat[i]);
	}

	free(S);

	return;
}

void destroy_seq_pair(seq_pair_t pair) {
	free(pair->a);
	free(pair->b);

	free(pair);

	return;
}

int smith_waterman(seq_pair_t problem) {
	unsigned int m = problem->alen + 1;
	unsigned int n = problem->blen + 1;
	matrix_t S = create_matrix(m, n);
	int result;
	unsigned int i, j, k, l;

	S->mat[0][0]->score = 0;
	S->mat[0][0]->prev[0] = 0;
	S->mat[0][0]->prev[1] = 0;
#pragma omp parallel for
	for (i = 1; i <= problem->alen; i++) {
		S->mat[i][0]->score = 0.0;
		S->mat[i][0]->prev[0] = i - 1;
		S->mat[i][0]->prev[1] = 0;
	}
#pragma omp parallel for
	for (j = 1; j <= problem->blen; j++) {
		S->mat[0][j]->score = 0.0;
		S->mat[0][j]->prev[0] = 0;
		S->mat[0][j]->prev[1] = j - 1;
	}
#pragma omp parallel for
	for (i = 1; i <= problem->alen; i++) {
		for (j = 1; j <= problem->blen; j++) {
			int nw_score =(*(problem->a + (i - 1))== *(problem->b + (j - 1))?MATCH : MISMATCH);

			S->mat[i][j]->score = 0;
			S->mat[i][j]->prev[0] = 0;
			S->mat[i][j]->prev[1] = 0;

			for (k = 0; k <= 1; k++) {
				for (l = 0; l <= 1; l++) {
					int val;

					if (k == 0 && l == 0) {
						continue;
					} else if (k > 0 && l > 0) {
						val = nw_score;
					} else if (k > 0 || l > 0) {
						if ((i == problem->alen && k == 0)
								|| (j == problem->blen && l == 0))
							val = 0.0;
						else
							val = GAP;
					} else {
						// do nothing..
					}
					val += S->mat[i - k][j - l]->score;

					if (val > S->mat[i][j]->score) {
						S->mat[i][j]->score = val;
						S->mat[i][j]->prev[0] = i - k;
						S->mat[i][j]->prev[1] = j - l;

					}
				}
			}
			//printf("%d\t", (int) S->mat[i][j]->score);

		}
		//printf("\n");
	}

	result = get_score(problem, S);

	destroy_matrix(S);

	return result;
}

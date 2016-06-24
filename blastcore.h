#ifndef MICSEARCH_H_INCLUDED
#define MICSEARCH_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <immintrin.h>
#include "scoremat.h"
#include "sort.h"
#include "macros.h"



void smith_waterman_vectorized (Query* query,
		Chunked_Database * db,
		char * submat,
		int open_gap,
		int extend_gap,
		int num_threads,
		int * scores,
		double * workTime,
		uint16_t query_length_threshold);


#endif

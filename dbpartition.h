/*
 * dbpartition.h
 *
 *  Created on: 25/04/2016
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>
#include "macros.h"
#include "sort.h"
#include "preprocess.h"


void divide_db (char * sequences_filename,
		uint64_t max_chunk_size,
		uint16_t * sequences_db_max_length,
		int * max_title_length,
		Database ** non_chunked,
		Chunked_Database** chunked);
// Load DB headers
void load_database_headers (char * sequences_filename,
		uint64_t sequences_count,
		int max_title_length,
		char *** ptr_sequences_db_headers);

void load_query_sequences(char * queries_filename,
		struct  query_data ** query_ret);

#ifndef ARGUMENTS_H_INCLUDED
#define ARGUMENTS_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include "scoremat.h"
#include <stdint.h>

#define DEBUG 2
#undef DEBUG

#define MAX_CHUNK_SIZE 100000000 // 128Mb

#define PREPROCESS_THREADS 4

#define HIT_THRESHOLD 1
#define NUM_THREADS 200
#define OPEN_GAP 10
#define EXTEND_GAP 3
#define TOP 10
#define QUERY_LENGTH_THRESHOLD 567 // max len


#define BUFFER_SIZE 1000 // for IO operations

#define UNIFY 'Z'+1 //out of aminoacid range
#define UNIFY_NUMERIC 24



#define BLOCK_SIZE 256
#define INTRINSIC_LEN 368
#define VECTOR_LENGTH 16
#define UNROLL 6

//Macro for Ceiling operator
#define CEILING(x,y) (((x) + (y) - 1) / (y))
#define min(a, b) (((a) < (b)) ? (a) : (b))

#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)
#define ALLOC_AND_FREE alloc_if(1) free_if(1)


extern char blosum45[];

typedef struct query_data{
	char* q_seq;
	char** query_headers;
	uint16_t* query_seq_len;
	uint64_t query_sequences_count;
	uint64_t lin_len_total;
	uint32_t* query_seq_disp;
}Query;
typedef struct chunked_db_data{
	uint64_t vect_sequences_db_count;
	char ** chunk_b;
	uint32_t chunk_count;
	uint32_t * chunk_vect_sequences_db_count;
	uint16_t ** chunk_n;
	uint32_t ** chunk_b_disp;
	uint64_t * chunk_vD;
}Chunked_Database;

typedef struct db_data{
	uint64_t  sequences_count;
	uint64_t  D;
	uint64_t  vD;
	uint16_t* sequences_lengths;
	uint32_t * sequences_disp;
	char* real_seq;
}Database;

#endif

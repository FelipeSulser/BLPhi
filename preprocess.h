#ifndef DB_H_INCLUDED
#define DB_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>
#include "macros.h"
#include "sort.h"
#include <stdint.h>


void convert_db (char * input_filename, char * out_filename) ;


void create_shifted_copy(char* seq,
		uint64_t D,
		uint64_t seq_count_db,
		uint32_t* shift,
		uint16_t* len,
		char** ptr_shifted_seq,
		int** ptr_extension,
		int** ptr_shift,
		uint64_t *total_linear
		);

void create_shifted_copy_nopad(char* seq,
		uint64_t D,
		uint64_t seq_count_db,
		uint32_t* shift,
		uint16_t* len,
		char** ptr_shifted_seq,
		int** ptr_extension,
		int** ptr_shift,
		uint64_t *total_linear
		);

#endif

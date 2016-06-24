/*
 * heuristic.h
 *
 *  Created on: 25/04/2016
 *      Author: Felipe Sulser Larraz
 */

#ifndef HEURISTIC_H_
#define HEURISTIC_H_
#include <stdint.h>
#include "macros.h"
#include <immintrin.h>
#include "sort.h"

void filter_sequences_single_db (char * query_sequences,
		char** query_headers,
		char** seq_headers,
		uint16_t * query_sequences_lengths,
		uint32_t query_sequences_count,
		uint64_t Q,
		uint32_t * query_disp,
		char* seq_db,
		uint64_t D,
		uint16_t* seq_len_db,
		uint64_t seq_count_db,
		uint32_t* seq_disp_db,
		int ** ptr_assigned,
		double* elapsed_time,
		uint64_t* total_done
		);
__declspec(target(mic))int popcount_hamming(int i);
/**
 * Only for test purposes, comparing hits and doing analysis.
 * Speed is not wanted here,therefore it is not optimized
 */
void filter_sequences_single_db_8 (char * query_sequences,
		char** query_headers,
		char** seq_headers,
		uint16_t * query_sequences_lengths,
		uint32_t query_sequences_count,
		uint64_t Q,
		uint32_t * query_disp,
		char* seq_db,
		uint16_t* seq_len_db,
		uint64_t seq_count_db,
		uint32_t* seq_disp_db,
		int ** ptr_assigned
		);
void filter_sequences_multiple_db (char * query_sequences,
		char** query_headers,
		uint16_t * query_sequences_lengths,
		uint32_t query_sequences_count,
		uint64_t Q,
		uint32_t * query_disp,
		uint64_t vect_sequences_db_count,
		char ** chunk_b,
		uint32_t chunk_count,
		uint32_t * chunk_vect_sequences_db_count,
		uint16_t ** chunk_n,
		uint32_t ** chunk_b_disp,
		uint64_t * chunk_vD,
		//novelty
		int *** ptr_assigned
		);
void filter_sequences_single_db_vectorized (char * query_sequences,
		char** query_headers,
		char** seq_headers,
		uint16_t * query_sequences_lengths,
		uint32_t query_sequences_count,
		uint64_t Q,
		uint32_t * query_disp,
		char* seq_db,
		uint64_t D,
		uint16_t* seq_len_db,
		uint64_t seq_count_db,
		uint32_t* seq_disp_db,
		int ** ptr_assigned,
		double* elapsed_time,
		uint64_t* total_done
		);

void filter_sequences_single_db_reduced (
		char * red_query_sequences,
		char** query_headers,
		char** seq_headers,
		uint16_t * red_query_sequences_lengths,
		uint32_t red_query_sequences_count,
		uint64_t red_Q,
		uint32_t * red_query_disp,
		char* red_seq_db,
		uint64_t red_D,
		uint16_t* red_seq_len_db,
		uint64_t red_seq_count_db,
		uint32_t* red_seq_disp_db,
		int ** ptr_assigned,
		double* elapsed_time,
		uint64_t* total_done
		);

void filter_sequences_shifted (char * query_sequences,
		char** query_headers,
		char** seq_headers,
		uint16_t * query_sequences_lengths,
		uint32_t query_sequences_count,
		uint64_t Q,
		uint32_t * query_disp,
		char* seq_db,
		uint64_t D,
		uint16_t* seq_len_db,
		uint64_t seq_count_db,
		uint32_t* seq_disp_db,
		int ** ptr_assigned,
		double* elapsed_time,
		uint64_t* total_done,
		char* shifted_copy,
		int* extension,
		int* shifted_shift,
		uint64_t shifted_linear
		);

void filter_sequences_shifted_vectorized (char * query_sequences,
		char** query_headers,
		char** seq_headers,
		uint16_t * query_sequences_lengths,
		uint32_t query_sequences_count,
		uint64_t Q,
		uint32_t * query_disp,
		char* seq_db,
		uint64_t D,
		uint16_t* seq_len_db,
		uint64_t seq_count_db,
		uint32_t* seq_disp_db,
		int ** ptr_assigned,
		double* elapsed_time,
		uint64_t* total_done,
		char* shifted_copy,
		int* extension,
		int* shifted_shift,
		uint64_t shifted_linear
		);

void filter_sequences_shifted_threshold_vectorized (char * query_sequences,
		char** query_headers,
		char** seq_headers,
		uint16_t * query_sequences_lengths,
		uint32_t query_sequences_count,
		uint64_t Q,
		uint32_t * query_disp,
		char* seq_db,
		uint64_t D,
		uint16_t* seq_len_db,
		uint64_t seq_count_db,
		uint32_t* seq_disp_db,
		int ** ptr_assigned,
		double* elapsed_time,
		uint64_t* total_done,
		char* shifted_copy,
		int* extension,
		int* shifted_shift,
		uint64_t shifted_linear,
		int threshold
		);

void filter_sequences_shifted_threshold (char * query_sequences,
		char** query_headers,
		char** seq_headers,
		uint16_t * query_sequences_lengths,
		uint32_t query_sequences_count,
		uint64_t Q,
		uint32_t * query_disp,
		char* seq_db,
		uint64_t D,
		uint16_t* seq_len_db,
		uint64_t seq_count_db,
		uint32_t* seq_disp_db,
		int ** ptr_assigned,
		double* elapsed_time,
		uint64_t* total_done,
		char* shifted_copy,
		int* extension,
		int* shifted_shift,
		uint64_t shifted_linear,
		int threshold
		);
__declspec(target(mic)) int popcount_4(uint64_t x);

void free_aligned(void *raw_data);
void *malloc_aligned(size_t alignment, size_t bytes);
//__declspec(target(mic)) uint16_t concatenate8(uint8_t x, uint8_t y) ;
//__declspec(target(mic)) uint32_t concatenate16(uint16_t x, uint16_t y) ;
#endif /* HEURISTIC_H_ */

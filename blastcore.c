#include "blastcore.h"



//Using KNC intrinsics
void smith_waterman_vectorized (Query* query,
		Chunked_Database* db,
		char * submat,
		int open_gap,
		int extend_gap,
		int num_threads,
		int * scores,
		double * workTime,
		uint16_t query_length_threshold){

	long int iter=0, i, j=0, k=0, l;
	double tick;

	char *q_seqs, * db_seqs, *queryProfiles;
	uint16_t * m, *n, sequences_db_max_length, query_sequences_max_length;
	uint32_t * q_disp, * db_disp = NULL, offload_max_vect_sequences_db_count=0, qp_count, sp_count;
	uint32_t c, offload_vect_sequences_db_count;
	uint64_t offload_vD;
	uint64_t Q = query->lin_len_total;

	q_seqs = query->q_seq;
	db_seqs = db->chunk_b[0];
	m = query->query_seq_len;
	n = db->chunk_n[0];
	q_disp = query->query_seq_disp;
	db_disp = db->chunk_b_disp[0];

	query_sequences_max_length = query->query_seq_len[query->query_sequences_count-1];
	sequences_db_max_length = db->chunk_n[db->chunk_count -1][db->chunk_vect_sequences_db_count[db->chunk_count-1]-1];

	// build query profile's
	queryProfiles = (char *)_mm_malloc(Q*BLOSUM_COLS*sizeof(char), 64);
	for (iter=0; iter<Q ; iter++)
		memcpy(queryProfiles+iter*BLOSUM_COLS,submat+q_seqs[iter]*BLOSUM_COLS,BLOSUM_COLS*sizeof(char));

	// calculate number of query sequences that are processed with query and score profile
	iter = 0;
	while ((iter < query->query_sequences_count) && (query->query_seq_len[iter] < query_length_threshold))
		iter++;
	qp_count = iter;
	sp_count = query->query_sequences_count-qp_count;

	tick = getTimestamp();
			uint64_t query_sequences_count = query->query_sequences_count;
			offload_vD = db->chunk_vD[0];
			offload_vect_sequences_db_count = db->chunk_vect_sequences_db_count[0];
			db_seqs = db->chunk_b[0];
			n = db->chunk_n[0];
			db_disp = db->chunk_b_disp[0];
			//printf("Lengths: %u %u %u\n",Q,offload_vD,query_sequences_count);
			#pragma offload target(mic:0) in(submat: length(BLOSUM_ELEMS) )\
			in(queryProfiles: length(Q*BLOSUM_COLS) )\
			in(db_seqs: length(offload_vD) ) \
			in(n: length(offload_vect_sequences_db_count) )\
			in(db_disp: length(offload_vect_sequences_db_count) ) \
			in(q_seqs: length(Q)) in(m: length(query_sequences_count) )\
			in(q_disp: length(query_sequences_count) ) \
			out(scores: length(query_sequences_count*offload_vect_sequences_db_count*VECTOR_LENGTH))
				#pragma omp parallel shared(offload_vect_sequences_db_count, query_sequences_count, open_gap, extend_gap, query_sequences_max_length, sequences_db_max_length) num_threads(num_threads)
				{

					__m512i  *vH, *vE, *vF, *lastCol;
					int  * ptr_scores;
					char * ptr_a, * ptr_b, *ptr_b_block, * scoreProfile, *queryProfile, *ptr_scoreProfile;

					__m512i vzero  __attribute__ ((aligned(64)))= _mm512_setzero_epi32(), score, vPrev, vCur, vAux, aux2, aux3, aux4, auxLastCol;
					__m512i vextend_gap  __attribute__ ((aligned(64)))= _mm512_set1_epi32(extend_gap), vopen_extend_gap = _mm512_set1_epi32(open_gap+extend_gap);
					__m512i v16  __attribute__ ((aligned(64)))= _mm512_set1_epi32(16), submat_hi, submat_lo, b_values;
					__mmask16 mask;

					uint32_t tid, iter, j, i, l, k, disp_1, disp_2, disp_3, dim1, dim2, nbb;
					uint64_t t, s, q;

					// allocate memory
					vH = (__m512i *) _mm_malloc((BLOCK_SIZE+1)*sizeof(__m512i), 64);
					vE = (__m512i *) _mm_malloc((BLOCK_SIZE+1)*sizeof(__m512i), 64);
					vF = (__m512i *) _mm_malloc((query_sequences_max_length)*sizeof(__m512i), 64);
					lastCol = (__m512i *) _mm_malloc((query_sequences_max_length)*sizeof(__m512i), 64);
					if (query_sequences_max_length >= query_length_threshold)
						scoreProfile = (char *) _mm_malloc(INTRINSIC_LEN*sequences_db_max_length*sizeof(char),16);

					// calculate chunk alignments using query profile technique
					#pragma omp for schedule(dynamic) nowait
					for (t=0; t< qp_count*offload_vect_sequences_db_count; t++) {

						q = (qp_count-1) - (t % qp_count);
						s = (offload_vect_sequences_db_count-1) - (t / qp_count);

						queryProfile = queryProfiles + q_disp[q]*BLOSUM_COLS;
						ptr_b = db_seqs + db_disp[s];
						ptr_scores = scores + (q*offload_vect_sequences_db_count+s)*VECTOR_LENGTH;

						// init buffers

						for (iter=0; iter<m[q] ; iter++ ) vF[iter] = _mm512_setzero_epi32(); // index 0 is not used

						for (iter=0; iter<m[q] ; iter++ ) lastCol[iter] = _mm512_setzero_epi32();
						
						// set score to 0
						score = _mm512_setzero_epi32();

						// calculate number of blocks
						nbb = ceil( (double) n[s] / (double) BLOCK_SIZE);

						for (k=0; k < nbb; k++){

							// calculate dim1
							disp_1 = k*BLOCK_SIZE;
							dim1 = (BLOCK_SIZE < n[s]-disp_1 ? BLOCK_SIZE : n[s]-disp_1);

							// get db_seqs block
							ptr_b_block = ptr_b + disp_1*VECTOR_LENGTH;

							// init buffers

							for (iter=1; iter<dim1+1 ; iter++ ) vE[iter] = _mm512_setzero_epi32(); //index 0 is not used

							for (iter=0; iter<dim1 ; iter++ ) vH[iter] = _mm512_setzero_epi32();
							auxLastCol = _mm512_setzero_epi32();

							for( iter = 0; iter < m[q]; iter++){
						
								// vPrev must start in 0
								vPrev = _mm512_setzero_epi32();
								// update vH[0] with lastCol elements
								vH[0] = lastCol[iter];
								// load submat values corresponding to vCur q_seqs residue
								disp_1 = iter*BLOSUM_COLS;

								submat_lo = _mm512_extload_epi32(queryProfile+disp_1, _MM_UPCONV_EPI32_SINT8, _MM_BROADCAST32_NONE, 0);
								submat_hi = _mm512_extload_epi32(queryProfile+disp_1+16, _MM_UPCONV_EPI32_SINT8, _MM_BROADCAST32_NONE, 0);

								// calculate dim2
								dim2 = dim1 / UNROLL;

								for (i=0; i<dim2 ; i++) {

									#pragma unroll(UNROLL)
									for( j=i*UNROLL+1, l=0; l < UNROLL;  l++, j++) {
										//calcuate the diagonal value

										b_values = _mm512_extload_epi32(ptr_b_block+(j-1)*VECTOR_LENGTH, _MM_UPCONV_EPI32_SINT8, _MM_BROADCAST32_NONE, 0);

										mask = _mm512_cmpge_epi32_mask(b_values,v16);
										vAux = _mm512_permutevar_epi32(b_values, submat_lo);
										vAux = _mm512_mask_permutevar_epi32(vAux, mask, b_values, submat_hi);
										vCur = _mm512_add_epi32(vH[j-1], vAux);
										// calculate vCur max value
										vCur = _mm512_max_epi32(vCur, vF[iter]);
										vCur = _mm512_max_epi32(vCur, vE[j]);
										vCur = _mm512_max_epi32(vCur, vzero);
										// update vF and vE
										vF[iter] = _mm512_sub_epi32(vF[iter], vextend_gap);
										vE[j] = _mm512_sub_epi32(vE[j], vextend_gap);
										vAux = _mm512_sub_epi32(vCur, vopen_extend_gap);
										vF[iter] = _mm512_max_epi32(vF[iter], vAux);
										vE[j] =  _mm512_max_epi32(vE[j], vAux);
										// update vH buffer
										vH[j-1] = vPrev;
										vPrev = vCur;
										// update max score
										score = _mm512_max_epi32(score,vCur);
									}
								}
								#pragma unroll
								for( j = dim2*UNROLL+1; j < dim1+1; j++) {
									//calcuate the diagonal value

									b_values = _mm512_extload_epi32(ptr_b_block+(j-1)*VECTOR_LENGTH, _MM_UPCONV_EPI32_SINT8, _MM_BROADCAST32_NONE, 0);

									mask = _mm512_cmpge_epi32_mask(b_values,v16);
									vAux = _mm512_permutevar_epi32(b_values, submat_lo);
									vAux = _mm512_mask_permutevar_epi32(vAux, mask, b_values, submat_hi);
									vCur = _mm512_add_epi32(vH[j-1], vAux);
									// calculate vCur max value
									vCur = _mm512_max_epi32(vCur, vF[iter]);
									vCur = _mm512_max_epi32(vCur, vE[j]);
									vCur = _mm512_max_epi32(vCur, vzero);
									// update vF and vE
									vF[iter] = _mm512_sub_epi32(vF[iter], vextend_gap);
									vE[j] = _mm512_sub_epi32(vE[j], vextend_gap);
									vAux = _mm512_sub_epi32(vCur, vopen_extend_gap);
									vF[iter] = _mm512_max_epi32(vF[iter], vAux);
									vE[j] =  _mm512_max_epi32(vE[j], vAux);
									// update vH buffer
									vH[j-1] = vPrev;
									vPrev = vCur;
									// update max score
									score = _mm512_max_epi32(score,vCur);
								}
								// update lastCol
								lastCol[iter] = auxLastCol;
								auxLastCol = vCur;

							}

						}

						// store max value
						
						_mm512_store_epi32(ptr_scores, score);
						

					}

					// calculate chunk alignments using score profile technique
					#pragma omp for schedule(dynamic) nowait
					for (t=0; t< sp_count*offload_vect_sequences_db_count; t++) {

						q = qp_count + (sp_count-1) - (t % sp_count);
						s = (offload_vect_sequences_db_count-1) - (t / sp_count);

						ptr_a = q_seqs + q_disp[q];
						ptr_b = db_seqs + db_disp[s];
						ptr_scores = scores + (q*offload_vect_sequences_db_count+s)*VECTOR_LENGTH;

						// build score profile
						disp_1 = n[s]*VECTOR_LENGTH;
						for (iter=0; iter<n[s] ;iter++ ) {

							vAux = _mm512_extload_epi32(ptr_b+iter*VECTOR_LENGTH, _MM_UPCONV_EPI32_SINT8, _MM_BROADCAST32_NONE, 0);

							disp_2 = iter*VECTOR_LENGTH;
							#pragma unroll(UNROLL)
							for (j=0; j< BLOSUM_ROWS; j++) {

								aux2 = _mm512_i32extgather_epi32(vAux, submat + j*BLOSUM_COLS, _MM_UPCONV_EPI32_SINT8  , 1, 0);
								_mm512_extstore_epi32(scoreProfile+disp_2+j*disp_1, aux2, _MM_DOWNCONV_EPI32_SINT8 , _MM_HINT_NONE );


							}
						}

						// init buffers

						for (iter=0; iter<m[q] ; iter++ ) vF[iter] = _mm512_setzero_epi32(); // index 0 is not used

						for (iter=0; iter<m[q] ; iter++ ) lastCol[iter] = _mm512_setzero_epi32();
						
						// set score to 0
						score = _mm512_setzero_epi32();

						// calculate number of blocks
						nbb = ceil( (double) n[s] / (double) BLOCK_SIZE);

						for (k=0; k < nbb; k++){

							// calculate dim1
							disp_2 = k*BLOCK_SIZE;
							dim1 = (BLOCK_SIZE < n[s]-disp_2 ? BLOCK_SIZE : n[s]-disp_2);

							// init buffers

							for (iter=1; iter<dim1+1 ; iter++ ) vE[iter] = _mm512_setzero_epi32(); //index 0 is not used

							for (iter=0; iter<dim1 ; iter++ ) vH[iter] = _mm512_setzero_epi32();
							auxLastCol = _mm512_setzero_epi32();

							for( iter = 0; iter < m[q]; iter++){
						
								// vPrev initialized
								vPrev = _mm512_setzero_epi32();

								// vh[0] contains lastcol elements
								vH[0] = lastCol[iter];

								// which score are we calculating? Add
								ptr_scoreProfile = scoreProfile + ((int)(ptr_a[iter]))*disp_1 + disp_2*VECTOR_LENGTH;

								//This is the main loop--> see tfg page 46 for the main loop description
								#pragma unroll(UNROLL)
								for( l=1; l < dim1+1; l++) {
									//vCur is
									vCur = _mm512_add_epi32(vH[l-1], _mm512_extload_epi32(ptr_scoreProfile+(l-1)*VECTOR_LENGTH, _MM_UPCONV_EPI32_SINT8, _MM_BROADCAST32_NONE, 0));
									// calculate vCur max value
									vCur = _mm512_max_epi32(vCur, vF[iter]);
									vCur = _mm512_max_epi32(vCur, vE[l]);
									vCur = _mm512_max_epi32(vCur, vzero);
									// update vF and vE
									vF[iter] = _mm512_sub_epi32(vF[iter], vextend_gap);
									vE[l] = _mm512_sub_epi32(vE[l], vextend_gap);
									vAux = _mm512_sub_epi32(vCur, vopen_extend_gap);
									vF[iter] = _mm512_max_epi32(vF[iter], vAux);
									vE[l] =  _mm512_max_epi32(vE[l], vAux);
									// update vH buffer
									vH[l-1] = vPrev;
									vPrev = vCur;
									// update max score
									score = _mm512_max_epi32(score,vCur);
								}
								// update lastCol
								lastCol[iter] = auxLastCol;
								auxLastCol = vCur;

							}

						}

						// store max value
						
						_mm512_store_epi32(ptr_scores, score);
						

					}

					_mm_free(vH);
					_mm_free(vE);
					_mm_free(vF);
					_mm_free(lastCol);
 					if (query_sequences_max_length > query_length_threshold) _mm_free(scoreProfile);
				}

	*workTime = getTimestamp()-tick;
	_mm_free(queryProfiles);
	
}


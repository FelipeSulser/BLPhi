/*
 * heuristic.c
 *
 *  Created on: 25/04/2016
 *      Author: Felipe Sulser Larraz
 */


#include "heuristic.h"

/**
 * __mmask64 _mm512_cmpeq_epi32_mask (__m512i a, __m512i b)
 */

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
		){
		/**
		 * Assigned will contain the following information:
		 *
		 *	assigned is a matrix
		 *
		 *	x x x x x x
		 *	x x x x x x
		 *	x x x x x x
		 *	x x x x x x
		 *	x x x x x x
		 *	x x x x x x
		 *
		 *	cols = query_sequences_count
		 *	row = db_count
		 *
		 * */

		uint32_t filtered_Q;
		//uint8_t first,second,third,fourth;
		//uint16_t firsthalf,secondhalf;
		double tick = 0;
		//1 since only one chunk

		//lets say the db are the columns
		//int assigned[query_sequences_count*seq_count_db+seq_count_db];
		int* assigned = (int*)malloc(sizeof(int)*(query_sequences_count*seq_count_db+seq_count_db));
		memset(assigned,0,sizeof(int)*query_sequences_count*seq_count_db+seq_count_db);
		//linear structure for the win

		int mic_no;
		int vector_length = VECTOR_LENGTH;
		//only one chunk
		//Due to offloads structure, cannot offload double references
		char* a = seq_db; //only 1
		uint32_t* a_disp =seq_disp_db;
		uint16_t* a_len = seq_len_db;
		uint32_t offload_db_count = seq_count_db;
		mic_no = omp_get_thread_num();
		uint32_t mem_ctr = 0;

		tick = getTimestamp();
		#pragma offload target(mic:mic_no)\
		inout(assigned:length(query_sequences_count*seq_count_db) )\
		in(a:length(D) )\
		in(a_disp:length(offload_db_count))\
		in(a_len:length(offload_db_count))\
		in(query_disp:length(query_sequences_count))\
		in(query_sequences_lengths:length(query_sequences_count))\
		in(query_sequences:length(Q))\
		in(query_sequences_count)\
		//out(offload_assigned:length(offload_db_count)ALLOC)
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			char q1, q2, q3, q4, db1, db2, db3, db4;
			uint32_t i;
			for(i = 0; i < query_sequences_count; i++){

				uint32_t len;
				len = query_sequences_lengths[i];
				uint32_t iter_seq;
				#pragma omp for ordered schedule(dynamic) nowait
				for(iter_seq = 0; iter_seq < len-3;iter_seq++){
					//int tid =  omp_get_thread_num();
					//printf("Thread:%u working on iter_seq:%u i:%u \n",tid,iter_seq,i);
					int diff;
					int diff2;
					diff = query_sequences_lengths[i] - iter_seq;
					if(diff >= 4){
						q1 = *(query_sequences+query_disp[i]+iter_seq);
						q2 = *(query_sequences+query_disp[i]+iter_seq+1);
						q3 = *(query_sequences+query_disp[i]+iter_seq+2);
						q4 = *(query_sequences+query_disp[i]+iter_seq+3);
						//now for all sequences in db
						uint32_t j;
						for(j = 0; j < seq_count_db; j++){
							//printf("This one is:%u\n",a_len[j]);
							uint32_t iter_db;
							for(iter_db = 0; iter_db+3 < a_len[j];iter_db++){
								//printf("%u ",a[iter_db+a_disp[j]]);
								diff2 = a_len[j] - iter_db;
								if(diff2 >= 4){
									db1 = *(a+iter_db+a_disp[j]);
									db2 = *(a+iter_db+a_disp[j]+1);
									db3 = *(a+iter_db+a_disp[j]+2);
									db4 = *(a+iter_db+a_disp[j]+3);
									if(db1 == q1 && db2 == q2 && db3 == q3 && db4 == q4){
										//printf("k:%u r:%u\n",i,j);
										assigned[seq_count_db*i+j] = 1;
										break; //yes, a break
									}
								}
							}
						}
					}
				}

				//eligible for hit

			}

		}
		*elapsed_time =  getTimestamp() - tick;
		int r= 0;
		int k = 0;
		uint64_t qs = 0;
				for(k = 0; k < query_sequences_count; k++){
					for(r = 0; r < seq_count_db;r++){
				//printf("%u\t",assigned[query_sequences_count*k+r]);
				if(assigned[seq_count_db*k+r]>0){
					//printf("Db:%u q:%u assigned:%u\n",r,k,assigned[query_sequences_count*k+r]);
					qs++;
				}
			}
			//printf("\n");
		}
		printf("Non-Vectorized filter:\t\tTotal comp:%lu Done: %lu (%f %)\n",(uint64_t)query_sequences_count*seq_count_db,
						qs,((1.0*qs)/(query_sequences_count*seq_count_db))*100);
		uint64_t ctr_q = 0;
		qs = 0;

		*ptr_assigned = assigned;

		*total_done = 0;
}


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
		){
		/**
		 * Assigned will contain the following information:
		 *
		 *	assigned is a matrix
		 *
		 *	x x x x x x
		 *	x x x x x x
		 *	x x x x x x
		 *	x x x x x x
		 *	x x x x x x
		 *	x x x x x x
		 *
		 *	cols = query_sequences_count
		 *	row = db_count
		 *
		 * */

		uint32_t filtered_Q;
		//uint8_t first,second,third,fourth;
		//uint16_t firsthalf,secondhalf;
		double tick = 0;
		//1 since only one chunk

		//lets say query is the columns
		int assigned[red_query_sequences_count*red_seq_count_db+red_seq_count_db];
		memset(assigned,0,sizeof(int)*red_query_sequences_count*red_seq_count_db+red_seq_count_db);
		//linear structure for the win



		int mic_no;
		int vector_length = VECTOR_LENGTH;
		//only one chunk
		//Due to offloads structure, cannot offload double references
		char* a = red_seq_db; //only 1
		uint32_t* a_disp =red_seq_disp_db;
		uint16_t* a_len = red_seq_len_db;
		uint32_t offload_db_count = red_seq_count_db;
		mic_no = omp_get_thread_num();
		uint32_t mem_ctr = 0;

		tick = getTimestamp();
		#pragma offload target(mic:mic_no)\
		inout(assigned:length(red_query_sequences_count*red_seq_count_db) ALLOC)\
		in(a:length(red_D) ALLOC)\
		in(a_disp:length(offload_db_count)ALLOC)\
		in(a_len:length(offload_db_count)ALLOC)\
		in(red_query_disp:length(red_query_sequences_count)ALLOC)\
		in(red_query_sequences_lengths:length(red_query_sequences_count)ALLOC)\
		in(red_query_sequences:length(red_Q)ALLOC)\
		in(red_query_sequences_count)
		#pragma omp parallel num_threads(NUM_THREADS)
		{

			uint32_t block,seq_num=0,query_num=0;
			uint32_t i, j, iter_seq, chunk_it, iter_db;
			char q1, q2, q3, q4, db1, db2, db3, db4;
			uint32_t val_db,val_q;
			uint32_t len;
			int hasHits = 0;
			for(i = 0; i < red_query_sequences_count; i++){
				hasHits = 0;
				len = red_query_sequences_lengths[i];
				#pragma omp for ordered schedule(dynamic) nowait
				for(iter_seq = 0; iter_seq < len-3;iter_seq++){
					//int tid =  omp_get_thread_num();
					//printf("Thread:%u working on iter_seq:%u i:%u \n",tid,iter_seq,i);
					int diff;
					int diff2;
					diff = red_query_sequences_lengths[i] - iter_seq;
					if(diff >= 4){
						q1 = *(red_query_sequences+red_query_disp[i]+iter_seq);
						q2 = *(red_query_sequences+red_query_disp[i]+iter_seq+1);
						q3 = *(red_query_sequences+red_query_disp[i]+iter_seq+2);
						q4 = *(red_query_sequences+red_query_disp[i]+iter_seq+3);
						//now for all sequences in db
						int gotHit = 0;
						for(j = 0; j < red_seq_count_db; j++){
							//printf("This one is:%u\n",a_len[j]);
							for(iter_db = 0; iter_db+3 < a_len[j];iter_db++){
								//printf("%u ",a[iter_db+a_disp[j]]);
								diff2 = a_len[j] - iter_db;
								if(diff2 >= 4){
									db1 = *(a+iter_db+a_disp[j]);
									db2 = *(a+iter_db+a_disp[j]+1);
									db3 = *(a+iter_db+a_disp[j]+2);
									db4 = *(a+iter_db+a_disp[j]+3);
									if(db1 == q1 && db2 == q2 && db3 == q3 && db4 == q4){
										assigned[red_seq_count_db*i+j]=1;
										break; //yes, a break
									}
								}
							}
						}
					}
				}

				//eligible for hit

			}

		}

		int r= 0;
		int k = 0;
		uint64_t qs = 0;

				for(k = 0; k < red_query_sequences_count; k++){
					for(r = 0; r < red_seq_count_db;r++){
				//printf("%u\t",assigned[query_sequences_count*k+r]);
				if(assigned[red_seq_count_db*k+r]>0){
					//printf("Db:%u q:%u assigned:%u\n",r,k,assigned[query_sequences_count*k+r]);
					qs++;
					//printf("%s\n",seq_headers[r]);
				}
			}
			//printf("\n");
		}
		printf("Non-Vectorized filter:\t\tTotal comp:%lu Done: %lu (%f %)\n",(uint64_t)red_query_sequences_count*red_seq_count_db,
						qs,((1.0*qs)/(red_query_sequences_count*red_seq_count_db))*100);



		*ptr_assigned = assigned;
		*elapsed_time =  getTimestamp() - tick;
		*total_done = qs;
}
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
		){
		uint32_t filtered_Q;
		//uint8_t first,second,third,fourth;
		//uint16_t firsthalf,secondhalf;

		//1 since only one chunk

		//lets say query is the columns
		uint32_t assigned[query_sequences_count*seq_count_db];
		memset(assigned,0,sizeof(uint32_t)*query_sequences_count*seq_count_db);
		//linear structure for the win
		printf("Eflag on\n");

		int mic_no;
		int vector_length = VECTOR_LENGTH;
		//only one chunk
		//Due to offloads structure, cannot offload double references
		char* a = seq_db; //only 1
		uint32_t* a_disp = seq_disp_db;
		uint16_t* a_len = seq_len_db;
		uint32_t offload_db_count = seq_count_db;
		mic_no = omp_get_thread_num();
		uint32_t mem_ctr = 0;

		#pragma omp parallel num_threads(4)
		{

			uint32_t block,seq_num=0,query_num=0;
			uint32_t i, j, iter_seq, chunk_it, iter_db;
			char q1, q2, q3, q4,q5,q6,q7,q8, db1, db2, db3, db4,db5,db6,db7,db8;
			uint32_t val_db,val_q;
			uint32_t len;
			int hasHits = 0;
			for(i = 0; i < query_sequences_count; i++){
				hasHits = 0;
				len = query_sequences_lengths[i];
				#pragma omp for ordered schedule(dynamic) nowait
				for(iter_seq = 0; iter_seq < len-7;iter_seq++){
					//int tid =  omp_get_thread_num();
					//printf("Thread:%u working on iter_seq:%u i:%u \n",tid,iter_seq,i);
					int diff;
					int diff2;
					diff = query_sequences_lengths[i] - iter_seq;

						q1 = *(query_sequences+query_disp[i]+iter_seq);
						q2 = *(query_sequences+query_disp[i]+iter_seq+1);
						q3 = *(query_sequences+query_disp[i]+iter_seq+2);
						q4 = *(query_sequences+query_disp[i]+iter_seq+3);
						q5 = *(query_sequences+query_disp[i]+iter_seq+4);
						q6 = *(query_sequences+query_disp[i]+iter_seq+5);
						q7 = *(query_sequences+query_disp[i]+iter_seq+6);
						q8 = *(query_sequences+query_disp[i]+iter_seq+7);
						//now for all sequences in db
						int gotHit = 0;
						for(j = 0; j < seq_count_db; j++){
							//printf("This one is:%u\n",a_len[j]);
							for(iter_db = 0; iter_db+7 < a_len[j];iter_db++){
								//printf("%u ",a[iter_db+a_disp[j]]);
								diff2 = a_len[j] - iter_db;

									db1 = *(a+iter_db+a_disp[j]);
									db2 = *(a+iter_db+a_disp[j]+1);
									db3 = *(a+iter_db+a_disp[j]+2);
									db4 = *(a+iter_db+a_disp[j]+3);
									db5 = *(a+iter_db+a_disp[j]+4);
									db6 = *(a+iter_db+a_disp[j]+5);
									db7 = *(a+iter_db+a_disp[j]+6);
									db8 = *(a+iter_db+a_disp[j]+7);
									if(db1 == q1 && db2 == q2 && db3 == q3 && db4 == q4 && db5 == q5 && db6 == q6 && db7 == q7 && db8 == q8){
										//assigned[query_sequences_count*i+j] = 1;
										//printf("already %u\n",assigned[query_sequences_count*i+j]);
										//printf("HERE\n");
										gotHit=1;

								}
							}
							if(gotHit > 0)
								assigned[seq_count_db*i+j]=1;
							gotHit = 0;
						}
						//printf("Max_j:%u total:%u\n",max_j,seq_count_db);

				}

				//eligible for hit

			}

		}

		int r= 0;
		int k = 0;
		uint32_t qs = 0;

				for(k = 0; k < query_sequences_count; k++){
					for(r = 0; r < seq_count_db;r++){
				//printf("%u\t",assigned[query_sequences_count*k+r]);
				if(assigned[query_sequences_count*k+r]>0){
					//printf("Db:%u q:%u assigned :%u value:\n",r,k,assigned[query_sequences_count*k+r]);
					qs++;
				}
			}
			//printf("\n");
		}
		printf("Non-Vectorized filter:\t\tTotal comp:%lu Done: %lu (%f %)\n",(uint64_t)query_sequences_count*seq_count_db,
						qs,((1.0*qs)/(query_sequences_count*seq_count_db))*100);



		*ptr_assigned = assigned;


}

//This is wrong, misses some subsequences!!
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
		){
		double tick;
		uint32_t filtered_Q;
		int assigned[query_sequences_count*seq_count_db+seq_count_db];
		//initialize to zero or else we wont know if it was a hit or not
		memset(assigned,0,sizeof(int)*query_sequences_count*seq_count_db+seq_count_db);
		int mic_no;
		int vector_length = VECTOR_LENGTH;
		//Due to offloads structure, cannot offload double references
		char* a __attribute__ ((aligned(64)))=seq_db; //only 1
		uint32_t* a_disp = seq_disp_db;
		uint16_t* a_len = seq_len_db;
		uint32_t offload_db_count = seq_count_db;
		mic_no = omp_get_thread_num();
		uint32_t mem_ctr = 0;
		uint32_t d_i;
		tick = getTimestamp();
		#pragma offload target(mic:mic_no)\
		in(a[0:D]:align(64)  )\
		in(a_disp:length(offload_db_count)  )\
		in(a_len:length(offload_db_count)  )\
		in(seq_count_db)\
		in(query_disp:length(query_sequences_count) )\
		in(query_sequences_lengths:length(query_sequences_count) )\
		in(query_sequences[0:Q]:align(64) )\
		in(query_sequences_count)\
		inout(assigned:length(query_sequences_count*offload_db_count) )
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			char* ptr_a_it __attribute__ ((aligned(64)));
			char* ptr_q_it __attribute__ ((aligned(64)));
			uint32_t block,nhits=0,seq_num=0,query_num=0;
			uint32_t i=0, j=0, iter_seq=0, chunk_it=0, gotHit=0, iter_db=0;
			int diff;
			uint8_t q1, q2, q3, q4, db1, db2, db3, db4;
			uint32_t window_it = 0;
			uint32_t window_db_iter = 0;
			int it_db = 0;
			__mmask16 mask;
			__m512i dataq __attribute__ ((aligned(64)))= _mm512_setzero_epi32();
			__m512i datadb __attribute__ ((aligned(64)))= _mm512_setzero_epi32();

			int disc = 0;
			int res;
			char q[64] __attribute__ ((aligned(64)));
			char db[64] __attribute__ ((aligned(64)));

			for(i = 0; i < query_sequences_count; i++){
				gotHit = 0;
				//simd takes 64 and vectorizes that
				//if(query_sequences_lengths[i] > 63){
					#pragma omp for schedule(dynamic) nowait
					for(iter_seq = 0; iter_seq < query_sequences_lengths[i] - 63;iter_seq++){
							//now construct the m512i from the 16 32 bits
							memcpy(q,&query_sequences[query_disp[i]+iter_seq],sizeof(char)*64);
							dataq = _mm512_load_epi32(q);
							//now for all sequences in db
							for(j = 0; j < seq_count_db; j++){
								if(a_len[j] > 63){
									//printf("alen:%u\n",a_len[j]);
									for(iter_db = 0; iter_db + 63 < a_len[j];iter_db++){
										memcpy(db,&a[iter_db+a_disp[j]],sizeof(char)*64);
										datadb = _mm512_load_epi32(db);
										mask = _mm512_cmpeq_epi32_mask(dataq,datadb);
										res = _mm512_mask2int(mask);
										if(res !=0 ){
											assigned[seq_count_db*i+j] = 1;
											break;
										}
									 }
									//case we have less than 64 char remaining in the sequence
										if(iter_db < a_len[j]){
											int info_len = a_len[j] - iter_db;
											//the rest will be filled with garbage
											for(it_db = 0; it_db < 64; it_db++){
												if(it_db < info_len)
													db[it_db]= *(a+iter_db+a_disp[j]+it_db);
												else
													db[it_db]= 0;
											}
											datadb = _mm512_load_epi32(db);
											mask = _mm512_cmpeq_epi32_mask(dataq,datadb);
											res = _mm512_mask2int(mask);
											int count_shifted = 16 - info_len/4;
											if((res << (16 - count_shifted)) != 0){
												assigned[seq_count_db*i+j] = 1;
											}
										}
									//case the sequence itself has less than 64
									}else{
										//Case db sequence non multiple of 64
										int info_len = a_len[j];
										int it_db;
										//the rest will be filled with garbage
										for(it_db = 0; it_db < info_len; it_db++){
											db[it_db]= *(a+iter_db+a_disp[j]+it_db);
											//printf("%u ",db[it_db]);
										}
										datadb = _mm512_load_epi32(db);
										mask = _mm512_cmpeq_epi32_mask(dataq,datadb);
										res = _mm512_mask2int(mask);
										int count_shifted = 16 - info_len/4;
										if(res << (16 - count_shifted) != 0){
											assigned[seq_count_db*i+j] = 1;
										}

								}

							}

						}
					//la query tiene menos de 64 caracteres, no podemos vectorizar
				//}
				//caso que no sea multiplo de 64
				if(iter_seq < query_sequences_lengths[i]){
					int len = query_sequences_lengths[i];
					#pragma omp for  schedule(dynamic) nowait
					for(iter_seq = 0; iter_seq < len-3;iter_seq++){
						int diff;
						int diff2;
						diff = query_sequences_lengths[i] - iter_seq;
						if(diff >= 4){
							q1 = *(query_sequences+query_disp[i]+iter_seq);
							q2 = *(query_sequences+query_disp[i]+iter_seq+1);
							q3 = *(query_sequences+query_disp[i]+iter_seq+2);
							q4 = *(query_sequences+query_disp[i]+iter_seq+3);
							//now for all sequences in db
							int gotHit = 0;
							for(j = 0; j < seq_count_db; j++){
								//printf("This one is:%u\n",a_len[j]);
								for(iter_db = 0; iter_db+3 < a_len[j];iter_db++){
									//printf("%u ",a[iter_db+a_disp[j]]);
									diff2 = a_len[j] - iter_db;
									if(diff2 >= 4){
										db1 = *(a+iter_db+a_disp[j]);
										db2 = *(a+iter_db+a_disp[j]+1);
										db3 = *(a+iter_db+a_disp[j]+2);
										db4 = *(a+iter_db+a_disp[j]+3);
										if(db1 == q1 && db2 == q2 && db3 == q3 && db4 == q4){
											//assigned[query_sequences_count*i+j] = 1;
											//printf("already %u\n",assigned[query_sequences_count*i+j]);
											gotHit=1;
											break;
										}
									}
								}
								if(gotHit > 0)
									assigned[seq_count_db*i+j]=1;
								gotHit = 0;
							}
						}
					}
				}
			}
		}
		*elapsed_time = getTimestamp()-tick;
		int r= 0;
		int k = 0;
		uint64_t qs = 0;

				for(k = 0; k < query_sequences_count; k++){
					for(r = 0; r < seq_count_db;r++){
				//printf("%u\t",assigned[query_sequences_count*k+r]);
				if(assigned[seq_count_db*k+r]>0){
					qs++;
				}
			}
			//printf("\n");
		}
		printf("Non-Vectorized filter:\t\tTotal comp:%lu Done: %lu (%f %)\n",(uint64_t)query_sequences_count*seq_count_db,qs,((1.0*qs)/(query_sequences_count*seq_count_db))*100);
		*ptr_assigned = assigned;

		*total_done = qs;

}

__declspec(target(mic)) int popcount_4(uint64_t x) {
    int count;
    for (count=0; x; count++)
        x &= x-1;
    return count;
}



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
		){
		double tick;
		uint32_t filtered_Q;
		int* assigned = (int*)malloc(sizeof(int)*(query_sequences_count*seq_count_db+seq_count_db));
		//int assigned[query_sequences_count*seq_count_db+seq_count_db];
		//initialize to zero or else we wont know if it was a hit or not
		memset(assigned,0,sizeof(int)*query_sequences_count*seq_count_db+seq_count_db);
		int mic_no;
		int vector_length = VECTOR_LENGTH;
		//Due to offloads structure, cannot offload double references
		char* a __attribute__ ((aligned(64)))=seq_db; //only 1
		uint32_t* a_disp = seq_disp_db;
		uint16_t* a_len = seq_len_db;
		uint32_t offload_db_count = seq_count_db;
		mic_no = omp_get_thread_num();
		uint32_t mem_ctr = 0;
		uint32_t d_i;
		tick = getTimestamp();
		#pragma offload target(mic:mic_no)\
		in(a[0:D]:align(64) ALLOC )\
		in(a_disp:length(offload_db_count) ALLOC )\
		in(a_len:length(offload_db_count) ALLOC )\
		in(seq_count_db)\
		in(query_disp:length(query_sequences_count) ALLOC)\
		in(query_sequences_lengths:length(query_sequences_count) ALLOC)\
		in(query_sequences[0:Q]:align(64) ALLOC)\
		in(query_sequences_count)\
		in(extension:length(seq_count_db)ALLOC)\
		in(shifted_shift:length(seq_count_db) ALLOC)\
		inout(assigned:length(query_sequences_count*offload_db_count) ALLOC)\
		in(shifted_copy:length(shifted_linear)ALLOC)
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			char* ptr_a_it __attribute__ ((aligned(64)));
			char* ptr_q_it __attribute__ ((aligned(64)));
			uint32_t block,nhits=0,seq_num=0,query_num=0;
			uint32_t i=0, j=0, iter_seq=0, chunk_it=0, gotHit=0, iter_db=0;
			int diff;
			uint8_t q1, q2, q3, q4, db1, db2, db3, db4;
			uint32_t window_it = 0;
			uint32_t window_db_iter = 0;
			int it_db = 0;
			__mmask16 mask;
			__m512i dataq __attribute__ ((aligned(64)))= _mm512_setzero_epi32();
			__m512i datadb __attribute__ ((aligned(64)))= _mm512_setzero_epi32();

			int disc = 0;
			int res;
			char q[64] __attribute__ ((aligned(64)));
			char db[64] __attribute__ ((aligned(64)));

			for(i = 0; i < query_sequences_count; i++){
				gotHit = 0;
					#pragma omp for schedule(dynamic) nowait
					for(iter_seq = 0; iter_seq  < query_sequences_lengths[i] -3  ;iter_seq++){

							q1 = *(query_sequences+query_disp[i]+iter_seq);
							q2 = *(query_sequences+query_disp[i]+iter_seq+1);
							q3 = *(query_sequences+query_disp[i]+iter_seq+2);
							q4 = *(query_sequences+query_disp[i]+iter_seq+3);
							//printf("q1:%u q2:%u q3:%u q4:%u\n",q1,q2,q3,q4);
							for(j = 0; j < seq_count_db; j++){

								int group_it = 0;
								if(a_len[j] >= 4){
									for(iter_db = 0; iter_db < a_len[j];iter_db++){

										//tantos grupos de 4 como longitud de cadena
										db1 = shifted_copy[shifted_shift[j]+group_it];
										db2 = shifted_copy[shifted_shift[j]+group_it+1];
										db3 = shifted_copy[shifted_shift[j]+group_it+2];
										db4 = shifted_copy[shifted_shift[j]+group_it+3];
										//printf("db1:%u db2:%u db3:%u db4:%u\n",db1,db2,db3,db4);
										group_it+=4;

										if(db1 == q1 && db2 == q2 && db3 == q3 && db4 == q4){
										//printf("we got a hit\n");
										//printf("k:%u r:%u\n",i,j);
											//printf("Hit between:%u and %u\n",i,j);
										assigned[seq_count_db*i+j] = 1;
										break; //yes, a break
									}
							}
						}

					}

				}
			}
		}
		int r= 0;
		int k = 0;
		uint64_t qs = 0;

				for(k = 0; k < query_sequences_count; k++){
					for(r = 0; r < seq_count_db;r++){
				//printf("%u\t",assigned[query_sequences_count*k+r]);
				if(assigned[seq_count_db*k+r]>0){
					qs++;
				}
			}
			//printf("\n");
		}
		printf("Non-Vectorized filter:\t\tTotal comp:%lu Done: %lu (%f %)\n",(uint64_t)query_sequences_count*seq_count_db,qs,((1.0*qs)/(query_sequences_count*seq_count_db))*100);
		*ptr_assigned = assigned;
		printf("%lf\n",tick);
		*elapsed_time = getTimestamp()-tick;
		*total_done = qs;

}
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
		){
		double tick;
		uint32_t filtered_Q;
		int* assigned = (int*)malloc(sizeof(int)*(query_sequences_count*seq_count_db+seq_count_db));
		//int assigned[query_sequences_count*seq_count_db+seq_count_db];
		//initialize to zero or else we wont know if it was a hit or not
		memset(assigned,0,sizeof(int)*query_sequences_count*seq_count_db+seq_count_db);
		int mic_no;
		int vector_length = VECTOR_LENGTH;
		//Due to offloads structure, cannot offload double references
		char* a __attribute__ ((aligned(64)))=seq_db; //only 1
		uint32_t* a_disp = seq_disp_db;
		uint16_t* a_len = seq_len_db;
		uint32_t offload_db_count = seq_count_db;
		mic_no = omp_get_thread_num();
		uint32_t mem_ctr = 0;
		uint32_t d_i;
		tick = getTimestamp();
		#pragma offload target(mic:mic_no)\
		in(a[0:D]:align(64) ALLOC )\
		in(a_disp:length(offload_db_count) ALLOC )\
		in(a_len:length(offload_db_count) ALLOC )\
		in(seq_count_db)\
		in(query_disp:length(query_sequences_count) ALLOC)\
		in(query_sequences_lengths:length(query_sequences_count) ALLOC)\
		in(query_sequences[0:Q]:align(64) ALLOC)\
		in(query_sequences_count)\
		in(extension:length(seq_count_db)ALLOC)\
		in(shifted_shift:length(seq_count_db) ALLOC)\
		inout(assigned:length(query_sequences_count*offload_db_count) ALLOC)\
		in(shifted_copy:length(shifted_linear)ALLOC)
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			char* ptr_a_it __attribute__ ((aligned(64)));
			char* ptr_q_it __attribute__ ((aligned(64)));
			uint32_t block,nhits=0,seq_num=0,query_num=0;
			uint32_t i=0, j=0, iter_seq=0, chunk_it=0, gotHit=0, iter_db=0;
			int diff;
			uint8_t q1, q2, q3, q4, db1, db2, db3, db4;
			uint32_t window_it = 0;
			uint32_t window_db_iter = 0;
			int it_db = 0;
			__mmask16 mask;
			__m512i dataq __attribute__ ((aligned(64)))= _mm512_setzero_epi32();
			__m512i datadb __attribute__ ((aligned(64)))= _mm512_setzero_epi32();

			int disc = 0;
			int res;
			char q[64] __attribute__ ((aligned(64)));
			char db[64] __attribute__ ((aligned(64)));

			for(i = 0; i < query_sequences_count; i++){
				gotHit = 0;
					#pragma omp for schedule(dynamic) nowait
					for(iter_seq = 0; iter_seq  < query_sequences_lengths[i] -3  ;iter_seq++){

							q1 = *(query_sequences+query_disp[i]+iter_seq);
							q2 = *(query_sequences+query_disp[i]+iter_seq+1);
							q3 = *(query_sequences+query_disp[i]+iter_seq+2);
							q4 = *(query_sequences+query_disp[i]+iter_seq+3);
							//printf("q1:%u q2:%u q3:%u q4:%u\n",q1,q2,q3,q4);
							for(j = 0; j < seq_count_db; j++){

								int group_it = 0;
								if(a_len[j] >= 4){
									for(iter_db = 0; iter_db < a_len[j];iter_db++){

										//tantos grupos de 4 como longitud de cadena
										db1 = shifted_copy[shifted_shift[j]+group_it];
										db2 = shifted_copy[shifted_shift[j]+group_it+1];
										db3 = shifted_copy[shifted_shift[j]+group_it+2];
										db4 = shifted_copy[shifted_shift[j]+group_it+3];
										//printf("db1:%u db2:%u db3:%u db4:%u\n",db1,db2,db3,db4);
										group_it+=4;

										if(db1 == q1 && db2 == q2 && db3 == q3 && db4 == q4){
										//printf("we got a hit\n");
										//printf("k:%u r:%u\n",i,j);
											//printf("Hit between:%u and %u\n",i,j);
										assigned[seq_count_db*i+j]++;
										//break; //yes, a break
									}
							}
						}

					}

				}
			}
		}
		int r= 0;
		int k = 0;
		uint64_t qs = 0;

				for(k = 0; k < query_sequences_count; k++){
					for(r = 0; r < seq_count_db;r++){
				//printf("%u\t",assigned[query_sequences_count*k+r]);
				if(assigned[seq_count_db*k+r]>=threshold){
					qs++;
				}
			}
			//printf("\n");
		}
		printf("Non-Vectorized filter:\t\tTotal comp:%lu Done: %lu (%f %)\n",(uint64_t)query_sequences_count*seq_count_db,qs,((1.0*qs)/(query_sequences_count*seq_count_db))*100);
		*ptr_assigned = assigned;
		*elapsed_time = getTimestamp()-tick;
		*total_done = qs;

}

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
		){
		double tick;
		uint32_t filtered_Q;
		//int assigned[query_sequences_count*seq_count_db+seq_count_db];
		int* assigned = (int*)malloc(sizeof(int)*(query_sequences_count*seq_count_db+seq_count_db));
		//initialize to zero or else we wont know if it was a hit or not
		memset(assigned,0,sizeof(int)*query_sequences_count*seq_count_db+seq_count_db);
		int mic_no;
		int vector_length = VECTOR_LENGTH;
		//Due to offloads structure, cannot offload double references
		char* a __attribute__ ((aligned(64)))=seq_db; //only 1
		uint32_t* a_disp = seq_disp_db;
		uint16_t* a_len = seq_len_db;
		uint32_t offload_db_count = seq_count_db;
		mic_no = omp_get_thread_num();
		uint32_t mem_ctr = 0;
		uint32_t d_i;
		tick = getTimestamp();
		#pragma offload target(mic:mic_no)\
		in(a:length(D) )\
		in(a_disp:length(offload_db_count)  )\
		in(a_len:length(offload_db_count)  )\
		in(seq_count_db)\
		in(query_disp:length(query_sequences_count) )\
		in(query_sequences_lengths:length(query_sequences_count) )\
		in(query_sequences[0:Q]:align(64) )\
		in(query_sequences_count)\
		in(extension:length(seq_count_db))\
		in(shifted_shift:length(seq_count_db) )\
		inout(assigned:length(query_sequences_count*seq_count_db) )\
		in(shifted_copy:length(shifted_linear))
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			char* ptr_a_it __attribute__ ((aligned(64)));
			char* ptr_q_it __attribute__ ((aligned(64)));
			uint32_t block,nhits=0,seq_num=0,query_num=0;
			uint32_t i=0, j=0, chunk_it=0, gotHit=0, iter_db=0;
			int diff;
			uint8_t q1, q2, q3, q4, db1, db2, db3, db4;
			uint32_t window_it = 0;
			uint32_t window_db_iter = 0;
			int it_db = 0;
			__mmask16 mask;
			__m512i dataq __attribute__ ((aligned(64)))= _mm512_setzero_epi32();
			__m512i datadb __attribute__ ((aligned(64)))= _mm512_setzero_epi32();
			int disc = 0;
			int res;
			char q[64] __attribute__ ((aligned(64)));
			char db[64] __attribute__ ((aligned(64)));


			for(i = 0; i < query_sequences_count; i++){
				gotHit = 0;
				int real_qsize = 0;
				int more_blocks = (query_sequences_lengths[i] < 4?0:1); //minimum size is 4
				int ltblock = (query_sequences_lengths[i] < 64?1:0);
				int qsize_group =0 ;
				int next_limit;
				int iter_seq = 0;
				int q_it;
				uint32_t val_q;
				int vec_size = 0;
				int n_iter = 0;
				#pragma omp for schedule(dynamic) nowait
				for(iter_seq = 0; iter_seq < query_sequences_lengths[i]-3;iter_seq++){
					real_qsize = 64;
					qsize_group = 16;
					/*int it_b =  0;
					for(it_b = 0; it_b < 16; it_b++){
						memcpy(&q[it_b*4],&query_sequences[query_disp[i]+iter_seq],sizeof(char)*4);
					}
					dataq = _mm512_load_epi32(q);
					*/
					 val_q =query_sequences[query_disp[i]+iter_seq] + (query_sequences[query_disp[i]+iter_seq+1] << 8) + (query_sequences[query_disp[i]+iter_seq+2] << 16) + (query_sequences[query_disp[i]+iter_seq+3] << 24);
					dataq = _mm512_set4_epi32(val_q,val_q,val_q,val_q);

					for(j = 0; j < seq_count_db; j++){
						if(a_len[j] >= 4){
							 vec_size = extension[j];
							 n_iter = vec_size/64;
							for(iter_db = 0; iter_db < n_iter;iter_db++){
								//memcpy(db,&shifted_copy[shifted_shift[j]+(iter_db*64)],64);
								datadb = _mm512_load_epi32(&shifted_copy[shifted_shift[j]+(iter_db*64)]);
								mask = _mm512_cmpeq_epi32_mask(dataq,datadb);
								res = _mm512_mask2int(mask);

							//if(res != 0){
									//printf("Hit between:%u and %u\n",i,j);
									/*pri1ntf("qsize_group:%u\n",qsize_group);
									printf("res:%u\n",res);
									printf("Comparing:\n");
									int ter = 0;
									printf("QUERY\n");
									for(ter = 0; ter < 64 ; ter++){
										printf("%u\t",q[ter]);
									}
									printf("\n\nDB\n");
									for(ter = 0; ter < 64; ter++){
										printf("%u\t",db[ter]);
									}
									printf("\n");*/
								//}
								if(res != 0){
									assigned[seq_count_db*i+j] = 1;
									break;
								}
							}
						}

					}
					//chheck more_block conditionç

				}
			}
		}
		*elapsed_time = getTimestamp()-tick;
		int r= 0;
		int k = 0;
		uint64_t qs = 0;

				for(k = 0; k < query_sequences_count; k++){
					for(r = 0; r < seq_count_db;r++){
				if(assigned[seq_count_db*k+r]>0){
					qs++;
				}
			}
			//printf("\n");
		}
		printf("Non-Vectorized filter:\t\tTotal comp:%lu Done: %lu (%f %)\n",(uint64_t)query_sequences_count*seq_count_db,qs,((1.0*qs)/(query_sequences_count*seq_count_db))*100);
		*ptr_assigned = assigned;
		*total_done = 0; //TO-DO implement this

}

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
		){
		double tick;
		uint32_t filtered_Q;
		//int assigned[query_sequences_count*seq_count_db+seq_count_db];
		int* assigned = (int*)malloc(sizeof(int)*(query_sequences_count*seq_count_db+seq_count_db));
		//initialize to zero or else we wont know if it was a hit or not
		memset(assigned,0,sizeof(int)*query_sequences_count*seq_count_db+seq_count_db);
		int mic_no;
		int vector_length = VECTOR_LENGTH;
		//Due to offloads structure, cannot offload double references
		char* a __attribute__ ((aligned(64)))=seq_db; //only 1
		uint32_t* a_disp = seq_disp_db;
		uint16_t* a_len = seq_len_db;
		uint32_t offload_db_count = seq_count_db;
		mic_no = omp_get_thread_num();
		uint32_t mem_ctr = 0;
		uint32_t d_i;
		tick = getTimestamp();
		#pragma offload target(mic:mic_no)\
		in(a:length(D) )\
		in(a_disp:length(offload_db_count)  )\
		in(a_len:length(offload_db_count)  )\
		in(seq_count_db)\
		in(query_disp:length(query_sequences_count) )\
		in(query_sequences_lengths:length(query_sequences_count) )\
		in(query_sequences[0:Q]:align(64) )\
		in(query_sequences_count)\
		in(extension:length(seq_count_db))\
		in(shifted_shift:length(seq_count_db) )\
		inout(assigned:length(query_sequences_count*seq_count_db) )\
		in(shifted_copy:length(shifted_linear))
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			char* ptr_a_it __attribute__ ((aligned(64)));
			char* ptr_q_it __attribute__ ((aligned(64)));
			uint32_t block,nhits=0,seq_num=0,query_num=0;
			uint32_t i=0, j=0, chunk_it=0, gotHit=0, iter_db=0;
			int diff;
			uint8_t q1, q2, q3, q4, db1, db2, db3, db4;
			uint32_t window_it = 0;
			uint32_t window_db_iter = 0;
			int it_db = 0;
			__mmask16 mask;
			__m512i dataq __attribute__ ((aligned(64)))= _mm512_setzero_epi32();
			__m512i datadb __attribute__ ((aligned(64)))= _mm512_setzero_epi32();
			int disc = 0;
			int res;
			char q[64] __attribute__ ((aligned(64)));
			char db[64] __attribute__ ((aligned(64)));


			for(i = 0; i < query_sequences_count; i++){
				gotHit = 0;
				int real_qsize = 0;
				int more_blocks = (query_sequences_lengths[i] < 4?0:1); //minimum size is 4
				int ltblock = (query_sequences_lengths[i] < 64?1:0);
				int qsize_group =0 ;
				int next_limit;
				int iter_seq = 0;
				int q_it;
				uint32_t val_q;
				int vec_size = 0;
				int n_iter = 0;
				#pragma omp for schedule(dynamic) nowait
				for(iter_seq = 0; iter_seq < query_sequences_lengths[i]-3;iter_seq++){
					real_qsize = 64;
					qsize_group = 16;
					/*int it_b =  0;
					for(it_b = 0; it_b < 16; it_b++){
						memcpy(&q[it_b*4],&query_sequences[query_disp[i]+iter_seq],sizeof(char)*4);
					}
					dataq = _mm512_load_epi32(q);
					*/
					 val_q =query_sequences[query_disp[i]+iter_seq] + (query_sequences[query_disp[i]+iter_seq+1] << 8) + (query_sequences[query_disp[i]+iter_seq+2] << 16) + (query_sequences[query_disp[i]+iter_seq+3] << 24);
					dataq = _mm512_set4_epi32(val_q,val_q,val_q,val_q);

					for(j = 0; j < seq_count_db; j++){
						if(a_len[j] >= 4){
							 vec_size = extension[j];
							 n_iter = vec_size/64;
							for(iter_db = 0; iter_db < n_iter;iter_db++){
								int popcount = 0;
								//memcpy(db,&shifted_copy[shifted_shift[j]+(iter_db*64)],64);
								datadb = _mm512_load_epi32(&shifted_copy[shifted_shift[j]+(iter_db*64)]);
								mask = _mm512_cmpeq_epi32_mask(dataq,datadb);
								res = _mm512_mask2int(mask);


							//if(res != 0){
									//printf("Hit between:%u and %u\n",i,j);
									/*pri1ntf("qsize_group:%u\n",qsize_group);
									printf("res:%u\n",res);
									printf("Comparing:\n");
									int ter = 0;
									printf("QUERY\n");
									for(ter = 0; ter < 64 ; ter++){
										printf("%u\t",q[ter]);
									}
									printf("\n\nDB\n");
									for(ter = 0; ter < 64; ter++){
										printf("%u\t",db[ter]);
									}
									printf("\n");*/
								//}
								if(res != 0){
									popcount = popcount_hamming(res);
									assigned[seq_count_db*i+j]+=popcount;
									//break;
								}
							}
						}

					}
					//chheck more_block conditionç

				}
			}
		}
		*elapsed_time = getTimestamp()-tick;
		int r= 0;
		int k = 0;
		uint64_t qs = 0;

				for(k = 0; k < query_sequences_count; k++){
					for(r = 0; r < seq_count_db;r++){
				if(assigned[seq_count_db*k+r]>= threshold){
					qs++;
				}
			}
			//printf("\n");
		}
		printf("Non-Vectorized filter:\t\tTotal comp:%lu Done: %lu (%f %)\n",(uint64_t)query_sequences_count*seq_count_db,qs,((1.0*qs)/(query_sequences_count*seq_count_db))*100);
		*ptr_assigned = assigned;
		*total_done = 0; //TO-DO implement this

}

__declspec(target(mic))int popcount_hamming(int i)
{
     // Java: use >>> instead of >>
     // C or C++: use uint32_t
     i = i - ((i >> 1) & 0x55555555);
     i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
     return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}
//target mic moves function to xeon phi memory
/*__declspec(target(mic))uint16_t concatenate8(uint8_t x, uint8_t y) {
    unsigned pow = 10;
    while(y >= pow)
        pow *= 10;
    return x * pow + y;
}
__declspec(target(mic)) uint32_t concatenate16(uint16_t x, uint16_t y) {
    unsigned pow = 10;
    while(y >= pow)
        pow *= 10;
    return x * pow + y;
}*/
//Deprecated
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
		int *** ptr_assigned //per chunk as well now
		){

		uint32_t filtered_Q;
		uint8_t first,second,third,fourth;
		uint16_t firsthalf,secondhalf;
		uint32_t block,nseq=0,nhits=0,seq_num=0,query_num=0;

		int pos;
		int mic_no;

		int vector_length = VECTOR_LENGTH;
		int* chunk_assi[chunk_count];

		//this contains the info to which queries we'll actually execute
		//Initial offload_transfer values
		char* a ;
		char* b;

		uint32_t* a_disp ;
		uint32_t* b_disp ;
		uint16_t * m, *n;
		int* ptr_chunk_assi;
		b = chunk_b[0];
		n = chunk_n[0];
		b_disp = chunk_b_disp[0];

		// calculate maximum chunk size

		mic_no = omp_get_thread_num();
		//now for all sequences in db

		uint64_t offload_vD;
		uint32_t * ptr_chunk_b_disp, offload_vect_sequences_db_count;
		uint16_t * ptr_chunk_n;
		char * ptr_chunk_b;
		uint32_t i,j;
		printf("we have %u chunks\n",chunk_count);
		for(i = 0; i < chunk_count; i++){
			ptr_chunk_assi = (int*) _mm_malloc(sizeof(int)*chunk_vect_sequences_db_count[i]*query_sequences_count,64);
			offload_vD = chunk_vD[i];
			offload_vect_sequences_db_count = chunk_vect_sequences_db_count[i];
			ptr_chunk_b = chunk_b[i];
			ptr_chunk_n = chunk_n[i];
			ptr_chunk_b_disp = chunk_b_disp[i];


			//ideally we could dedicate a xeon phi to each block, instead if do it iteratively which is a bit resource heavy
			#pragma offload target(mic:mic_no) \
			inout(ptr_chunk_assi:length(offload_vect_sequences_db_count*query_sequences_count)ALLOC align(64))\
			in(ptr_chunk_b:length(offload_vD) ALLOC align(64)) \
			in(ptr_chunk_n:length(offload_vect_sequences_db_count) ALLOC align(64)) \
			in(ptr_chunk_b_disp:length(offload_vect_sequences_db_count) ALLOC align(64)) \
			in(query_sequences: length(Q) ALLOC align(64))\
			in(query_disp: length(query_sequences_count) ALLOC align(64))\
			in(offload_vect_sequences_db_count,offload_vD,query_sequences_count )\
			in(query_sequences_lengths:length(query_sequences_count)ALLOC align(64))\
			in(i)\
			inout(nhits)
			#pragma omp parallel num_threads(200)
			{
				//int tid = omp_get_thread_num();
				//printf("Im thread:%u and im doing chunk %u\n",tid,i);
				uint32_t j,iter_seq,seq_i;
				int gotHit = 0;
				int diff;
				char db1,db2,db3,db4,q1,q2,q3,q4;
				#pragma omp for
				for(j = 0; j < offload_vect_sequences_db_count; j++){
					uint32_t iter_db;

					for(iter_db = 0; iter_db+3 < ptr_chunk_n[j];iter_db++){
						//printf("%u ",ptr_chunk_b[iter_db+ptr_chunk_b_disp[j]]);
						diff = ptr_chunk_n[j] - iter_db;
						if(diff >= 4){
							db1 = *(ptr_chunk_b+iter_db+ptr_chunk_b_disp[j]);
							db2 = *(ptr_chunk_b+iter_db+ptr_chunk_b_disp[j]+1);
							db3 = *(ptr_chunk_b+iter_db+ptr_chunk_b_disp[j]+2);
							db4 = *(ptr_chunk_b+iter_db+ptr_chunk_b_disp[j]+3);

							for(seq_i = 0; seq_i < query_sequences_count; seq_i++){
								gotHit = 0;
								for(iter_seq = 0; iter_seq+3 < query_sequences_lengths[seq_i];iter_seq+=1){
									diff = query_sequences_lengths[seq_i] - iter_seq;
									if(diff >= 4){
										q1 = *(query_sequences+query_disp[seq_i]+iter_seq);
										q2 = *(query_sequences+query_disp[seq_i]+iter_seq+1);
										q3 = *(query_sequences+query_disp[seq_i]+iter_seq+2);
										q4 = *(query_sequences+query_disp[seq_i]+iter_seq+3);
										if(db1 == q1 && db2 == q2 && db3 == q3 && db4 == q4){
											gotHit = 1;
										}

									}
								}
								if(gotHit == 1){
									//printf("Chunk:%u DB:%u Q:%u\n",i,j,seq_i);
									ptr_chunk_assi[query_sequences_count*seq_i+j] = 1;
#pragma omp critical
									nhits++;

								}else{
									ptr_chunk_assi[query_sequences_count*seq_i+j] = 0;

								}
								gotHit = 0;


							}

						}

					}

				}


			}

			chunk_assi[i] = ptr_chunk_assi;
		}

		int c = 0;
		int r= 0;
		int k = 0;
		uint32_t qs = 0;
		uint32_t tc = 0;
		printf("Nhits:%u\n",nhits);
		for(c = 0; c < chunk_count; c++){
			tc += chunk_vect_sequences_db_count[c];
			printf("Seq:%u\n",chunk_vect_sequences_db_count[c]);
			for(r = 0; r < chunk_vect_sequences_db_count[c];r++){
				for(k = 0; k < query_sequences_count; k++){
					int* ass = chunk_assi[c];
					if(ass[query_sequences_count*k+r]==1){
						//printf("Db:%u q:%u assigned :%u\n",r,k,assigned[query_sequences_count*k+r]);
						qs++;
					}
				}
			}
		}

		printf("Chunked Non-Vectorized filter:\tTotal comp:%u Done: %lu (%f %)\n",query_sequences_count*tc,
						qs,((1.0*qs)/(query_sequences_count*chunk_vect_sequences_db_count[0]))*100);

		*ptr_assigned = chunk_assi;

}




void *malloc_aligned(size_t alignment, size_t bytes)
{
    // we need to allocate enough storage for the requested bytes, some
    // book-keeping (to store the location returned by malloc) and some extra
    // padding to allow us to find an aligned byte.  im not entirely sure if
    // 2 * alignment is enough here, its just a guess.
    const size_t total_size = bytes + (2 * alignment) + sizeof(size_t);

    // use malloc to allocate the memory.
    char *data = malloc(sizeof(char) * total_size);

    if (data)
    {
        // store the original start of the malloc'd data.
        const void * const data_start = data;

        // dedicate enough space to the book-keeping.
        data += sizeof(size_t);

        // find a memory location with correct alignment.  the alignment minus
        // the remainder of this mod operation is how many bytes forward we need
        // to move to find an aligned byte.
        const size_t offset = alignment - (((size_t)data) % alignment);

        // set data to the aligned memory.
        data += offset;

        // write the book-keeping.
        size_t *book_keeping = (size_t*)(data - sizeof(size_t));
        *book_keeping = (size_t)data_start;
    }

    return data;
}

void free_aligned(void *raw_data)
{
    if (raw_data)
    {
        char *data = raw_data;

        // we have to assume this memory was allocated with malloc_aligned.
        // this means the sizeof(size_t) bytes before data are the book-keeping
        // which points to the location we need to pass to free.
        data -= sizeof(size_t);

        // set data to the location stored in book-keeping.
        data = (char*)(*((size_t*)data));

        // free the memory.
        free(data);
    }
}


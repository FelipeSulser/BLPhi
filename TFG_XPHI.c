#include "TFG_XPHI.h"

// Global options
char *db_file=NULL;
char *query_file=NULL;
char *output_file=NULL;
char *scoremat=blosum45;
char prof=0;

int vector_length=VECTOR_LENGTH;
int num_threads= NUM_THREADS;
int open_gap=OPEN_GAP;
int extend_gap=EXTEND_GAP;

int hit_threshold = HIT_THRESHOLD;
uint16_t query_length_threshold = QUERY_LENGTH_THRESHOLD;

uint64_t max_db_size=MAX_CHUNK_SIZE, top=TOP;

char *action = NULL;

int main(int argc, char *argv[]) {

	Query * query_s = malloc(sizeof * query_s);
	Chunked_Database* ch_database_s = malloc(sizeof * ch_database_s);
	Database* database_s = malloc(sizeof* database_s);


	uint64_t sequences_count =0;
	uint64_t db_linear_size;
	uint64_t vect_sequences_db_count;
	uint64_t vD;
	uint64_t * block_inter_disp;
	uint64_t * vect_sequences_db_disp;
	uint64_t query_sequences_count;
	uint64_t query_linear_size;
	uint64_t total_done = 0;
	uint64_t shifted_total_linear;

	uint32_t num_blocks;
	uint32_t * chunk_vect_sequences_db_count;
	uint32_t** chunk_vect_sequences_db_disp;
	uint32_t* query_sequences_disp;

	int max_title_length;

	int *scores;
	int *ptr_assigned_single;
	int **ptr_assigned_multiple;
	int* extension;
	int* newdata_shift;

	unsigned char* query_seq_red;
	uint64_t Q_red;
	uint32_t * query_disp_red;
	uint16_t* query_lengths_red;

	unsigned char* db_seq_red;
	uint64_t D_red;
	uint32_t* db_disp_red;
	uint16_t* db_lengths_red;

	uint16_t ** block_db_len __attribute__ ((aligned(64)));
	uint16_t* vect_sequences_db_lengths __attribute__ ((aligned(64)));
	uint16_t sequences_db_max_length __attribute__ ((aligned(64)));
	uint16_t* query_sequences_lengths __attribute__ ((aligned(64)));
	uint16_t* real_seq_len  __attribute__ ((aligned(64)));
	uint32_t * real_seq_disp __attribute__ ((aligned(64)));

	char* real_seq;
	//The new datastructure to do the vectorized approach
	char* newdatastruc __attribute__ ((aligned(64)));

	char ** chunk_vect_sequences_db __attribute__ ((aligned(64)));
	char *vect_sequences_db __attribute__ ((aligned(64)));
	char *query_sequences __attribute__ ((aligned(64)));
	char **query_headers __attribute__ ((aligned(64)));
	char **sequence_db_headers __attribute__ ((aligned(64)));
	char **tmp_sequence_db_headers __attribute__ ((aligned(64)));
    time_t time_init = time(NULL);

	double elapsed_time, tick,filter_time,start_encoding_time,end_encoding_time;
	//argument parsing variables
		int c;
		extern char* optarg;
		extern int optind;
		int mflag = 0;
		int fflag = 0;
		int tflag = 0;
		int qflag = 0;
		int oflag = 0;
		int vflag = 0;
		int red_flag = 0;
		int bflag = 0;
		int eflag = 0;
		int sflag = 0;
		int xflag = 0;
		int zflag = 0;
		int lflag = 0;
		while((c = getopt(argc,argv,"m:f:q:o:b:rvesxztl:")) != -1){
				switch(c){
				//modus of operation
				//can be -m search
				//can be -m format
				case 'm':
					mflag = 1;
					action = optarg;
					break;
				case 'f':
					fflag = 1;
					db_file = optarg;
					break;
				case 'q':
					qflag = 1;
					query_file = optarg;
					break;
				case 'r':
					red_flag = 1;
					break;
				case 'o':
					oflag = 1;
					output_file = optarg;
					break;
				case 'b':
					bflag = 1;
					max_db_size =  atoi(optarg)*1024; //in kB
					break;
				case 'v':
					vflag = 1;
					break;
				case 'e':
					eflag = 1;
					break;
				case 's':
					sflag = 1;
					break;
				case 'x':
					xflag = 1;
					break;
				case 'z':
					zflag = 1;
					break;
				case 't':
					tflag = 1;
					break;
				case 'l':
					lflag = 1;
					hit_threshold = atoi(optarg);
					break;
				default:
					printf("Error\nUsage: ./SW_X ...\n");
					exit(1);
					break;
					}
				}

	//preprocess the fasta database
	//in our case the swissprot
	if(strcmp(action,"makedb") == 0)
		convert_db (db_file,output_file);
	else {
		if(xflag){

			double startTime = getTimestamp();
			load_query_sequences(query_file,&query_s);

			query_sequences = query_s->q_seq;
			query_sequences_lengths = query_s->query_seq_len;
			query_sequences_count = query_s->query_sequences_count;
			query_headers = query_s->query_headers;
			query_linear_size = query_s->lin_len_total;
			query_sequences_disp = query_s->query_seq_disp;


			divide_db (db_file,
					max_db_size,
					&sequences_db_max_length,
					&max_title_length,
					&database_s,
					&ch_database_s);

				sequences_count = database_s->sequences_count;
				db_linear_size = database_s->D;
				real_seq_len = database_s->sequences_lengths;
				real_seq_disp = database_s->sequences_disp;
				real_seq = database_s->real_seq;
				vD = database_s->vD;

				vect_sequences_db_count = ch_database_s->vect_sequences_db_count;
				chunk_vect_sequences_db = ch_database_s->chunk_b;
				num_blocks = ch_database_s->chunk_count;
				chunk_vect_sequences_db_count = ch_database_s->chunk_vect_sequences_db_count;
				block_inter_disp = ch_database_s->chunk_vD;
				block_db_len = ch_database_s->chunk_n;
				chunk_vect_sequences_db_disp = ch_database_s->chunk_b_disp;

			if(zflag){
			printf("Creating the shifted copy\n");
			create_shifted_copy_nopad(real_seq,
							db_linear_size,
							sequences_count,
							real_seq_disp,
							real_seq_len,
							&newdatastruc,
							&extension,
							&newdata_shift,
							&shifted_total_linear);

			printf("Finished creating shifted copy\n");
			}

			printf("NUM BLOCKS:\t\t\t%u\n",num_blocks);

			Assignment assi;
			assi = (assignment_tuple*)malloc(sizeof(assignment_tuple));
			assi -> orig[0] = 'A';
			assi -> new[0] = 'A';
			assi -> orig[1] = 'R';
			assi -> new[1] = 'R';
			assi -> orig[2] = 'N';
			assi -> new[2] = 'N';
			assi -> orig[3] = 'D';
			assi -> new[3] = 'D';
			assi -> orig[4] = 'C';
			assi -> new[4] = 'Y';
			assi -> orig[5] = 'Q';
			assi -> new[5] = 'Q';
			assi -> orig[6] = 'E';
			assi -> new[6] = 'E';
			assi -> orig[7] = 'G';
			assi -> new[7] = 'G';
			assi -> orig[8] = 'H';
			assi -> new[8] = 'M';
			assi -> orig[9] = 'I';
			assi -> new[9] = 'I';
			assi -> orig[10] = 'L';
			assi -> new[10] = 'L';
			assi -> orig[11] = 'K';
			assi -> new[11] = 'K';
			assi -> orig[12] = 'M';
			assi -> new[12] = 'M';
			assi -> orig[13] = 'F';
			assi -> new[13] = 'Q';
			assi -> orig[14] = 'P';
			assi -> new[14] = 'P';
			assi -> orig[15] = 'S';
			assi -> new[15] = 'S';
			assi -> orig[16] = 'T';
			assi -> new[16] = 'T';
			assi -> orig[17] = 'W';
			assi -> new[17] = 'Y';
			assi -> orig[18] = 'Y';
			assi -> new[18] = 'Y';
			assi -> orig[19] = 'V';
			assi -> new[19] = 'V';
			assi -> orig[20] = 'B';
			assi -> new[20] = 'Y';
			assi -> orig[21] = 'Z';
			assi -> new[21] = 'Y';
			assi -> orig[22] = 'X';
			assi -> new[22] = 'Y';
			assi->reduced_set[0] = 'A';
			assi->reduced_set[1] = 'R';
			assi->reduced_set[2] = 'N';
			assi->reduced_set[3] = 'D';
			assi->reduced_set[4] = 'Y';
			assi->reduced_set[5] = 'Q';
			assi->reduced_set[6] = 'E';
			assi->reduced_set[7] = 'G';
			assi->reduced_set[8] = 'M';
			assi->reduced_set[9] = 'I';
			assi->reduced_set[10] = 'L';
			assi->reduced_set[11] = 'K';
			assi->reduced_set[12] = 'P';
			assi->reduced_set[13] = 'S';
			assi->reduced_set[14] = 'T';
			assi->reduced_set[15] = 'V';

		if(red_flag){
			start_encoding_time = getTimestamp();
			change_encoding_query(
					query_sequences,
					query_sequences_count,
					query_sequences_disp,
					query_linear_size,
					query_sequences_lengths,
					assi,
					&query_seq_red,
					&Q_red,
					&query_disp_red,
					&query_lengths_red
				);
			change_encoding_db(
					real_seq,
					sequences_count,
					real_seq_disp,
					db_linear_size,
					real_seq_len,
					assi,
					&db_seq_red,
					&D_red,
					&db_disp_red,
					&db_lengths_red
					);

			end_encoding_time = getTimestamp();
			printf("Encoding done in:\t\t%f seconds\n",end_encoding_time - start_encoding_time);
		}

		load_database_headers (db_file, sequences_count, max_title_length, &sequence_db_headers);
		if(num_blocks < 2){
			if(vflag){
				printf("Executing vectorized version\n");

				if(red_flag){

				}else{

				filter_sequences_single_db_vectorized (query_sequences,
						query_headers,
						sequence_db_headers,
						query_sequences_lengths,
						query_sequences_count,
						query_linear_size,
						query_sequences_disp,
						real_seq,
						db_linear_size,
						real_seq_len,
						sequences_count,
						real_seq_disp,
						&ptr_assigned_single,
						&filter_time,
						&total_done
						);
				}
			}else{
				if(eflag)
				filter_sequences_single_db_8(query_sequences,
									query_headers,
									sequence_db_headers,
									query_sequences_lengths,
									query_sequences_count,
									query_linear_size,
									query_sequences_disp,
									real_seq,
									real_seq_len,
									sequences_count,
									real_seq_disp,
									&ptr_assigned_single
							);
				else{
					if(red_flag){
						printf("Executing reduced filter\n");
						filter_sequences_single_db_reduced(
						query_seq_red,
						query_headers,
						sequence_db_headers,
						query_lengths_red,
						query_sequences_count,
						Q_red,
						query_disp_red,
						db_seq_red,
						D_red,
						db_lengths_red,
						sequences_count,
						db_disp_red,
						&ptr_assigned_single,
						&filter_time,
						&total_done);
					}else{
						if(zflag){
							printf("Filter mode:\t\t\tVectorized with threshold\n");
							filter_sequences_shifted_threshold_vectorized (
									query_sequences,
									query_headers,
									sequence_db_headers,
									query_sequences_lengths,
									query_sequences_count,
									query_linear_size,
									query_sequences_disp,
									real_seq,
									db_linear_size,
									real_seq_len,
									sequences_count,
									real_seq_disp,
									&ptr_assigned_single,
									&filter_time,
									&total_done,
									newdatastruc,
									extension,
									newdata_shift,
									shifted_total_linear,
									hit_threshold
									);
						}else{
							printf("Filter mode:\t\t\tNon-Vectorized\n");
							filter_sequences_single_db(query_sequences,
							query_headers,
							sequence_db_headers,
							query_sequences_lengths,
							query_sequences_count,
							query_linear_size,
							query_sequences_disp,
							real_seq,
							db_linear_size,
							real_seq_len,
							sequences_count,
							real_seq_disp,
							&ptr_assigned_single,
							&filter_time,
							&total_done);
						}



					}


				}

			}

		}
		else{
			printf("Filtering multi sequences\n");
			filter_sequences_multiple_db (query_sequences,
							query_headers,
							query_sequences_lengths,
							query_sequences_count,
							query_linear_size,
							query_sequences_disp,
							vect_sequences_db_count,
							chunk_vect_sequences_db,
							num_blocks,
							chunk_vect_sequences_db_count,
							block_db_len,
							chunk_vect_sequences_db_disp,
							block_inter_disp,
							&ptr_assigned_multiple
					);
		}





		// Print database
		printf("Filtering time:\t\t\t%lf seconds\n",filter_time);
		printf("Database size:\t\t\t%ld sequences (%ld residues) \n",sequences_count,db_linear_size);
		printf("Longest database sequence: \t%d residues\n",sequences_db_max_length);
		printf("Gap open penalty:\t\t%d\n",open_gap);
		printf("Gap extend penalty:\t\t%d\n",extend_gap);
		printf("Query filename:\t\t\t%s\n",query_file);
		query_length_threshold = query_sequences_lengths[query_sequences_count-1]+1;
		if(sflag){

			// maximum number of shown hits
			top = (sequences_count < top ? sequences_count : top);
			scores = (int*) _mm_malloc(query_sequences_count*(vect_sequences_db_count*vector_length)*sizeof(int), 64);
			tmp_sequence_db_headers = (char**) malloc(sequences_count*sizeof(char *));

			double init_time,time_slow_smith;
			int* result;
			int* index;
			seq_pair problem;
			uint32_t iterq,iterdb;
			printf("Performing SLOW SW\n");
			init_time = getTimestamp();
			for(iterq = 0; iterq < query_sequences_count; iterq++){

				result = (int*) malloc(sizeof(int)*sequences_count);
				index = (int*)malloc(sizeof(int)*sequences_count);
				//memcpy(buffq,&query_sequences[query_sequences_disp[iterq]],query_sequences_lengths[iterq]*sizeof(char));
				#pragma omp parallel for private(problem) num_threads(8) schedule(dynamic)
				for(iterdb = 0; iterdb < sequences_count; iterdb++){
					if(ptr_assigned_single[sequences_count*iterq+iterdb] >= hit_threshold){

						//memcpy(buffdb,&real_seq[real_seq_disp[iterdb]],sizeof(char)*real_seq_len[iterdb]);
						problem.a = &query_sequences[query_sequences_disp[iterq]];
						problem.alen = query_sequences_lengths[iterq];
						problem.b = &real_seq[real_seq_disp[iterdb]];
						problem.blen = real_seq_len[iterdb];
						result[iterdb] = smith_waterman(&problem);
						index[iterdb] = iterdb;
						//printf("Score:%u\n",result[iterdb]);
					}else{
						result[iterdb] = 0;
						index[iterdb] = iterdb;
					}
				}
				memcpy(tmp_sequence_db_headers,sequence_db_headers,sequences_count*sizeof(char *));
				sort_scores(result,tmp_sequence_db_headers,sequences_count);
				int score_j;
				//Important, show only those we have computed.
				top = min(top,total_done);
				printf("\nScore\tSequence description\n");
				for (score_j=0; score_j<top; score_j++){
					printf("%u\t%s",result[score_j],tmp_sequence_db_headers[score_j]);
				}
				free(result);
				free(index);
			}
			time_slow_smith = getTimestamp();
			printf("\nSearch date:\t\t\t%s",ctime(&time_init));
			printf("Search time:\t\t\t%lf seconds\n",time_slow_smith - init_time);
			printf("Total time:\t\t\t%lf seconds\n",time_slow_smith - startTime);

			// Free allocated memory
			//_mm_free(query_sequences_lengths);
			uint64_t i;
			for (i=0; i<query_sequences_count ; i++ )
				free(query_headers[i]);
			free(query_headers);
			for (i=0; i<sequences_count ; i++ )
				free(sequence_db_headers[i]);
			free(sequence_db_headers);
			free(tmp_sequence_db_headers);
		}else{
			//free the datastructure used in filtering
			if(zflag){
				_mm_free(newdatastruc);
				free(extension);
				free(newdata_shift);
				}
			double ini_sw = getTimestamp();
			//here, feed the algorithm the comparisons
			uint32_t vect_seq_count = chunk_vect_sequences_db_count[0];
			uint64_t q_i = 0;
			uint64_t db_i = 0;
			uint64_t iter_seq = 0;
			uint64_t k;
			uint64_t tmp_db_D =0;
			uint64_t tmp_total = 0;
			uint64_t lin_ctr = 0;
			uint64_t tmp_vD =0;
			uint64_t tmp_total_vec = 0;
			uint64_t lin_size = 0;
			uint64_t shifter = 0;
			int tmp_count;
			uint64_t* original_index;
			char** seqs;
			uint64_t tmp_Q = 0;

			char* tmp_q_seq;

			//tflag is done for testing purposes,sets assigned to 1 so that it does all comparisons
			//used to compare times
			if(tflag){
				printf("t-flag enabled\n");
				uint64_t it_ass = 0;
				for(it_ass = 0; it_ass < sequences_count*query_sequences_count+sequences_count;it_ass++){
					ptr_assigned_single[it_ass] = 1;
				}
			}

			for(q_i = 0; q_i < query_sequences_count; q_i++){
				double time_start_sw = getTimestamp();
				tmp_q_seq = (char*) malloc(sizeof(char)*query_sequences_lengths[q_i]);
				tmp_Q = query_sequences_lengths[q_i]*sizeof(char);

				memcpy(tmp_q_seq,query_sequences+query_sequences_disp[q_i],query_sequences_lengths[q_i]*sizeof(char));
				tmp_total = 0;
				lin_size = 0;
				shifter = 0;
				tmp_vD = 0;
				shifter = 0;
				//printf("TMP_TOTAL= %u\n",tmp_total);


				//fetches the hits of the row of the matrix that is assigned to the current q_i
				for(db_i = 0; db_i < sequences_count; db_i++){
					if(ptr_assigned_single[sequences_count*q_i+db_i] >= hit_threshold){
						//printf("%u %u\n",q_i,db_i);
						tmp_db_D += block_db_len[0][db_i];
						tmp_total++;
					}
				}

				original_index = (uint64_t*) malloc(sizeof(uint64_t)*tmp_total);
				tmp_total_vec = ceil( (double) (tmp_total) / (double) vector_length);
				uint16_t *tmp_seq_len = (uint16_t*) malloc(tmp_total*sizeof(uint16_t));
				uint16_t* tmp_vect_sequences_len = (uint16_t*) malloc(sizeof(uint16_t)*(tmp_total_vec));
				uint32_t* tmp_vect_sequences_disp = (uint32_t*) malloc(sizeof(uint32_t)*(tmp_total_vec+1));
				lin_ctr = 0;

				for(db_i = 0; db_i < sequences_count; db_i++){
					if(ptr_assigned_single[sequences_count*q_i+db_i] >= hit_threshold){
						original_index[lin_ctr] = db_i; //surjective
						//printf("%s\n",sequence_db_headers[db_i]);
						tmp_seq_len[lin_ctr] = real_seq_len[db_i];
						lin_ctr++;
					}
				}


				for(db_i = 0; db_i < tmp_total_vec - 1; db_i++){
					tmp_vect_sequences_len[db_i] =tmp_seq_len[(db_i+1)*vector_length-1];

				}
				tmp_vect_sequences_len[tmp_total_vec-1] = tmp_seq_len[tmp_total-1];
				for(db_i = 0; db_i < tmp_total_vec; db_i++){
					tmp_vect_sequences_len[db_i] = ceil((double)tmp_vect_sequences_len[db_i] /(double) 4)*4;


				}

				tmp_vect_sequences_disp[0] = 0;
				for(db_i = 1; db_i < tmp_total_vec+1; db_i++){
					tmp_vect_sequences_disp[db_i] = tmp_vect_sequences_disp[db_i-1] + (vector_length*tmp_vect_sequences_len[db_i-1]);

				}
				lin_size = 0;
				#pragma omp parallel for reduction(+:lin_size) num_threads(4)
				for(db_i = 0; db_i < tmp_total; db_i++)
					lin_size = lin_size + tmp_seq_len[db_i];

				tmp_vD = 0;
				#pragma omp parallel for reduction(+:tmp_vD) num_threads(4)
				for(db_i = 0; db_i < tmp_total_vec; db_i++){
					tmp_vD = tmp_vD + tmp_vect_sequences_len[db_i]*vector_length;

				}


				char* tmp_seq_db __attribute__ ((aligned(64)))= (char*) _mm_malloc(sizeof(char)*lin_size,64);

				//printf("LIN_SIZE:%lu\n",lin_size);
				lin_ctr = 0;
				shifter = 0;
				//here we copy the normal db
				for(db_i = 0; db_i < sequences_count; db_i++){
					if(ptr_assigned_single[sequences_count*q_i+db_i] >= hit_threshold){
						memcpy(tmp_seq_db+shifter,real_seq+real_seq_disp[db_i],real_seq_len[db_i]*sizeof(char));
						shifter+=real_seq_len[db_i];
					}
				}
				shifter = 0;
				seqs = (char**) malloc(sizeof(char*)*tmp_total);
				seqs[0] = tmp_seq_db;
				for(db_i = 1; db_i < tmp_total; db_i++){
					seqs[db_i] = seqs[db_i-1]+tmp_seq_len[db_i-1];
				}

				char* tmp_vec_seq_db __attribute__ ((aligned(64)))= (char*) _mm_malloc(sizeof(char)*tmp_vD,16);
				for(db_i = 0; db_i < tmp_total_vec-1; db_i++){
					for(iter_seq = 0; iter_seq < tmp_vect_sequences_len[db_i];iter_seq++){
						for(k = 0; k < vector_length; k++){
							if(iter_seq < tmp_seq_len[db_i*vector_length+k])
								*(tmp_vec_seq_db+tmp_vect_sequences_disp[db_i]+(iter_seq*vector_length)+k) = seqs[db_i*vector_length+k][iter_seq];
							else
								*(tmp_vec_seq_db+tmp_vect_sequences_disp[db_i]+(iter_seq*vector_length)+k) = UNIFY_NUMERIC;
						}
					}
				}

			for(db_i = tmp_total_vec-1,iter_seq = 0;iter_seq < tmp_vect_sequences_len[db_i];iter_seq++){
				for(k = 0; k < vector_length; k++){
					if(db_i*vector_length+k < tmp_total){
						if(iter_seq < tmp_seq_len[db_i*vector_length+k])
							*(tmp_vec_seq_db+tmp_vect_sequences_disp[db_i]+(iter_seq*vector_length)+k) = seqs[db_i*vector_length+k][iter_seq];
						else
							*(tmp_vec_seq_db+tmp_vect_sequences_disp[db_i]+(iter_seq*vector_length)+k) = UNIFY_NUMERIC;
					}else
						*(tmp_vec_seq_db+tmp_vect_sequences_disp[db_i]+(iter_seq*vector_length)+k) = UNIFY_NUMERIC;
				}
			}
			uint32_t* tmp_chunk_count = (uint32_t*) _mm_malloc(sizeof(uint32_t),64);
			uint32_t* tmp_chunk_disp = (uint32_t*) _mm_malloc(sizeof(uint32_t),64);
			uint16_t* tmp_query_seq_len __attribute__ ((aligned(64)))= (uint16_t*) _mm_malloc(sizeof(uint16_t),64);
			tmp_query_seq_len[0] = query_sequences_lengths[q_i];
			tmp_chunk_disp[0] = 0;
			tmp_chunk_count[0] = tmp_total_vec;

			tmp_count = 1;
			uint32_t* tmp_chunk_vect_sequences_count = (uint32_t*) malloc((tmp_count)*sizeof(uint32_t));

			db_i = 0;
			k = 0;
			uint64_t accum;
			uint32_t chunk_size= 0;
			while(db_i < tmp_total_vec){
				iter_seq = 0;
				chunk_size = 0;
				accum = tmp_vect_sequences_len[db_i]*vector_length*sizeof(char)+sizeof(uint16_t)+ sizeof(uint32_t);
				while((db_i < tmp_total_vec)){
					chunk_size += accum;
					iter_seq++;
					db_i++;
					if(db_i < tmp_total_vec)
						accum = tmp_vect_sequences_len[db_i]*vector_length*sizeof(char)+sizeof(uint16_t)+sizeof(uint32_t);
				}

				tmp_chunk_vect_sequences_count[k] = iter_seq;
				tmp_count++;
				k++;
				tmp_chunk_vect_sequences_count = (uint32_t*) realloc(tmp_chunk_vect_sequences_count,tmp_count*sizeof(uint32_t));

			}
			tmp_count--;

			uint16_t** tmp_chunk_vect_sequences_len = (uint16_t**)_mm_malloc(tmp_count*sizeof(uint16_t*),64);
			uint64_t offset = 0;
			for(db_i = 0; db_i < tmp_count; db_i++){
				tmp_chunk_vect_sequences_len[db_i] = (uint16_t*) _mm_malloc(tmp_chunk_vect_sequences_count[db_i]*sizeof(uint16_t),64);
				memcpy(tmp_chunk_vect_sequences_len[db_i],tmp_vect_sequences_len+offset,(tmp_chunk_vect_sequences_count[db_i]*sizeof(uint16_t)));
				offset += tmp_chunk_vect_sequences_count[db_i];

			}

			accum = 0;
			uint32_t **tmp_chunk_vect_sequences_disp = (uint32_t**)_mm_malloc(tmp_count*sizeof(uint32_t*),64);
			for(db_i = 0; db_i < tmp_count; db_i++){
				tmp_chunk_vect_sequences_disp[db_i] = (uint32_t*) _mm_malloc(tmp_chunk_vect_sequences_count[db_i]*sizeof(uint32_t),64);
				offset = tmp_vect_sequences_disp[accum];
				for(iter_seq = 0, k = accum; iter_seq < tmp_chunk_vect_sequences_count[db_i];iter_seq++,k++){
					tmp_chunk_vect_sequences_disp[db_i][iter_seq] = (uint32_t)(tmp_vect_sequences_disp[k] - offset);
				}
				accum += tmp_chunk_vect_sequences_count[db_i];
			}

			uint64_t* tmp_chunk_vD = (uint64_t*) malloc(tmp_count*sizeof(uint64_t));

			offset = 0;
			uint64_t ii;
			for(db_i = 0; db_i < tmp_count; db_i++){
				ii = offset + tmp_chunk_vect_sequences_count[db_i];
				tmp_chunk_vD[db_i] = tmp_vect_sequences_disp[ii] - tmp_vect_sequences_disp[offset];
				offset = ii;
			}

			char** tmp_chunk_vect_sequences_db = (char**) _mm_malloc(sizeof(char*)*tmp_count,64);

			offset = 0;
			for(db_i = 0; db_i < tmp_count; db_i++){
				tmp_chunk_vect_sequences_db[db_i] = tmp_vec_seq_db+offset;
				offset += tmp_chunk_vD[db_i];
			}
			scores = (int*) _mm_malloc(query_sequences_count*(tmp_total_vec*vector_length)*sizeof(int), 64);
			Query* tmp_q_s = malloc(sizeof * tmp_q_s);
			tmp_q_s->q_seq = tmp_q_seq;
			tmp_q_s->lin_len_total = tmp_Q;
			tmp_q_s->query_seq_disp = query_sequences_disp;
			tmp_q_s->query_seq_len = tmp_query_seq_len;
			tmp_q_s->query_sequences_count = 1;

			Chunked_Database * tmp_ch = malloc(sizeof * tmp_ch);
			tmp_ch->chunk_b = tmp_chunk_vect_sequences_db;
			tmp_ch->chunk_b_disp = tmp_chunk_vect_sequences_disp;
			tmp_ch->chunk_count = tmp_count;
			tmp_ch->chunk_n = tmp_chunk_vect_sequences_len;
			tmp_ch->chunk_vD = tmp_chunk_vD;
			tmp_ch->chunk_vect_sequences_db_count = tmp_chunk_vect_sequences_count;


			double difft = getTimestamp() - time_start_sw;
			smith_waterman_vectorized (
					tmp_q_s,
					tmp_ch,
					scoremat,
					open_gap,
					extend_gap,
					num_threads,
					scores,
					&elapsed_time,
					query_length_threshold);

			top = (tmp_total < top ? tmp_total : top);
			tmp_sequence_db_headers = (char**) malloc(sizeof(char*)* tmp_total);
			//copy to buffer, order them and show them on console
			for(db_i = 0; db_i < tmp_total; db_i++){
				tmp_sequence_db_headers[db_i] = (char*) malloc(sizeof(char*)*strlen(sequence_db_headers[original_index[db_i]])+1);
				strcpy(tmp_sequence_db_headers[db_i],sequence_db_headers[original_index[db_i]]);
				//printf("%s\n",tmp_sequence_db_headers[db_i]);
			}

			sort_scores(scores,tmp_sequence_db_headers,tmp_total);
			printf("\nQuery no.\t\t\t%d\n",q_i+1);
			printf("Query description: \t\t%s\n",query_headers[q_i]);
			printf("Query length:\t\t\t%d residues\n",query_sequences_lengths[q_i]);
			printf("\nScore\tSequence description\n");
			uint32_t i,j;
			for (j=0; j<top; j++)
				printf("%d\t%s",scores[j],tmp_sequence_db_headers[j]+1);


			printf("\nPreprocess time:\t\t%lf seconds\n",difft);
			printf("Search time:\t\t\t%lf seconds\n",elapsed_time);
			for(i = 0; i < tmp_total;i++){
				free(tmp_sequence_db_headers[i]);
			}
			free(tmp_sequence_db_headers);
			_mm_free(scores);
			_mm_free(tmp_chunk_disp);
			_mm_free(tmp_chunk_count);
			_mm_free(tmp_seq_db);
			free(tmp_vect_sequences_disp);
			free(tmp_vect_sequences_len);
			free(tmp_seq_len);
			free(original_index);
			free(tmp_q_seq);
			}
			printf("Total SW execution time:\t%lf\n",getTimestamp()-ini_sw);
			printf("Search date:\t\t\t%s",ctime(&time_init));
			printf("TOTAL TIME:\t%lf\n",getTimestamp()-startTime);

		}

	}else{
		//Execute SW without any filters
		double startTime = getTimestamp();
				load_query_sequences(query_file,
						&query_s);

				query_sequences = query_s->q_seq;
				query_sequences_lengths = query_s->query_seq_len;
				query_sequences_count = query_s->query_sequences_count;
				query_headers = query_s->query_headers;
				query_linear_size = query_s->lin_len_total;
				query_sequences_disp = query_s->query_seq_disp;

				/*
				divide_db (db_file,
									max_db_size,
									&sequences_count,
									&db_linear_size,
									&sequences_db_max_length,
									&max_title_length,
									&vect_sequences_db_count,
									&vD,
									&chunk_vect_sequences_db,
									&num_blocks,
									&chunk_vect_sequences_db_count,
									&block_inter_disp,
									&block_db_len,
									&chunk_vect_sequences_db_disp,
									&real_seq_len,
									&real_seq_disp,
									&real_seq);*/
				divide_db (db_file,max_db_size,&sequences_db_max_length,&max_title_length,&database_s,&ch_database_s);

								sequences_count = database_s->sequences_count;
								db_linear_size = database_s->D;
								real_seq_len = database_s->sequences_lengths;
								real_seq_disp = database_s->sequences_disp;
								real_seq = database_s->real_seq;
								vD = database_s->vD;

								vect_sequences_db_count = ch_database_s->vect_sequences_db_count;
								chunk_vect_sequences_db = ch_database_s->chunk_b;
								num_blocks = ch_database_s->chunk_count;
								chunk_vect_sequences_db_count = ch_database_s->chunk_vect_sequences_db_count;
								block_inter_disp = ch_database_s->chunk_vD;
								block_db_len = ch_database_s->chunk_n;
								chunk_vect_sequences_db_disp = ch_database_s->chunk_b_disp;

				printf("NUM BLOCKS:\t\t\t%u\n",num_blocks);
				// maximum number of shown hits
				// maximum number of shown hits
				top = (sequences_count < top ? sequences_count : top);


				// Print database
				printf("Database size:\t\t\t%ld sequences (%ld residues) \n",sequences_count,db_linear_size);
				printf("Longest database sequence: \t%d residues\n",sequences_db_max_length);
				printf("Gap open penalty:\t\t%d\n",open_gap);
				printf("Gap extend penalty:\t\t%d\n",extend_gap);
				printf("Query filename:\t\t\t%s\n",query_file);

				query_length_threshold = query_sequences_lengths[query_sequences_count-1]+1;

				//divide query
				char* tmp_q_seq;
				uint64_t tmp_Q;
				uint64_t q_i = 0;
				uint64_t db_i = 0;
				for(q_i = 0; q_i < query_sequences_count; q_i ++){
					elapsed_time = getTimestamp();
					double time_start_sw = getTimestamp();
					tmp_q_seq = (char*) malloc(sizeof(char)*query_sequences_lengths[q_i]);
					tmp_Q = query_sequences_lengths[q_i]*sizeof(char);
					scores = (int*) _mm_malloc(query_sequences_count*(vect_sequences_db_count*vector_length)*sizeof(int), 64);
					tmp_sequence_db_headers = (char**) malloc(sequences_count*sizeof(char *));
					memcpy(tmp_q_seq,query_sequences+query_sequences_disp[q_i],query_sequences_lengths[q_i]*sizeof(char));
					uint16_t* tmp_query_seq_len __attribute__ ((aligned(64)))= (uint16_t*) _mm_malloc(sizeof(uint16_t),64);
					tmp_query_seq_len[0] = query_sequences_lengths[q_i];

					Query* tmp_q_s = malloc(sizeof * tmp_q_s);

					tmp_q_s ->lin_len_total = tmp_Q;
					tmp_q_s->q_seq = tmp_q_seq;
					tmp_q_s->query_seq_disp = query_sequences_disp;
					tmp_q_s->query_seq_len = tmp_query_seq_len;
					tmp_q_s->query_sequences_count = 1;

					smith_waterman_vectorized (tmp_q_s,
					ch_database_s,
					scoremat,
					open_gap,
					extend_gap,
					num_threads,
					scores,
					&elapsed_time,
					query_length_threshold);

					//Database headers
					uint64_t i;
					load_database_headers (db_file, sequences_count, max_title_length, &sequence_db_headers);
					// Print the sequences that scored most
					memcpy(tmp_sequence_db_headers,sequence_db_headers,sequences_count*sizeof(char *));
					sort_scores(scores,tmp_sequence_db_headers,sequences_count);
					printf("\nQuery no.\t\t\t%d\n",i+1);
					printf("Query description: \t\t%s\n",query_headers[q_i]+1);
					printf("Query length:\t\t\t%d residues\n",query_sequences_lengths[q_i]);
					printf("\nScore\tSequence description\n");

					uint32_t j;
					for (j=0; j<top; j++)
						printf("%d\t%s",scores[j],tmp_sequence_db_headers[j]+1);



					printf("\nSearch date:\t\t\t%s",ctime(&time_init));
					printf("Search time:\t\t\t%lf seconds\n",elapsed_time);


				}

			_mm_free(query_sequences);
			_mm_free(query_sequences_disp);

			_mm_free(chunk_vect_sequences_db[0]);
			_mm_free(chunk_vect_sequences_db);
			uint64_t i;
			for (i=0; i< num_blocks ; i++ )
				_mm_free(block_db_len[i]);
			_mm_free(block_db_len);
			for (i=0; i< num_blocks ; i++ )
				_mm_free(chunk_vect_sequences_db_disp[i]);
			_mm_free(chunk_vect_sequences_db_disp);
			free(chunk_vect_sequences_db_count);
			free(block_inter_disp);
			_mm_free(query_sequences_lengths);
					_mm_free(scores);
					for (i=0; i<query_sequences_count ; i++ )
						free(query_headers[i]);
					free(query_headers);
					for (i=0; i<sequences_count ; i++ )
						free(sequence_db_headers[i]);
					free(sequence_db_headers);
					free(tmp_sequence_db_headers);
		}
	}


	return 0;
}

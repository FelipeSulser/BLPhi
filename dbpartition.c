/*
 * dbpartition.c
 *
 *  Created on: 25/04/2016
 *      Author: Felipe Sulser Larraz
 */
#include "dbpartition.h"


void load_query_sequences(char * queries_filename,
		Query ** query_ret) {

	struct query_data* query = malloc(sizeof * query);

	long int i, j, k;
	uint64_t sequences_count=0, Q=0, disp, accum, chunk_size;
	uint32_t * sequences_disp;
	uint16_t *sequences_lengths, * title_lengths, *tmp, length=0, tmp_length, ok;
	char ** sequences=NULL, **titles, buffer[BUFFER_SIZE], filename[BUFFER_SIZE], * bin_filename, * res, *tmp_seq, *a, diff, new_line='\n';
	FILE * sequences_file;

	// open query sequence filename
	if((sequences_file = fopen(queries_filename,"r")) == NULL){
		fprintf(stderr,"Error while opening sequence file\n");
		exit(1);
	}
	// Allocate memory for sequences_lengths array
	sequences_lengths = (uint16_t *) malloc (BUFFER_SIZE*sizeof(uint16_t));
	title_lengths = (uint16_t *) malloc (BUFFER_SIZE*sizeof(uint16_t));

	// Calculate number of sequences in database and its lengths
	sequences_count=0;

	res = fgets(buffer,BUFFER_SIZE,sequences_file);
	while (res != NULL) {
		length = 0;
		// read title
		while (strrchr(buffer,new_line) == NULL) {
			length += strlen(buffer);
			res = fgets(buffer,BUFFER_SIZE,sequences_file);
		}
		title_lengths[sequences_count] = length + strlen(buffer) + 1;
		// read sequence
		length = 0;
		res = fgets(buffer,BUFFER_SIZE,sequences_file);
		while ((res != NULL) && (buffer[0] != '>')) {
			length += strlen(buffer)-1;
			res = fgets(buffer,BUFFER_SIZE,sequences_file);
		}
		sequences_lengths[sequences_count] = length;
		(sequences_count)++;
		if ((sequences_count) % BUFFER_SIZE == 0) {
			sequences_lengths = (uint16_t *) realloc(sequences_lengths,((sequences_count)+BUFFER_SIZE)*sizeof(uint16_t));
			title_lengths = (uint16_t *) realloc(title_lengths,((sequences_count)+BUFFER_SIZE)*sizeof(uint16_t));
		}
	}

	// copy lengths to aligned buffer
	tmp = sequences_lengths;
	sequences_lengths = (uint16_t *) _mm_malloc (sequences_count*sizeof(uint16_t), 64);
	memcpy(sequences_lengths,tmp,sequences_count*sizeof(uint16_t));
	free(tmp);

	// Allocate memory for sequences array
	sequences = (char **) malloc(sequences_count*sizeof(char *));

	for (i=0; i<sequences_count; i++ ) {
		sequences[i] = (char *) malloc(sequences_lengths[i]*sizeof(char));

	}

	// Rewind sequences database file
	rewind(sequences_file);

	// Read sequences from the database file and load them in sequences array
	i = 0;
	res = fgets(buffer,BUFFER_SIZE,sequences_file);
	while (res != NULL) {
		// read title
		while (strrchr(buffer,new_line) == NULL)
			res = fgets(buffer,BUFFER_SIZE,sequences_file);
		// read sequence
		length = 1;
		res = fgets(buffer,BUFFER_SIZE,sequences_file);
		while ((res != NULL) && (buffer[0] != '>')) {
			//printf("%s %d\n",buffer,strlen(buffer));
			strncpy(sequences[i]+(length-1),buffer,strlen(buffer)-1);
			length += strlen(buffer)-1;
			res = fgets(buffer,BUFFER_SIZE,sequences_file);
		}
		i++;
	}

	// Rewind sequences database file
	rewind(sequences_file);

	// Allocate memory for titles array
	titles = (char **) malloc(sequences_count*sizeof(char *));
	for (i=0; i<sequences_count; i++ ) {
		titles[i] = (char *) malloc(title_lengths[i]*sizeof(char));
	}

	i = 0;
	res = fgets(buffer,BUFFER_SIZE,sequences_file);
	while (res != NULL) {
		// discard sequences
		while ((res != NULL) && (buffer[0] != '>'))
			res = fgets(buffer,BUFFER_SIZE,sequences_file);
		if (res != NULL){
			// read header
			length = 1;
			do{
				strncpy(titles[i]+(length-1),buffer,strlen(buffer)-1);
				length += strlen(buffer)-1;
				res = fgets(buffer,BUFFER_SIZE,sequences_file);
			} while (strrchr(buffer,new_line) == NULL);
			titles[i][length] = '\0';
			i++;
		}
	}

	// Close sequences database file
	fclose(sequences_file);

	// Sort sequence array by length
	sort_sequences(sequences,titles,sequences_lengths, sequences_count);

	// calculate total number of residues
	#pragma omp parallel for reduction(+:Q) num_threads(PREPROCESS_THREADS)
	for (i=0; i< sequences_count; i++ )
		Q = Q +  sequences_lengths[i];


	a = (char *) _mm_malloc(Q*sizeof(char), 64);

	disp = 0;
	for (i=0; i< sequences_count; i++ ) {
		memcpy(a+disp,sequences[i],sequences_lengths[i]);
		disp += sequences_lengths[i];
	}

	// process vect sequences DB
	#pragma omp parallel for private(diff) num_threads(PREPROCESS_THREADS) schedule(dynamic)
	for (i=0; i< Q; i++) {
		a[i] = ((a[i] == 'J') ? UNIFY : a[i]);
		a[i] = ((a[i] == 'O') ? UNIFY : a[i]);
		a[i] = ((a[i] == 'U') ? UNIFY : a[i]);
		diff = 'A';
		diff = (a[i] > 'J' ? diff+1 : diff);
		diff = (a[i] > 'O' ? diff+1 : diff);
		diff = (a[i] > 'U' ? diff+1 : diff);
		a[i] -= diff;
	}

	// Calculate displacement for current sequences db
	sequences_disp = (uint32_t *) _mm_malloc((sequences_count+1)*sizeof(uint32_t), 64);

	sequences_disp[0] = 0;
	for (i=1; i < sequences_count+1; i++)
		sequences_disp[i] = sequences_disp[i-1] + sequences_lengths[i-1];
	query->q_seq = a;
	query->lin_len_total = Q;
	query->query_headers = titles;
	query->query_seq_disp = sequences_disp;
	query->query_seq_len = sequences_lengths;
	query->query_sequences_count = sequences_count;

	*query_ret = query;

	// Free memory
	for (i=0; i< sequences_count; i++ )
		free(sequences[i]);
	free(sequences);
	free(title_lengths);
}

/**
 * 	uint64_t vect_sequences_db_count;
	char ** chunk_b;
	uint32_t chunk_count;
	uint32_t * chunk_vect_sequences_db_count;
	uint16_t ** chunk_n;
	uint32_t ** chunk_b_disp;
	uint64_t * chunk_vD;
 */
void divide_db (char * sequences_filename,
		uint64_t max_chunk_size,
		uint16_t * sequences_db_max_length,
		int * max_title_length,
		Database ** non_chunked,
		Chunked_Database** chunked){


	Chunked_Database* ret_db = malloc(sizeof * ret_db);

	Database * database = malloc(sizeof * database);

	char ** sequences, *s __attribute__ ((aligned(64))), **chunk_vect_sequences_db __attribute__ ((aligned(64))), filename[200], * header, *b;
	uint16_t ** chunk_vect_sequences_db_lengths, ** chunk_n, * sequences_lengths, * vect_sequences_lengths;
	uint64_t i, ii, j, jj, k, * chunk_vD, accum, aux_vD=0, offset, chunk_size, * vect_sequences_disp;
	uint32_t * chunk_vect_sequences_count, **chunk_vect_sequences_disp, * tmp_chunk_vect_sequences_disp, c;
	FILE * sequences_file, * info_file;

	uint32_t*seq_disp;

	uint64_t * D;
	uint64_t * sequences_count;

	// Open info file
	sprintf(filename,"%s.stats",sequences_filename);


	if((info_file = fopen(filename,"r")) == NULL){
		fprintf(stderr,"Error while opening info file\n");
		exit(1);
	}
	fscanf(info_file,"%ld %ld %d",&database->sequences_count,&database->D,max_title_length);
	D = &database->D;
	sequences_count = &database->sequences_count;


	fclose(info_file);
	// Open sequences file
	sprintf(filename,"%s.bin",sequences_filename);

	if((sequences_file = fopen(filename,"r")) == NULL){
		fprintf(stderr,"Error while opening seq file\n");
		exit(1);
	}

	// Read sequences lengths
	sequences_lengths = (uint16_t *) malloc((*sequences_count)*sizeof(uint16_t));
	fread(sequences_lengths,sizeof(uint16_t),*sequences_count,sequences_file);

	seq_disp = (uint32_t*) malloc((*sequences_count)*sizeof(uint32_t));


	/**
	 * vect_sequences_disp[0] = 0;
	for (k=1; k < ret_db->vect_sequences_db_count+1; k++)
		vect_sequences_disp[k] = vect_sequences_disp[k-1] + (vect_sequences_lengths[k-1]*VECTOR_LENGTH);
	 */
	uint32_t iter_n;
	seq_disp[0] = 0;
	for(iter_n = 1; iter_n < (*sequences_count); iter_n++){
		seq_disp[iter_n] = seq_disp[iter_n-1]+(sequences_lengths[iter_n-1]);
	}
	database->sequences_disp = seq_disp;

	// Read sequences
	s = (char *) _mm_malloc((*D)*sizeof(char),64);
	fread(s,sizeof(char),*D,sequences_file);


	fclose(sequences_file);


	sequences = (char **) malloc((*sequences_count)*sizeof(char *));

	sequences[0] = s;
	for (i=1; i<*sequences_count ; i++)
		sequences[i] = sequences[i-1] + sequences_lengths[i-1];
	database->real_seq = s;
	// calculate vect_sequences_count
	ret_db->vect_sequences_db_count = ceil( (double) (*sequences_count) / (double) VECTOR_LENGTH);
	//Allocate memory for vect_sequences_lengths
	vect_sequences_lengths = (uint16_t *) malloc((ret_db->vect_sequences_db_count)*sizeof(uint16_t));

	vect_sequences_disp = (uint64_t *) malloc((ret_db->vect_sequences_db_count+1)*sizeof(uint64_t));


	// calculate values for vect_sequences_lengths array
	for (i=0; i< ret_db->vect_sequences_db_count - 1; i++ )
		vect_sequences_lengths[i] = sequences_lengths[(i+1)*VECTOR_LENGTH-1];
	vect_sequences_lengths[ret_db->vect_sequences_db_count-1] = sequences_lengths[*sequences_count-1];

	// make length multiple of 4 to allow 32/64 bytes aligned data
	for (i=0; i< ret_db->vect_sequences_db_count; i++ ){
		vect_sequences_lengths[i] = ceil( (double) vect_sequences_lengths[i] / (double) 4) * 4;
		//printf("%u %u\n",vect_sequences_lengths[i],sequences_lengths[i]);
	}

	// Calculate displacement for current sequences db
	vect_sequences_disp[0] = 0;
	for (k=1; k < ret_db->vect_sequences_db_count+1; k++)
		vect_sequences_disp[k] = vect_sequences_disp[k-1] + (vect_sequences_lengths[k-1]*VECTOR_LENGTH);

	#pragma omp parallel for reduction(+:aux_vD) num_threads(PREPROCESS_THREADS)
	for (i=0; i< ret_db->vect_sequences_db_count; i++ )
		aux_vD = aux_vD + vect_sequences_lengths[i]*VECTOR_LENGTH;
	database->vD = aux_vD;

	b = (char *) _mm_malloc((database->vD)*sizeof(char),16);

	// Copy sequences db to host buffers reordering elements to get better locality when computing alignments
	for (i=0; i < ret_db->vect_sequences_db_count-1; i++) {
		for (j=0; j< vect_sequences_lengths[i]; j++ ) {
			for (k=0;k< VECTOR_LENGTH; k++)
				if (j < sequences_lengths[i*VECTOR_LENGTH+k])
					*(b+vect_sequences_disp[i]+(j*VECTOR_LENGTH)+k) = sequences[i*VECTOR_LENGTH+k][j];
				else
					*(b+vect_sequences_disp[i]+(j*VECTOR_LENGTH)+k) = UNIFY_NUMERIC;
		}
	}
	//rest = sequences_count % VECTOR_LENGTH;
	for (i=ret_db->vect_sequences_db_count-1, j=0; j< vect_sequences_lengths[i]; j++ ) {
		for (k=0;k< VECTOR_LENGTH; k++)
			if (i*VECTOR_LENGTH+k < *sequences_count){
				if (j < sequences_lengths[i*VECTOR_LENGTH+k])
					*(b+vect_sequences_disp[i]+(j*VECTOR_LENGTH)+k) = sequences[i*VECTOR_LENGTH+k][j];
				else
					*(b+vect_sequences_disp[i]+(j*VECTOR_LENGTH)+k) = UNIFY_NUMERIC;
			} else
				*(b+vect_sequences_disp[i]+(j*VECTOR_LENGTH)+k) = UNIFY_NUMERIC;
	}

	// calculate chunks
	ret_db->chunk_count = 1;


	chunk_vect_sequences_count = (uint32_t *) malloc((ret_db->chunk_count)*sizeof(uint32_t));


	i = 0;
	c = 0;
	while (i< ret_db->vect_sequences_db_count) {
		// group sequences till reach max chunk size
		j = 0;
		chunk_size = 0;
		accum = vect_sequences_lengths[i]*VECTOR_LENGTH*sizeof(char) + sizeof(uint16_t) + sizeof(uint32_t); // secuencias + longitud + desplazamiento
		while ((i< ret_db->vect_sequences_db_count) && (chunk_size <= max_chunk_size)) {
			chunk_size += accum;
			j++;
			i++;
			if (i < ret_db->vect_sequences_db_count)
				accum = vect_sequences_lengths[i]*VECTOR_LENGTH*sizeof(char) + sizeof(uint16_t) + sizeof(uint32_t); // secuencias + longitud + desplazamiento
		}
		// number of sequences in chunk
		chunk_vect_sequences_count[c] = j;
		// increment chunk_count
		(ret_db->chunk_count)++;
		c++;
		chunk_vect_sequences_count = (uint32_t *) realloc(chunk_vect_sequences_count,(ret_db->chunk_count)*sizeof(uint32_t));


	}
	// update chunk count
	(ret_db->chunk_count)--;

	// calculate chunk_vect_sequences_db_lengths
	chunk_vect_sequences_db_lengths = (uint16_t **) _mm_malloc((ret_db->chunk_count)*sizeof(uint16_t *),64);


	offset = 0;
	for (i=0; i< ret_db->chunk_count ; i++) {
		chunk_vect_sequences_db_lengths[i] = (uint16_t *) _mm_malloc((chunk_vect_sequences_count[i])*sizeof(uint16_t),64);
		memcpy(chunk_vect_sequences_db_lengths[i],vect_sequences_lengths+offset,(chunk_vect_sequences_count[i])*sizeof(uint16_t));
		offset += chunk_vect_sequences_count[i];
	}

	// calculate chunk_vect_sequences_db_disp
	accum = 0;
	chunk_vect_sequences_disp = (uint32_t **) _mm_malloc((ret_db->chunk_count)*sizeof(uint32_t *),64);
	for (i=0; i< ret_db->chunk_count ; i++){
		chunk_vect_sequences_disp[i] = (uint32_t *) _mm_malloc(chunk_vect_sequences_count[i]*sizeof(uint32_t),64);

		// adapt sequence displacements to chunk
		offset = vect_sequences_disp[accum];
		for ( j=0, jj=accum; j<chunk_vect_sequences_count[i] ; j++, jj++)
			chunk_vect_sequences_disp[i][j] = (uint32_t)(vect_sequences_disp[jj] - offset);
		accum += chunk_vect_sequences_count[i];
	}

	// calculate chunk_vD
	chunk_vD = (uint64_t *) malloc((ret_db->chunk_count)*sizeof(uint64_t));


	offset = 0;
	for (i=0; i< ret_db->chunk_count; i++){
		ii = offset + chunk_vect_sequences_count[i];
		chunk_vD[i] = vect_sequences_disp[ii] - vect_sequences_disp[offset];
		offset = ii;
	}

	// calculate chunk_vect_sequences_db
	chunk_vect_sequences_db = (char **) _mm_malloc((ret_db->chunk_count)*sizeof(char *), 64);

	offset = 0;
	for (i=0; i< ret_db->chunk_count ; i++) {
		chunk_vect_sequences_db[i] = b + offset;
		offset += chunk_vD[i];
	}
	ret_db->chunk_b= chunk_vect_sequences_db;
	ret_db->chunk_n = chunk_vect_sequences_db_lengths;
	ret_db->chunk_vect_sequences_db_count = chunk_vect_sequences_count;
	ret_db->chunk_b_disp = chunk_vect_sequences_disp;
	ret_db->chunk_vD = chunk_vD;

	*chunked = ret_db;
	*sequences_db_max_length = sequences_lengths[*sequences_count-1];
	database->sequences_lengths = sequences_lengths;

	*non_chunked = database;

	//_mm_free(s);

	free(vect_sequences_lengths);
	free(vect_sequences_disp);
}

void load_database_headers (char * sequences_filename, uint64_t sequences_count, int max_title_length, char *** ptr_sequences_db_headers) {

	char ** sequences_db_headers, filename[200], * header;
	FILE * header_file;
	uint64_t i;

	// Load sequence headers

	// Open header file
	sprintf(filename,"%s.title",sequences_filename);
	if((header_file = fopen(filename,"r")) == NULL){
		fprintf(stderr,"Erorr while opening title file\n");
		exit(1);
	}

	// Read sequences lengths
	sequences_db_headers = (char **) malloc(sequences_count*sizeof(char *));
	header = (char *) malloc((max_title_length+1)*sizeof(char));

	for (i=0; i<sequences_count; i++){
		fgets(header,max_title_length,header_file);
		sequences_db_headers[i] = (char *) malloc((strlen(header)+1)*sizeof(char));
		strcpy(sequences_db_headers[i],header);
	}

	fclose(header_file);
	free(header);

	*ptr_sequences_db_headers = sequences_db_headers;
}

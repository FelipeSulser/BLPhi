#include "preprocess.h"



void convert_db (char * input_filename, char * out_filename) {

	uint64_t sequences_count=0;
	uint64_t linear_len=0;
	uint64_t disp;
	uint64_t accum;
	uint64_t chunk_size;
	uint64_t i;
	uint64_t j;
	uint64_t k;

	uint16_t *sequences_lengths=NULL;
	uint16_t* title_lengths=NULL;
	uint16_t length=0;
	uint16_t tmp_length;

	char ** sequences=NULL;
	char **titles=NULL;
	char buffer[BUFFER_SIZE];
	char filename[BUFFER_SIZE];
	char* bin_filename;
	char * res;
	char *tmp_seq;
	char *linear_struc __attribute__ ((aligned(64)))=NULL;
	char diff;

	int max_title_length;
	double epoch_ini= getTimestamp();

	FILE * sequences_file;
	FILE *titles_file;
	FILE *info_file;
	FILE * bin_file;

	if((sequences_file = fopen(input_filename,"r")) == NULL){
		fprintf(stderr,"Error while opening sequence file\n");
		exit(1);
	}

	// Allocate memory for sequences_lengths array
	sequences_lengths = (uint16_t *) malloc (BUFFER_SIZE*sizeof(uint16_t));
	title_lengths = (uint16_t *) malloc (BUFFER_SIZE*sizeof(uint16_t));

	// Calculate number of sequences in database and its lengths
	sequences_count=0;

	res = fgets(buffer,BUFFER_SIZE,sequences_file);
	while (res ) {
		length = 0;
		// read title
		while (strrchr(buffer,'\n') == NULL) {
			length += strlen(buffer);
			res = fgets(buffer,BUFFER_SIZE,sequences_file);
		}
		title_lengths[sequences_count] = length + strlen(buffer) + 1;
		// read sequence
		length = 0;
		res = fgets(buffer,BUFFER_SIZE,sequences_file);
		while ((res) && (buffer[0] != '>')) {
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

	// Allocate memory for sequences array
	if((sequences = (char **) malloc(sequences_count*sizeof(char *))) == NULL){
		fprintf(stderr,"Error while allocating seq memory\n");
		exit(1);
	}
	for (i=0; i<sequences_count; i++ ) {
		if((sequences[i] = (char *) malloc(sequences_lengths[i]*sizeof(char))) == NULL){
			fprintf(stderr,"Error while allocating memory for individual aminoacids");
			exit(1);
		}
	}


	rewind(sequences_file);


	i = 0;
	res = fgets(buffer,BUFFER_SIZE,sequences_file);
	while (res) {
		//  title
		while (strrchr(buffer,'\n') == NULL)
			res = fgets(buffer,BUFFER_SIZE,sequences_file);
		//  sequence
		length = 1;
		res = fgets(buffer,BUFFER_SIZE,sequences_file);
		while ((res) && (buffer[0] != '>')) {

			strncpy(sequences[i]+(length-1),buffer,strlen(buffer)-1);
			length += strlen(buffer)-1;
			res = fgets(buffer,BUFFER_SIZE,sequences_file);
		}
		i++;
	}

	// Rewind sequences
	rewind(sequences_file);

	// Allocate memory for titles array
	if((titles = (char **) malloc(sequences_count*sizeof(char *))) == NULL){
		fprintf(stderr,"Error allocating memory for titles array\n");
		exit(1);
	}
	for (i=0; i<sequences_count; i++ ) {
		if((titles[i] = (char *) malloc(title_lengths[i]*sizeof(char))) == NULL){
			fprintf(stderr,"Error allocating memory for single title\n");
			exit(1);
		}
	}

	// calculate max title length
	max_title_length = 0;
	for (i=0; i<sequences_count ; i++)
		max_title_length = (max_title_length > title_lengths[i] ? max_title_length : title_lengths[i]);


	free(title_lengths);

	// read sequence headers
	i = 0;
	res = fgets(buffer,BUFFER_SIZE,sequences_file);
	while (res ) {
		// discard
		while ((res ) && (buffer[0] != '>'))
			res = fgets(buffer,BUFFER_SIZE,sequences_file);
		if (res ){

			length = 1;
			do{
				strncpy(titles[i]+(length-1),buffer,strlen(buffer)-1);
				length += strlen(buffer)-1;
				res = fgets(buffer,BUFFER_SIZE,sequences_file);
			} while (strrchr(buffer,'\n') == NULL);
			titles[i][length] = '\0';
			i++;
		}
	}



	fclose(sequences_file);

	sort_sequences(sequences,titles,sequences_lengths, sequences_count);


	sprintf(filename,"%s.title",out_filename);
	if((titles_file = fopen(filename,"w")) == NULL){
		fprintf(stderr,"Error opening title file to write, maybe already existed and no privilege to write?\n");
		exit(1);
	}


	// write titles
	for (i=0; i<sequences_count ; i++)
		fprintf(titles_file,"%s\n",titles[i]);

	// close titles file


	for (i=0; i< sequences_count; i++ )
		linear_len = linear_len +  sequences_lengths[i];

	// transform bidimensional sequence array to a unidimensional one
	if((linear_struc = (char *) _mm_malloc(linear_len*sizeof(char),64)) == NULL){
		fprintf(stderr,"Error allocating memory for linear sequence datastruct\n");
		exit(1);
	}

	disp = 0;
	for (i=0; i< sequences_count; i++ ) {
		memcpy(linear_struc+disp,sequences[i],sequences_lengths[i]);
		disp += sequences_lengths[i];
	}

	// Free memory
	for (i=0; i< sequences_count; i++ )
		free(sequences[i]);
	free(sequences);
	// preprocess vect sequences DB
		//Original protein alphabet 24 chars --> treat JOU as the same
		#pragma omp parallel for private(diff) num_threads(PREPROCESS_THREADS) schedule(dynamic)
		for (i=0; i< linear_len; i++) {
			linear_struc[i] = ((linear_struc[i] == 'J') ? UNIFY : linear_struc[i]);
			linear_struc[i] = ((linear_struc[i] == 'O') ? UNIFY : linear_struc[i]);
			linear_struc[i] = ((linear_struc[i] == 'U') ? UNIFY : linear_struc[i]);
			diff = 'A';
			diff = (linear_struc[i] > 'J' ? diff+1 : diff);
			diff = (linear_struc[i] > 'O' ? diff+1 : diff);
			diff = (linear_struc[i] > 'U' ? diff+1 : diff);
			linear_struc[i] -= diff;
		}


	// Create info file: this file contains sequences count, number of residues and the maximum title length
	sprintf(filename,"%s.stats",out_filename);

	if((info_file = fopen(filename,"w")) == NULL){
		fprintf(stderr,"Error opening file to write stats\n");
		exit(1);
	}

	// Write info
	fprintf(info_file,"%ld %ld %d",sequences_count,linear_len,max_title_length);
	// close info file
		fclose(info_file);

	// Create sequences binary file: this file contains first the sequences lengths and then the preprocessed sequences residues
	sprintf(filename,"%s.bin",out_filename);

	if((bin_file = fopen(filename,"wb")) == NULL){
		fprintf(stderr,"Erorr opening binary file\n");
		exit(1);
	}
	// Write vectorized sequences lengths
	fwrite(sequences_lengths,sizeof(uint16_t),sequences_count,bin_file);

	//Write sequences
	fwrite(linear_struc,sizeof(char),linear_len,bin_file);


	fclose(titles_file);

	// Close bin file
	fclose(bin_file);

	// free memory
	free(sequences_lengths);
	_mm_free(linear_struc);

	printf("Database input:\t\t\t%s\n",input_filename);
	printf("Protein db size:\t\t\t%ld sequences and %lu aminoacids \n",sequences_count,linear_len);
	printf("Output files:\t\t\t%s\n",out_filename);
	printf("Elapsed time:\t\t\t%lf seconds\n\n",getTimestamp()-epoch_ini);

}

void reduce_encoding_db(){}

void reduce_encoding_query(char * query_sequences,
		uint16_t *query_sequences_lengths,
		uint64_t  query_sequences_count,
		uint64_t  Q,
		uint32_t * query_sequences_disp){

	uint32_t i,iter_seq;
	char a,b;
	char encoded;

	for(i = 0; i < query_sequences_count; i++){
		//were going to fit in a single char two old chars
		for(iter_seq = 0; iter_seq+1 < query_sequences_lengths[i]; iter_seq+=2){

		}
		if(iter_seq < query_sequences_lengths[i]){
			//no era multiplo de 2, insertar el ultimo con padding
		}
	}

}

void create_shifted_copy_nopad(char* seq,
		uint64_t D,
		uint64_t seq_count_db,
		uint32_t* shift,
		uint16_t* len,
		char** ptr_shifted_seq,
		int** ptr_extension,
		int** ptr_shift,
		uint64_t *total_linear
		){

	uint32_t j;
	int* extension = (int*) malloc(sizeof(int)*seq_count_db);
	int* new_shift = (int*) malloc(sizeof(int)*seq_count_db);

	//lets calculate total size to alloc memory
	int sub_len = 0;
	uint64_t added_total_len = 0;
	for(j = 0; j < seq_count_db; j++){
		sub_len = 0;
		if(len[j] >= 4){
			int needs_padding;
			int pad_len;
			uint32_t it = 0;
			uint32_t it_s;
			int pos_new = 0;

			int truelen = 0;

			for(it = 0; it < 4; it++){
				int size = len[j]-it;
				it_s = size;
				truelen += size;
				pos_new+=size;
				sub_len += size;
				//printf("it_s:%u pos_new:%u\n",it_s,pos_new);

				if( it_s % 4 != 0){
					int red_size = it_s %4;
					pos_new-= red_size;
					truelen -= red_size;
					sub_len -=red_size;
				}
			}

			needs_padding = (truelen % 64 == 0?0:1);
			pad_len = 64 - (truelen %64);


			if(needs_padding == 1){


				for(it = 0; it < pad_len; it++){
					sub_len++;
					//printf("added ad pos %u\n",pos_new);
					pos_new++;
				}
			}

		added_total_len += sub_len;
		extension[j] = sub_len;

		}else{
			extension[j] = 0;
		}
	}





	new_shift[0] = 0;
	for (j=1; j < seq_count_db;j++){
		new_shift[j] = new_shift[j-1] + extension[j-1];
		//printf("%u\n",new_shift[j]);
	}

	char* new_datastruc = (char*)_mm_malloc(sizeof(char)*added_total_len,64);
	for(j = 0; j < seq_count_db; j++){
		if(len[j] >= 4){
			//create the new data structure
			//first copy it four times-> best case


			int needs_padding = (len[j] % 16 == 0?0:1);

			int pad_len = 64 - (len[j] % 16)*4; //total pad size to add
			int total_len;
			if(needs_padding == 1)
				total_len = len[j]*4+pad_len;
			else
				total_len = len[j]*4;

			uint32_t it = 0;
			uint32_t it_s;
			int pos_new = 0;

			int truelen = 0;

			for(it = 0; it < 4; it++){
				int size = len[j]-it;
				memcpy(&new_datastruc[new_shift[j]+pos_new],&seq[shift[j]+it],size*sizeof(char));
				it_s = size;
				truelen += size;
				pos_new+=size;
				//printf("it_s:%u pos_new:%u\n",it_s,pos_new);

				if( it_s % 4 != 0){
					int red_size = it_s %4;
					pos_new-= red_size;
					truelen -= red_size;
				}
			}

			needs_padding = (truelen % 64 == 0?0:1);
			pad_len = 64 - (truelen%64);
			/*printf("truelen:%u\n",truelen);
			printf("pad_len:%u\n",pad_len);
			printf("pos_new:%u\n",pos_new);*/

			if(needs_padding == 1){
				//printf("Padding len:%u\n",pad_len);
				//printf("total_size %u\n",total_len);
				for(it = 0; it < pad_len; it++){
					new_datastruc[new_shift[j]+pos_new] = 127;
					pos_new++;
				}
			}


			}else{


			}
		}


	/*int it;
	for(j = 0; j < seq_count_db; j++){
		printf("\nORIGINAL:\n");
		for(it = 0; it < len[j]; it++){
			printf("%c\t",seq[shift[j]+it]);
		}

		printf("\n\nNEW\n");

		for(it = 0; it <extension[j]; it++){

			printf("%c\t",new_datastruc[it+new_shift[j]]);
		}
		printf("\n\n\n");
	}

	printf("LINEAR STRUCTURE is\n");
	for(it = 0; it < added_total_len; it++){

		printf("%c\t",new_datastruc[it]);


	}
	printf("\n");*/
	*ptr_shifted_seq = new_datastruc;
	*ptr_extension = extension;
	*ptr_shift = new_shift;
	*total_linear = added_total_len;
}

void create_shifted_copy(char* seq,
		uint64_t D,
		uint64_t seq_count_db,
		uint32_t* shift,
		uint16_t* len,
		char** ptr_shifted_seq,
		int** ptr_extension,
		int** ptr_shift,
		uint64_t *total_linear
		){

	uint32_t j;
	int* extension = (int*) malloc(sizeof(int)*seq_count_db);
	int* new_shift = (int*) malloc(sizeof(int)*seq_count_db);

	//lets calculate total size to alloc memory
	uint64_t added_total_len = 0;
	for(j = 0; j < seq_count_db; j++){
		if(len[j] >= 4){
		int needs_padding = (len[j] % 16 == 0?0:1);

		int pad_len = 64 - (len[j] % 16)*4; //total pad size to add
		int total_len;
		if(needs_padding == 1)
			total_len = len[j]*4+pad_len;
		else
			total_len = len[j]*4;

		added_total_len += total_len;
		extension[j] = total_len;
		}else{
			extension[j] = 0;

		}
	}
	new_shift[0] = 0;
	for (j=1; j < seq_count_db;j++){
		new_shift[j] = new_shift[j-1] + extension[j-1];
		//printf("%u\n",new_shift[j]);
	}


	char* new_datastruc = (char*)_mm_malloc(sizeof(char)*added_total_len,64);
	for(j = 0; j < seq_count_db; j++){
		if(len[j] >= 4){
			//create the new data structure
			//first copy it four times-> best case


			int needs_padding = (len[j] % 16 == 0?0:1);

			int pad_len = 64 - (len[j] % 16)*4; //total pad size to add
			int total_len;
			if(needs_padding == 1)
				total_len = len[j]*4+pad_len;
			else
				total_len = len[j]*4;

			uint32_t it = 0;
			uint32_t it_s;
			int pos_new = 0;
			for(it = 0; it < 4; it++){
				int size = len[j]-it;
				memcpy(&new_datastruc[new_shift[j]+pos_new],&seq[shift[j]+it],size*sizeof(char));
				it_s = size;
				pos_new+=size;
				//printf("it_s:%u pos_new:%u\n",it_s,pos_new);

				if( it_s % 4 != 0){
					//needs a padding
					int pad_size = 4 - it_s % 4;
					//printf("pad_size:%u\n",pad_size);
					for(it_s = 0; it_s < pad_size; it_s++){
						new_datastruc[new_shift[j]+pos_new] = 127;
						pos_new++;

					}
				}
			}

			if(needs_padding == 1){
				//printf("Padding len:%u\n",pad_len);
				//printf("total_size %u\n",total_len);
				for(it = 0; it < pad_len; it++){
					new_datastruc[new_shift[j]+pos_new] = 127;
					//printf("added ad pos %u\n",pos_new);
					pos_new++;
				}
			}


			}else{


			}
		}

/*
	int it;
	for(j = 0; j < seq_count_db; j++){
		printf("\nORIGINAL:\n");
		for(it = 0; it < len[j]; it++){
			printf("%c\t",seq[shift[j]+it]);
		}

		printf("\n\nNEW\n");

		for(it = 0; it <extension[j]; it++){

			printf("%c\t",new_datastruc[it+new_shift[j]]);
		}
		printf("\n\n\n");
	}

	printf("LINEAR STRUCTURE is\n");
	for(it = 0; it < added_total_len; it++){

		printf("%c\t",new_datastruc[it]);


	}
	printf("\n");
	*/
	*ptr_shifted_seq = new_datastruc;
	*ptr_extension = extension;
	*ptr_shift = new_shift;
	*total_linear = added_total_len;
}



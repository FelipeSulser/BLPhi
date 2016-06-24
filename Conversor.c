/*
 * Conversor.c
 *
 *  Created on: 23/02/2016
 *      Author: Felipe Sulser
 */

#include "Conversor.h"

// Our file type is FILE and operation to read is read


int readDatabase(char* file_name,db_t * dbt)
{

	db_t  db = NULL;
	char *seq;
	char *name;
	int  length;
	int l;
	uint32_t pos;
	int i = 0;
	fasta_file_t * fp = open_fasta(file_name);

	if(fp == NULL){
		perror("Open error\n");
	}
	 while (read_fasta(fp, &seq, &name, &length,&pos))
	    {
		 i++;
		 addNode(&db,name,seq,length,pos);
	     // printf(">%s\n", name);
	     // printf("%s\n",  seq);

	     // free(seq);
	     // free(name);
	    }
	close_fasta(fp);

	//traverse_list(db);
	*dbt = db;
	return i;
}


void change_encoding_db(
		char* db_seq,
		uint64_t sequences_count,
		uint32_t* seq_disp,
		uint64_t D,
		uint16_t* sequences_lengths,
		Assignment a,
		unsigned char** ptr_db_seq_red,
		uint64_t* ptr_D_red,
		uint32_t** ptr_db_disp_red,
		uint16_t** ptr_db_lengths_red
		){
		uint32_t i,j,k;
		//calculate new total len
		uint64_t new_size = 0;

		for(i = 0 ; i < sequences_count; i++){
				new_size += sequences_lengths[i]%2 == 0? (sequences_lengths[i]/2):(sequences_lengths[i]/2+1);
			}
		//new_size = D%2 == 0? D/2:D/2+1;
		unsigned char b1,b2;
		unsigned char red1,red2;
		unsigned char red_buffer = 0;
		uint8_t mask1 = 0xf0;
		uint8_t mask2 = 0x0f;
		int enc;

		unsigned char* db_seq_red = (unsigned char*) _mm_malloc(sizeof(unsigned char)*new_size,64);

		//create set with new alphabet
		char* my_set ;
		char* my_red;
		create_set(a,&my_red,&my_set);

		//create new query_disp
		uint16_t red_len;
		uint16_t reduced_len_curr;
		unsigned  int* db_disp_red = (unsigned  int*) _mm_malloc(sizeof(unsigned  int)*sequences_count+1,64);
		uint16_t* db_lengths_red = (uint16_t*)_mm_malloc(sizeof(unsigned short)*sequences_count+1,64);
		db_disp_red[0] =0;
		for(i =1; i < sequences_count; i++){
			red_len = (sequences_lengths[i-1] %2 == 0? sequences_lengths[i-1]/2:(sequences_lengths[i-1]/2)+1);
			db_disp_red[i] = db_disp_red[i-1]+ red_len;
		}
		for(i = 0; i < sequences_count; i++){
			reduced_len_curr = (sequences_lengths[i] %2 == 0? sequences_lengths[i]/2:(sequences_lengths[i]/2)+1);
			//printf("Red_len:%u\n",reduced_len_curr);
			db_lengths_red[i] = reduced_len_curr;
		}
		//query len is the same but per nibble! Caution
		//	least significant nibble is curr
		//	most significant nibble is next
		//  therefore [   i+1  ,  i    ]
		//			    4bit    4 bit

		uint32_t iter_lin = 0;
		int diff;
		for(i = 0; i < sequences_count; i++){

			//De dos en dos
			for(j = 0; j+1  < sequences_lengths[i]; j+=2){
				red_buffer = 0;
				b1 = db_seq[seq_disp[i]+j];

				b2 = db_seq[seq_disp[i]+j+1];
				//printf("%u %u\n",b1,b2);
				red1 = my_set[b1];
				red2 = my_set[b2];
				red_buffer = (red_buffer & mask1) | red1;
				red_buffer = (red_buffer & mask2) | (red2<<4);

				//printf("RED BUFFER:%x\n",red_buffer);
				db_seq_red[db_disp_red[i]+iter_lin] = red_buffer;
				iter_lin++;
			}
			if(j < sequences_lengths[i]){
				//Impar
				b1 = db_seq[seq_disp[i]+j];
				red1 = my_set[b1];
				red_buffer = (red_buffer & mask1)|red1;
				db_seq_red[db_disp_red[i]+iter_lin] = red_buffer;
				iter_lin++;

			}
			iter_lin = 0;
		}
	*ptr_db_seq_red = db_seq_red;
	*ptr_db_disp_red  =db_disp_red;
	*ptr_D_red = new_size;
	*ptr_db_lengths_red = db_lengths_red;
}


/**
 * In this method, we are going to change the encoding to 4 bits per base.
 * The encoding is established by the Assignment datastructure
 */
void change_encoding_query(
	char* query_seq,
	uint64_t query_count,
	uint32_t* query_disp,
	uint64_t Q,
	uint16_t* query_lengths,
	Assignment a,
	unsigned char** ptr_query_seq_red,
	uint64_t* ptr_Q_red,
	uint32_t ** ptr_query_disp_red,
	uint16_t** ptr_query_lengths_red
){
	uint32_t i,j,k;
	//calculate new total len
	uint64_t new_size = 0;
	//calculate new size
	for(i = 0 ; i < query_count; i++){
		new_size += query_lengths[i]%2 == 0? (query_lengths[i]/2):(query_lengths[i]/2+1);
	}
	//new_size = Q%2 == 0? Q/2:Q/2+1;
	unsigned char b1,b2;
	unsigned char red1,red2;
	unsigned char red_buffer = 0;
	uint8_t mask1 = 0xf0;
	uint8_t mask2 = 0x0f;
	int enc;
	unsigned char* query_seq_red = (unsigned char*) _mm_malloc(sizeof(unsigned char)*new_size,64);

	//create set with new alphabet
	char* my_set ;
	char* my_red;
	create_set(a,&my_red,&my_set);

	//create new query_disp
	uint16_t red_len;
	uint16_t reduced_len_curr;
	unsigned  int* query_disp_red = (unsigned  int*) _mm_malloc(sizeof(unsigned  int)*query_count+1,64);
	uint16_t* query_lengths_red = (uint16_t*)_mm_malloc(sizeof(unsigned short)*query_count+1,64);
	query_disp_red[0] =0;
	for(i =1; i < query_count; i++){
		red_len = (query_lengths[i-1] %2 == 0? query_lengths[i-1]/2:(query_lengths[i-1]/2)+1);
		query_disp_red[i] = query_disp_red[i-1]+ red_len;
	}
	for(i = 0; i < query_count; i++){
		reduced_len_curr = (query_lengths[i] %2 == 0? query_lengths[i]/2:(query_lengths[i]/2)+1);
		query_lengths_red[i] = reduced_len_curr;
	}
	//query len is the same but per nibble! Caution
	//	least significant nibble is curr
	//	most significant nibble is next
	//  therefore [   i+1  ,  i    ]
	//			    4bit    4 bit

	uint32_t iter_lin = 0;
	int diff;
	for(i = 0; i < query_count; i++){

		//De dos en dos
		for(j = 0; j+1  < query_lengths[i]; j+=2){
			red_buffer = 0;
			b1 = query_seq[query_disp[i]+j];

			b2 = query_seq[query_disp[i]+j+1];

			//printf("%u %u\n",b1,b2);
			red1 = my_set[b1];
			red2 = my_set[b2];
			red_buffer = (red_buffer & mask1) | red1;
			red_buffer = (red_buffer & mask2) | (red2<<4);
			//printf("RED_BUFFER:%x\n",red_buffer);
			query_seq_red[query_disp_red[i]+iter_lin] = red_buffer;
			iter_lin++;
		}
		if(j < query_lengths[i]){
			//Impar
			b1 = query_seq[query_disp[i]+j];

			red1 = my_set[b1];
			red_buffer = (red_buffer & mask1)|red1;

			query_seq_red[query_disp_red[i]+iter_lin] = red_buffer;
			iter_lin++;

		}
		iter_lin = 0;
	}
*ptr_query_seq_red = query_seq_red;
*ptr_query_disp_red  =query_disp_red;
*ptr_Q_red = new_size;
*ptr_query_lengths_red = query_lengths_red;


}



void destroy_db(node_t ** head){
	node_t * iter = *head;
	while(iter != NULL){
		free(iter->name);
		free(iter->seq);
		iter = iter->next;

	}
}

void destroy_reduced_db(reduced_node_t ** head){
	reduced_node_t * iter = *head;
	int i =0;
	while(iter != NULL){
		printf("At sequence %d\n",iter->location);
		printf("Length is %u\n",iter->length);
		for(i = 0; i < iter->length/2; i++){
					uint8_t twonib = iter->seq[i].my_byte;
					printf("%x,",twonib);
				}
				printf("\n");
		iter = iter->next;
	}
}

//Debug function to check all prot sequences added
void traverse_list(node_t* head){
	node_t * iter = head;
	while(iter != NULL){
		printf("name: %s\n", iter->name);
		printf("seq len %d\n",iter->seq_len);
		printf("at line %d\n",iter->pos_seq);
		printf("seq: %s\n", iter->seq);
		iter = iter->next;

	}
}



// As of now, format does not align sequences in blocks of 61 characters
// Can be improved
int storeDatabase(node_t ** head, char* file_name){


	FILE* f = fopen(file_name,"w");
	if(f == NULL) return -1;

	node_t * iter = *head;

	while(iter != NULL){
		if(fprintf(f,">%s\n",iter->name)  < 0) return -1;
		if(fprintf(f,"%s\n",iter->seq) < 0) return -1;
		iter = iter->next;

	}

	return 0;
}

//in binary format
int storeReducedDatabase(reduced_node_t **head, char * file_name){
	FILE* f = fopen(file_name,"w+");
	if(f == NULL) return -1;

	reduced_node_t* iter = *head;

	while(iter != NULL){
		fprintf(f,"%u\n",iter->location);
		fprintf(f,"%u\n",iter->length); // in nibbles
		int i = 0;
		for(i = 0; i < iter->length; i++){
			nibble_t* my_n = malloc(sizeof(nibble_t));
			my_n = &iter->seq[i];
			fprintf(f,"%c",my_n->my_byte);
		}
		fprintf(f,"\n");
		fwrite(iter->seq,1,iter->length*2,f);
		iter = iter->next;
	}

	fclose(f);
	return 0;
}

reduced_node_t* readReducedDatabase(char * file_name){
	FILE * f = fopen(file_name,"r");
	reduced_node_t * head;


	return head;
}

void addNode(node_t ** head, char* name, char* seq, int len, int pos) {
    node_t * new_node ;

    new_node = malloc(sizeof(node_t));
    new_node->name = name;
    new_node->seq_len = len;
    new_node->seq = seq;
    new_node->pos_seq = pos;

    new_node->next = *head;

    *head = new_node;
}


void addReducedNode(reduced_node_t ** head, uint32_t location, uint16_t len, nibble_t* seq){
	reduced_node_t * new_node;

	new_node = malloc(sizeof(reduced_node_t));
	new_node->length = len;
	new_node -> location = location;
	new_node ->seq = seq;
	new_node->next = *head;
	*head = new_node;



}

//Debug function to check all prot sequences added
void traverse_reduced_node(reduced_node_t* head){
	reduced_node_t * iter = head;
	while(iter != NULL){
		printf("loc: %u\n", iter->location);
		printf("seq len in nibbles %u\n",iter->length);
		int i =0;
		for(i = 0; i < iter->length/2; i++){
			uint8_t twonib = iter->seq[i].my_byte;
			printf("%x,",twonib);
		}
		printf("\n");
		iter = iter->next;

	}
}








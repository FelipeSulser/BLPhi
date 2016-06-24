/*
 * Conversor.h
 *
 *  Created on: 23/02/2016
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Set.h"
#include "fasta.h"


#ifndef CONVERSOR_H_
#define CONVERSOR_H_


typedef struct node {
    int16_t seq_len;
    uint32_t pos_seq;
    char* name;
    char* seq ;
    struct node * next;
    //starts counting with 0

} node_t;


typedef  node_t* db_t;


//encoded aminoacids in 4 bits each

typedef  struct nibble {
	uint8_t my_byte;
} nibble_t;

typedef  struct reduced_node {
	//pos starts at 0
	uint32_t location ;
	uint16_t length ; //nibble length--> if size is 4 means 2 bytes
	nibble_t * seq ;
	struct reduced_node* next;
}reduced_node_t;

typedef  reduced_node_t* reduced_db;

/**
 *	Changes the encoding of the current protein database to a 16 bit based encoding.
 *
 * Args:
 *           db		   - the database that has been opened previously with "readDatabase"
 *           a       - the assignment policy we want to follow (statistical or biological reasons)
 *           r_db      - RETURN: The encoded reduced database
 *           numEntries    - Establishes how many node the database has (not necessary but included because if we do so we can use a for loop that can be SIMD'd eventually)
 *
 *
 * Returns:
 *
 *           changeEncoding() "always succeeds" and returns void.
 *
 *
 *
 *	EXAMPLE:
 *
 *	For example we have the following change for the first entry of swissprot:
 *
 *		>sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) GN=FV3-001R PE=4 SV=1
		MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS
		EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD
		AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL
		EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD
		SFRKIYTDLGWKFTPL



		changed to (omit the commas please)

		loc: 0 (where it starts on the file)
		seq len in nibbles 256 (per 4 bits)
		(the actual sequence, does not have line breaks)
		8,ec,60,3,ba,46,13,11,81,6,aa,ea,4a,d4,32,b1,aa,43,6b,e0,dd
		,1,5,6,bd,d0,60,20,d2,ed,b6,a7,9,7,ec,97,4b,7b,b6,50,e0,6,3
		,b2,8,0,e0,cb,3b,80,11,54,97,f5,b0,d9,b7,a0,3e,3a,b0,b9,40,
		a2,f,76,60,c7,10,e4,1,bf,5,0,c0,ba,a6,1,be,45,26,20,9a,40,a9
		,3f,1b,30,59,a0,b6,a3,b0,c3,b,a,e6,0,81,51,7,98,2,4b,a9,54,aa,
		bb,70,70,3d,d7,93,fa,b0,7f,be,7,4a,33,ce,b1,49,3f,7a,b0,fc,ad,


		as we can see:
		 8 = 0000 1000 --> A = 0000, M = 1000 (the 8th new character established in our encoding policy)
		 ec = 1100 1110 --> S = 1100 F = 1110
		 ...
		 and so on


 *
 */
void changeEncoding(node_t * db, Assignment a, reduced_db* r_db, int numEntries);

//Used to load the database into memory
int readDatabase(char* file_name,db_t *dbt);


/**
 * Method that shows the content of protein database
 */
void traverse_list(node_t* head);


/**
 * Method used to add a node on the protein database with all the parameter that a protein sequence has
 *
 */
void addNode(node_t** head, char* name, char* seq, int len, int pos);
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
);

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
		);
/**
 * Destroy the database, free the memory
 *
 */
void destroy_db(node_t ** head);


void destroy_reduced_db(reduced_node_t ** head);



void addReducedNode(reduced_node_t ** head, uint32_t location, uint16_t len, nibble_t*seq);


/**
 * Method that shows the content of the reduced encoded database
 */
void traverse_reduced_node(reduced_node_t* head);


int storeReducedDatabase(reduced_node_t **head, char * file_name);

#endif /* CONVERSOR_H_ */

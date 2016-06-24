/*
 * Set.c
 *
 *  Created on: 25/02/2016
 *      Author: galvez
 */

#include "Set.h"


void create_set(Assignment a, char** ptr_set,char** ptr_enc_red){
	//given the assignment, get all chars of new encoding
	int index = 0;
	int diff;
	int i = 0;
	//find biggest
	char biggest = 0;
	char *enc_red;
	char big_red = 0;
	//waste of memory, but assured not to be greater than the whole ascii table so it doesnt really matter
	//Not a bijective function nor injective
	for(i = 0; i < 23; i++){
		biggest = a->orig[i]>biggest?a->orig[i]:biggest;
		big_red = a->new[i]>big_red?a->new[i]:big_red;

	}

	enc_red = (char*) malloc(sizeof(char)*big_red);
	for( i = 0; i < 16; i++){
		diff = 'A';
		diff = (a->reduced_set[i] >'J'?diff+1:diff);
		diff = (a->reduced_set[i] > 'O'?diff+1:diff);
		diff = (a->reduced_set[i] > 'U'? diff+1:diff);
		enc_red[a->reduced_set[i] - diff] = i;

		//enc_red[a->reduced_set[i]] = i;
	}

	char* set = (char*) malloc(sizeof(char)*biggest);
	int ctr = 0;
	for(i = 0; i < 23; i++){
		diff = 'A';
		diff = (a->reduced_set[i] >'J'?diff+1:diff);
		diff = (a->reduced_set[i] > 'O'?diff+1:diff);
		diff = (a->reduced_set[i] > 'U'? diff+1:diff);
		enc_red[a->reduced_set[i] - diff] = i;
		set[a->orig[i]] = enc_red[a->new[i] - diff];
	}

	*ptr_set = set;
	*ptr_enc_red = enc_red;

}

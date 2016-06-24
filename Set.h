/*
 * Set.h
 *
 *  Created on: 25/02/2016
 *      Author: galvez
 */

#ifndef SET_H_
#define SET_H_



#include <stdio.h>
#include <stdlib.h>
#include "Assignment.h"

//Obsolete
typedef struct set {
	char letters[23];
}red_set_t;


void create_set(Assignment a, char** ptr_set,char** ptr_enc_red);
#endif /* SET_H_ */

/*
 * Assignment.h
 *
 *  Created on: 01/03/2016
 *      Author: galvez
 */

#ifndef ASSIGNMENT_H_
#define ASSIGNMENT_H_

#define PROTEIN_STD_ENCODING 23
//the assignment will be 1 to many
typedef __attribute__ ((aligned(64))) struct{
    char orig[PROTEIN_STD_ENCODING];
    char new[PROTEIN_STD_ENCODING];

    char reduced_set[16];


}assignment_tuple; // 2 seqs to be compared

typedef __attribute__ ((aligned(64))) assignment_tuple* Assignment;


#endif /* ASSIGNMENT_H_ */

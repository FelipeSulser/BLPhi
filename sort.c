#include "sort.h"


void merge_score_simple(int*score,int*index,unsigned long int size){
	unsigned long int i1 = 0;
	unsigned long int i2 = size / 2;
	unsigned long int it = 0;

	int* tmp1 =(int*) malloc(sizeof(int)*size);
	int*tmp2 = (int*) malloc(sizeof(int)*size);

	while(i1 < size/2 && i2 < size){
		if(score[i1] > score[i2]){
			tmp1[it] = score[i1];
			tmp2[it] = index[i1];
			i1++;
		}
		else{
			tmp1[it] = score[i2];
			tmp2[it] = index[i2];
			i2++;
		}
		it++;
	}
}
void merge_scores(int * scores, char ** titles, unsigned long int size) {
	unsigned long int i1 = 0;
	unsigned long int i2 = size / 2;
	unsigned long int it = 0;
	// allocate memory for temporary buffers
	char ** tmp2 = (char **) malloc(size*sizeof(char *));
	int * tmp3 = (int *) malloc (size*sizeof(int));

	while(i1 < size/2 && i2 < size) {
		if (scores[i1] > scores[i2]) {
			tmp2[it] = titles[i1];
			tmp3[it] = scores[i1];
			i1++;
		}
		else {
			tmp2[it] = titles[i2];
			tmp3[it] = scores[i2];
			i2 ++;
		}
		it ++;
	}

	while (i1 < size/2) {
		tmp2[it] = titles[i1];
		tmp3[it] = scores[i1];
	    i1++;
	    it++;
	}
	while (i2 < size) {
		tmp2[it] = titles[i2];
		tmp3[it] = scores[i2];
	    i2++;
	    it++;
	}

	memcpy(titles, tmp2, size*sizeof(char *));
	memcpy(scores, tmp3, size*sizeof(int));

	free(tmp2);
	free(tmp3);

}

void sort_score_simple(int* score, int* index, unsigned long int size){
	int tmp_score;
	int tmp_index;
	if(size == 2){
		if(score[0] <= score[1]){
			tmp_score = score[0];
			tmp_index = index[0];

			score[0] = score[1];
			score[1] = tmp_score;

			index[0] = index[1];
			index[1] = tmp_index;
			return;
		}
	}else{
		if(size > 2){
			sort_score_simple(score,index,size/2);
			sort_score_simple(score + size/2,index+size/2,size - size/2);
			merge_score_simple(score,index,size);
		}
	}
}
void sort_scores(int * scores, char ** titles, unsigned long int size) {
	int tmp_score;
	char * tmp_seq;

	if (size == 2) { 
		if (scores[0] <= scores[1]) {
			// swap scores
			tmp_score = scores[0];
			scores[0] = scores[1];
			scores[1] = tmp_score;
			// swap titles
			tmp_seq = titles[0];
			titles[0] = titles[1];
			titles[1] = tmp_seq;
		}
	} else {
		if (size > 2){
			sort_scores(scores, titles, size/2);
			sort_scores(scores + size/2, titles + size/2, size - size/2);
			merge_scores(scores, titles, size);
		}
	}
}

// Wall time
double getTimestamp()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}




void merge_sequences(char ** sequences, char ** titles, unsigned short int * sequences_lengths, unsigned long int size) {
	unsigned long int i1 = 0;
	unsigned long int i2 = size / 2;
	unsigned long int it = 0;
	// allocate memory for temporary buffers
	char ** tmp1 = (char **) malloc(size*sizeof(char *));
	char ** tmp2 = (char **) malloc(size*sizeof(char *));
	unsigned short int * tmp3 = (unsigned short int *) malloc (size*sizeof(unsigned short int));

	while(i1 < size/2 && i2 < size) {
		if (sequences_lengths[i1] <= sequences_lengths[i2]) {
			tmp1[it] = sequences[i1];
			tmp2[it] = titles[i1];
			tmp3[it] = sequences_lengths[i1];
			i1++;
		}
		else {
			tmp1[it] = sequences[i2];
			tmp2[it] = titles[i2];
			tmp3[it] = sequences_lengths[i2];
			i2 ++;
		}
		it ++;
	}

	while (i1 < size/2) {
		tmp1[it] = sequences[i1];
		tmp2[it] = titles[i1];
		tmp3[it] = sequences_lengths[i1];
	    i1++;
	    it++;
	}
	while (i2 < size) {
		tmp1[it] = sequences[i2];
		tmp2[it] = titles[i2];
		tmp3[it] = sequences_lengths[i2];
	    i2++;
	    it++;
	}

	memcpy(sequences, tmp1, size*sizeof(char *));
	memcpy(titles, tmp2, size*sizeof(char *));
	memcpy(sequences_lengths, tmp3, size*sizeof(unsigned short int));

	free(tmp1);
	free(tmp2);
	free(tmp3);

}


void sort_sequences (char ** sequences, char ** titles, unsigned short int * sequences_lengths, unsigned long int size) {
	char * tmp_seq;
	unsigned short int tmp_seq_len;

	if (size == 2) {
		if (sequences_lengths[0] > sequences_lengths[1]) {
			// swap sequences
			tmp_seq = sequences[0];
			sequences[0] = sequences[1];
			sequences[1] = tmp_seq;
			// swap titles
			tmp_seq = titles[0];
			titles[0] = titles[1];
			titles[1] = tmp_seq;
			// swap sequences lengths
			tmp_seq_len = sequences_lengths[0];
			sequences_lengths[0] = sequences_lengths[1];
			sequences_lengths[1] = tmp_seq_len;
			return;
		}
	} else {
		if (size > 2){
			sort_sequences(sequences, titles, sequences_lengths, size/2);
			sort_sequences(sequences + size/2, titles + size/2, sequences_lengths + size/2, size - size/2);
			merge_sequences(sequences, titles, sequences_lengths, size);
		}
	}
}

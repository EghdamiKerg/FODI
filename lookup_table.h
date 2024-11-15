#ifndef LOOKUP_TABLE_H
#define LOOKUP_TABLE_H

#include "split_database.h"

typedef struct
{
	int*  kmer_counts;
	int** kmer_starts;
	int*  kmer_offsets;
} ref_index;

typedef struct supporting_node
{
    int loc1;
    //int loc2;
    supporting_node* next_node;
};

ref_index*
destroy_ref_index(ref_index* ridx);

ref_index*
create_ref_index(volume_t* v, int kmer_size, const int num_threads);

ref_index*
create_second_index(volume_t* v, int kmer_size, int second_kmer_size, const int num_threads, ref_index* first_index);
#endif // LOOKUP_TABLE_H

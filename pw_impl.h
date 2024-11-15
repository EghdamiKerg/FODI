#ifndef PW_IMPL_H
#define PW_IMPL_H

#include <iostream>
#include<math.h>

#include "alignment.h"
#include "packed_db.h"
#include "lookup_table.h"

#define RM 			        100000
#define DN 			        500
#define BC 			        10
#define SM 			        3
#define SI 			        4
#define MAX_GAP             1000
#define NUM_RANDOMIZED_KMERS 5
#define CHUNK_SIZE 	        500
#define BLOCK_SIZE 	        256
#define initial_arr_size    50000
#define READ_KMER_SIZE      8
#define POSITIONS_PER_KMER  5
#define MUL_BLOCK_SIZE(a) 	((a)*BLOCK_SIZE)
#define DIV_BLOCK_SIZE(a) 	((a)/BLOCK_SIZE)
#define MOD_BLOCK_SIZE(a) 	((a)%BLOCK_SIZE)

struct candidate_save
{
    int loc1,loc2,left1,left2,right1,right2,score,num1,num2,readno,readstart;
    char chain;
};

typedef candidate_save Candidate;


struct PWThreadData
{
	options_t*				options;
	int 					used_thread_id;
	pthread_mutex_t 		id_lock;
	volume_t* 				reference;
	volume_t* 				reads;
	ref_index* 				first_index;
	ref_index*              second_index;
	bool                    self_detection;
 	std::ostream*			out;
	M4Record** 				m4_results;
	ExtensionCandidate**	ec_results;
	static const int		kResultListSize = 10000;
	pthread_mutex_t			result_write_lock;
	int						next_processed_id;
	pthread_mutex_t			read_retrieve_lock;

	PWThreadData(options_t* opt, volume_t* ref, volume_t* rd, ref_index* first_index, ref_index* second_index, std::ostream* o);
	~PWThreadData();
};

struct temp_arr_element
{
    int score;
    int index;
};

/*struct marked_blocks
{
    int row_index;
    int column_index;
    short block_score;
    int chaining_score;
    int prev;
    bool visited;
};*/

struct Short_seed
{
    int seed_loc_on_ref;
    int seed_loc_on_read;//not exactly location
    int chaining_score;
    int prev;
    bool visited;

};

struct Long_seed{
    int seed_loc_on_ref;
    int seed_loc_on_read;
};

struct candidate_region
{
    int start_loc_on_ref;
    int end_loc_on_ref;
    //int start_intercept;
    //int end_intercept;
    int num_seeds;
	Short_seed seeds[10000];
};

struct SeedingBK
{
	//marked_blocks* marked_blocks_list;
	//int** database;
	int num_candidate_regions;
	candidate_region* candidate_regions;
	uint32_t* kmer_ids;
	//SeedingBK();];
    int num_long_seeds;
    Long_seed* long_seeds;
    int size_of_long_seeds_array=100000;
    SeedingBK();
	~SeedingBK();
};

void
process_one_volume(options_t* options, const int svid, const int evid, volume_names_t* vn, std::ostream* out);

#endif // PW_IMPL_H

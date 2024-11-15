#include "split_database.h"
#include "pw_options.h"
#include "diff_gapalign.h"
#include "xdrop_gapalign.h"
#include "packed_db.h"
#include "lookup_table.h"
#include "second_lookup_table.h"
#include "pw_impl.h"

#include <algorithm>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <tgmath.h>

#define MSS MAX_SEQ_SIZE

static int MAXC = 100;
static int output_gapped_start_point = 1;
static int kmer_size = 13;
static int second_kmer_size = 20;
static const double pacbio_error_rate = 0.10;
static const double nanopore_error_rate = 0.15;
static double error_rate=pacbio_error_rate;
int index_separator;
//static const double ddfs_cutoff_pacbio = 0.25;
//static const double ddfs_cutoff_nanopore = 0.25;
//static double ddfs_cutoff = ddfs_cutoff_pacbio;
static int min_align_size = 0;
//static int min_kmer_match = 0;
static int min_kmer_dist = 0;

using namespace std;

PWThreadData::PWThreadData(options_t* opt, volume_t* ref, volume_t* rd, ref_index* fidx, ref_index* sidx, std::ostream* o)
	: options(opt), used_thread_id(0), reference(ref), reads(rd), first_index(fidx), second_index(sidx), out(o), m4_results(NULL), ec_results(NULL), next_processed_id(0), self_detection(false)
{
	pthread_mutex_init(&id_lock, NULL);
	if (options->task == TASK_SEED)
	{
		safe_malloc(ec_results, ExtensionCandidate*, options->num_threads);
		for (int i = 0; i < options->num_threads; ++i)
			safe_malloc(ec_results[i], ExtensionCandidate, kResultListSize);
	}
	else if (options->task == TASK_ALN)
	{
		safe_malloc(m4_results, M4Record*, options->num_threads);
		for (int i = 0; i < options->num_threads; ++i)
			safe_malloc(m4_results[i], M4Record, kResultListSize);
	}
	else
	{
		LOG(stderr, "Task must be either %d or %d, not %d!", TASK_SEED, TASK_ALN, options->task);
		abort();
	}
	pthread_mutex_init(&result_write_lock, NULL);
	pthread_mutex_init(&read_retrieve_lock, NULL);
}

PWThreadData::~PWThreadData()
{
	if (ec_results)
	{
		for (int i = 0; i < options->num_threads; ++i) safe_free(ec_results[i]);
		safe_free(ec_results);
	}
	if (m4_results)
	{
		for (int i = 0; i < options->num_threads; ++i) safe_free(m4_results[i]);
		safe_free(m4_results);
	}
}

void
reverse_complement(char* dst, const char* src, const int size)
{
	const uint8_t* rc_table = get_dna_complement_table();
	int i;
	for (i = 0; i < size; ++i)
	{
		uint8_t c = src[i];
		assert(c >= 0 && c < 4);
		c = rc_table[c];
		dst[size - 1 - i] = c;
	}
}

int
extract_kmers(const char* s, const int sstart, const int ssize, uint32_t* kmer_ids)
{
	uint32_t max_id = 1 << (2 * kmer_size);
	//int num_kmers = ssize / kmer_size;
	int read_end = ssize - 1;
	int i=0, j;
    uint32_t eit = 0;
    //uint32_t leftnum = 34 - 2 * kmer_size;
    int p=kmer_size-1;
    for (j = 0; j < ssize; ++j)
    {
        eit = (eit << 2) | s[j];
			//assert(eit < index_count);

        if (j==p)
        {
            assert(eit < max_id);
            kmer_ids[i] = eit;
            eit= 0;
            ++i;
            p=(i+1)*kmer_size-1;
            //p+=(kmer_size/2);
        }
        //if (j>=kmer_size-1){
            //eit <<= leftnum;
            //eit >>= leftnum;
        //}
    }
	return i-1;
}

int
extract_kmers2(const char* s, const int ssize, uint32_t* kmer_ids)
{
	uint32_t max_id = 1 << (2 * kmer_size);
	//int num_kmers = ssize / kmer_size;
	//int read_end = ssize - 1;
	int i=0, j;
    uint32_t eit = 0;
    uint32_t leftnum = 34 - 2 * kmer_size;
    //int p=kmer_size-1;
    for (j = 0; j < ssize; ++j)
    {
        eit = (eit << 2) | s[j];
			//assert(eit < index_count);

        if (j>=kmer_size-1)
        {
            assert(eit < max_id);
            kmer_ids[i] = eit;
            //eit= 0;
            ++i;
            //p=(i+1)*kmer_size-1;
            //p+=(kmer_size/2);
        //}
        //if (j>=kmer_size-1){
            eit <<= leftnum;
            eit >>= leftnum;
        }
    }
	return i;
}

void
extract_randomized_kmers(volume_t* v, const int start_pos, int* kmer_ids, uint32_t index_count)
{
	//int max_id = 1 << (2 * kmer_size);
	//int num_kmers = ssize / kmer_size;
	int i,j;
	bool bad_pos=false;
	//int end_pos = ssize - 1;
    uint32_t eit = 0;
    srand((unsigned) time(NULL));
    int random_pos;
    for (j = 0; j < NUM_RANDOMIZED_KMERS; ++j)
    {
        random_pos = start_pos + (rand() % (BLOCK_SIZE-READ_KMER_SIZE));
        for (i=random_pos; i<random_pos+READ_KMER_SIZE; ++i){
            uint8_t c = PackedDB::get_char(v->data, i);
            if (c>3 || c<0){
                bad_pos=true;
                break;
            }
			eit = (eit << 2) | c;
			assert(eit < index_count);//***************************
        }
        if(bad_pos){
            j--;
            eit=0;
            bad_pos=false;
            continue;
        }
        kmer_ids[j] = eit;
        eit= 0;
    }
}

SeedingBK::SeedingBK()
{
	num_candidate_regions=0;
	num_long_seeds=0;
	size_of_long_seeds_array=initial_arr_size;
	safe_malloc(long_seeds, Long_seed, initial_arr_size);
	safe_malloc(candidate_regions, candidate_region,20000);
	safe_malloc(kmer_ids, uint32_t, MAX_SEQ_SIZE);///kmer_size+kmer_size);
}

SeedingBK::~SeedingBK()
{
	safe_free(kmer_ids);
	safe_free(long_seeds);
	safe_free(candidate_regions);
}

bool compare_long_seeds(Long_seed i, Long_seed j){
    if(i.seed_loc_on_ref==j.seed_loc_on_ref)
        return i.seed_loc_on_read<j.seed_loc_on_read;
    else
        return i.seed_loc_on_ref<j.seed_loc_on_ref;

}

bool compare_short_seeds (Short_seed i, Short_seed j){
    if(i.seed_loc_on_ref==j.seed_loc_on_ref)
        return i.seed_loc_on_read<j.seed_loc_on_read;
    else
        return i.seed_loc_on_ref<j.seed_loc_on_ref;
}


Long_seed* resize_array(Long_seed* arr, int& size){
    Long_seed* new_arr;
    safe_malloc(new_arr, Long_seed, size*2);
    for(int i=0; i<size; i++){
        new_arr[i]=arr[i];
    }
    size*=2;
    safe_free(arr);
    return new_arr;
}

void
seeding(const char* read, const int read_start, const int read_size, const int chain, volume_t* v, ref_index* first_index, ref_index* second_index, SeedingBK* sbk)
{
	uint32_t* kmer_ids = sbk->kmer_ids;
	sbk->num_candidate_regions=0;
    candidate_region* candidate_regions = sbk->candidate_regions;
    candidate_region* candidate_regions_spr;
	sbk->num_long_seeds=0;
	Long_seed* long_seeds = sbk->long_seeds;
	Long_seed* long_seeds_spr=long_seeds;
	int num_kmers = extract_kmers2(read, read_size, kmer_ids);
	int km;
	int loc_on_ref;
	int loc_on_read=read_start;

    int ind;
	for (km = 0; km < num_kmers; ++km){
        int num_seeds=second_index->kmer_counts[kmer_ids[km]];
        int* seed_arr = second_index->kmer_starts[kmer_ids[km]];
        ind=0;
        //if (chain==REV){
            for (int j=0;j<num_seeds;j++){
                if (seed_arr[j]==index_separator){
                    continue;
                }
                int t1=km+kmer_size;
                int t2=seed_arr[j]+kmer_size;
                int sim=0;
                for (int k=0;k<second_kmer_size-kmer_size;k++){
                    char ch=PackedDB::get_char(v->data, t2+k);
                    if (read[t1+k]==ch)
                        sim++;
                    else
                        break;
                }
                if (sim==second_kmer_size-kmer_size){
                    for (int m=j; seed_arr[m]!=index_separator; m++,j++){
                        //if(seed_arr[m]>=read_start) continue;
                        long_seeds_spr->seed_loc_on_read=loc_on_read-read_start;
                        long_seeds_spr->seed_loc_on_ref=seed_arr[m];
                        ++sbk->num_long_seeds;
                        if (sbk->num_long_seeds >= sbk->size_of_long_seeds_array){
                            sbk->long_seeds = resize_array(sbk->long_seeds,sbk->size_of_long_seeds_array);
                            long_seeds = sbk->long_seeds;
                            long_seeds_spr = long_seeds + sbk->num_long_seeds;
                        }
                        else ++long_seeds_spr;
                    }
                    break;
                }
                else{
                    while(j<num_seeds && seed_arr[j]!=index_separator){
                        j++;
                    }
                }

            }
        //}
        //loc_on_read+=kmer_size;
		loc_on_read++;
    }


	if (sbk->num_long_seeds>600000)
        std::cout<<"ooops"<<std::endl;
	sort(long_seeds, long_seeds+sbk->num_long_seeds,compare_long_seeds);
	int it=0, start, end;
	double error_re;
	candidate_regions_spr=candidate_regions;
	long_seeds_spr=long_seeds;
	while(it<sbk->num_long_seeds){
        sbk->num_candidate_regions++;
        if (sbk->num_candidate_regions>20000){
            std::cout<<"oooops"<<std::endl;
        }
        start = long_seeds_spr->seed_loc_on_ref - long_seeds_spr->seed_loc_on_read;
        error_re = long_seeds_spr->seed_loc_on_read*error_rate;
        start-=error_re;
        candidate_regions_spr->start_loc_on_ref=start;
        end=long_seeds_spr->seed_loc_on_ref + (read_size-long_seeds_spr->seed_loc_on_read);
        error_re = (read_size-long_seeds_spr->seed_loc_on_read)*error_rate;
        end+=error_re;
        candidate_regions_spr->end_loc_on_ref=end;
        candidate_regions_spr->num_seeds=0;
        long_seeds_spr++;
        it++;
        while (it<sbk->num_long_seeds){
            if (long_seeds_spr->seed_loc_on_ref<=candidate_regions_spr->end_loc_on_ref){//(long_seeds_spr->seed_loc_on_ref>=start && long_seeds_spr->seed_loc_on_ref<=end){
                start = long_seeds_spr->seed_loc_on_ref - long_seeds_spr->seed_loc_on_read;
                error_re = long_seeds_spr->seed_loc_on_read*error_rate;
                start-=error_re;
                candidate_regions_spr->start_loc_on_ref= min(start,candidate_regions_spr->start_loc_on_ref);
                end=long_seeds_spr->seed_loc_on_ref + (read_size-long_seeds_spr->seed_loc_on_read);
                error_re = (read_size-long_seeds_spr->seed_loc_on_read)*error_rate;
                end+=error_re;
                candidate_regions_spr->end_loc_on_ref=max(end,candidate_regions_spr->end_loc_on_ref);
            }
            else{
                break;
            }
            long_seeds_spr++;
            it++;
        }
        candidate_regions_spr++;
	}
    loc_on_read=0;
    int L, R, mid;
	for (km = 0; km < num_kmers; km+=kmer_size)
	{
		int num_seeds = first_index->kmer_counts[kmer_ids[km]];
		int* seed_arr = first_index->kmer_starts[kmer_ids[km]];

		//int* seed_arr= NULL;
		//int num_seeds = find_relevant_kmers(kmer_ids[km] ,kmer_gc_percents[km], ridx, &seed_arr);
            for (int sid = 0; sid < num_seeds; ++sid)
            {
                L=0;
                R = sbk->num_candidate_regions-1;
                while(L<=R){
                    mid = (L+R)/2;
                    candidate_regions_spr=candidate_regions+mid;
                    if (seed_arr[sid]>=candidate_regions_spr->start_loc_on_ref && seed_arr[sid]<=candidate_regions_spr->end_loc_on_ref){
                        int nn=candidate_regions_spr->num_seeds;
                        if (nn<10000){
                            candidate_regions_spr->seeds[nn].seed_loc_on_ref=seed_arr[sid];
                            candidate_regions_spr->seeds[nn].seed_loc_on_read=loc_on_read;
                            candidate_regions_spr->num_seeds++;
                        }else{
                            std::cout<<"oh"<<std::endl;
                        }
                        break;
                    }
                    else if (seed_arr[sid]<candidate_regions_spr->start_loc_on_ref){
                        R=mid-1;
                    }
                    else{
                        L=mid+1;
                    }
                }
            }
        loc_on_read+=kmer_size;
	}
}

void
self_seeding(const char* read, const int read_start, const int read_size, const int chain, volume_t* v, ref_index* first_index, ref_index* second_index, SeedingBK* sbk)
{
	uint32_t* kmer_ids = sbk->kmer_ids;
	sbk->num_candidate_regions=0;
    candidate_region* candidate_regions = sbk->candidate_regions;
    candidate_region* candidate_regions_spr;
	sbk->num_long_seeds=0;
	Long_seed* long_seeds = sbk->long_seeds;
	Long_seed* long_seeds_spr=long_seeds;
	int num_kmers = extract_kmers2(read, read_size, kmer_ids);
	int km;
	int loc_on_ref;
	int loc_on_read=read_start;
    int ind;
	for (km = 0; km < num_kmers; ++km){
        int num_seeds=second_index->kmer_counts[kmer_ids[km]];
        int* seed_arr = second_index->kmer_starts[kmer_ids[km]];
        ind=0;
        if (chain==REV){
            for (int j=0;j<num_seeds;j++){
                if (seed_arr[j]==index_separator){
                    continue;
                }
                int t1=km+kmer_size;
                int t2=seed_arr[j]+kmer_size;
                int sim=0;
                for (int k=0;k<second_kmer_size-kmer_size;k++){
                    char ch=PackedDB::get_char(v->data, t2+k);
                    if (read[t1+k]==ch)
                        sim++;
                    else
                        break;
                }
                if (sim==second_kmer_size-kmer_size){
                    for (int m=j; seed_arr[m]!=index_separator; m++,j++){
                        if(seed_arr[m]>=read_start) continue;
                        long_seeds_spr->seed_loc_on_read=loc_on_read-read_start;
                        long_seeds_spr->seed_loc_on_ref=seed_arr[m];
                        ++sbk->num_long_seeds;
                        if (sbk->num_long_seeds >= sbk->size_of_long_seeds_array){
                            sbk->long_seeds = resize_array(sbk->long_seeds,sbk->size_of_long_seeds_array);
                            long_seeds = sbk->long_seeds;
                            long_seeds_spr = long_seeds + sbk->num_long_seeds;
                        }
                        else ++long_seeds_spr;
                    }
                    break;;
                }
                else{
                    while(j<num_seeds && seed_arr[j]!=index_separator){
                        j++;
                    }
                }

            }
        }
        else {
            for (int j=0;j<num_seeds;j++){
                if (seed_arr[j]==index_separator){
                    ind=j+1;
                    continue;
                }
                if (seed_arr[j]==loc_on_read){
                    for (int m=ind; seed_arr[m]!=index_separator; m++){
                        if(seed_arr[m]>=read_start) continue;
                        long_seeds_spr->seed_loc_on_read=loc_on_read-read_start;
                        long_seeds_spr->seed_loc_on_ref=seed_arr[m];
                        ++sbk->num_long_seeds;
                        if (sbk->num_long_seeds >= sbk->size_of_long_seeds_array){
                            sbk->long_seeds = resize_array(sbk->long_seeds,sbk->size_of_long_seeds_array);
                            long_seeds = sbk->long_seeds;
                            long_seeds_spr = long_seeds + sbk->num_long_seeds;
                        }
                        else ++long_seeds_spr;
                    }
                }
            }
        }
        //loc_on_read+=kmer_size;
		loc_on_read++;
    }


	if (sbk->num_long_seeds>50000)
        std::cout<<"ooops"<<std::endl;
	sort(long_seeds, long_seeds+sbk->num_long_seeds,compare_long_seeds);
	int it=0, start, end;
	double error_re;
	candidate_regions_spr=candidate_regions;
	long_seeds_spr=long_seeds;
	while(it<sbk->num_long_seeds){
        sbk->num_candidate_regions++;
        if (sbk->num_candidate_regions>20000){
            std::cout<<"oooops"<<std::endl;
        }
        start = long_seeds_spr->seed_loc_on_ref - long_seeds_spr->seed_loc_on_read;
        error_re = long_seeds_spr->seed_loc_on_read*error_rate;
        start-=error_re;
        candidate_regions_spr->start_loc_on_ref=start;
        end=long_seeds_spr->seed_loc_on_ref + (read_size-long_seeds_spr->seed_loc_on_read);
        error_re = (read_size-long_seeds_spr->seed_loc_on_read)*error_rate;
        end+=error_re;
        candidate_regions_spr->end_loc_on_ref=end;
        candidate_regions_spr->num_seeds=0;
        long_seeds_spr++;
        it++;
        while (it<sbk->num_long_seeds){
            if (long_seeds_spr->seed_loc_on_ref<=candidate_regions_spr->end_loc_on_ref){//(long_seeds_spr->seed_loc_on_ref>=start && long_seeds_spr->seed_loc_on_ref<=end){
                start = long_seeds_spr->seed_loc_on_ref - long_seeds_spr->seed_loc_on_read;
                error_re = long_seeds_spr->seed_loc_on_read*error_rate;
                start-=error_re;
                candidate_regions_spr->start_loc_on_ref= min(start,candidate_regions_spr->start_loc_on_ref);
                end=long_seeds_spr->seed_loc_on_ref + (read_size-long_seeds_spr->seed_loc_on_read);
                error_re = (read_size-long_seeds_spr->seed_loc_on_read)*error_rate;
                end+=error_re;
                candidate_regions_spr->end_loc_on_ref=max(end,candidate_regions_spr->end_loc_on_ref);
            }
            else{
                break;
            }
            long_seeds_spr++;
            it++;
        }
        candidate_regions_spr++;
	}
    loc_on_read=0;
    int L, R, mid;
	for (km = 0; km < num_kmers; km+=kmer_size)
	{
		int num_seeds = first_index->kmer_counts[kmer_ids[km]];
		int* seed_arr = first_index->kmer_starts[kmer_ids[km]];

		//int* seed_arr= NULL;
		//int num_seeds = find_relevant_kmers(kmer_ids[km] ,kmer_gc_percents[km], ridx, &seed_arr);
            for (int sid = 0; sid < num_seeds; ++sid)
            {
		    //************************May change*****************************
                if(seed_arr[sid] >= read_start){
                    continue;
                }
            //**************************************************************
                /*candidate_regions_spr=candidate_regions;
                for (int reg_id=0; reg_id<sbk->num_candidate_regions;reg_id++, candidate_regions_spr++){
                    if (seed_arr[sid]>=candidate_regions_spr->start_loc_on_ref && seed_arr[sid]<=candidate_regions_spr->end_loc_on_ref){
                        int nn=candidate_regions_spr->num_seeds;
                        if (nn<10000){
                            candidate_regions_spr->seeds[nn].seed_loc_on_ref=seed_arr[sid];
                            candidate_regions_spr->seeds[nn].seed_loc_on_read=loc_on_read;
                            candidate_regions_spr->num_seeds++;
                        }else{
                            std::cout<<"oh"<<std::endl;
                        }
                        break;
                    }
                }*/
                L=0;
                R = sbk->num_candidate_regions-1;
                while(L<=R){
                    mid = (L+R)/2;
                    candidate_regions_spr=candidate_regions+mid;
                    if (seed_arr[sid]>=candidate_regions_spr->start_loc_on_ref && seed_arr[sid]<=candidate_regions_spr->end_loc_on_ref){
                        int nn=candidate_regions_spr->num_seeds;
                        if (nn<10000){
                            candidate_regions_spr->seeds[nn].seed_loc_on_ref=seed_arr[sid];
                            candidate_regions_spr->seeds[nn].seed_loc_on_read=loc_on_read;
                            candidate_regions_spr->num_seeds++;
                        }else{
                            std::cout<<"oh"<<std::endl;
                        }
                        break;
                    }
                    else if (seed_arr[sid]<candidate_regions_spr->start_loc_on_ref){
                        R=mid-1;
                    }
                    else{
                        L=mid+1;
                    }
                }
            }
        loc_on_read+=kmer_size;
	}
	//return used_segs;
}


bool compare_temp_arr_elements(temp_arr_element i, temp_arr_element j){
    return i.score>j.score;
}


int chaining (Short_seed* seeds_list, int num_seeds, int read_size){
    Short_seed* seeds_list_spr1 = seeds_list;
    Short_seed* seeds_list_spr2;
    int i, j;
    double temp_score;
    int temp_gap_penalty;
    double max_score=-1;
    int max_index=-1;
    for (i=0; i<num_seeds; ++i, ++seeds_list_spr1){
        seeds_list_spr1->visited=false;
        seeds_list_spr1->chaining_score=kmer_size;
        seeds_list_spr1->prev=-1;
        for (j=0 , seeds_list_spr2 = seeds_list_spr1-1;j<i && j<20; ++j , --seeds_list_spr2){
            if(seeds_list_spr1->seed_loc_on_read<=seeds_list_spr2->seed_loc_on_read ||
                seeds_list_spr1->seed_loc_on_ref==seeds_list_spr2->seed_loc_on_ref){
                continue;
            }
            temp_gap_penalty = max(seeds_list_spr1->seed_loc_on_read - seeds_list_spr2->seed_loc_on_read,
                                   seeds_list_spr1->seed_loc_on_ref - seeds_list_spr2->seed_loc_on_ref);
            if(temp_gap_penalty>MAX_GAP)
                continue;
            temp_gap_penalty = abs((seeds_list_spr1->seed_loc_on_read - seeds_list_spr2->seed_loc_on_read)-
                                   (seeds_list_spr1->seed_loc_on_ref - seeds_list_spr2->seed_loc_on_ref));
            temp_score=seeds_list_spr2->chaining_score + min(seeds_list_spr1->seed_loc_on_read - seeds_list_spr2->seed_loc_on_read,
                                   seeds_list_spr1->seed_loc_on_ref - seeds_list_spr2->seed_loc_on_ref);
            temp_score-=(temp_gap_penalty*0.1);//????????????????????????????????
            if (temp_score > seeds_list_spr1->chaining_score){
                seeds_list_spr1->chaining_score=temp_score;
                seeds_list_spr1->prev=i-j-1;
            }
        }
        if (seeds_list_spr1->chaining_score>max_score){
            max_score=seeds_list_spr1->chaining_score;
            max_index=i;
        }
    }
    return max_index;
}


int
get_candidates(volume_t* ref,
			   SeedingBK* sbk,
			   const char* read,
			   const int read_id,
			   const int read_size,
			   const char chain,
			   candidate_save* candidates,
			   int candidatenum)
{
	candidate_region* candidate_regions=sbk->candidate_regions;
	candidate_region* candidate_regions_spr=candidate_regions;

	int index;
	float chain_score;
	candidate_save *candidate_loc = candidates, candidate_temp;
	int location_loc[4];
	int i, j, k;
	bool sw=false;
    if(read_id==3092){
        int y=0;
    }
    for (int reg_id=0;reg_id<sbk->num_candidate_regions;reg_id++, candidate_regions_spr++){
        int nn=candidate_regions_spr->num_seeds;
        if(nn== 0) continue;
        sort(candidate_regions_spr->seeds, candidate_regions_spr->seeds+nn,compare_short_seeds);
        //temp_arr_element* temp_arr;
        //safe_malloc(temp_arr,temp_arr_element,nn);
        index=chaining (candidate_regions_spr->seeds, nn, read_size);
        chain_score=candidate_regions_spr->seeds[index].chaining_score;
        location_loc[2] = candidate_regions_spr->seeds[index].seed_loc_on_ref;
        location_loc[3] = candidate_regions_spr->seeds[index].seed_loc_on_read;//??????????????
        while (candidate_regions_spr->seeds[index].prev!=-1){
                index=candidate_regions_spr->seeds[index].prev;
        }

        location_loc[0] = candidate_regions_spr->seeds[index].seed_loc_on_ref;
        location_loc[1] = candidate_regions_spr->seeds[index].seed_loc_on_read;//???????????????????????????
        if (location_loc[2]-location_loc[0]<2000)
            continue;
        candidate_temp.chain = chain;
        candidate_temp.score=chain_score;
			//int loc_list = location_loc[0];
        int sid = get_read_id_from_offset_list(ref->offset_list, location_loc[0]);
        int sstart = ref->offset_list->offset_list[sid].offset;
        int ssize = ref->offset_list->offset_list[sid].size;
        int send = sstart + ssize + 1;/////////////why +1
        sid += ref->start_read_id;//////important
        candidate_temp.readno = sid;
        candidate_temp.readstart = sstart;
        int left_length1 = location_loc[0] - sstart + kmer_size - 1;
        int right_length1 = send - location_loc[0];
        int left_length2 = location_loc[1] + kmer_size - 1;
        int right_length2 = read_size - location_loc[1];
        int num1 = (left_length1 > left_length2) ? left_length2 : left_length1;////////////////what
        int num2 = (right_length1 > right_length2) ? right_length2 : right_length1;
        if (num1 + num2 < min_kmer_dist) {
            continue;/////*****************why
        }
        candidate_temp.loc1=location_loc[0] - sstart;
        candidate_temp.num1=num1;
        candidate_temp.loc2=location_loc[1];
        candidate_temp.num2=num2;
        candidate_temp.left1=left_length1;
        candidate_temp.left2=left_length2;
        candidate_temp.right1=right_length1;
        candidate_temp.right2=right_length2;
            //insert canidate position or delete this position
        int low=0;
        int high=candidatenum-1;
        int mid;
        while(low<=high)
        {
            mid=(low+high)/2;
            if(mid>=candidatenum||candidate_loc[mid].score<candidate_temp.score) high=mid-1;
            else low=mid+1;
        }
        if(candidatenum<MAXC)for(j=candidatenum-1; j>high; j--)candidate_loc[j+1]=candidate_loc[j];
        else for(j=candidatenum-2; j>high; j--)candidate_loc[j+1]=candidate_loc[j];
        if(high+1<MAXC)candidate_loc[high+1]=candidate_temp;
        if(candidatenum<MAXC)candidatenum++;
        else candidatenum=MAXC;
    }
	sbk->num_candidate_regions=0;
	sbk->num_long_seeds=0;
	if (sbk->size_of_long_seeds_array!=initial_arr_size){
        safe_free(sbk->long_seeds);
        safe_malloc(sbk->long_seeds, Long_seed, initial_arr_size);
        sbk->size_of_long_seeds_array=initial_arr_size;
	}

	return candidatenum;
}



void
fill_m4record(GapAligner* aligner, const int qid, const int sid,
			  const char qchain, int qsize, int ssize,
			  int qstart, int sstart, int vscore, M4Record* m)
{
	if (qchain == 'F')
	{
		m->qid = sid;
		m->sid = qid;
		m->ident = aligner->calc_ident();
		m->vscore = vscore;
		m->qdir = 0;
		m->qoff = aligner->target_start();
		m->qend = aligner->target_end();
		m->qsize = ssize;
		m->sdir = 0;
		m->soff = aligner->query_start();
		m->send = aligner->query_end();
		m->ssize = qsize;
		m->qext = sstart;
		m->sext = qstart;
	}
	else
	{
		m->qid = sid;
		m->sid = qid;
		m->ident = aligner->calc_ident();
		m->vscore = vscore;
		m->qdir = 0;
		m->qoff = aligner->target_start();
		m->qend = aligner->target_end();
		m->qsize = ssize;
		m->sdir = 1;
		m->soff = qsize - aligner->query_end();
		m->send = qsize - aligner->query_start();
		m->ssize = qsize;
		m->qext = sstart;
		m->sext = qsize - 1 - qstart;
	}
}


void output_m4record(ostream& out, const M4Record& m4)
{
	const char sep = '\t';

	out << m4qid(m4)    << sep
	    << m4sid(m4)    << sep
	    << m4ident(m4) << sep
	    << m4vscore(m4)   << sep
	    << m4qdir(m4)   << sep
	    << m4qoff(m4)   << sep
	    << m4qend(m4)   << sep
	    << m4qsize(m4)  << sep
	    << m4sdir(m4)   << sep
	    << m4soff(m4)   << sep
	    << m4send(m4)   << sep
	    << m4ssize(m4);

	if (output_gapped_start_point)
		out << sep
			<< m4qext(m4)	<< sep
			<< m4sext(m4);
	out << "\n";
}

void
print_m4record_list(ostream* out, M4Record* m4_list, int num_m4)
{
	for (int i = 0; i < num_m4; ++i) output_m4record(*out, m4_list[i]);
}

struct CmpM4RecordByQidAndOvlpSize
{
	bool operator()(const M4Record& a, const M4Record& b)
	{
		if (m4qid(a) != m4qid(b)) return m4qid(a) < m4qid(b);
		int o1 = M4RecordOverlapSize(a);
		int o2 = M4RecordOverlapSize(b);
		return o1 > o2;
	}
};

void
check_records_containment(M4Record* m4v, int s, int e, int* valid)
{
	const int soft = 100;

	for (int i = s; i < e; ++i)
	{
		if (!valid[i]) continue;
		int qb1 = m4qoff(m4v[i]);
		int qe1 = m4qend(m4v[i]);
		int sb1 = m4soff(m4v[i]);
		int se1 = m4send(m4v[i]);
		for (int j = i + 1; j < e; ++j)
		{
			if (!valid[j]) continue;
			if (m4sdir(m4v[i]) != m4sdir(m4v[j])) continue;
			int qb2 = m4qoff(m4v[j]);
			int qe2 = m4qend(m4v[j]);
			int sb2 = m4soff(m4v[j]);
			int se2 = m4send(m4v[j]);

			if (qb2 + soft >= qb1 && qe2 - soft <= qe1 && sb2 + soft >= sb1 && se2 - soft <= se1) valid[j] = 0;
		}
	}
}

void
append_m4v(M4Record* glist, int* glist_size,
		   M4Record* llist, int* llist_size,
		   ostream* out, pthread_mutex_t* results_write_lock)
{
	sort(llist, llist + *llist_size, CmpM4RecordByQidAndOvlpSize());
	int i = 0, j;
	int valid[*llist_size];
	fill(valid, valid + *llist_size, 1);
	while (i < *llist_size)
	{
		idx_t qid = m4qid(llist[i]);
		j = i + 1;
		while (j < *llist_size && m4qid(llist[j]) == qid) ++j;
		if (j - i > 1) check_records_containment(llist, i, j, valid);
		i = j;
	}

	if ((*glist_size) + (*llist_size) > PWThreadData::kResultListSize)
	{
		pthread_mutex_lock(results_write_lock);
		print_m4record_list(out, glist, *glist_size);
		*glist_size = 0;
		pthread_mutex_unlock(results_write_lock);
	}

	for (i = 0; i < *llist_size; ++i)
		if (valid[i])
		{
			glist[*glist_size] = llist[i];
			++(*glist_size);
		}

	*llist_size = 0;
}

inline void
get_next_chunk_reads(PWThreadData* data, int& Lid, int& Rid)
{
	pthread_mutex_lock(&data->read_retrieve_lock);
	Lid = data->next_processed_id;
	Rid = Lid + CHUNK_SIZE;
	if (Rid > data->reads->num_reads) Rid = data->reads->num_reads;
	data->next_processed_id += CHUNK_SIZE;
	pthread_mutex_unlock(&data->read_retrieve_lock);
}

/*void
pairwise_mapping(PWThreadData* data, int tid)
{
}*/

void
pairwise_mapping(PWThreadData* data, int tid)
{
	char *read, *read1, *read2, *subject;
	safe_malloc(read1, char, MSS);
	safe_malloc(read2, char, MSS);
	safe_malloc(subject, char, MSS);
	SeedingBK* sbk = new SeedingBK();
	candidate_save candidates[MAXC];//candidate
	int num_candidates = 0;
	r_assert(data->m4_results);
	M4Record* m4_list = data->m4_results[tid];
	int m4_list_size = 0;
	M4Record* m4v = new M4Record[MAXC];
	bool self_detection = data->self_detection;
	int num_m4 = 0;
	GapAligner* aligner = NULL;
	if (data->options->tech == TECH_PACBIO) {
		aligner = new DiffAligner(0);
	} else if (data->options->tech == TECH_NANOPORE) {
		aligner = new XdropAligner(0);
	} else {
		ERROR("TECH must be either %d or %d", TECH_PACBIO, TECH_NANOPORE);
	}

	int rid, Lid, Rid;
	while (1)
	{
		get_next_chunk_reads(data, Lid, Rid);
		if (Lid >= data->reads->num_reads) break;
		for (rid = Lid; rid < Rid; ++rid)
		{
			int rsize = data->reads->offset_list->offset_list[rid].size;
			if (rsize >= MAX_SEQ_SIZE) {
                cout << "rsize = " << rsize << "\t" << MAX_SEQ_SIZE << endl;
                abort();
            }
			extract_one_seq(data->reads, rid, read1);
			reverse_complement(read2, read1, rsize);
			int s;
			int chain;
			num_candidates = 0;
			for (s = 0; s < 2; ++s)
			{
				if (s%2) { chain = REV; read = read2; }
				else { chain = FWD; read = read1; }
				int num_segs;
				int rstart= data->reads->offset_list->offset_list[rid].offset;

				// = seeding(read, rsize, data->ridx, sbk);
				if (self_detection)
                    self_seeding(read, rstart, rsize, chain, data->reference, data->first_index, data->second_index, sbk);
                else
                    seeding(read, rstart, rsize, chain, data->reference, data->first_index, data->second_index, sbk);
				if (sbk->num_candidate_regions>0)
                    num_candidates = get_candidates(data->reference,
											sbk,
											read,
											rid + data->reads->start_read_id,
											rsize,
											chain,
											candidates,
											num_candidates);
			}

			for (s = 0; s < num_candidates; ++s)
			{
				if (candidates[s].chain == 'F') read = read1;
				else read = read2;
				extract_one_seq(data->reference, candidates[s].readno - data->reference->start_read_id, subject);
				int sstart = candidates[s].loc1;
				int qstart = candidates[s].loc2;
				if (qstart && sstart)
				{
					qstart += kmer_size / 2;
					sstart += kmer_size / 2;
				}
				int ssize = data->reference->offset_list->offset_list[candidates[s].readno - data->reference->start_read_id].size;

				int flag = aligner->go(read, qstart, rsize, subject, sstart, ssize, min_align_size);

				if (flag)
				{
					fill_m4record(aligner, rid + data->reads->start_read_id,
								  candidates[s].readno, candidates[s].chain,
								  rsize, ssize, qstart, sstart, candidates[s].score,
								  m4v + num_m4);
					++num_m4;
				}
			}

			append_m4v(m4_list, &m4_list_size, m4v, &num_m4, data->out, &data->result_write_lock);
		}
	}

		if (m4_list_size)
		{
			pthread_mutex_lock(&data->result_write_lock);
			print_m4record_list(data->out, m4_list, m4_list_size);
			m4_list_size = 0;
			pthread_mutex_unlock(&data->result_write_lock);
		}

		safe_free(read1);
		safe_free(read2);
		safe_free(subject);
		delete sbk;
		delete aligner;
		delete[] m4v;
}

void
candidate_detect(PWThreadData* data, int tid)
{
	char *read, *read1, *read2, *subject;
	safe_malloc(read1, char, MAX_SEQ_SIZE);
	safe_malloc(read2, char, MAX_SEQ_SIZE);
	safe_malloc(subject, char, MAX_SEQ_SIZE);
	SeedingBK* sbk = new SeedingBK();//sh uld be canged, new can be removed
	//safe_calloc(sbk,SeedingBK,1);
	Candidate candidates[MAXC];
	int num_candidates = 0;
	r_assert(data->ec_results);
	ExtensionCandidate* eclist = data->ec_results[tid];
	bool self_detection = data->self_detection;
	int nec = 0;
	ExtensionCandidate ec;

	int rid, Lid, Rid;
	while (1)
	{
		get_next_chunk_reads(data, Lid, Rid);
		if (Lid >= data->reads->num_reads) break;
		//if (Lid >=1000) break;
        for (rid = Lid; rid < Rid; ++rid)
        {
            int rsize = data->reads->offset_list->offset_list[rid].size;
            //int rstart= data->reads->offset_list->offset_list[rid].offset;
            if (rsize >= MAX_SEQ_SIZE) {
                cout << "rsize = " << rsize << "\t" << MAX_SEQ_SIZE << endl;
                abort();
            }
            extract_one_seq(data->reads, rid, read1);
            reverse_complement(read2, read1, rsize);
            int s;
            int chain;
            num_candidates = 0;
            for (s = 0; s < 2; ++s)
            {
                if (s%2) { chain = REV; read = read2; }
                else { chain = FWD; read = read1; }
                int num_segs;
                //std::cout<<rid<<std::endl;
                int rstart= data->reads->offset_list->offset_list[rid].offset;

                if (self_detection)
                    self_seeding(read, rstart, rsize, chain, data->reference, data->first_index, data->second_index, sbk);
                else
                    seeding(read, rstart, rsize, chain, data->reference, data->first_index, data->second_index, sbk);
                if (sbk->num_candidate_regions>0)
                    num_candidates = get_candidates(data->reference,
											sbk,
											read,
											rid + data->reads->start_read_id,
											rsize,
											chain,
											candidates,
											num_candidates);
            }

            for (s = 0; s < num_candidates; ++s)
            {
                int qstart = candidates[s].loc2;
                int sstart = candidates[s].loc1;
                if (qstart && sstart)
                {
                    qstart += kmer_size / 2;
                    sstart += kmer_size / 2;
                }
                int qdir = candidates[s].chain;
                int sdir = FWD;
                int qid = rid + data->reads->start_read_id;
                int sid = candidates[s].readno;
                int score = candidates[s].score;

                ec.qid = qid;
                ec.qdir = qdir;
                ec.qext = qstart;
                ec.sid = sid;
                ec.sdir = sdir;
                ec.sext = sstart;
                ec.score = score;
                ec.qsize = rsize;
                ec.ssize = data->reference->offset_list->offset_list[sid - data->reference->start_read_id].size;
                if (ec.qdir == REV) ec.qext = ec.qsize - 1 - ec.qext;
                 if (ec.sdir == REV) ec.sext = ec.ssize - 1 - ec.sext;
                eclist[nec] = ec;
                ++nec;
                if (nec == PWThreadData::kResultListSize)
                {
                    pthread_mutex_lock(&data->result_write_lock);
                    for (int i = 0; i < nec; ++i) (*data->out) << eclist[i];
                    nec = 0;
                    pthread_mutex_unlock(&data->result_write_lock);
                }
            }
        }
	}

	if (nec)
	{
		pthread_mutex_lock(&data->result_write_lock);
		for (int i = 0; i < nec; ++i) (*data->out) << eclist[i];
		nec = 0;
		pthread_mutex_unlock(&data->result_write_lock);
	}

	safe_free(read1);
	safe_free(read2);
	safe_free(subject);
	delete sbk;
	//destroy_read_index(read_index, read_index_count, POSITIONS_PER_KMER);
}


void*
multi_thread_func(void* p)
{
	PWThreadData* data = (PWThreadData*)p;
	int t = 0;
	pthread_mutex_lock(&data->id_lock);
	t = data->used_thread_id;
	++data->used_thread_id;
	pthread_mutex_unlock(&data->id_lock);
	r_assert(data->options->task == TASK_SEED || data->options->task == TASK_ALN);
	if (data->options->task == TASK_SEED) candidate_detect(data, t);
	else pairwise_mapping(data, t);
	return NULL;
}

void
process_one_volume(options_t* options, const int svid, const int evid, volume_names_t* vn, ostream* out)
{
	MAXC = options->num_candidates;
	output_gapped_start_point = options->output_gapped_start_point;
	min_align_size = options->min_align_size;
	//min_kmer_match = options->min_kmer_match;

	if (options->tech == TECH_PACBIO) {
		error_rate=pacbio_error_rate;
		 min_kmer_dist= 1800;
	} else if (options->tech == TECH_NANOPORE) {
		error_rate=nanopore_error_rate;
		min_kmer_dist = 400;
	} else {
		ERROR("TECH must be either %d or %d", TECH_PACBIO, TECH_NANOPORE);
	}

	const char* ref_name = get_vol_name(vn, svid);
	volume_t* ref = load_volume(ref_name);
	//ref_index* first_index = (ref_index*)malloc(sizeof(ref_index));
	//ref_index* second_index = (ref_index*)malloc(sizeof(ref_index));
    ref_index* first_index = create_ref_index(ref, kmer_size, options->num_threads);
	ref_index* second_index = create_second_index(ref, kmer_size, second_kmer_size, options->num_threads, first_index);
	index_separator=ref->curr+1;
	pthread_t tids[options->num_threads];
	char volume_process_info[1024];;
	int vid, tid;
	for(vid = svid; vid < evid; ++vid)
	{
		sprintf(volume_process_info, "process volume %d", vid);
		DynamicTimer dtimer(volume_process_info);
		const char* read_name = get_vol_name(vn, vid);
		LOG(stderr, "processing %s\n", read_name);
		volume_t* read = load_volume(read_name);
		PWThreadData* data = new PWThreadData(options, ref, read, first_index, second_index, out);
		if (vid == svid)
            data->self_detection=true;
        else data->self_detection=false;
		for (tid = 0; tid < options->num_threads; ++tid)
		{
			int err_code = pthread_create(tids + tid, NULL, multi_thread_func, (void*)data);
			if (err_code)
			{
				LOG(stderr, "Error: return code is %d\n", err_code);
				abort();
			}
		}
		for (tid = 0; tid < options->num_threads; ++tid) pthread_join(tids[tid], NULL);
		read = delete_volume_t(read);
		delete data;
	}
	ref = delete_volume_t(ref);
	first_index = destroy_ref_index(first_index);
	second_index = destroy_ref_index(second_index);
}

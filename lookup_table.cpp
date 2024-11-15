#include "lookup_table.h"
#include "packed_db.h"

#include <cstdio>

ref_index*
destroy_ref_index(ref_index* ridx)
{
	safe_free(ridx->kmer_counts);
	safe_free(ridx->kmer_starts);
	safe_free(ridx->kmer_offsets);
	safe_free(ridx);
	return NULL;
}

typedef struct
{
	ref_index* ridx;
	uint32_t min_key;
	uint32_t max_key;
	volume_t* v;
	int kmer_size;
} ref_index_thread_info;

typedef struct
{
	ref_index* first_ridx;
	uint32_t min_key;
	uint32_t max_key;
	volume_t* v;
	int kmer_size;
	int second_kmer_size;
	ref_index* second_ridx;
	supporting_node** output_kmer_lists;
	int* all_long_kmers_per_thread;
} create_supporting_list_for_second_index_thread_info;

typedef struct
{
	ref_index* second_index;
	uint32_t min_key;
	uint32_t max_key;
	supporting_node** kmer_lists;
	int offset_start_point;
} fill_second_index_thread_info;


int quick_partition(float *a,int *b,int start_index,int end_index)
{
    float pivot=a[end_index];

    int P_index=start_index;
    int i,t1;
    float t2;

    for(i=start_index;i<end_index;i++)
    {
    	if(a[i]<=pivot)
        {
            t2=a[i];
            a[i]=a[P_index];
            a[P_index]=t2;
            t1=b[i];
            b[i]=b[P_index];
            b[P_index]=t1;
            P_index++;
        }
     }
     //Now exchanging value of
     //pivot and P-index
      t2=a[end_index];
      a[end_index]=a[P_index];
      a[P_index]=t2;
      t1=b[end_index];
      b[end_index]=b[P_index];
      b[P_index]=t1;

     //at last returning the pivot value index
     return P_index;
 }
 void Quicksort(float *a, int *b, int start_index,int end_index)
 {
    if(start_index<end_index)
    {
         int P_index=quick_partition(a,b,start_index,end_index);
             Quicksort(a,b,start_index,P_index-1);
             Quicksort(a,b,P_index+1,end_index);
    }
}


void*
fill_ref_index_offsets_func(void* arg)
{
	ref_index_thread_info* riti = (ref_index_thread_info*)(arg);
	volume_t* v = riti->v;
	int num_reads = v->num_reads;
	ref_index* index = riti->ridx;
	int kmer_size = riti->kmer_size;
	uint32_t index_count = 1 << (kmer_size * 2);
	uint32_t leftnum = 34 - 2 * kmer_size;
	int i, j;
	for (i = 0; i < num_reads; ++i)
	{
		int read_start = v->offset_list->offset_list[i].offset;
		int read_size = v->offset_list->offset_list[i].size;
		uint32_t eit = 0;
		for (j = 0; j < read_size; ++j)
		{
			int k = read_start + j;
			uint8_t c = PackedDB::get_char(v->data, k);
			eit = (eit << 2) | c;
			assert(eit < index_count);
			if (j >= kmer_size - 1)
			{
				if (index->kmer_starts[eit] && eit >= riti->min_key && eit <= riti->max_key)
				{
					index->kmer_starts[eit][index->kmer_counts[eit]] = k + 1 - kmer_size;
					++index->kmer_counts[eit];
				}
				eit <<= leftnum;
				eit >>= leftnum;
			}
		}
	}

	return NULL;
}

ref_index*
create_ref_index(volume_t* v, int kmer_size, int num_threads)
{
	DynamicTimer dtimer(__func__);
	uint32_t index_count = 1 << (kmer_size * 2);
	uint32_t leftnum = 34 - 2 * kmer_size;
	ref_index* index = (ref_index*)malloc(sizeof(ref_index));
	safe_calloc(index->kmer_counts, int, index_count);
	int num_reads = v->num_reads;
	for (uint32_t i = 0; i != index_count; ++i) assert(index->kmer_counts[i] == 0);
	for (int i = 0; i != num_reads; ++i)
	{
		int read_start = v->offset_list->offset_list[i].offset;
		int read_size = v->offset_list->offset_list[i].size;
		uint32_t eit = 0;
		for (int j = 0; j < read_size; ++j)
		{
			int k = read_start + j;
			uint8_t c = PackedDB::get_char(v->data, k);
			assert(c>= 0 && c < 4);
			eit = (eit << 2) | c;
			if (j >= kmer_size - 1)
			{
				assert(eit < index_count);
				++index->kmer_counts[eit];
				eit = eit << leftnum;
				eit = eit >> leftnum;
			}
		}
	}

	int num_kmers = 0;
	for (uint32_t i = 0; i != index_count; ++i)
	{
		if (index->kmer_counts[i] > 128) index->kmer_counts[i] = 0;
		num_kmers += index->kmer_counts[i];
	}
	printf("number of kmers: %d\n", num_kmers);
	safe_malloc(index->kmer_offsets, int, num_kmers);
	safe_malloc(index->kmer_starts, int*, index_count);

	if (v->curr < 10 * 1000000) num_threads = 1;
	int kmers_per_thread = (num_kmers + num_threads - 1) / num_threads;
	fprintf(stderr, "%d threads are used for filling offset lists.\n", num_threads);
	uint32_t hash_boundaries[2 * num_threads];
	uint32_t L = 0;
	num_kmers = 0;
	int kmer_cnt = 0;
	int tid = 0;
	for (uint32_t i = 0; i != index_count; ++i)
	{
		if (index->kmer_counts[i])
		{
			index->kmer_starts[i] = index->kmer_offsets + num_kmers;
			num_kmers += index->kmer_counts[i];
			kmer_cnt += index->kmer_counts[i];
			index->kmer_counts[i] = 0;

			if (kmer_cnt >= kmers_per_thread)
			{
				printf("thread %d: %d\t%d\n", tid, L, i);
				hash_boundaries[2 * tid] = L;
				hash_boundaries[2 * tid + 1] = i;
				++tid;
				L = i + 1;
				kmer_cnt = 0;
			}
		}
		else
		{
			index->kmer_starts[i] = NULL;
		}
	}
	if (kmer_cnt)
	{
		printf("thread %d: %d\t%d\n", tid, L, index_count - 1);
		hash_boundaries[2 * tid] = L;
		hash_boundaries[2 * tid + 1] = index_count - 1;
	}

	ref_index_thread_info ritis[num_threads];
	for (int i = 0; i != num_threads; ++i)
	{
		ritis[i].ridx = index;
		ritis[i].min_key = hash_boundaries[2 * i];
		ritis[i].max_key = hash_boundaries[2 * i + 1];
		ritis[i].v = v;
		ritis[i].kmer_size = kmer_size;
	}

	pthread_t tids[num_threads];
	for (int j = 0; j < num_threads; ++j)
		pthread_create(tids + j, NULL, fill_ref_index_offsets_func, (void*)(ritis + j));
	for (int j = 0; j < num_threads; ++j)
		pthread_join(tids[j], NULL);

	return index;
}

void* create_supporing_lists_for_second_index_func(void* arg){
    create_supporting_list_for_second_index_thread_info* csl=(create_supporting_list_for_second_index_thread_info*)(arg);
    volume_t* v = csl->v;
	//int num_reads = v->num_reads;
	ref_index* first_index = csl->first_ridx;
	ref_index* second_index=csl->second_ridx;
	int kmer_size = csl->kmer_size;
	int second_kmer_size=csl->second_kmer_size;
	int t=second_kmer_size-kmer_size;
	//int end_th=v->curr-kmer_size+1;
	int* num_long_kmers = csl->all_long_kmers_per_thread;
	(*num_long_kmers)=0;
	supporting_node** long_kmer_lists = csl->output_kmer_lists;
	//uint32_t index_count = 1 << (kmer_size * 2);
	//uint32_t leftnum = 34 - 2 * kmer_size;
	bool visited[128];///////////////////////////////////////////////////////////////////*****************/
	for (uint32_t i=csl->min_key; i<=csl->max_key; i++){
        //supporting_node* curr_list_head=long_kmer_lists[i];
        //curr_list_head=NULL;
        if(i==16776865){
            int y=0;
        }
        long_kmer_lists[i]=NULL;
        if (first_index->kmer_counts[i]<2) continue;
        //safe_calloc(visited,bool,first_index->kmer_counts[i]);
        for (int j=0; j<first_index->kmer_counts[i];j++)
            visited[j]=false;
        for (int j=0; j<first_index->kmer_counts[i]-1;j++){
            if (visited[j]) {
                continue;
            }
            visited[j]=true;
            bool first_seed=true;
            int loc1=first_index->kmer_starts[i][j];
            /*(*num_long_kmers)++;
            second_index->kmer_counts[i]++;
            supporting_node* new_node = (supporting_node*)malloc(sizeof(supporting_node));
            new_node->loc1=loc1;
            new_node->next_node=long_kmer_lists[i];
            long_kmer_lists[i]=new_node;*/
            int nextLoc1=loc1+kmer_size;
            //if (nextLoc1>end_th) continue;
            for (int k=j+1; k<first_index->kmer_counts[i];k++){
                if (visited[k]){
                    continue;
                }
                int loc2=first_index->kmer_starts[i][k];
                int nextLoc2=loc2+kmer_size;
                //if (nextLoc2>end_th) continue;
                int num_sim=0;
                for (int m=0;m<t;m++){
                    char ch1 = PackedDB::get_char(v->data, nextLoc1+m);
                    char ch2 = PackedDB::get_char(v->data, nextLoc2+m);
                    if(ch1 == ch2)
                        num_sim++;
                    else
                        break;
                }
                //uint16_t ch1 = PackedDB::get_8_next_nuc(v->data, nextLoc1);
                //uint16_t ch2 = PackedDB::get_8_next_nuc(v->data, nextLoc2);

                //if (ch1==ch2){
                if (num_sim==t){
                    visited[k]=true;
                    if (first_seed){
                            first_seed=false;
                            //add separator
                            (*num_long_kmers)++;
                            second_index->kmer_counts[i]++;
                            supporting_node* new_node1 = (supporting_node*)malloc(sizeof(supporting_node));
                            new_node1->loc1=v->curr+1;
                            new_node1->next_node=long_kmer_lists[i];
                            long_kmer_lists[i]=new_node1;
                            //add first and second seeds
                            (*num_long_kmers)++;
                            second_index->kmer_counts[i]++;
                            supporting_node* new_node2 = (supporting_node*)malloc(sizeof(supporting_node));
                            new_node2->loc1=loc1;
                            new_node2->next_node=long_kmer_lists[i];
                            long_kmer_lists[i]=new_node2;

                            (*num_long_kmers)++;
                            second_index->kmer_counts[i]++;
                            supporting_node* new_node3 = (supporting_node*)malloc(sizeof(supporting_node));
                            new_node3->loc1=loc2;
                            new_node3->next_node=long_kmer_lists[i];
                            long_kmer_lists[i]=new_node3;
                    }
                    else{
                        (*num_long_kmers)++;
                        second_index->kmer_counts[i]++;
                        supporting_node* new_node3 = (supporting_node*)malloc(sizeof(supporting_node));
                        new_node3->loc1=loc2;
                        new_node3->next_node=long_kmer_lists[i];
                        long_kmer_lists[i]=new_node3;
                    }
                }
            }
        }
        //safe_free(visited);
    }
    return NULL;
}

void* fill_second_index_offsets_func (void* arg){
    fill_second_index_thread_info* siti = (fill_second_index_thread_info*)(arg);
    ref_index* second_index = siti->second_index;
    int offset_point = siti->offset_start_point;
    supporting_node** kmer_lists = siti->kmer_lists;
    supporting_node* curr_kmer_list;
    supporting_node* next_node;
    for (uint32_t i=siti->min_key; i<=siti->max_key;++i){
        if (second_index->kmer_counts[i]==0){
            second_index->kmer_starts[i]=NULL;
            continue;
        }
        if(i==16776865){
            int y=0;
        }
        second_index->kmer_starts[i]=second_index->kmer_offsets+offset_point;
        curr_kmer_list=kmer_lists[i];
        for (int j=0; j<second_index->kmer_counts[i];j++){
            second_index->kmer_offsets[offset_point++]=curr_kmer_list->loc1;
            //second_index->kmer_offsets[offset_point++]=curr_kmer_list->loc2;
            if (curr_kmer_list->loc1==8894425){
                int y=0;
            }
            next_node=curr_kmer_list->next_node;

            safe_free(curr_kmer_list);//???????????????????????????
            curr_kmer_list=next_node;
        }
    }
}

ref_index*
create_second_index(volume_t* v, int kmer_size, int second_kmer_size, const int num_threads, ref_index* first_index){
    DynamicTimer dtimer(__func__);
	uint32_t index_count = 1 << (kmer_size * 2);
	uint32_t leftnum = 34 - 2 * kmer_size;
	ref_index* second_index = (ref_index*)malloc(sizeof(ref_index));
	safe_calloc(second_index->kmer_counts, int, index_count);
	safe_malloc(second_index->kmer_starts, int*, index_count);
	//int num_reads = v->num_reads;
	for (uint32_t i = 0; i != index_count; ++i) {
        assert(second_index->kmer_counts [i] == 0);
	}

	uint32_t hash_boundaries[2 * num_threads];
	uint32_t H, L = 0;
	int rows_per_thread= index_count/num_threads;
	for(int i=0;i<num_threads;i++){
        hash_boundaries[i*2]=L;
        H=L+rows_per_thread-1;
        hash_boundaries[i*2+1]=H;
        L=H+1;
	}
	hash_boundaries[2*num_threads-1]=index_count-1;
	create_supporting_list_for_second_index_thread_info cslsis[num_threads];
	int long_kmers_per_thread[num_threads];
	supporting_node** kmer_lists_per_thread;
	safe_malloc(kmer_lists_per_thread,supporting_node*, index_count);

	for (int i = 0; i != num_threads; ++i)
	{
		cslsis[i].first_ridx = first_index;
		cslsis[i].min_key = hash_boundaries[2 * i];
		cslsis[i].max_key = hash_boundaries[2 * i + 1];
		cslsis[i].v = v;
		cslsis[i].kmer_size = kmer_size;
		cslsis[i].second_kmer_size = second_kmer_size;
		cslsis[i].second_ridx=second_index;
		cslsis[i].all_long_kmers_per_thread=long_kmers_per_thread+i;
		cslsis[i].output_kmer_lists=kmer_lists_per_thread;
	}

	pthread_t tids[num_threads];
	for (int j = 0; j < num_threads; ++j)
		pthread_create(tids + j, NULL, create_supporing_lists_for_second_index_func, (void*)(cslsis + j));
	for (int j = 0; j < num_threads; ++j)
		pthread_join(tids[j], NULL);
    ///***********complete second index***************
    int all_long_kmers=0;
    int offset_start_point[num_threads];
    for (int j=0; j<num_threads; ++j){
        offset_start_point[j]=all_long_kmers;
        all_long_kmers+=long_kmers_per_thread[j];
    }
    safe_malloc(second_index->kmer_offsets, int, all_long_kmers);
    fill_second_index_thread_info sitis[num_threads];
    for (int i=0; i!=num_threads; ++i){
        sitis[i].kmer_lists=kmer_lists_per_thread;
        sitis[i].min_key=hash_boundaries[2 * i];
        sitis[i].max_key=hash_boundaries[2 * i + 1];
        sitis[i].offset_start_point=offset_start_point[i];
        sitis[i].second_index=second_index;
    }
    for (int j = 0; j < num_threads; ++j)
		pthread_create(tids + j, NULL, fill_second_index_offsets_func, (void*)(sitis + j));
	for (int j = 0; j < num_threads; ++j)
		pthread_join(tids[j], NULL);

    safe_free(kmer_lists_per_thread);
    return second_index;
}

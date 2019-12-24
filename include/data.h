#ifndef _DATA_H
#define _DATA_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct data_t{
	int data_cnt;
	int* alloc_len;
	int* busy_len;
	int** data;
} data_t;

//create dynamicaly data object with "data_size" sections with size "defolt_size" in every section
data_t data_t_construct(int data_size, unsigned int defolt_size);

void data_t_destruct(data_t data); //deep destructor

//add information - "elem_id" into section "coord"
//and add "addition_size" of space into this section if it doesn't have enough place
void add_elem_to_data_ext(data_t data, int elem_id, int coord, int addition_size);

//add information - "elem_id" into section "coord"
//cares about memory automatically
void add_elem_to_data(data_t data, int elem_id, int coord);

//debug print of data_t object
void data_t_dump(data_t data);

//print used size of every section of "data"
void data_t_busy_len_dump(data_t data);

#ifdef __cplusplus
}
#endif


#endif

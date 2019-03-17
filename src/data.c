#include <stdlib.h>
#include <stdio.h>
#include "data.h"

data_t data_t_construct(int data_size, unsigned int defolt_size){
	data_t data;
	data.data = (int**)calloc(data_size, sizeof(int*));
	data.alloc_len = (int*)malloc(data_size * sizeof(int));
	data.busy_len = (int*)calloc(data_size, sizeof(int));
	data.data_cnt = data_size;
	int i;
	for (i = 0; i < data_size; i++){
		data.data[i] = (int*)malloc(defolt_size * sizeof(int));
		data.alloc_len[i] = defolt_size;
	}
	
	return data;
}

void data_t_destruct(data_t data){
	int data_cnt = data.data_cnt;
	for (int i = 0; i < data_cnt; i++) 
		if (data.data[i] != NULL) free(data.data[i]);
	free(data.data);
	free(data.alloc_len);
	free(data.busy_len);
	data.data_cnt = -1;
}

void add_elem_to_data_ext(data_t data, int elem_id, int coord, int addition_size){
	//data_t data = box.data;
	if (data.busy_len[coord] >= data.alloc_len[coord]){
		data.alloc_len[coord] += addition_size;
		data.data[coord] = (int*) realloc(data.data[coord], data.alloc_len[coord] * sizeof(int));
	}
	
	data.data[coord][data.busy_len[coord]++] = elem_id;
}

void add_elem_to_data(data_t data, int elem_id, int coord){
	add_elem_to_data_ext(data, elem_id, coord, 3);
}

void data_t_dump(data_t data){
	printf("data:{\n");
	printf("data_cnt = %d\n", data.data_cnt);
	printf("alloc_len: %p\nbusy_len: %p\ndata: %p\n", data.alloc_len, data.busy_len, data.data);
	printf("}\n");
}

void data_t_busy_len_dump(data_t data){
	for (int i = 0; i < data.data_cnt; i++)
		printf("%d ", data.busy_len[i]);
}

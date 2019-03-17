#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include "data.h"
#include "format_in.h"
//############formatted_input###########################################
unsigned int get_state_from_input(unsigned int state){
	int res = state;
	if (state == 2) res = 0;
	if (state == 3) res = 2;
	if (state == 4) res = 3;
	if (state == 1) res = 1;
	return res;
}

net_t get_net_from_file(FILE* input){
	unsigned int vrtx_cnt = 0, elems_cnt = 0, i;
	if (fscanf (input, "%u %u", &vrtx_cnt, &elems_cnt) == EOF)
        perror("Can't read vertex count and elems count "), exit(-1);
	vrtx_t vrtx = vrtx_t_construct(vrtx_cnt);
	for (i = 0; i < vrtx_cnt; i++){
		double x, y, z, h = 0.5;
		unsigned int state = 0;
		if (fscanf (input, "%lg %lg %lg %u", &x, &y, &z, &state)== EOF)
            perror("Can't read node data "), exit(-1);
		state = get_state_from_input(state);
		vrtx.nodes[i] = node_t_construct(x, y, z, h, state, i);
	}
	elems_t elems = elems_t_construct(elems_cnt);
	for (i = 0; i < elems_cnt; i++){
		int x, y, z;
		double n1, n2, n3;
		if (fscanf (input, "%d %d %d %lg %lg %lg", &x, &y, &z, &n1, &n2, &n3)== EOF)
            perror("Can't read element data "), exit(-1);
		elems.elems[i] = elem_t_construct(vrtx.nodes[x - 1], vrtx.nodes[y - 1], vrtx.nodes[z - 1], i);
		//(*((*elems).elems[i])).id = i;
		if (DOT((*(elems).elems[i]).or_area, *(point_t_construct(n1, n2, n3))) < 0) {
			(*(elems).elems[i]).coef = -1;
		}
	}
	springs_t sprs = springs_t_construct(3 * elems_cnt);
	net_t net = net_t_get(vrtx, elems, sprs);

	return net;
}

nets_t formated_in(char* file_name){
	FILE* input = fopen (file_name, "r");
	unsigned int count_nets = 0, i;
	if (fscanf (input, "%u", &count_nets) == EOF)
         perror("Can't read count of nets "), exit(-1);;
	nets_t nets = nets_t_get_net(count_nets);
	for (i = 0; i < count_nets; i++){
		nets.nets[i] = get_net_from_file(input);
		net_t_set_springs(&nets.nets[i]);
	}

	fclose(input);
	init_nets(nets);

	return nets;
}

int get_f_point_t_id(f_point_t point, f_point_t* p, int p_len){
	int id = -1;
	for (int i = 0; i < p_len; i++)
		if (f_point_t_eq(p[i], point)) return i;
	return id;
}

#define READDATA(TYPE, NAME)				\
NAME = *((TYPE*) &buffer[offset]);			\
	offset += sizeof(TYPE)

net_t read_net_from_stl_buf(char* buffer, int len){
	const int header_size = 80;
	buffer = buffer + header_size; //ignore header
	int offset = 0;
	uint32_t e_cnt;
	READDATA(uint32_t, e_cnt);
	printf("e_cnt = %u\n", e_cnt);
	f_point_t* pelems = (f_point_t*)calloc(e_cnt, sizeof(f_point_t));
	f_point_t* pvrtx = (f_point_t*)calloc(3 * e_cnt, sizeof(f_point_t));
	for (unsigned int i = 0; i < e_cnt; i++){
		READDATA(f_point_t, pelems[i]);
		READDATA(f_point_t, pvrtx[3 * i]);
		READDATA(f_point_t, pvrtx[3 * i + 1]);
		READDATA(f_point_t, pvrtx[3 * i + 2]);
		const int attr_byte_cnt_size = 2;
		offset += attr_byte_cnt_size;
	}

	printf("First\n");
	f_point_t* nvrtx = (f_point_t*)calloc(3 * e_cnt, sizeof(f_point_t));
	data_t data = data_t_construct(e_cnt, 7);
	int nv_len = 0;
	for (unsigned int i = 0; i < 3 * e_cnt; i++){
		int flag = 0;
		for (int j = 0; j < nv_len && !flag; j++){
			flag += f_point_t_eq(pvrtx[i], nvrtx[j]);
			if (flag) add_elem_to_data(data, j, i / 3);
		}
		if (flag) continue;
		add_elem_to_data(data, nv_len, i / 3);
		nvrtx[nv_len++] = pvrtx[i];
		if (nv_len % 1000 == 0) printf("nv_len = %d, i = %d\n", nv_len, i);
	}

	printf("Second\n");
	vrtx_t vrtx = vrtx_t_construct(nv_len);
	for (int i = 0; i < nv_len; i++)
		vrtx.nodes[i] = node_t_construct(nvrtx[i].coord[0], nvrtx[i].coord[1], nvrtx[i].coord[2], 0.5, 1, i);
	printf("Three\n");
	elems_t elems = elems_t_construct(e_cnt);
	for (unsigned int i = 0; i < e_cnt; i++){
		node_t* node[3];
		node[0] = vrtx.nodes[data.data[i][0]];
		node[1] = vrtx.nodes[data.data[i][1]];
		node[2] = vrtx.nodes[data.data[i][2]];
		elems.elems[i] = elem_t_construct(node[0], node[1], node[2], i);
	}
	printf("Four\n");
	springs_t sprs = springs_t_construct(3 * e_cnt);

	net_t net = net_t_get(vrtx, elems, sprs);

	data_t_destruct(data);
	free(nvrtx);
	free(pvrtx);
	free(pelems);

	return net;
}
#undef READDATA

int file_len (FILE* file){
    fseek (file, 0, SEEK_END);
	int len = ftell (file);
	rewind (file);

	if (!len) return 0;

	return len;
}

net_t read_net_from_stl(char* file_name){
	if (strcmp(file_name + strlen(file_name) - strlen(".stl"), ".stl")){
		perror("Unknown file type\n");
		net_t net = {0};
		return net;
	}

	FILE* fp = fopen (file_name, "rb");
	int len = file_len(fp);
	char* buffer = (char*)calloc(len, sizeof(char));
	if (fread (buffer, len, sizeof (char), fp) < len)
        perror("Can't read stl file into buffer "), exit(-1);
	fclose(fp);

	net_t net = read_net_from_stl_buf(buffer, len);
	net_t_set_springs(&net);
	init_net(&net);

	free(buffer);

	return net;
}

bnds_t read_bnds(char* f_name){
	if (strcmp(f_name + strlen(f_name) - strlen(".bnd"), ".bnd")){
		perror("Unknown file type\n");
		bnds_t bnds = {0};
		return bnds;
	}

	FILE* input = fopen (f_name, "r");
	int cnt_bnds = 0;
	if (fscanf (input, "%u", &cnt_bnds) == EOF)
        perror("Can't read count of boundary nodes "), exit(-1);
	bnds_t bnds = {cnt_bnds, (line_t*)calloc(cnt_bnds, sizeof(line_t))};
	for (int i = 0; i < cnt_bnds; i++){
		int pnt_cnt = 0;
		if (fscanf (input, "%u", &pnt_cnt) == EOF)
            perror("Can't read count of points in .bnd "), exit(-1);
		line_t bnd = {pnt_cnt, (point_t*)calloc(pnt_cnt, sizeof(point_t))};
		for (int j = 0; j < pnt_cnt; j++){
			if (fscanf (input, "%lg %lg %lg", &bnd.line[j].coord[0], &bnd.line[j].coord[1], &bnd.line[j].coord[2]) == EOF)
                perror("Can't read point in .bnd "), exit(-1);
		}
		bnds.bnds[i] = bnd;
	}

	fclose(input);

	return bnds;
}
//############formatted_input###########################################


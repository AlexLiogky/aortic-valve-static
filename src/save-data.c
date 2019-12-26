#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "save-data.h"

//############save_coord################################################
void save_net_coord(FILE* fp, net_t net){
	int cnt = net.vrtx.count, i;
	for (i = 0; i < cnt; i++){
	point_t point = net.vrtx.nodes[i][0].coord;
	fprintf(fp, "%lg %lg %lg\n", point.coord[0], point.coord[1], point.coord[2]);
	}
}

void save_coord(nets_t nets, char* file_name){
	FILE* fp = fopen (file_name, "w");
	unsigned int i;
	for (i = 0; i < nets.count; i++)
		save_net_coord(fp, nets.nets[i]);

	fclose(fp);
	return;

}
//############save_coord################################################
//############plot_output###############################################
void point_t_plot(point_t point, FILE* fp){
	fprintf(fp, "%.6f %.6f %.6f", point.coord[0], point.coord[1], point.coord[2]);
}


void plot_net(FILE* fp, net_t net){
	int cnt = net.springs.count, i;
	for (i = 0; i < cnt; i++){
		spring_t spring = *(net.springs.springs[i]);
		point_t_plot((*(spring.ends[0])).coord, fp);
		fprintf(fp, " ");
		point_t_plot((*(spring.ends[1])).coord, fp);
		fprintf(fp, " ");
		int coopt = ((*(spring.ends[0])).coaptative > -1) * ((*(spring.ends[1])).coaptative > -1);
		fprintf(fp, "%d ", (coopt > 0));
		fprintf(fp, "\n");
	}
}

void plot(nets_t nets, char* file_name){
	FILE* fp = fopen (file_name, "w");

	unsigned int i;
	for (i = 0; i < nets.count; i++)
		plot_net(fp, nets.nets[i]);

	fclose(fp);
	return;
}
//############plot_output###############################################
void print_stl_header(int fd){
	int header_size = 80;
	char buffer[80] = {0};
	write(fd, buffer, header_size);
}

void print_stl_vec(int fp, f_point_t point){
	write(fp, (char*)point.coord, DIM * sizeof(float));
}

void print_stl_block(int fp, elem_t* elem){
	point_t norm = (*elem).or_area;
	point_t_coef_mul(1 / point_t_length(norm), &norm);
	f_point_t normal = f_point_t_convert(norm);
	f_point_t vrts[3];
	int i;
	for (i = 0; i < 3; i++)
		vrts[i] = f_point_t_convert((*elem).vrts[i][0].coord);
	print_stl_vec(fp, normal);
	if ((*elem).coef < 0)
		for (i = 0; i < 3; i++)
			print_stl_vec(fp, vrts[i]);
	else
	for (i = 2; i >= 0; i--)
			print_stl_vec(fp, vrts[i]);
	uint16_t attr = 0;
	write(fp, (char*)&attr, sizeof(uint16_t));
}

void to_stl1(net_t net, const char* file_name){
    char max_file_name[1024] = {};
    strcat(max_file_name, file_name);
    strcat(max_file_name + strlen(file_name), ".stl");
    int fd = open (max_file_name, O_TRUNC | O_CREAT | O_RDWR, S_IRWXU);
    print_stl_header(fd);
    unsigned int i, cnt = net.elems.count;
    write(fd, (char*)&cnt, sizeof(unsigned int));
    int elem_cnt = net.elems.count, j;
    for (j = 0; j < elem_cnt; j++)
        print_stl_block(fd, net.elems.elems[j]);

    close(fd);
    return;
}

void to_stl(nets_t nets, const char* file_name){
	char max_file_name[1024] = {};
	strcat(max_file_name, file_name);
	strcat(max_file_name + strlen(file_name), ".stl");
	int fd = open (max_file_name, O_TRUNC | O_CREAT | O_RDWR, S_IRWXU);
	print_stl_header(fd);
	unsigned int i, cnt = 0;
	for (i = 0; i < nets.count; i++)
		cnt += nets.nets[i].elems.count;
	write(fd, (char*)&cnt, sizeof(unsigned int));
	for (i = 0; i < nets.count; i++){
		int elem_cnt = nets.nets[i].elems.count, j;
		for (j = 0; j < elem_cnt; j++)
			print_stl_block(fd, nets.nets[i].elems.elems[j]);
	}

	close(fd);
	return;
}

//############save_nets#################################################
void save_nets_to_file(nets_t nets, const char* file_name){
	char max_file_name[1024] = {};
	strcat(max_file_name, file_name);
	strcat(max_file_name + strlen(file_name), ".nts");
	int fd = open (max_file_name, O_TRUNC | O_CREAT | O_RDWR, S_IRWXU);
	int buflen = 0;
	char* buf_nets = serialize_nets(nets, &buflen);
	write(fd, buf_nets, buflen);

	free(buf_nets);
	close(fd);
	return;
}

int _file_len (FILE* file){
    fseek (file, 0, SEEK_END);
	int len = ftell (file);
	rewind (file);

	if (!len) return 0;

	return len;
}

nets_t download_nets_from_file(const char* file_name){
	if (strcmp(file_name + strlen(file_name) - strlen(".nts"), ".nts")){
		perror("Unknown file type ");
		return nets_t_get_net(0);
	}

	FILE* fp = fopen (file_name, "rb");
	if (!fp)
		return perror("Can't open file "), nets_t_get_net(0);
	int len = _file_len(fp);
	char* buffer = (char*)calloc(len, sizeof(char));
	fread (buffer, len, sizeof (char), fp);

	nets_t nets = deserialize_nets(buffer);

	free(buffer);
	fclose(fp);

	return nets;
}


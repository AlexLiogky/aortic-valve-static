#ifndef _FORMAT_IN_H
#define _FORMAT_IN_H

#include "nets.h"
#include <stdio.h>

//get single net from txt file where single net saved in specific form
net_t get_net_from_file(FILE* input);

//get nets from txt file where nets saved in specific format
nets_t formated_in(const char* file_name);

//read net from ".stl" file
net_t read_net_from_stl(char* file_name);

typedef struct bnds_t{
	int cnt;
	line_t* bnds;
}bnds_t;

//read bound range from file where it saved in special format
//ATTENTION: file must have type - ".bnd"
bnds_t read_bnds(char* f_name);
//############formatted_input###########################################

#endif

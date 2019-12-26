#ifndef _FORMAT_OUT_H
#define _FORMAT_OUT_H

#include "nets.h"

#ifdef __cplusplus
extern "C" {
#endif

//save list of coords into file "file_name"
void save_coord(nets_t nets, char* file_name);

//save data in format for plotting
void plot(nets_t nets, char* file_name);

//save nets in stl format into file "file_name"
void to_stl1(net_t net, const char* file_name);
void to_stl(nets_t nets, const char* file_name);

//save serialized nets container into file "file_name"
//add to file type ".nts"
void save_nets_to_file(nets_t nets, const char* file_name);

//download serialized nets container from file "file_name"
//ATTENTION: file must have type - ".nts"
nets_t download_nets_from_file(const char* file_name);

#ifdef __cplusplus
}
#endif


#endif

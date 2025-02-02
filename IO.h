#ifndef IO_H
#define IO_H

#include "classes.h"

Params read_param(std::string fname);
std::vector<Cell> read_ic(Params& p);
void write_output(std::vector<Cell> c, Params p, Vars& v);

#endif // IO_H
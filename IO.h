#ifndef IO_H
#define IO_H

#include "classes.h"

Params read_param(std::string fname);
std::vector<Cell> read_ic(Params& p, Grid& g);
void write_output(const std::vector<Cell> &c, Params p, Vars& v, const Grid& g);
std::vector<Cell> read_ic_hdf5(Params& p, Grid& g);
void write_output_hdf5(const std::vector<Cell>& c, Params& p, Vars& v, const Grid& g);

#endif // IO_H
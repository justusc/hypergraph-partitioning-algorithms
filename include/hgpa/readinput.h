#ifndef HGPA_READINPUT_INCLUDED
#define HGPA_READINPUT_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */

#include "hgpa/defs.h"

/* read input hypergraph and construct cells, nets, cnets, & ncells arrays */
void read_hgraph(char fname[],
                 hypergraph_t* hgraph);

#endif

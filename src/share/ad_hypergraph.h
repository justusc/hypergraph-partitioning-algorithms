#ifndef AD_HYPERGRAPH_H__INCLUDED
#define AD_HYPERGRAPH_H__INCLUDED

#include "ad_defs.h"

/* initialize cells array */
void init_cells(int nocells, cells_t cells[]);

/* initialize nets array */
void init_nets(int nonets,
               nets_t nets[]);

/* initialize netlist pointers */
void init_netlist(hypergraph_t* hgraph);

/* Allocate hypergraph memory and initialize variables */
void init_hypergraph(int nocells,
                     int nonets,
                     int nopins,
                     hypergraph_t* hgraph);

/* Free hypergraph memory */
void free_hypergraph(hypergraph_t* hgraph);

#endif /* AD_HYPERGRAPH_H__INCLUDED */

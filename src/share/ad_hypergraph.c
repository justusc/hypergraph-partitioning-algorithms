#include "ad_hypergraph.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/* initialize cells array */
void init_cells(int nocells, cells_t cells[])
{
    for (int i = 0; i < nocells; i++) {
        cells[i].cno_nets = 0;
        cells[i].cno_inets = 0;
        cells[i].netlist = NIL;
    }   /* for */
}   /* init_cells */

/* initialize nets array */
void init_nets(int nonets,
               nets_t nets[])
{
    for (int i = 0; i < nonets; i++) {
        nets[i].nno_cells = 0;
        nets[i].celllist = NIL;
        nets[i].npartdeg = NULL;
    }   /* for i */
}   /* init_nets */

/* initialize netlist pointers */
void init_netlist(hypergraph_t* hgraph)
{
    int nonets     = hgraph->nonets;
    cells_t* cells = hgraph->cells;
    nets_t* nets   = hgraph->nets;
    corn_t* cnets  = hgraph->cnets;
    corn_t* ncells = hgraph->ncells;

    for (int i = 0; i < nonets; i++) {
        for (int j = 0; j < nets[i].nno_cells; j++) {
            int nets_i_celllist_j = nets[i].celllist + j;
            int cell_no = ncells[nets_i_celllist_j].corn_no;
            int cnets_inx = cells[cell_no].netlist + cells[cell_no].cno_inets;
            cnets[cnets_inx].corn_no = i;
            cells[cell_no].cno_inets++;
            if (cells[cell_no].cno_inets > cells[cell_no].cno_nets) {
                printf("Error: Inconsistency in cell_%d degrees.\n", j);
                exit(1);
            }   /* if */
        }   /* for j */
    }   /* for i */
}   /* init_netlist */

/* Allocate hypergraph memory and initialize variables */
void init_hypergraph(int nocells,
                     int nonets,
                     int nopins,
                     hypergraph_t* hgraph)
{
    cells_t *cells = (cells_t *) calloc(nocells, sizeof(cells_t));
    assert(cells != NULL);

    nets_t *nets = (nets_t *) calloc(nonets, sizeof(nets_t));
    assert(nets != NULL);

    /* cells of nets */
    corn_t *cnets = (corn_t *) calloc(nopins, sizeof(corn_t));
    assert(cnets != NULL);
    /* nets of cells */
    corn_t *ncells = (corn_t *) calloc(nopins, sizeof(corn_t));
    assert(ncells != NULL);

    /* initialize cells & nets arrays */
    init_cells(nocells, cells);
    init_nets(nonets, nets);

    /* initialize the hgraph data */
    hgraph->nocells     = nocells;
    hgraph->nonets      = nonets;
    hgraph->nopins      = nopins;
    hgraph->totcellsize = 0;
    hgraph->totnetsize  = 0;
    hgraph->max_cdeg    = -1;
    hgraph->max_ndeg    = -1;
    hgraph->max_cweight = -1;
    hgraph->max_nweight = -1;
    hgraph->cells       = cells;
    hgraph->nets        = nets;
    hgraph->ncells      = ncells;
    hgraph->cnets       = cnets;
}

/* Free hypergraph memory */
void free_hypergraph(hypergraph_t* hgraph) {
    /* free memory for all data structures */
    free(hgraph->cnets);
    free(hgraph->ncells);

    for (int i = 0; i < hgraph->nonets; i++) {
        free(hgraph->nets[i].npartdeg);
    }
    free(hgraph->nets);

    free(hgraph->cells);
}


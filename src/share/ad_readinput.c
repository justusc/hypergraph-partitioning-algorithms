
/* COPYRIGHT C 1991- Ali Dasdan */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include "hpga/ad_defs.h"
#include "hpga/ad_fileio.h"
#include "hpga/ad_hypergraph.h"
#include "hpga/ad_readinput.h"


/* read input hypergraph and construct cells, nets, cnets, & ncells arrays */
void read_hgraph(char fname[],
                 hypergraph_t* hgraph)
{
    int nocells;           /* number of cells */
    int nonets;            /* number of nets */
    int nopins;            /* number of pins */
    int totcellsize;       /* total cell weight of the partition */
    int totnetsize;        /* total net weight of the partition */
    int max_cdeg;          /* max density of a cell */
    int max_ndeg;          /* max density of a net */
    int max_cweight;       /* max cell weight */
    int max_nweight;       /* max net weight */

    FILE *fp;
    open_file(&fp, fname, "r");

    /* Read hypergraph size */
    if (fscanf(fp, "%d%d%d", &nocells, &nonets, &nopins) == EOF) {
        printf("Error: Cannot read from %s: errno= %d error= %s\n", fname, errno, strerror(errno));
        close_file(&fp);
        exit(1);
    }

    if ((nocells < 0) || (nonets < 0) || (nopins < 0)) {
        printf("Error: Invalid attributes of graph.\n");
        close_file(&fp);
        exit(1);
    }   /* if */

    init_hypergraph(nocells, nonets, nopins, hgraph);

    /* initialize local variables */
    totcellsize = hgraph->totcellsize;
    totnetsize  = hgraph->totnetsize;
    max_cdeg    = hgraph->max_cdeg;
    max_ndeg    = hgraph->max_ndeg;
    max_cweight = hgraph->max_cweight;
    max_nweight = hgraph->max_nweight;
    cells_t * restrict cells = hgraph->cells;
    nets_t * restrict nets   = hgraph->nets;
    corn_t * restrict ncells = hgraph->ncells;
    corn_t * restrict cnets  = hgraph->cnets;


    /* read nets */
    int ncells_inx = 0;
    for (int i = 0; i < nonets; i++) {
        if (fscanf(fp, "%d%d", &nets[i].nweight, &nets[i].nno_cells) == EOF) {
            printf("Error: Cannot read from %s: errno= %d error= %s\n", fname, errno, strerror(errno));
            close_file(&fp);
            exit(1);
        }
        totnetsize += nets[i].nweight;
        if (nets[i].nweight > max_nweight) {
            max_nweight = nets[i].nweight;
        }
        if (nets[i].nno_cells > max_ndeg) {
            max_ndeg = nets[i].nno_cells;
        }
        nets[i].celllist = ncells_inx;

        for (int j = 0; j < nets[i].nno_cells; j++) {
            int cell_no;
            if (fscanf(fp, "%d", &cell_no) == EOF) {
                printf("Error: Cannot read from %s: errno= %d error= %s\n", fname, errno, strerror(errno));
                close_file(&fp);
                exit(1);
            }
            ncells[ncells_inx].corn_no = cell_no;
            ncells_inx++;
            cells[cell_no].cno_nets++;
        }   /* for j */
    }   /* for i */

    /* read  cell weights */
    int cnets_inx = 0;
    for (int i = 0; i < nocells; i++) {
        if (fscanf(fp, "%d", &cells[i].cweight) == EOF) {
            printf("Error: Cannot read from %s: errno= %d error= %s\n", fname, errno, strerror(errno));
            close_file(&fp);
            exit(1);
        }
        totcellsize += cells[i].cweight;
        if (cells[i].cweight > max_cweight)
            max_cweight = cells[i].cweight;
        if (cells[i].cno_nets > max_cdeg)
            max_cdeg = cells[i].cno_nets;
        cells[i].netlist = cnets_inx;
        cnets_inx += cells[i].cno_nets;
    }   /* for i */

    close_file(&fp);

    /* create netlists */
    init_netlist(hgraph);

    /* initialize the hgraph data */
    hgraph->totcellsize = totcellsize;
    hgraph->totnetsize  = totnetsize;
    hgraph->max_cdeg    = max_cdeg;
    hgraph->max_ndeg    = max_ndeg;
    hgraph->max_cweight = max_cweight;
    hgraph->max_nweight = max_nweight;
}   /* read_hgraph */

/* EOF */

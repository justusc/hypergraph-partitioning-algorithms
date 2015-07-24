
/* COPYRIGHT C 1991- Ali Dasdan */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "ad_defs.h"
#include "ad_random.h"
#include "ad_fileio.h"
#include "ad_readinput.h"
#include "ad_partition.h"
#include "ad_print.h"
#include "ad_bucketio.h"
#include "ad_lib.h"
#include "ad_lib_fms.h"

/* FOR SANCHIS' VERSION OF MULTI-WAY PARTITIONING */
/* Also mentioned as the SN algorithm */
/* Direct multi-way partitioning.
   Locking is used.
   Cells are moved wrt their gains
*/

int main(int argc, char *argv[])
{
    /* the hypergraph to partition */
    hypergraph_t hgraph;

    /* partitioning variables */
    int noparts;           /* number of partitions */
    int cutsize;           /* cutsize of the partition */
    int max_gain;          /* max gain of a cell */

    /* FMS algorithm variables */
    int bucketsize;        /* max size of a bucket array */
    int msize;             /* index to mcells */

    if (argc < 3) {
        printf("\nUsage: %s InputFileName NoParts [Seed]\n", argv[0]);
        exit(1);
    }  /* if */

    char fname[STR_SIZE];
    sprintf(fname, "%s", argv[1]);

    noparts = atoi(argv[2]);

    long seed;
    if (argc > 3) {
        seed = (long) atoi(argv[3]);
    } else {
        seed = -1;
    }
    seed = randomize((long) seed);
    printf("SEED = %ld fname = %s\n", seed, fname);

    read_hgraph(fname, &hgraph);

    cells_info_t *cells_info = (cells_info_t *) calloc(hgraph.nocells, sizeof(cells_info_t));
    assert(cells_info != NULL);
    for (int i = 0; i < hgraph.nocells; i++) {
        cells_info[i].mgain = (int *) calloc(noparts, sizeof(int));
        assert(cells_info[i].mgain != NULL);
        cells_info[i].partb_ptr = (bnode_ptr_t *) calloc(noparts - 1, sizeof(bnode_ptr_t));
        assert(cells_info[i].partb_ptr != NULL);
        cells_info[i].partb_gain_inx = (int *) calloc(noparts - 1, sizeof(int));
        assert(cells_info[i].partb_gain_inx != NULL);
    }

    nets_info_t *nets_info = (nets_info_t *) calloc(hgraph.nonets, sizeof(nets_info_t));
    assert(nets_info != NULL);
    for (int i = 0; i < hgraph.nonets; i++) {
        nets_info[i].npartdeg = (int *) calloc(noparts, sizeof(int));
        assert(nets_info[i].npartdeg != NULL);
        hgraph.nets[i].npartdeg = (int *) calloc(noparts, sizeof(int));
        assert(hgraph.nets[i].npartdeg != NULL);
        for (int j = 0; j < noparts; j++) {
            hgraph.nets[i].npartdeg[j] = 0;
        }
    }

    /* partition buckets */
    partb_t partb[noparts][noparts - 1];
    parts_info_t parts_info[noparts];

    /* population (w/ one individual!) */
    ind_t pop[MAX_POP];
    for (int i = 0; i < MAX_POP; i++) {
        pop[i].chrom = (allele *) calloc(hgraph.nocells, sizeof(allele));
        assert(pop[i].chrom != NULL);
        pop[i].parts = (parts_t *) calloc(noparts, sizeof(parts_t));
        assert(pop[i].parts != NULL);
    }

    /* selected cell */
    selected_cell_t scell[1];

    /* moved cells */
    mcells_t *mcells = (mcells_t *) calloc(hgraph.nocells, sizeof(mcells_t));
    assert(mcells != NULL);

    max_gain = hgraph.max_cdeg * hgraph.max_nweight;
    bucketsize = 2 * max_gain + 1;

    /* alloc memory (statically if possible) */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < noparts - 1; ++j) {
            partb[i][j].bnode_ptr = (bnode_ptr_t *) calloc(bucketsize, sizeof(bnode_ptr_t));
        }
    }

    float off_ratio = (float) 0.1;
    create_partition(hgraph.nocells, noparts, hgraph.totcellsize, hgraph.max_cweight, &off_ratio,
                     hgraph.cells, hgraph.nets, hgraph.cnets, &pop[0]);

#ifdef DEBUG
    printf("off=%f\n", off_ratio);
    printf("Initial : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("II %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0]. parts[i].pmax_size);
    }
#endif

    init_buckets(noparts, bucketsize, partb);
    cutsize = find_cut_size(hgraph.nonets, noparts, hgraph.totnetsize, hgraph.nets, &pop[0]);

#ifdef DEBUG
    printf("Totalsize = %d Initial cutsize = %d\n", hgraph.totnetsize, cutsize);
#endif

    int gain_sum;
    int no_iter = 0;
    int glob_inx = 0;
    do {

        copy_partition(noparts, parts_info, &pop[0]);

        copy_nets_info(hgraph.nonets, noparts, hgraph.nets, nets_info);

        compute_gains(hgraph.nocells, noparts, hgraph.cells, hgraph.nets, hgraph.cnets,
                      cells_info, pop[0].chrom);

        create_buckets(hgraph.nocells, noparts, max_gain, pop[0].chrom,
                       partb, cells_info);

        msize = 0;

        int nlocked = 0;
        do {
            int move_possible = select_cell(noparts, scell, parts_info, hgraph.cells,
                                            partb, cells_info);

            delete_partb_nodes_of_cell(noparts, scell[0].mov_cell_no,
                                       scell[0].from_part, partb, cells_info);

            /* lock cell */
            cells_info[scell[0].mov_cell_no].locked = True;
            if (move_possible == True) {
                move_cell(mcells, msize, scell);
                msize++;
                update_gains(noparts, max_gain, scell,
                             hgraph.cells, hgraph.nets, hgraph.cnets, hgraph.ncells, nets_info,
                             partb, cells_info, pop[0].chrom);
            }   /* if */
            nlocked++;

        } while (nlocked < hgraph.nocells);

        int max_mcells_inx;
        gain_sum = find_move_set(mcells, msize, &max_mcells_inx);

#ifdef DEBUG
        printf("gain_sum=%d max_mcells_inx=%d msize = %d\n",
               gain_sum, max_mcells_inx, msize);
#endif

        if (gain_sum > 0) {
            int cut_gain = move_cells(False, hgraph.nocells, msize, mcells, max_mcells_inx,
                                      cutsize, &glob_inx, &pop[0], hgraph.cells, hgraph.nets, hgraph.cnets);
            cutsize -= cut_gain;
        }   /* if */
        no_iter++;

#ifdef DEBUG
        printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n",
               no_iter, cutsize, find_cut_size(hgraph.nonets, noparts, hgraph.totnetsize, hgraph.nets, &pop[0]));
#endif

    } while ((gain_sum > 0) && (cutsize > 0) && (no_iter < NO_ITERATIONS));

    printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n", no_iter,
           cutsize, find_cut_size(hgraph.nonets, noparts, hgraph.totnetsize, hgraph.nets, &pop[0]));

    free_nodes(noparts, bucketsize, partb);

#ifdef DEBUG
    printf("Final : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("FF %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif

    /* free memory for all data structures */
    for (int i = 0; i < hgraph.nocells; i++) {
        free(cells_info[i].mgain);
        free(cells_info[i].partb_ptr);
        free(cells_info[i].partb_gain_inx);
    }
    free(cells_info);

    for (int i = 0; i < hgraph.nonets; i++) {
        free(nets_info[i].npartdeg);
    }
    free(nets_info);

    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < noparts - 1; ++j) {
            free(partb[i][j].bnode_ptr);
        }
    }

    for (int i = 0; i < MAX_POP; i++) {
        free(pop[i].chrom);
        free(pop[i].parts);
    }

    free(mcells);

    free_hypergraph(&hgraph);

    return (0);
}   /* main-fms */

/* EOF */


/* COPYRIGHT C 1991- Ali Dasdan */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "hpga/ad_defs.h"
#include "hpga/ad_random.h"
#include "hpga/ad_fileio.h"
#include "hpga/ad_hypergraph.h"
#include "hpga/ad_readinput.h"
#include "hpga/ad_partition.h"
#include "hpga/ad_print.h"
#include "hpga/ad_bucketio.h"
#include "hpga/ad_lib.h"
#include "hpga/ad_lib_pfm.h"

/* PARTITIONING BY FREE MOVES */

int main(int argc, char *argv[])
{
    /* the hypergraph to partition */
    hypergraph_t hgraph;

    /* partitioning variables */
    int noparts;           /* number of partitions */
    int cutsize;           /* cutsize of the partition */
    int max_gain;          /* max gain of a cell */

    if (argc < 5) {
        printf("\nUsage: %s InputFileName NoParts Version BucketSizeFactor [Seed]\n", argv[0]);
        exit(1);
    }  /* if */

    char fname[STR_SIZE];
    sprintf(fname, "%s", argv[1]);

    noparts = atoi(argv[2]);
    int version = atoi(argv[3]);

    int bucketsize_factor = atoi(argv[4]);

    long seed;
    if (argc > 5) {
        seed = (long) atoi(argv[5]);
    } else {
        seed = (long) -1;
    }
    seed = randomize((long)  seed);
    printf("SEED = %ld fname = %s\n", seed, fname);

    /* in read_graph, change tnoparts -> noparts & delete the line (D) */
    read_hgraph(fname, &hgraph);

    /* determine max_noiter based on pfm version */
    /* pfm1: size=max_cells;
       pfm2: size=max_cells * max_parts;
       pfm3: size=max_cells * max_parts^2 */
    int max_noiter = hgraph.nocells;
    switch (version) {
    case 1 : break;
    case 2 : max_noiter *= noparts; break;
    case 3 : max_noiter *= noparts * noparts; break;
    default : break;
    }

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
    selected_cell_t prev_scell[1];

    /* moved cells */
    mcells_t *mcells = (mcells_t *) calloc(max_noiter, sizeof(mcells_t));
    assert(mcells != NULL);

    /* temp chrom */
    allele *tchrom = (allele *) calloc(hgraph.nocells, sizeof(allele));
    assert(tchrom != NULL);

    float off_ratio = (float) 0.1; /* alpha in initial partitioning */
    create_partition(hgraph.nocells, noparts, hgraph.totcellsize, hgraph.max_cweight, &off_ratio,
                     hgraph.cells, hgraph.nets, hgraph.cnets, &pop[0]);

#ifndef NDEBUG
    printf("Initial : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("II %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif

    /* bucketsize has impact on cutsize and runtime */
    float K;
    max_gain = hgraph.max_cdeg * hgraph.max_nweight;
    int bucketsize = 2 * max_gain + 1;
    if (bucketsize_factor > 0) {
        bucketsize *= bucketsize_factor;
    }

    /* cache to speed up math heavy function evals */
    eval_t *eval = (eval_t *) calloc(2 * max_gain + 1, sizeof(eval_t));
    calculate_scale(hgraph.nocells, noparts, max_gain, &K);
    fill_eval(max_gain, K, eval);

    /* alloc memory (statically if possible) */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < noparts - 1; ++j) {
            partb[i][j].bnode_ptr = (bnode_ptr_t *) calloc(bucketsize, sizeof(bnode_ptr_t));
        }
    }

    init_buckets(noparts, bucketsize, partb);
    cutsize = find_cut_size(hgraph.nonets, noparts, hgraph.totnetsize, hgraph.nets, &pop[0]);

#ifndef NDEBUG
    printf("BB bucketsize=%d\n", bucketsize);
    printf("Totalsize = %d Initial cutsize = %d\n", hgraph.totnetsize, cutsize);
#endif

    int gain_sum;
    int glob_inx = 0;
    int pass_no = 0;
    do {

        copy_partition(noparts, parts_info, &pop[0]);

        copy_nets_info(hgraph.nonets, noparts, hgraph.nets, nets_info);

        for (int i = 0; i < hgraph.nocells; i++) {
            tchrom[i] = pop[0].chrom[i];
        }

        compute_gains(hgraph.nocells, noparts, hgraph.cells, hgraph.nets, hgraph.cnets,
                      cells_info, pop[0].chrom);

        create_buckets(bucketsize, hgraph.nocells, noparts, max_gain, eval, pop[0].chrom,
                       partb, cells_info);

        prev_scell[0].mov_cell_no = -1;

        int msize = 0; /* index to mcells */
        int mov_count = 0; /* count of total moves */
        while (mov_count < max_noiter) {

            select_cell(noparts, scell, parts_info, hgraph.cells, partb, cells_info);

            move_cell(mcells, msize, scell, tchrom, cells_info);

            msize++;

            /* insert previous cell */
            create_partb_nodes_of_cell(bucketsize, noparts, max_gain,
                                       prev_scell[0].mov_cell_no,
                                       prev_scell[0].to_part,
                                       eval, partb, cells_info);

            update_gains(bucketsize, noparts, max_gain, eval, scell,
                         hgraph.cells, hgraph.nets, hgraph.cnets, hgraph.ncells, nets_info,
                         partb, cells_info, tchrom);

            /* currently moved cell becomes previous cell */
            bcopy((void *) scell, (void *) prev_scell, sizeof(scell[0]));

            /* delete current cell */
            delete_partb_nodes_of_cell(noparts, scell[0].mov_cell_no,
                                       scell[0].from_part, partb, cells_info);

            mov_count++;
        }   /* while */

        int max_mcells_inx;
        gain_sum = find_move_set(mcells, msize, &max_mcells_inx);

#ifndef NDEBUG
        printf("gain_sum=%d max_mcells_inx=%d msize = %d\n",
               gain_sum, max_mcells_inx, msize);
#endif

        if (gain_sum > 0) {
            int cut_gain = move_cells(False, hgraph.nocells, msize, mcells, max_mcells_inx,
                                      cutsize, &glob_inx, &pop[0], hgraph.cells, hgraph.nets, hgraph.cnets);
            cutsize -= cut_gain;
        }   /* if */
        pass_no++;

#ifndef NDEBUG
        printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n",
               pass_no, cutsize,
               find_cut_size (hgraph.nonets, noparts, hgraph.totnetsize, hgraph.nets, &pop[0]));
#endif

        free_nodes(noparts, bucketsize, partb);

    }  while ((gain_sum > 0) && (cutsize > 0)); /* while pass_no */
    /* or (pass_no < noparts) */

    printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n", pass_no,
           cutsize, find_cut_size(hgraph.nonets, noparts, hgraph.totnetsize, hgraph.nets, &pop[0]));

#ifndef NDEBUG
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

    free(tchrom);

    free(eval);

    free_hypergraph(&hgraph);

    return (0);
}  /* main_pfm */

/* EOF */

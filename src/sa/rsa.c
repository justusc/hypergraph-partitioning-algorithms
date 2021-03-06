
/* COPYRIGHT C 1991- Ali Dasdan */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "hgpa/defs.h"
#include "hgpa/fileio.h"
#include "hgpa/hgpa.h"
#include "hgpa/hypergraph.h"
#include "hgpa/partition.h"
#include "hgpa/readinput.h"
#include "hgpa/sa.h"
#include "hgpa/random.h"

/* Simulated Annealing for Multiple-way Hypergraph Partitioning -
   based on the paper: Johnson et al., Optimization by Simulated
   Annealing, Operations Research, 37(6), 1989. */

/* Uses the ratio cut metric (cutsize divided by the product of the
   part sizes) to relax the balance constraint to a softer
   constraint. This implementation works for 2-way partitioning
   only. */

int main(int argc, char *argv[])
{
    /* the hypergraph to partition */
    hypergraph_t hgraph;   /* hypergraph structure */

    /* partitioning variables */
    int noparts;           /* number of partitions */
    int cutsize;           /* cutsize of the partition */
    int max_gain;          /* max gain of a cell */

    /* SA algorithm variables */
    float temperature;     /* temperature in SA alg. */
    float tempfactor;      /* cooling ratio */
    float minpercent;      /* percentage of accepted moves */
    float cutoff;          /* a speedup option - see the paper */
    int delta;             /* cost difference between the last two states */
    int freezelimit;       /* when freezecount reaches this limit, the system is frozen */
    int freezecount;
    int nochanges;         /* number of accepted states */
    int changeslimit;      /* limit on max number of accepted changes */
    int notrials;          /* number of tried states */
    int trialslimit;       /* limit on max number of tried states */
    int sizefactor;        /* propartionality constant for temperature length */
    int neigh_size;        /* neighborhood size */
    int selected;          /* set if a feasible move is found */
    int changed;           /* set if a change occurs in best_cutsize */
    int same;              /* count the number of times cutsize remains the same */
    int samecount;         /* limit on the counter "same" */
    int prev_cutsize;      /* previous cutsize value - used with "same" */
    int pass_no;           /* pass number */

    /* additional variables for ratio cut partitioning */
    float ratiosize;       /* current ratio size */
    float rdelta;          /* decrease in ratiosize due to a move */
    float bprod;           /* product of sizes of from and to parts before the move */
    float aprod;           /* product of sizes of from and to parts after the move */
    float prev_ratiosize;  /* previous ratiosize */
    float best_ratiosize;  /* best ratiosize so far */

    if (argc < 3) {
        printf("\nUsage: %s InputFileName NoParts [Seed]\n", argv[0]);
        exit(1);
    }  /* if */

    char fname[STR_SIZE];
    sprintf(fname, "%s", argv[1]);

    noparts = atoi(argv[2]);
    if (noparts > 2) {
        printf("Warning: SA with ratio cut works for 2-way partitioning only.\n");
        noparts = 2;
    }

    long seed;
    if (argc > 3) {
        seed = (long) atoi(argv[3]);
    } else {
        seed = -1;
    }
    seed = randomize((long) seed);
    printf("SEED = %ld fname = %s\n", seed, fname);

    read_hgraph(fname, &hgraph);

    /* alloc memory for all data structures */
    cells_info_t *cells_info = (cells_info_t *) calloc(hgraph.nocells, sizeof(cells_info_t));
    assert(cells_info != NULL);
    for (int i = 0; i < hgraph.nocells; i++) {
        cells_info[i].mgain = (int *) calloc(noparts, sizeof(int));
        assert(cells_info[i].mgain != NULL);
        cells_info[i].partb_ptr = NULL;
        cells_info[i].partb_gain_inx = NULL;
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

    /* properties of the best solution */
    int best_cutsize;

    nets_t *best_nets = (nets_t *) calloc(hgraph.nonets, sizeof(nets_t));
    assert(best_nets != NULL);
    for (int i = 0; i < hgraph.nonets; i++) {
        best_nets[i].npartdeg = (int *) calloc(noparts, sizeof(int));
        assert(best_nets[i].npartdeg != NULL);
    }

    ind_t best_pop[MAX_POP];
    for (int i = 0; i < MAX_POP; i++) {
        best_pop[i].chrom = (allele *) calloc(hgraph.nocells, sizeof(allele));
        assert(best_pop[i].chrom != NULL);
        best_pop[i].parts = (parts_t *) calloc(noparts, sizeof(parts_t));
        assert(best_pop[i].parts != NULL);
    }

    /* create the initial partition or solution */
    float off_ratio = (float) 0.05;   /* alpha in initial partitioning */
    create_partition(hgraph.nocells, noparts, hgraph.totcellsize, hgraph.max_cweight, &off_ratio,
                     hgraph.cells, hgraph.nets, hgraph.cnets, &pop[0]);

#ifdef AD_TRACE
    printf("off=%f\n", off_ratio);
    printf("Initial : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("II %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif

    max_gain = hgraph.max_cdeg * hgraph.max_nweight;

    /* find initial cutsize */
    cutsize = find_cut_size(hgraph.nonets, noparts, hgraph.totnetsize, hgraph.nets, &pop[0]);

    /* find initial ratiosize */
    ratiosize = (float) cutsize;
    for (int i = 0; i < noparts; i++) {
        ratiosize /= (float) pop[0].parts[i].pcurr_size;
    }

#ifdef AD_TRACE
    printf("Totalsize = %d Initial cutsize = %d\n", hgraph.totnetsize, cutsize);
#endif

    /* compute costs (gains) of cells in hypergraph */
    compute_gains(hgraph.nocells, noparts, hgraph.cells, hgraph.nets, hgraph.cnets, cells_info, pop[0].chrom);

    /* initialize the current champion solution */
    best_ratiosize = ratiosize;
    best_cutsize = cutsize;
    copy_pop(hgraph.nocells, noparts, pop, best_pop);
    copy_nets(hgraph.nonets, noparts, hgraph.nets, best_nets);

    /* initialize variables of SA algorithm */
    changed = False;
    neigh_size = hgraph.nocells * (noparts - 1);
    sizefactor = 16;      /* 16 in paper */
    minpercent = 0.005;   /* 0.02 in paper */
    cutoff = 1.0;         /* 1.0 in paper */
    temperature = 10.0;   /* 0.6 - not given in paper */
    tempfactor = 0.95;    /* 0.95 in paper */
    freezelimit = 5;      /* 5 in paper */
    freezecount = 0;
    pass_no = 0;
    same = 0;
    samecount = hgraph.nocells;
    prev_ratiosize = -1.0;
    prev_cutsize = -1;
    trialslimit = sizefactor * neigh_size;
    changeslimit = (int) (cutoff * sizefactor * neigh_size);

    /* while not yet frozen */
    while (freezecount < freezelimit)  {

        /* perform the following exploration trialslimit times */
        nochanges = notrials = 0;
        while ((notrials < trialslimit) && (nochanges < changeslimit)) {
            notrials++;

            /* randomly select a cell to move to a randomly selected part */
            selected = select_cell(hgraph.nocells, noparts, scell, pop[0].chrom,
                                   hgraph.cells, pop[0].parts, cells_info);
            if (! selected) {
                printf("Error: Cannot find a move to select.\n");
                exit(1);
            }   /* if */

            /* delta is the change in the cutsize based cost function */
            delta = -scell[0].mov_gain;

            /* rdelta is the change in the ratiocut based cost function */
            int from_part_size = pop[0].parts[scell[0].from_part].pcurr_size;
            int to_part_size = pop[0].parts[scell[0].to_part].pcurr_size;
            bprod = (float) (from_part_size * to_part_size);
            aprod = bprod + (float) hgraph.cells[scell[0].from_part].cweight *
                (float) (from_part_size - to_part_size);
            rdelta = -(ratiosize * ((float) 1.0 - (bprod / aprod)) +
                       (float) scell[0].mov_gain / aprod);

            /* if delta is negative, this change is a downhill move so
               accept the new solution */
            /* if delta is positive, this change is an uphill move so
               accept with an ever decreasing probability that depends
               on delta and the current temperature */
            if (((float) rdelta <= 0.0) ||
                (((float) rdelta > 0.0) &&
                 (rand01() <= (float) exp((double) -rdelta / (double) temperature)))) {

#ifdef AD_TRACE
                printf("cell=%d from=%d to=%d gain=%d ratio=%f\n",
                       scell[0].mov_cell_no, scell[0].from_part,
                       scell[0].to_part, scell[0].mov_gain, rdelta);
#endif

                nochanges++;

                /* move the selected cell */
                move_cell(scell, pop[0].chrom, hgraph.cells, pop[0].parts);

                /* update the costs of the neighbor cells */
                update_gains(scell, hgraph.cells, hgraph.nets, hgraph.cnets, hgraph.ncells,
                             cells_info, pop[0].chrom);
                cutsize += delta;
                ratiosize += rdelta;

                /* update the current champion solution */
                if (ratiosize < best_ratiosize) {
                    changed = True;
                    best_cutsize = cutsize;
                    best_ratiosize = ratiosize;
                    copy_pop(hgraph.nocells, noparts, pop, best_pop);
                    copy_nets(hgraph.nonets, noparts, hgraph.nets, best_nets);
                }   /* if */

            }   /* if the selected move is accepted */
        }   /* while */

#ifdef AD_TRACE
        printf("changes=%d trials=%d accept=%f temperature=%f\n",
               nochanges, notrials,
               (100.0 * (float) nochanges / notrials), temperature);
#endif

        /* reduce the temperature and update the other variables of SA
           algorithm */
        temperature = tempfactor * temperature;
        if (changed) {
            freezecount = 0;
        }
        if (((float) nochanges / notrials) < minpercent) {
            freezecount++;
        }
        pass_no++;
        changed = False;

        float diff = ratiosize - prev_ratiosize;
        if (diff < 0.0) diff = -diff;
        if (diff <= 0.0001 * prev_ratiosize) {
            /* if the change is less than 1 in 10000, it is practially
               the same */
            same++;
        } else {
            same = 0;
            prev_ratiosize = ratiosize;
        }
        /* if has seen the same solution enough number of times, it is
           time to quit */
        if (same >= samecount) {
            freezecount = freezelimit;  /* exit */
        }

#ifdef AD_TRACE
        printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n",
               pass_no, best_cutsize,
               find_cut_size(hgraph.nonets, noparts, hgraph.totnetsize,
                             best_nets, &best_pop[0]));
#endif

    }   /* while */

#ifdef AD_TRACE
    printf("Why : same=%d accept rate=%f\n", same, (float) nochanges / notrials);
#endif

    printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n",
           pass_no, best_cutsize,
           find_cut_size(hgraph.nonets, noparts, hgraph.totnetsize, best_nets, &best_pop[0]));

#ifdef AD_TRACE
    printf("Final : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("FF %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif

    /* free memory for all data structures */
    for (int i = 0; i < hgraph.nocells; i++) {
        free(cells_info[i].mgain);
    }
    free(cells_info);

    for (int i = 0; i < hgraph.nonets; i++) {
        free(nets_info[i].npartdeg);
        free(best_nets[i].npartdeg);
    }
    free(nets_info);
    free(best_nets);

    for (int i = 0; i < MAX_POP; i++) {
        free(pop[i].chrom);
        free(pop[i].parts);
        free(best_pop[i].chrom);
        free(best_pop[i].parts);
    }

    free_hypergraph(&hgraph);

    return (0);
}   /* main-sa1 */

/* EOF */

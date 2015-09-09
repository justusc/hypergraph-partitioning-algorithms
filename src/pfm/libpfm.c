
/* COPYRIGHT C 1991- Ali Dasdan */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hgpa/bucketio.h"
#include "hgpa/defs.h"
#include "hgpa/hgpa.h"
#include "hgpa/pfm.h"

/*
  ALSO TRY THESE:
  mov_value = ((float) (mov_gain + max_gain)) /
  ((float) mcount * 2.0 * max_gain);

  mov_value = ((float) (mov_gain + max_gain)) /
  ((float) log((double) (mcount + 3.0)) * 2.0 * max_gain);
*/

/* map a given mov_gain into index to a bucket array */
int map_gain(int bucketsize,
             int mov_gain,
             int mov_count,
             int max_gain,
             eval_t eval[])
{
    int mcount = mov_count;
    if (mcount == 0) {
        mcount = 1;
    }
    float mov_value = ((float) (bucketsize) /
                       (1.0 + (float) sqrt ((double) (mcount)) *
                        (float) eval[max_gain - mov_gain].val));
    return ((int) mov_value);
}   /* map_gain */

/* fill all bucket arrays */
void create_buckets(int bucketsize,
                    int nocells,
                    int noparts,
                    int max_gain,
                    eval_t eval[],
                    allele chrom[],
                    partb_t partb[][noparts - 1],
                    cells_info_t cells_info[])
{
    /* read cell gains from cells_info and insert them into buckets */
    for (int cell_no = 0; cell_no < nocells; cell_no++) {

        int part_no = chrom[cell_no];

        /* for each possible move direction */
        for (int mov_part_no = 0; mov_part_no < noparts; mov_part_no++) {
            if (mov_part_no != part_no) {

                int mov_gain = calculate_gain(cell_no, part_no, mov_part_no, cells_info);

                /* find mapped_after calculating gain & max */
                int mapped_part_no = map_part_no(mov_part_no, part_no);
                int gain_inx = map_gain(bucketsize, mov_gain, cells_info[cell_no].mcount,
                                        max_gain, eval);

                /* create a partb node */
                create_partb_node(noparts, cell_no, part_no, mapped_part_no,
                                  gain_inx, partb, cells_info);

            }   /* if mapped_part_no not equal part_no */
        }   /* for mapped_part_no */
    }   /* for part_no */
}   /* create_buckets */

/* select a cell to move */
void select_cell(int noparts,
                 selected_cell_t scell[],
                 parts_info_t parts_info[],
                 cells_t cells[],
                 partb_t partb[][noparts - 1],
                 cells_info_t cells_info[])
{
    int max_mov_value = -1;
    scell[0].mov_cell_no = -1;
    for (int from = 0; from < noparts; from++) {
        int fpcurr_size = parts_info[from].pcurr_size;
        int fpmin_size = parts_info[from].pmin_size;
        if (fpcurr_size > fpmin_size) {
            for (int to = 0; to < noparts; to++) {
                int tpmax_size = parts_info[to].pmax_size;
                int tpcurr_size = parts_info[to].pcurr_size;
                if ((from != to) && (tpcurr_size < tpmax_size)) {
                    int dest_part = map_part_no(to, from);
                    int max_inx = partb[from][dest_part].max_inx;
                    if (max_inx < 0) {
                        printf("Error: max_inx cannot be negative.\n");
                        exit(1);
                    }   /* if */
                    int cell_no = partb[from][dest_part].bnode_ptr[max_inx]->cell_no;
                    if ((max_mov_value < max_inx) &&
                        (fpcurr_size >= (fpmin_size + cells[cell_no].cweight)) &&
                        ((tpcurr_size + cells[cell_no].cweight) <= tpmax_size)) {
                        max_mov_value = max_inx;
                        scell[0].mov_cell_no = cell_no;
                        scell[0].from_part = from;
                        scell[0].to_part = to;
                    }   /* if many conditions */
                }   /* if */
            }   /* for to */
        }   /* if curr > min */
    }   /* for from */

    if (scell[0].mov_cell_no == -1) {
        printf("Error: Cannot find a node to move.\n");
        exit(1);
    }   /* if */

    scell[0].mov_gain = calculate_gain(scell[0].mov_cell_no,
                                       scell[0].from_part,
                                       scell[0].to_part,
                                       cells_info);
    parts_info[scell[0].from_part].pmax_cells--;
    parts_info[scell[0].to_part].pmax_cells++;
    parts_info[scell[0].from_part].pcurr_size -=
        cells[scell[0].mov_cell_no].cweight;
    parts_info[scell[0].to_part].pcurr_size +=
        cells[scell[0].mov_cell_no].cweight;
}   /* select_cell */

/* move selected cell, and save the move in a file */
void move_cell(mcells_t mcells[],
               int msize,
               selected_cell_t scell[],
               allele tchrom[],
               cells_info_t cells_info[])
{
    tchrom[scell[0].mov_cell_no] = scell[0].to_part;
    cells_info[scell[0].mov_cell_no].mcount++;
    mcells[msize].cell_no = scell[0].mov_cell_no;
    mcells[msize].from = scell[0].from_part;
    mcells[msize].to = scell[0].to_part;
    mcells[msize].mgain = scell[0].mov_gain;
}   /* move_cell */

/* update gains after a move */
void update_gains(int bucketsize,
                  int noparts,
                  int max_gain,
                  eval_t eval[],
                  selected_cell_t scell[],
                  cells_t cells[],
                  nets_t nets[],
                  corn_t cnets[],
                  corn_t ncells[],
                  nets_info_t nets_info[],
                  partb_t partb[][noparts - 1],
                  cells_info_t cells_info[],
                  allele tchrom[])
{
    int net_ptr = cells[scell[0].mov_cell_no].netlist;
    int mov_cell_no = scell[0].mov_cell_no;
    int from_part = scell[0].from_part;
    int to_part = scell[0].to_part;

    /* for each neighor net */
    for (int i = 0; i < cells[mov_cell_no].cno_nets; i++) {

        int net_no = cnets[net_ptr + i].corn_no;
        int net_weight = nets[net_no].nweight;
        int cell_ptr = nets[net_no].celllist;

        /* do operations before the move */
        if (nets_info[net_no].npartdeg[from_part] == nets[net_no].nno_cells) {

            update1(False, bucketsize, noparts, max_gain, from_part, mov_cell_no,
                    cell_ptr, net_no, net_weight,
                    eval, nets, ncells, partb, cells_info, tchrom);

        } else if (nets_info[net_no].npartdeg[from_part] == (nets[net_no].nno_cells - 1)) {

            update2(False, bucketsize, noparts, max_gain, from_part, mov_cell_no,
                    cell_ptr, net_no, net_weight,
                    eval, nets, ncells, partb, cells_info, tchrom);

        }   /* else */

        /* update net info */
        nets_info[net_no].npartdeg[from_part]--;
        nets_info[net_no].npartdeg[to_part]++;

        /* do operations after the move */
        if (nets_info[net_no].npartdeg[to_part] == nets[net_no].nno_cells) {

            update1(True, bucketsize, noparts, max_gain, to_part, mov_cell_no,
                    cell_ptr, net_no, net_weight,
                    eval, nets, ncells, partb, cells_info, tchrom);

        } else if (nets_info[net_no].npartdeg[to_part] == (nets[net_no].nno_cells - 1)) {

            update2(True, bucketsize, noparts, max_gain, to_part, mov_cell_no,
                    cell_ptr, net_no, net_weight,
                    eval, nets, ncells, partb, cells_info, tchrom);

        }   /* else */

    }   /* for i */
}   /* update_gains */

/* update gain of a cell */
void update1(int flag,
             int bucketsize,
             int noparts,
             int max_gain,
             int dest_part,
             int mov_cell_no,
             int cell_ptr,
             int net_no,
             int net_weight,
             eval_t eval[],
             nets_t nets[],
             corn_t ncells[],
             partb_t partb[][noparts - 1],
             cells_info_t cells_info[],
             allele tchrom[])
{
    int i = 0;
    while (i < nets[net_no].nno_cells) {

        int other_cell = ncells[cell_ptr + i].corn_no;
        int other_part_no = tchrom[other_cell];
        if (other_cell != mov_cell_no) {

            if (flag == False) {
                cells_info[other_cell].mgain[dest_part] -= net_weight;
            } else {
                cells_info[other_cell].mgain[dest_part] += net_weight;
            }

            /* since internal gains change, update all gains */
            delete_partb_nodes_of_cell(noparts, other_cell, other_part_no,
                                       partb, cells_info);
            create_partb_nodes_of_cell(bucketsize, noparts, max_gain,
                                       other_cell, other_part_no,
                                       eval, partb, cells_info);
        }   /* if */

        i++;
    }   /* while */
}   /* update1 */

/* update gain of a cell */
void update2(int flag,
             int bucketsize,
             int noparts,
             int max_gain,
             int dest_part,
             int mov_cell_no,
             int cell_ptr,
             int net_no,
             int net_weight,
             eval_t eval[],
             nets_t nets[],
             corn_t ncells[],
             partb_t partb[][noparts - 1],
             cells_info_t cells_info[],
             allele tchrom[])
{
    int other_cell, other_part_no;

    int found = False;
    int i = 0;
    while ((i < nets[net_no].nno_cells) && (!found)) {

        other_cell = ncells[cell_ptr + i].corn_no;
        other_part_no = tchrom[other_cell];
        if ((other_cell != mov_cell_no) &&
            (other_part_no != dest_part)) {
            found = True;
        }   /* if */

        i++;
    }   /* while */

    if (!found) return;

    if (flag == False) {
        cells_info[other_cell].mgain[dest_part] -= net_weight;
    } else {
        cells_info[other_cell].mgain[dest_part] += net_weight;
    }

    int mov_gain = calculate_gain(other_cell, other_part_no,
                                  dest_part, cells_info);
    int gain_inx = map_gain(bucketsize, mov_gain, cells_info[other_cell].mcount,
                            max_gain, eval);
    int mapped_part_no = map_part_no(dest_part, other_part_no);
    bnode_ptr_t tnode_ptr = delete_partb_node(False, mapped_part_no,
                                              &partb[other_part_no][mapped_part_no],
                                              &cells_info[other_cell]);
    insert_partb_node(tnode_ptr, mapped_part_no, gain_inx,
                      &partb[other_part_no][mapped_part_no],
                      &cells_info[other_cell]);
}   /* update2 */

/* fill eval array */
void fill_eval(int max_gain,
               float K,
               eval_t eval[])
{
    for (int i = 0; i <= (2 * max_gain); i++)
        eval[i].val = (float) exp((double) (-K * (max_gain - i)));
}   /* fill_eval */

void create_partb_nodes_of_cell(int bucketsize,
                                int noparts,
                                int max_gain,
                                int cell_no,
                                int part_no,
                                eval_t eval[],
                                partb_t partb[][noparts - 1],
                                cells_info_t cells_info[])
{
    if (cell_no == -1) {
        return;
    }

    /* for each possible move direction */
    for (int mov_part_no = 0; mov_part_no < noparts; mov_part_no++) {

        if (mov_part_no != part_no) {  /* part_no is home_part of cell_no */

            int mov_gain = calculate_gain(cell_no, part_no, mov_part_no, cells_info);

            /* find mapped_after calculating gain & max */
            int mapped_part_no = map_part_no(mov_part_no, part_no);
            int gain_inx = map_gain(bucketsize, mov_gain, cells_info[cell_no].mcount,
                                    max_gain, eval);

            /* create a partb node */
            create_partb_node(noparts, cell_no, part_no, mapped_part_no,
                              gain_inx, partb, cells_info);

        }   /* if mapped_part_no not equal part_no */

    }   /* for mapped_part_no */
}   /* create_partb_nodes_of_cell */

/* calculate K and scale factor */
void calculate_scale(int nocells,
                     int noparts,
                     int max_gain,
                     float *K)
{
    double ratio = log((double) ((1.0 - EPSILON) / EPSILON));
    *K = (1.0 / max_gain) * (float) ratio;
    double scale = (double) pow(ratio, (double) 1.0 / max_gain);
    scale = (scale + 1.0 + ratio * noparts * sqrt((double) nocells)) / (scale - 1.0);
    scale = (double) pow(ratio, (double) 1.0 / max_gain);
    scale = (scale + 1.0 + ratio * noparts) / (scale - 1.0);
}   /* calculate_scale */

/* EOF */

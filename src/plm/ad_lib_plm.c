
/* COPYRIGHT C 1991- Ali Dasdan */

#include <stdio.h>
#include <stdlib.h>
#include "hpga/ad_defs.h"
#include "hpga/ad_bucketio.h"
#include "hpga/ad_lib.h"
#include "hpga/ad_lib_plm.h"

/* map a given mov_gain into index to a bucket array */
int map_gain(int mov_gain, int max_gain)
{
    return (mov_gain + max_gain);
}   /* map_gain */

/* fill all bucket arrays */
void create_buckets(int nocells,
                    int noparts,
                    int max_gain,
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
                int gain_inx = map_gain(mov_gain, max_gain);

                /* create a partb node */
                create_partb_node(noparts, cell_no, part_no, mapped_part_no,
                                  gain_inx, partb, cells_info);

            }   /* if mapped_part_no not equal part_no */
        }   /* for mapped_part_no */
    }   /* for part_no */
}   /* create_buckets */

/* select a cell to move */
int select_cell(int noparts,
                selected_cell_t scell[],
                parts_info_t parts_info[],
                cells_t cells[],
                partb_t partb[][noparts - 1],
                cells_info_t cells_info[])
{
    /*
      add to the 1st "if"   && (parts_info[from].pmax_cells > 0)
    */
    int move_possible = True;
    int max_mov_value = -1;
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
                    }
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

    if (max_mov_value != -1) {

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

    } else {

        move_possible = False;
        max_mov_value = -1;
        for (int from = 0; from < noparts; from++) {
            // int fpcurr_size = parts_info[from].pcurr_size;
            // int fpmin_size = parts_info[from].pmin_size;
            for (int to = 0; to < noparts; to++) {
                // int tpmax_size = parts_info[to].pmax_size;
                // int tpcurr_size = parts_info[to].pcurr_size;
                if (from != to) {
                    int dest_part = map_part_no(to, from);
                    int max_inx = partb[from][dest_part].max_inx;
                    if (max_inx > 0) {
                        int cell_no = partb[from][dest_part].bnode_ptr[max_inx]->cell_no;
                        if (max_mov_value < max_inx) {
                            max_mov_value = max_inx;
                            scell[0].mov_cell_no = cell_no;
                            scell[0].from_part = from;
                        }   /* if many conditions */
                    }   /* if */
                }   /* if */
            }   /* for to */
        }   /* for from */

    }   /* else */

    return move_possible;
}   /* select_cell */

/* move selected cell, and save the move in a file */
void move_cell(mcells_t mcells[],
               int msize,
               selected_cell_t scell[],
               allele tchrom[])
{
    tchrom[scell[0].mov_cell_no] = scell[0].to_part;
    mcells[msize].cell_no = scell[0].mov_cell_no;
    mcells[msize].from = scell[0].from_part;
    mcells[msize].to = scell[0].to_part;
    mcells[msize].mgain = scell[0].mov_gain;
}   /* move_cell */

/* update gains after a move */
void update_gains(int noparts,
                  int max_gain,
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

            update1(False, noparts, max_gain, from_part, mov_cell_no,
                    cell_ptr, net_no, net_weight,
                    nets, ncells, partb, cells_info, tchrom);

        } else if (nets_info[net_no].npartdeg[from_part] == (nets[net_no].nno_cells - 1)) {

            update2(False, noparts, max_gain, from_part, mov_cell_no,
                    cell_ptr, net_no, net_weight,
                    nets, ncells, partb, cells_info, tchrom);

        }

        /* update net info */
        nets_info[net_no].npartdeg[from_part]--;
        nets_info[net_no].npartdeg[to_part]++;

        /* do operations after the move */
        if (nets_info[net_no].npartdeg[to_part] == nets[net_no].nno_cells) {

            update1(True, noparts, max_gain, to_part, mov_cell_no,
                    cell_ptr, net_no, net_weight,
                    nets, ncells, partb, cells_info, tchrom);

        } else if (nets_info[net_no].npartdeg[to_part] == (nets[net_no].nno_cells - 1)) {

            update2(True, noparts, max_gain, to_part, mov_cell_no,
                    cell_ptr, net_no, net_weight,
                    nets, ncells, partb, cells_info, tchrom);

        }   /* else */

    }   /* for i */
}   /* update_gains */

/* update gain of a cell */
void update1(int flag,
             int noparts,
             int max_gain,
             int dest_part,
             int mov_cell_no,
             int cell_ptr,
             int net_no,
             int net_weight,
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
        if ((other_cell != mov_cell_no) && (cells_info[other_cell].locked == False)) {

            if (flag == False) {
                cells_info[other_cell].mgain[dest_part] -= net_weight;
            } else {
                cells_info[other_cell].mgain[dest_part] += net_weight;
            }

            /* since internal gains change, update all gains */
            delete_partb_nodes_of_cell(noparts, other_cell, other_part_no,
                                       partb, cells_info);
            create_partb_nodes_of_cell(noparts, max_gain,
                                       other_cell, other_part_no,
                                       partb, cells_info);
        }   /* if */

        i++;
    }   /* while */
}   /* update1 */

/* update gain of a cell */
void update2(int flag,
             int noparts,
             int max_gain,
             int dest_part,
             int mov_cell_no,
             int cell_ptr,
             int net_no,
             int net_weight,
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
            (cells_info[other_cell].locked == False) &&
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
    int gain_inx = map_gain(mov_gain, max_gain);
    int mapped_part_no = map_part_no(dest_part, other_part_no);
    bnode_ptr_t tnode_ptr = delete_partb_node(False, mapped_part_no,
                                              &partb[other_part_no][mapped_part_no],
                                              &cells_info[other_cell]);
    insert_partb_node(tnode_ptr, mapped_part_no, gain_inx,
                      &partb[other_part_no][mapped_part_no],
                      &cells_info[other_cell]);
}   /* update2 */

void create_partb_nodes_of_cell(int noparts,
                                int max_gain,
                                int cell_no,
                                int part_no,
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
            int gain_inx = map_gain(mov_gain, max_gain);

            /* create a partb node */
            create_partb_node(noparts, cell_no, part_no, mapped_part_no,
                              gain_inx, partb, cells_info);

        }   /* if mapped_part_no not equal part_no */
    }   /* for mapped_part_no */
}   /* create_partb_nodes_of_cell */

/* EOF */

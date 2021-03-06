#ifndef HGPA_PARTITION_INCLUDED
#define HGPA_PARTITION_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */

#include "hgpa/partition.h"

#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "hgpa/fileio.h"
#include "hgpa/random.h"
#include "hgpa/defs.h"

/* create initial partition */
int create_partition(int nocells,
                     int noparts,
                     int totcellsize,
                     int max_cweight,
                     float *off_ratio,
                     cells_t cells[],
                     nets_t nets[],
                     corn_t cnets[],
                     ind_t *ind);

/* copy partition properties to parts_info for temporary use */
void copy_partition(int noparts,
                    parts_info_t parts_info[],
                    ind_t *ind);

/* read a partition prepared beforehand */
int read_partition(FILE *fp,
                   char *filename,
                   int noparts,
                   cells_t cells[],
                   nets_t nets[],
                   corn_t cnets[],
                   ind_t *ind);

/* write a partition */
void write_partition(FILE *fp,
                     char *filename,
                     int nocells,
                     int noparts,
                     ind_t *ind);

#endif

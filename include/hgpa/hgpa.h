#ifndef HGPA_LIB_INCLUDED
#define HGPA_LIB_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* initialize all bucket indices and pointers */
void init_buckets(int noparts, 
                  int bucketsize,
                  partb_t partb[][noparts - 1]);
 
/* map part no such that home_part is excluded */
int map_part_no(int dest_part, int home_part);
 
/* compute move gain from home_part to dest_part */
int calculate_gain(int cell_no, 
                   int home_part, 
                   int dest_part,
                   cells_info_t cells_info[]);

/* compute gains of all cells and place them into cells_info */
void compute_gains(int nocells, 
                   int noparts,
                   cells_t cells[],
                   nets_t nets[],
                   corn_t cnets[],
                   cells_info_t cells_info[],
                   allele tchrom[]);

/* compute gains of all cells and place them into cells_info */
void compute_gains2(int nocells, 
                    int noparts,
                    cells_t cells[],
                    nets_t nets[],
                    corn_t cnets[],
                    cells_info_t cells_info[],
                    allele tchrom[],
                    nets_info_t nets_info[]);

/* free all allocated nodes */
void free_nodes(int noparts, 
                int bucketsize,
                partb_t partb[][noparts - 1]);

/* count number of bucket nodes */
void number_nodes(int noparts, 
                  int bucketsize, 
                  int *npartb,
                  partb_t partb[][noparts - 1]);

/* find set of cells to be actually moved */
int find_move_set(mcells_t mcells[],
                  int msize,
                  int *max_mcells_inx);

/* move cells actually - permanently */
int move_cells(int wflag,
               int nocells, 
               int msize,
               mcells_t mcells[],
               int max_mcells_inx, 
               int cutsize, 
               int *glob_inx,
               ind_t *ind,
               cells_t cells[],
               nets_t nets[],
               corn_t cnets[]);

/* finds cut size of a given partition - used for control */
int find_cut_size(int nonets, 
                  int noparts, 
                  int totnetsize,
                  nets_t nets[],
                  ind_t *ind);

/* save nets info in nets_info */
void copy_nets_info(int nonets, 
                    int noparts,
                    nets_t nets[],
                    nets_info_t nets_info[]);

#endif

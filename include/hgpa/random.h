#ifndef HGPA_RANDOM_INCLUDED
#define HGPA_RANDOM_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */

/* These two lines of code are TERRIBLY necessary */
double drand48();
void srand48();

/* initializes random number generator with seed or */
/* with any value if seed = -1 */
long randomize(long seed);

/* generates a float random number in [0, 1) */
float rand01();

/* generates a int random number in [min, max] */
/* if int output is needed, prms are true min & max */
int irandom(int min, int max);

#endif

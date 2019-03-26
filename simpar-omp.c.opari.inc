#include "pomp_lib.h"

struct ompregdescr omp_rd_1 = {
  "parallelfor", "", 0, "simpar-omp.c", 55, 55, 59, 59
};

#define POMP_DLIST_00001 shared(omp_rd_1)

struct ompregdescr omp_rd_2 = {
  "parallelfor", "", 0, "simpar-omp.c", 75, 75, 83, 83
};

#define POMP_DLIST_00002 shared(omp_rd_2)


#include "pomp_lib.h"


extern struct ompregdescr omp_rd_1;

int POMP_MAX_ID = 2;

struct ompregdescr* pomp_rd_table[2] = {
  0,
  &omp_rd_1,
};

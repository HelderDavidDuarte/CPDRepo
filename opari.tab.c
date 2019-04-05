#include "pomp_lib.h"


extern struct ompregdescr omp_rd_1;
extern struct ompregdescr omp_rd_2;

int POMP_MAX_ID = 3;

struct ompregdescr* pomp_rd_table[3] = {
  0,
  &omp_rd_1,
  &omp_rd_2,
};

#ifndef main_headers
#define main_headers
#include "libs.hpp"
#endif

void ReadAAPPSeq(char* filename, char** aa_seq, void*** pp_seq, int* len_seq);
float CalcPPEuclidianDist(void* elem_1, void* elem_2);
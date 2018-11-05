#ifndef main_headers
#define main_headers
#include "libs.hpp"
#endif

int arg_max(float* theArray, float* pMaxValue, int len);
void NW_Align(void** seq_1, char* chr_seq_1, int len_seq_1, void** seq_2, char* chr_seq_2, int len_seq_2, float (*sim_func)(void*, void*), int gap_start, int gap_ext);
void ReverseString(char * theString);
 
void print_scores_mat(float** scoresMat, int len_seq_1, int len_seq_2);
void print_dir_mat(char** dirMat, int len_seq_1, int len_seq_2);
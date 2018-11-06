#ifndef main_headers
#define main_headers
#include "libs.hpp"
#endif

#ifndef univ_nw
#define univ_nw
#include "univ_nw.hpp"
#endif

#ifndef blosum
#define blosum
#include "blosum.hpp"
#endif

#ifndef phipsi
#define phipsi
#include "phipsi_calc.hpp"
#endif

using namespace std;

float relay_func(void*, void*);
Blosum_mat* b_mat = NULL;

int main() {
    int len_seq_1 = 0;    
    char* chr_seq_1 = NULL;
    void** seq_1 = NULL;
    
    int len_seq_2 = 0;    
    char* chr_seq_2 = NULL;
    void** seq_2 = NULL;

    char* filename_1 = "s_1.pre";
    char* filename_2 = "s_2.pre";

    ReadAAPPSeq("s_1.pre", &chr_seq_1, &seq_1, &len_seq_1);
    ReadAAPPSeq("s_2.pre", &chr_seq_2, &seq_2, &len_seq_2);

    NW_Align(seq_1, chr_seq_1, len_seq_1, seq_2, chr_seq_2, len_seq_2, relay_func, 10, 5);

    // b_mat = new Blosum_mat();

    // char* test_seq_1 = "IPGAWD";
    // char* test_seq_2 = "AAAAIPGAWDDDDD";

    // char* seq_1[6] = { NULL };
    // char* seq_2[14] = { NULL };

    // for (int i = 0; i < 6; i++) {
    //     seq_1[i] = test_seq_1 + i;
    // }
    
    // for (int i = 0; i < 14; i++) {
    //     seq_2[i] = test_seq_2 + i;
    // }

    // // NW_Align((void**)seq_2, "AAAAIPGAWDDDDD", 14, (void**)seq_1, "IPGAWD", 6, relay_func, 12, 7);
    // return 0;
}

float relay_func(void* elem_1, void* elem_2) {

    // return b_mat->BlosumSimFunc(elem_1, elem_2);
    return CalcPPEuclidianDist(elem_1, elem_2);
}
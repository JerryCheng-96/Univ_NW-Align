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

using namespace std;

float relay_func(void*, void*);
Blosum_mat* b_mat = NULL;

int main() {
    b_mat = new Blosum_mat();

    char* test_seq_1 = "IPGAWD";
    char* test_seq_2 = "AAAAIPGAWDDDDD";

    char* seq_1[6] = { NULL };
    char* seq_2[15] = { NULL };

    for (int i = 0; i < 6; i++) {
        seq_1[i] = test_seq_1 + i;
    }
    
    for (int i = 0; i < 15; i++) {
        seq_2[i] = test_seq_2 + i;
    }

    NW_Align((void**)seq_1, 6, (void**)seq_2, 15, relay_func, 12, 7);
    return 0;
}

float relay_func(void* elem_1, void* elem_2) {
    return b_mat->BlosumSimFunc(elem_1, elem_2);
}
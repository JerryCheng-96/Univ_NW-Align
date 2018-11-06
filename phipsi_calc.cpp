#ifndef phipsi
#define phipsi
#include "phipsi_calc.hpp"
#endif

using namespace std;

void ReadAAPPSeq(char* filename, char** aa_seq, void*** pp_seq, int* len_seq) {
    ifstream fin(filename);
    string buf;

    getline(fin, buf);

    *len_seq = buf.length(); 

    *aa_seq = (char*)malloc(sizeof(char) * (*len_seq + 1));
    buf.copy(*aa_seq, *len_seq, 0);

    getline(fin, buf);

    *pp_seq = (void**)malloc(sizeof(void*) * *len_seq);
    float* pp_val_array = (float*)malloc(sizeof(float) * 2 * *len_seq);

    int pp_seq_idx = 0;

    while (getline(fin, buf)) {
        char tmp[buf.length() + 1] = {0};
        buf.copy(tmp, buf.length(), 0);

        int idx = 0;
        for (idx = 0; tmp[idx] != ',' && idx <= buf.length() ; idx++) { }
        tmp[idx] = '\0';
        *(pp_val_array + 2 * pp_seq_idx) = atof(tmp);
        *(pp_val_array + 2 * pp_seq_idx + 1) = atof(tmp + idx + 1);
        (*pp_seq)[pp_seq_idx] = (void*)(pp_val_array + 2 * pp_seq_idx);

        pp_seq_idx++;
        printf("%s", tmp);
    }

    return;
}

float CalcPPEuclidianDist(void* elem_1, void* elem_2) {
    float* aa_1 = (float*)elem_1;
    float* aa_2 = (float*)elem_2;

    return sqrt(pow(aa_1[0] - aa_2[0], 2) + pow(aa_1[1] - aa_2[1], 2));
}
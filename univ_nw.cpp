#ifndef univ_nw
#define univ_nw
#include "univ_nw.hpp"
#endif

using namespace std;

int arg_max(float* theArray, float* pMaxValue, int len) {
    float maxValue = theArray[0];
    int maxIndex = 0;

    for(int i = 0; i < len; i++) {
        if (theArray[i] > maxValue) {
            maxValue = theArray[i];
            maxIndex = i;
        }
    }

    *pMaxValue = maxValue;
    return maxIndex;
}

void NW_Align(void** seq_1, int len_seq_1, 
              void** seq_2, int len_seq_2, 
              float (*sim_func)(void*, void*), 
              int gap_start, int gap_ext) {

    // Initializing the matrices
    float** scoresMat = (float**)malloc((len_seq_1 + 1) * sizeof(float*));
    char** dirMat = (char**)malloc((len_seq_1 + 1) * sizeof(char*));

    
    for(int i = 0; i <= len_seq_1 + 1; i++) {
        scoresMat[i] = (float*)malloc((len_seq_2 + 1) * sizeof(float));
        dirMat[i] = (char*)malloc((len_seq_2 + 1) * sizeof(char));
    }
    

    // Setting init values
    scoresMat[0][0] = 0;
    scoresMat[0][1] = -gap_start;
    scoresMat[1][0] = -gap_start;
    dirMat[0][0] = 'd';
    dirMat[0][1] = 'l';
    dirMat[1][0] = 't';

    for(int j = 2; j <= len_seq_2; j++) {
        scoresMat[0][j] = -(gap_ext * (j - 1)) - gap_start;
        dirMat[0][j] = 'l';
    }
    
    for(int i = 2; i <= len_seq_1; i++) {
        scoresMat[i][0] = -(gap_ext * (i - 1)) - gap_start;
        dirMat[i][0] = 't';
    }

    float scoreValues[] = {0,0,0};

    for(int i = 1; i <= len_seq_1; i++) {
        for(int j = 1; j <= len_seq_2; j++) {
            float sim_val = sim_func(seq_1[i - 1], seq_2[j - 1]);
            scoreValues[0] = scoresMat[i - 1][j - 1] + sim_func(seq_1[i - 1], seq_2[j - 1]);
            scoreValues[1] = scoresMat[i - 1][j] - (dirMat[i - 1][j] == 'd' ? gap_start : gap_ext);
            scoreValues[2] = scoresMat[i][j - 1] - (dirMat[i][j - 1] == 'd' ? gap_start : gap_ext);

            switch (arg_max(scoreValues, &scoresMat[i][j], 3)) {
                case 0:
                    dirMat[i][j] = 'd';
                    break;
                case 1: 
                    dirMat[i][j] = 't';
                    break;
                case 2: 
                    dirMat[i][j] = 'l';
                    break;
            }
        }
    }

    print_scores_mat(scoresMat, len_seq_1, len_seq_2);
    print_dir_mat(dirMat, len_seq_1, len_seq_2);

    return;
}

void print_scores_mat(float** scoresMat, int len_seq_1, int len_seq_2) {
    for(int i = 0; i <= len_seq_1; i++)
    {
        for(int j = 0; j <= len_seq_2; j++)
        {
            printf("%.2f\t", scoresMat[i][j]);
        }
        printf("\n");
    }
}

void print_dir_mat(char** dirMat, int len_seq_1, int len_seq_2) {
    for(int i = 0; i <= len_seq_1; i++)
    {
        for(int j = 0; j <= len_seq_2; j++)
        {
            printf("%c\t", dirMat[i][j]);
        }
        printf("\n");
    }
}
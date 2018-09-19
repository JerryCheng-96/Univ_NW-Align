#ifndef main_headers
#define main_headers
#include "libs.hpp"
#endif

#define BLOSUM_SIZE 26
#define MAX_LINE_LENGTH 1024
#define MAX_ARRAY_SIZE 50

class Blosum_mat
{
  private:
    void PrintBlosum(int *blosum);

	int blosum_mat_vals[BLOSUM_SIZE * BLOSUM_SIZE] = { 4 ,-2 ,0 ,-2 ,-1 ,-2 ,0 ,-2 ,-1 ,0 ,-1 ,-1 ,-1 ,-2 ,-4 ,-1 ,-1 ,-1 ,1 ,0 ,0 ,0 ,-3 ,0 ,-2 ,-1 ,-2 ,4 ,-3 ,4 ,1 ,-3 ,-1 ,0 ,-3 ,0 ,0 ,-4 ,-3 ,3 ,-4 ,-2 ,0 ,-1 ,0 ,-1 ,0 ,-3 ,-4 ,-1 ,-3 ,1 ,0 ,-3 ,9 ,-3 ,-4 ,-2 ,-3 ,-3 ,-1 ,0 ,-3 ,-1 ,-1 ,-2 ,-4 ,-3 ,-3 ,-3 ,-1 ,-1 ,0 ,-1 ,-2 ,-2 ,-2 ,-3 ,-2 ,4 ,-3 ,6 ,2 ,-3 ,-1 ,-1 ,-3 ,0 ,-1 ,-4 ,-3 ,0 ,-4 ,-1 ,0 ,-2 ,0 ,-1 ,0 ,-3 ,-4 ,-1 ,-3 ,1 ,-1 ,1 ,-4 ,2 ,5 ,-3 ,-2 ,0 ,-3 ,0 ,1 ,-3 ,-2 ,-2 ,-4 ,-1 ,2 ,0 ,0 ,-1 ,0 ,-2 ,-3 ,-1 ,-2 ,4 ,-2 ,-3 ,-2 ,-3 ,-3 ,6 ,-3 ,-1 ,0 ,0 ,-3 ,0 ,0 ,-3 ,-4 ,-4 ,-3 ,-3 ,-2 ,-2 ,0 ,-1 ,1 ,-1 ,3 ,-3 ,0 ,-1 ,-3 ,-1 ,-2 ,-3 ,6 ,-2 ,-4 ,0 ,-2 ,-4 ,-3 ,0 ,-4 ,-2 ,-2 ,-2 ,0 ,-2 ,0 ,-3 ,-2 ,-1 ,-3 ,-2 ,-2 ,0 ,-3 ,-1 ,0 ,-1 ,-2 ,8 ,-3 ,0 ,-1 ,-3 ,-2 ,-1 ,-4 ,-2 ,0 ,0 ,-1 ,-2 ,0 ,-3 ,-2 ,-1 ,2 ,0 ,-1 ,-3 ,-1 ,-3 ,-3 ,0 ,-4 ,-3 ,4 ,0 ,-3 ,2 ,1 ,-1 ,-4 ,-3 ,-3 ,-3 ,-2 ,-1 ,0 ,3 ,-3 ,-1 ,-1 ,-3 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,-1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,-1 ,0 ,-3 ,-1 ,1 ,-3 ,-2 ,-1 ,-3 ,0 ,5 ,-2 ,-1 ,0 ,-4 ,-1 ,1 ,2 ,0 ,-1 ,0 ,-2 ,-3 ,-1 ,-2 ,1 ,-1 ,-4 ,-1 ,-4 ,-3 ,0 ,-4 ,-3 ,2 ,0 ,-2 ,4 ,2 ,-4 ,-4 ,-3 ,-2 ,-2 ,-2 ,-1 ,0 ,1 ,-2 ,-1 ,-1 ,-3 ,-1 ,-3 ,-1 ,-3 ,-2 ,0 ,-3 ,-2 ,1 ,0 ,-1 ,2 ,5 ,-1 ,-4 ,-2 ,0 ,-1 ,-1 ,-1 ,0 ,1 ,-1 ,-1 ,-1 ,-1 ,-2 ,3 ,-3 ,1 ,0 ,-3 ,0 ,1 ,-3 ,0 ,0 ,-3 ,-2 ,-1 ,-4 ,-2 ,0 ,0 ,1 ,0 ,0 ,-3 ,-4 ,-1 ,-2 ,0 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,0 ,-4 ,-4 ,-4 ,-4 ,1 ,-4 ,-4 ,-4 ,-4 ,-4 ,0 ,-4 ,-4 ,-4 ,-4 ,-4 ,-1 ,-2 ,-3 ,-1 ,-1 ,-4 ,-2 ,-2 ,-3 ,0 ,-1 ,-3 ,-2 ,1 ,-4 ,7 ,-1 ,-2 ,-1 ,-1 ,0 ,-2 ,-4 ,-2 ,-3 ,-1 ,-1 ,0 ,-3 ,0 ,2 ,-3 ,-2 ,0 ,-3 ,0 ,1 ,-2 ,0 ,0 ,-4 ,-1 ,5 ,1 ,0 ,-1 ,0 ,-2 ,-2 ,-1 ,-1 ,3 ,-1 ,-1 ,-3 ,-2 ,0 ,-3 ,-2 ,0 ,-3 ,0 ,2 ,-2 ,-1 ,0 ,-4 ,-2 ,1 ,5 ,-1 ,-1 ,0 ,-3 ,-3 ,-1 ,-2 ,0 ,1 ,0 ,-1 ,0 ,0 ,-2 ,0 ,-1 ,-2 ,0 ,0 ,-2 ,-1 ,0 ,-4 ,-1 ,0 ,-1 ,4 ,1 ,0 ,-2 ,-3 ,0 ,-2 ,0 ,0 ,-1 ,-1 ,-1 ,-1 ,-2 ,-2 ,-2 ,-1 ,0 ,-1 ,-1 ,-1 ,-3 ,-4 ,-1 ,-1 ,-1 ,1 ,5 ,0 ,0 ,-2 ,0 ,-2 ,-1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,-3 ,-1 ,-3 ,-2 ,-1 ,-3 ,-3 ,3 ,0 ,-2 ,1 ,1 ,-3 ,-4 ,-2 ,-2 ,-3 ,-2 ,0 ,0 ,4 ,-3 ,-1 ,-1 ,-2 ,-3 ,-4 ,-2 ,-4 ,-3 ,1 ,-2 ,-2 ,-3 ,0 ,-3 ,-2 ,-1 ,-1 ,-4 ,-4 ,-2 ,-3 ,-3 ,-2 ,0 ,-3 ,11 ,-2 ,2 ,-3 ,0 ,-1 ,-2 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,0 ,-1 ,-1 ,-1 ,-1 ,-4 ,-2 ,-1 ,-1 ,0 ,0 ,0 ,-1 ,-2 ,-1 ,-1 ,-1 ,-2 ,-3 ,-2 ,-3 ,-2 ,3 ,-3 ,2 ,-1 ,0 ,-2 ,-1 ,-1 ,-2 ,-4 ,-3 ,-1 ,-2 ,-2 ,-2 ,0 ,-1 ,2 ,-1 ,7 ,-2 ,-1 ,1 ,-3 ,1 ,4 ,-3 ,-2 ,0 ,-3 ,0 ,1 ,-3 ,-1 ,0 ,-4 ,-1 ,3 ,0 ,0 ,-1 ,0 ,-2 ,-3 ,-1 ,-2 ,4 };

  public:
    Blosum_mat();
    ~Blosum_mat();

    float BlosumSimFunc(void* elem_1, void* elem_2);
};

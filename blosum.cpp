#ifndef blosum
#define blosum
#include "blosum.hpp"
#endif

#pragma region GettingBlosumMatrix

void Blosum_mat::PrintBlosum(int* blosum_mat)
{
	for (int i = 0; i < BLOSUM_SIZE; i++)
	{
		printf("%c\n", i + 'A');
		for (int j = 0; j < BLOSUM_SIZE; j++)
		{
			printf("%d\t", blosum_mat[i * BLOSUM_SIZE + j]);
		}
		printf("\n");
	}
}

float Blosum_mat::BlosumSimFunc(void* elem_1, void* elem_2) {
	char* aa_1 = (char*) elem_1;
	char* aa_2 = (char*) elem_2;

	return (float)blosum_mat_vals[(*aa_1 - 'A') * BLOSUM_SIZE + (*aa_2 - 'A')];
}

Blosum_mat::Blosum_mat() {
	// PrintBlosum(blosum_mat_vals);
}

#pragma endregion

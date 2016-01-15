
#include <stdio.h>
#include <vector>
#include "nautyinit.h"

using namespace std;

struct nauty_output {
	int label[MAXN], ptn[MAXN], orbits[MAXN];
	graph *canonical = (graph *)malloc(MAXN * sizeof(int) * 2);
};

vector<unsigned int> mkLibAdj();
void mkLibAdjHelper(vector<unsigned int> &, int, unsigned int, int);
unsigned int* convert_NR_to_Set_array(unsigned int);
nauty_output ComputeCanonAdj(unsigned int, unsigned int*);
unsigned int convert_R_Adjacency_to_Bitset(unsigned int*);
vector<int> Graphlet_Lib(unsigned int);
vector<int> Find_Orbits(unsigned int);
vector<vector<int>> mkDict();


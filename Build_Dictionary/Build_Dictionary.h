
#include <stdio.h>
#include <vector>
#include "nautyinit.h"

using namespace std;

vector<unsigned int> mkLibAdj();
void mkLibAdjHelper(vector<unsigned int> &, int, unsigned int, int);
unsigned int* convert_NR_to_Set_array(unsigned int);
unsigned int * ComputeCanonAdj(unsigned int, unsigned int*);
unsigned int convert_R_Adjacency_to_Bitset(unsigned int*);
vector<int> Graphlet_Lib(unsigned int);


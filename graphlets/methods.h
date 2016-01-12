// methods.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include "nautyinit.h"
#include "Build_Dictionary.h"

using namespace std;

const vector<vector<int>> enumerateSubgraphs(vector<int> &, vector<int> & , int);
void extendSubgraph(vector<int>, vector<int>, int, vector<int>&, vector<int>&, const int, vector<vector<int>> &);
void removeSmallerThanW(vector<int>&, int & );
vector<int> findNeighborhood(vector<int>&, vector<int>&, vector<int>&);
bool ifContains(vector<int>&, int);
vector<int> vertexNeigbors(vector<int>&, vector<int>&, int ); 
vector<int> adjacencyFromGraphlets(vector<int>&, vector<int>&, vector<vector<int>> &);
vector<vector<unsigned int>> adjacencyNSFromGraphlets(vector<int>&, vector<int>&, vector<vector<int>> &);
int ifConnected(int ,int ,vector<int>& ,vector<int>& );
unsigned int *  ComputeLabel(unsigned int n, unsigned int *adjacencyMatrix);
vector<unsigned int *> getCanons(unsigned int, vector<vector<unsigned int>>);
// methods.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include "nautyinit.h"

using namespace std;

struct network {
	vector<int> graphVector;		//graphVector is the finall vector used to represent the network
	vector<int> delimiterInfo;	//the array used to store the information on where in graphVector each array starts
	int Number_of_Edges;
};

struct nauty_output {
	int label[MAXN], ptn[MAXN], orbits[MAXN];
	graph *canonical = (graph *)malloc(MAXN * sizeof(int) * 2);
};

vector<vector<int>> sparse_network(network&);
void enumerateSubgraphs(network&, vector<vector<int>>&, const int, vector<vector<int>>&, vector<vector<int>>&);
void extendSubgraph(vector<int>&, vector<int>&, int, network&, vector<vector<int>>&, const int, vector<vector<int>>&, vector<vector<int>>&);
void removeSmallerThanW(vector<int>&, int & );
vector<int> findNeighborhood(vector<int>&, vector<int>&, vector<int>&);
bool ifContains(vector<int>&, int);
int ifConnected_fast(int, int, vector<vector<int>>&);
vector<int> vertexNeigbors(vector<int>&, vector<int>&, int ); 
vector<unsigned int> adjacencyNSFromGraphlets(vector<vector<int>>&, vector<int> &);
int ifConnected(int ,int ,vector<int>& ,vector<int>& );
nauty_output  ComputeCanon(unsigned int, unsigned int);
nauty_output getCanons(unsigned int, vector<unsigned int>);
network read_network_from_file(string);
vector<vector<int>> read_dictionary_from_file(string);
unsigned int convert_Redundant_Adjacency_to_Bitset(unsigned int*, int);
void update_GDD(vector<int>&, network&, vector<vector<int>>&, vector<vector<int>>&, vector<vector<int>>&);
int find_value(int*, int, int);
void write_in_file(vector<vector<int>>);
bool check_answers(vector<vector<int>>&);
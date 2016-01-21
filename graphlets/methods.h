
#include <vector>
#include <stdlib.h>
#include "nautyinit.h"

using namespace std;

struct network {
	//vector<int> graphVector;		//graphVector is the finall vector used to represent the network
	//vector<int> delimiterInfo;	//the array used to store the information on where in graphVector each array starts
	int* gVector = static_cast<int *>(malloc(350 * sizeof(int)*2));;		//graphVector is the finall vector used to represent the network
	int* dInfo = static_cast<int *>(malloc(150 * sizeof(int)));	//the array used to store the information on where in graphVector each array starts

	int Number_of_Edges;
	int Number_of_Vertices;
};

struct nauty_output {
	int label[MAXN], ptn[MAXN], orbits[MAXN];
	graph *canonical = static_cast<graph *>(malloc(MAXN * sizeof(int) * 2));
};

vector<vector<int>> sparse_network(network&);
void enumerateSubgraphs(network&, vector<vector<int>>&, const int, vector<vector<int>>&, vector<vector<int>>&);
void extendSubgraph(vector<int>&, vector<int>&, int, network&, vector<vector<int>>&, const int, vector<vector<int>>&, vector<vector<int>>&);
void removeSmallerThanW(vector<int>&, int &);
vector<int> findNeighborhood(network&, vector<int>&);
bool ifContains(vector<int>&, int);
int ifConnected_fast(int, int, vector<vector<int>>&);
vector<int> vertexNeigbors(int*, int*, int);
vector<unsigned int> adjacencyNSFromGraphlets(vector<vector<int>>&, vector<int> &);
int ifConnected(int, int, network&);
nauty_output  ComputeCanon(unsigned int, unsigned int*);
nauty_output getCanons(unsigned int, vector<unsigned int>);
network read_network_from_file(string);
vector<vector<int>> read_dictionary_from_file(string);
unsigned int convert_Redundant_Adjacency_to_Bitset(unsigned int*, int);
void update_GDD(vector<int>&, network&, vector<vector<int>>&, vector<vector<int>>&, vector<vector<int>>&);
int find_value(int*, int, int);
void write_in_file(vector<vector<int>>);
bool check_answers(vector<vector<int>>&);
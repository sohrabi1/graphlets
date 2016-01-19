// graphlets.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <ctime>
#include "methods.h"


using namespace std;

clock_t t_start;
clock_t t_end;
clock_t t_end_s;
clock_t t_start_s;
clock_t t_update_GDD=0;
clock_t t_extendSubgraph = 0;
clock_t t_adjacencyNSFromGraphlets=0;
clock_t t_getCanons=0;
clock_t t_update_loop=0;
clock_t t_findNeighborhood = 0;
clock_t t_addNeighbor = 0;

int main()
{
	clock_t startTime;
	clock_t endTime;
	startTime = clock();

	vector<vector<int>> dictionary = read_dictionary_from_file("dictionary.bin");
	network net = read_network_from_file("exmpl100.in");
	vector<vector<int>> network_adjacency = sparse_network(net);

	vector<vector<int>> gdd(net.delimiterInfo.size()-1, vector<int>(73,0));
	enumerateSubgraphs(net, network_adjacency, 5, gdd, dictionary);
	enumerateSubgraphs(net, network_adjacency, 4, gdd, dictionary);
	enumerateSubgraphs(net, network_adjacency, 3, gdd, dictionary);
	enumerateSubgraphs(net, network_adjacency, 2, gdd, dictionary);
	write_in_file(gdd);
	endTime = clock();
	double t= (endTime - startTime) / CLOCKS_PER_SEC;
	//double tsec_update_GDD = t_update_GDD / CLOCKS_PER_SEC;
	//double tsec_extendSubgraph = t_extendSubgraph / CLOCKS_PER_SEC;
	//double tsec_adjacencyNSFromGraphlets = t_adjacencyNSFromGraphlets / CLOCKS_PER_SEC;
	//double tsec_getCanons = t_getCanons / CLOCKS_PER_SEC;
	//double tsec_update_loop = t_update_loop / CLOCKS_PER_SEC;
	//double tsec_findNeighborhood = t_findNeighborhood / CLOCKS_PER_SEC;
	//double tsec_addNeighbor = t_addNeighbor / CLOCKS_PER_SEC;
	bool check = check_answers(gdd);
	return 0;
}


vector<vector<int>> sparse_network(network& net)
{
	int r = (net.delimiterInfo.size() - 1) % 32;
	int fraction = (net.delimiterInfo.size() - 1) / 32;
	int nEdges = net.delimiterInfo.size() - 1;
	vector<vector<int>> network_adjacency(nEdges, vector<int>((fraction+1),0));
	int if_connected;
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdges; j++)
		{
			if_connected = ifConnected(i, j, net.graphVector, net.delimiterInfo);
			network_adjacency[i][j / 32] = network_adjacency[i][j / 32] << 1;
			network_adjacency[i][j / 32] = network_adjacency[i][j / 32] | if_connected;
		}
	for (int i = 0; i < nEdges; i++)
		network_adjacency[i][fraction] = network_adjacency[i][fraction] << (32 - r);
	return network_adjacency ;
}

void enumerateSubgraphs(network& n, vector<vector<int>>& net_adjacency, const int k, vector<vector<int>>& GDD, vector<vector<int>>& orbit_dict)
{
	for (int i = 0; i < (n.delimiterInfo.size() - 1) ; i++)
	{
		vector<int> extensionV = vertexNeigbors(n.graphVector, n.delimiterInfo, i);	//stores neigbors for each vertex
		vector<int> subgraphV(1, i);
		extendSubgraph(subgraphV, extensionV, i, n, net_adjacency, k, GDD, orbit_dict);
	}

}


void extendSubgraph(vector<int>& subgraphVertices, vector<int>& extensionVertices, int v, network& n, vector<vector<int>>& net_adjacency,  const int k, vector<vector<int>>& GDD, vector<vector<int>>& orbit_dict)
{
	if (subgraphVertices.size() == k)
	{
		t_start = clock();
		
		update_GDD(subgraphVertices, n, net_adjacency, GDD, orbit_dict);

		t_end = clock();
		t_update_GDD += t_end - t_start;

		return;
	}

	int w = 0;

	for (int i = 0; i < extensionVertices.size(); i++)		//remove an arbitrary chosen vertex every time from extensionVertices called w
	{
		w = extensionVertices[i];
		vector<int> extVer = extensionVertices;

		vector<int> subVer = subgraphVertices;
		subVer.push_back(w);
		extVer.erase(extVer.begin() + i);		//remove w from extVer
		removeSmallerThanW(extVer, w);			//remove elements in extVer if they are smaller than w (avoid graphlet repitition)

		//vector<int> neigbors = findNeighborhood(n.graphVector, n.delimiterInfo, subgraphVertices);	//find the neighborhood of the remaining vertices



		//add neighboring vertices of w if they are bigger than v and in exclusive neighborhood of other vertices in subgraphVertices
		vector<bool> if_in_extVer (n.delimiterInfo.size(), FALSE);
		for (int j = 0; j < extVer.size(); j++)
			if_in_extVer[extVer[j]] = TRUE;

		for (int j = n.delimiterInfo[w]; j < n.delimiterInfo[w + 1]; j++)
		{
			if (n.graphVector[j] > v)
			{
				bool a = FALSE;
				//go over a loop to see if n.graphVector[j] is a neigbor of subgraphVertices
				for (int p = 0; p < subgraphVertices.size(); p++)
					for (int q = n.delimiterInfo[subgraphVertices[p]]; q < n.delimiterInfo[subgraphVertices[p] + 1]; q++)
						if (n.graphVector[j] == n.graphVector[q])
						{
							a = TRUE;
							break;
						}

				
				if(!a)
					if_in_extVer[n.graphVector[j]] = TRUE;
			}		
		}
		int count = 0;
		for (int j = 0; j < if_in_extVer.size(); j++)
			if (if_in_extVer[j] == TRUE)
				count++;
		vector<int> extVer2(count);
		count = 0;
		for (int j = 0; j < if_in_extVer.size(); j++)
			if (if_in_extVer[j] == TRUE){
				extVer2[count] = j;
				count++;
			}
	


		extendSubgraph(subVer, extVer2, v, n, net_adjacency, k, GDD, orbit_dict);
	}
	return;
}

//stores neigbors for each vertex
vector<int> vertexNeigbors(vector<int>& graph, vector<int>& delimiter, int vertex)
{
	vector<int> extension;
	for (int i = delimiter[vertex]; i < delimiter[vertex + 1]; i++)
	{
		if (graph[i] > vertex)
			extension.push_back(graph[i]);
	}
	return extension;
}

//remove elements in inputVector if they are smaller than w 
void removeSmallerThanW(vector<int>& inputVector, int & w)
{
	for (int i = 0; i < inputVector.size(); i++)
	{
		if (inputVector[i] < w)
		{
			inputVector.erase(inputVector.begin() + i);
			i--;
		}
	}
}


//find the neighborhood of the remaining vertices in ExtensionVertices(EV)
vector<int> findNeighborhood(vector<int>& graph, vector<int>& delimiter, vector<int>& EV)
{
	int max_size=0;
	for (int i = 0; i < EV.size(); i++)
		max_size += (delimiter[EV[i] + 1] - delimiter[EV[i]]);
	vector<int> neighborhood;
	neighborhood.reserve(max_size);

	for (int i = 0; i < EV.size(); i++)
		for (int j = delimiter[EV[i]]; j < delimiter[EV[i] + 1]; j++)
		{
			neighborhood.push_back(graph[j]);
		}
	return neighborhood;


}


//check if input vector contains an element
bool ifContains(vector<int>& input, int element)
{
	//input.push_back(-1);
	bool found = find(input.begin(), input.end(), element) != input.end();
	return found;
}

//outputs sparse adjacency matrix
vector<int> adjacencyFromGraphlets(vector<int>& graph, vector<int>& delimiter, vector<vector<int>> & graphlets)
{
	vector<int> graphletsAdjacency;
	for (int i = 0; i < graphlets.size(); i++)
	{
		int bit = 0;
		int ifEdge;
		for (int j = 1; j < graphlets[0].size();j++)
		{
			for (int k = 0; k < j; k++)
			{
				bit = bit << 1;
				ifEdge = ifConnected(graphlets[i][j], graphlets[i][k], graph, delimiter);
				bit = bit | ifEdge;
			}
		}
		graphletsAdjacency.push_back(bit);
	}
	return graphletsAdjacency;
}

//outputs nonsparse adjacency matrix
vector<unsigned int> adjacencyNSFromGraphlets(vector<vector<int>>& net_adjacency, vector<int> & graphlet)
{
	int s = graphlet.size();
	unsigned int a = 1 << (sizeof(unsigned int)*8 -1);
	vector<unsigned int> graphletsAdjacency(s);
		for (int j = 0; j < s; j++)
		{
			unsigned int bit = 0;
			int ifEdge;
			for (int k = s-1; k >=0 ; k--)
			{
				bit = bit >> 1;
				ifEdge = ifConnected_fast(graphlet[j], graphlet[k], net_adjacency);
				bit = bit | (ifEdge == 1 ? a : 0);
			}
			graphletsAdjacency[j] = bit;
		}
	return graphletsAdjacency;
}

//check if two vertices are connected in a graph as network
int ifConnected(int vertex1, int vertex2, vector<int>& graph, vector<int>& delimiter)
{
	if (find(graph.begin() + delimiter[vertex1], graph.begin() + delimiter[vertex1 + 1] , vertex2) != (graph.begin() + delimiter[vertex1 + 1] ))
		return 1;
	return 0;

}

//check if two vertices are connected in a graph as adjacency
int ifConnected_fast(int vertex1, int vertex2, vector<vector<int>>& net_adjacency)
{
	return (net_adjacency[vertex1][vertex2 / 32] >> (31 - vertex2 % 32)) & 1;
}

nauty_output  ComputeCanon(unsigned int n, unsigned int *NCAdjacency)
{
	nauty_output output;
	int m = 1;
	DEFAULTOPTIONS(options);
	statsblk(stats);
	setword workspace[160 * MAXM];

	options.writeautoms = FALSE;
	options.getcanon = TRUE;

	nauty(NCAdjacency, output.label, output.ptn, NULL, output.orbits, &options, &stats, workspace, 160 * MAXM, m, n, output.canonical); 

	return output;
}

//outputs canonical labling vector of an all canons of all graphlets
nauty_output getCanons(unsigned int n, vector<unsigned int> graphletAdj){
	nauty_output graphletCanon = ComputeCanon(n, &(graphletAdj[0]));
	return graphletCanon;
}

//gets the file name and reads the *txt file and outputs the readout as a network structure
network read_network_from_file(string file_name){

	network net;

	//read network from file
	ifstream readNetwork(file_name, ios::in);

	int v1, v2;	//v1 and v2 store the vertices for each edge entry
	net.Number_of_Edges = 0;			//NumberEdges counts the number of edges as input file streams
	vector<vector<int>> preVector;
	while (readNetwork >> v1)
	{
		readNetwork >> v2;
		net.Number_of_Edges++;
		if (preVector.size() <= max(v1, v2))
		{
			preVector.resize(max(v1, v2) + 1);
		}
		preVector[v1].push_back(v2);	//add corresponding info for an edge to the to vertices
		preVector[v2].push_back(v1);
	}
	readNetwork.close();
	//

	net.graphVector.resize(net.Number_of_Edges * 2);		//graphVector is the finall vector used to represent the network
	net.delimiterInfo.resize(preVector.size() + 1);			//the array used to store the information on where in graphVector each array starts
	net.delimiterInfo[0] = 0;
	int current = 0;								//to indicate where in the graphVector we are
	for (int i = 0; i < preVector.size(); i++)
	{
		net.delimiterInfo[i + 1] = net.delimiterInfo[i] + preVector[i].size();
		for (int j = 0; j < preVector[i].size(); j++)
		{
			net.graphVector[current] = preVector[i][j];
			current++;
		}
	}
	return net;
}

//reads the dictionary of orbits from a binary file
vector<vector<int>> read_dictionary_from_file(string dict_name){

	vector<vector<int>> dictionary (1024, vector<int>(5, -1));
	ifstream inputFileDict(dict_name, ios::in | ios::binary);
	for (int i = 0; i < 1024; i++)
		for (int j = 0; j < 5; j++)
			inputFileDict.read((char*)&dictionary[i][j], sizeof(1));
	inputFileDict.close();
	return dictionary;
}

//convert a reduntand array of sets to a nonreduntant adjacency bitset.
unsigned int convert_Redundant_Adjacency_to_Bitset(unsigned int* input, int k){
	unsigned int output = 0;  //a bit set representing the adjacency matrix
	const unsigned int a = 1 << (sizeof(unsigned int) * 8 - 1);
	unsigned int b;

	//fill adjacency matrix information into the bitset
	int count = 0;
	for (int i = 1; i < k; i++)
		for (int j = 0; j < i; j++){
			output = output << 1;
			count++;
			b = input[i] & a;
			output = output | (b != 0 ? 1 : 0);
			input[i] = input[i] << 1;
		}
	output = output << (32-count);
	return output;
}



//compute Graphlet Degree Distribution
void update_GDD(vector<int>& graphlet, network& n, vector<vector<int>>& net_adjacency, vector<vector<int>>& GDD, vector<vector<int>>& orbit_dict){
	
	int k = graphlet.size();

	t_start = clock();
	vector<unsigned int> graphletAdjacency = adjacencyNSFromGraphlets(net_adjacency, graphlet);
	t_end = clock();
	t_adjacencyNSFromGraphlets += t_end - t_start;

	t_start = clock();
	nauty_output graphlet_canon = getCanons(k, graphletAdjacency);
	t_end = clock();
	t_getCanons += t_end - t_start;

	t_start = clock();
	unsigned int canon_bitset = convert_Redundant_Adjacency_to_Bitset(graphlet_canon.canonical, k);
	canon_bitset = canon_bitset >> 22;
	int p, orb;
	for (int i = 0; i < k; i++){
		p = find_value(graphlet_canon.label, i, k);
		orb = orbit_dict[canon_bitset][p];
		GDD[graphlet[i]][orb]++;
	}
	t_end = clock();
	t_update_loop += t_end - t_start;

	free(graphlet_canon.canonical); 

}

int find_value(int* search, int number, int size){
	for (int i = 0; i < size; i++)
		if (search[i] == number)
			return i;
	return -1;
}

void write_in_file(vector<vector<int>> gdd){
	ofstream answers("answers.out", ios::out | ios::binary);
	for (int i = 0; i < gdd.size(); i++)
		for (int j = 0; j < 73; j++)
			answers.write(reinterpret_cast<const char*>(&gdd[i][j]), sizeof(gdd[i][j]));
	answers.close();
}

bool check_answers(vector<vector<int>>& gdd){
	vector<vector<int>> ans(gdd.size(), vector<int>(73, 0));
	ifstream inputFile("answers.out", ios::in | ios::binary);
	for (int i = 0; i < gdd.size(); i++)
		for (int j = 0; j < 73; j++)
			inputFile.read(reinterpret_cast<char*>(&ans[i][j]), sizeof(1));
	inputFile.close();
	for (int i = 0; i < gdd.size(); i++)
		if (gdd[i] != ans[i])
			return FALSE;
	return TRUE;
}
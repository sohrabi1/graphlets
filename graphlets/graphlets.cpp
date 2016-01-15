// graphlets.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include "methods.h"


using namespace std;



int main()
{

	vector<vector<int>> dictionary = read_dictionary_from_file("dictionary.bin");
	network net = read_network_from_file("network.txt");

	vector<vector<int>> gdd(net.delimiterInfo.size(), vector<int>(73,0));
	enumerateSubgraphs(net, 3, gdd, dictionary);

	return 0;
}



void enumerateSubgraphs(network& n, const int k, vector<vector<int>>& GDD, vector<vector<int>>& orbit_dict)
{
	const int nVertices = n.delimiterInfo.size() - 1;		//number of vertices in the graph
	vector<vector<int>> answer;					//stores graphlets info

	for (int i = 0; i < nVertices; i++)
	{
		vector<int> extensionV = vertexNeigbors(n.graphVector, n.delimiterInfo, i);	//stores neigbors for each vertex
		vector<int> subgraphV(1, i);
		extendSubgraph(subgraphV, extensionV, i, n, k, GDD, orbit_dict);
	}

}

void extendSubgraph(vector<int> subgraphVertices, vector<int> extensionVertices, int v, network& n, const int k, vector<vector<int>>& GDD, vector<vector<int>>& orbit_dict)
{
	if (subgraphVertices.size() == k)
	{
		update_GDD(subgraphVertices, n, GDD, orbit_dict);
		return;
	}

	int w = 0;

	for (int i = 0; i < extensionVertices.size(); i++)		//remove an arbitrary chosen vertex every time from extensionVertices called w
	{
		w = extensionVertices[i];
		vector<int> extVer;
		extVer = extensionVertices;

		vector<int> subVer;
		subVer = subgraphVertices;
		subVer.push_back(w);
		extVer.erase(extVer.begin() + i);		//remove w from extVer
		removeSmallerThanW(extVer, w);			//remove elements in extVer if they are smaller than w (avoid graphlet repitition)

		vector<int> neigbors = findNeighborhood(n.graphVector, n.delimiterInfo, extVer);	//find the neighborhood of the remaining vertices


		//add neighboring vertices of w if they are bigger than v and in exclusive neighborhood of other vertices in extensionVertices
		for (int j = n.delimiterInfo[w]; j < n.delimiterInfo[w + 1]; j++)
		{
			if ((n.graphVector[j] > w) && (!ifContains(neigbors, n.graphVector[j])) && (!ifContains(extVer, n.graphVector[j])))
				extVer.push_back(n.graphVector[j]);
		}

		extendSubgraph(subVer, extVer, v, n, k, GDD, orbit_dict);
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
	vector<int> neighborhood;

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
	return find(input.begin(), input.end(), element) != input.end();
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
vector<unsigned int> adjacencyNSFromGraphlets(network& n, vector<int> & graphlet)
{
	unsigned int a = 1 << (sizeof(unsigned int)*8 -1);
	vector<unsigned int> graphletsAdjacency;
		for (int j = 0; j < graphlet.size(); j++)
		{
			unsigned int bit = 0;
			int ifEdge;
			for (int k = graphlet.size()-1; k >=0 ; k--)
			{
				bit = bit >> 1;
				ifEdge = ifConnected(graphlet[j], graphlet[k], n.graphVector, n.delimiterInfo);
				bit = bit | (ifEdge == 1 ? a : 0);
			}
			graphletsAdjacency.push_back(bit);
		}
	return graphletsAdjacency;
}

//check if two vertices are connected in a graph
int ifConnected(int vertex1, int vertex2, vector<int>& graph, vector<int>& delimiter)
{
	if (find(graph.begin() + delimiter[vertex1], graph.begin() + delimiter[vertex1 + 1] , vertex2) != (graph.begin() + delimiter[vertex1 + 1] ))
		return 1;
	return 0;

}


nauty_output  ComputeCanon(unsigned int n, unsigned int *NCAdjacency)
{
	nauty_output output;
	int m = 1;
	DEFAULTOPTIONS(options);
	statsblk(stats);
	setword workspace[160 * MAXM];
	set *gv;

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
	unsigned int b = 0;

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
void update_GDD(vector<int> graphlet, network& n, vector<vector<int>>& GDD, vector<vector<int>>& orbit_dict){
	int k = graphlet.size();

	vector<unsigned int> graphletAdjacency = adjacencyNSFromGraphlets(n, graphlet);
	nauty_output graphlet_canon = getCanons(k, graphletAdjacency);
	unsigned int canon_bitset = convert_Redundant_Adjacency_to_Bitset(graphlet_canon.canonical, k);
		canon_bitset = canon_bitset >> 22;
		int p, orb;
		for (int i = 0; i < k; i++){
			p = find_value(graphlet_canon.label, i, k);
			orb = orbit_dict[canon_bitset][p];
			GDD[graphlet[i]][orb]++;
		}
}

int find_value(int* search, int number, int size){
	for (int i = 0; i < size; i++)
		if (search[i] == number)
			return i;
	return -1;
}
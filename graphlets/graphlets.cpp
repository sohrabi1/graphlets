// graphlets.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "methods.h"


using namespace std;



int main()
{
	vector<vector<int>> dictionary(1024, vector<int>(5, -1));
	ifstream inputFileDict("dictionary.bin", ios::in | ios::binary);
	for (int i = 0; i < 1024; i++)
		for (int j = 0; j < 5; j++)
			inputFileDict.read((char*)&dictionary[i][j], sizeof(1));
	inputFileDict.close();

	vector<vector<int>> preVector;

	//read network from file
	ifstream readNetwork("network.txt", ios::in);
	int v1, v2, NumberEdges;	//v1 and v2 store the vertices for each edge entry
	NumberEdges = 0;			//NumberEdges counts the number of edges as input file streams

	while (readNetwork >> v1)
	{
		readNetwork >> v2;
		NumberEdges++;
		if (preVector.size() <= max(v1, v2))
		{
			preVector.resize(max(v1, v2)+1);
		}
		preVector[v1].push_back(v2);	//add corresponding info for an edge to the to vertices
		preVector[v2].push_back(v1);
	}
	readNetwork.close();
	//

	vector<int> graphVector(NumberEdges * 2);		//graphVector is the finall vector used to represent the network
	vector<int> delimiterInfo(preVector.size()+1);	//the array used to store the information on where in graphVector each array starts
	delimiterInfo[0] = 0;
	int current = 0;								//to indicate where in the graphVector we are
	for (int i = 0; i < preVector.size(); i++)
	{
		delimiterInfo[i + 1] = delimiterInfo[i] + preVector[i].size();
		for (int j = 0; j < preVector[i].size(); j++)
		{
			graphVector[current] = preVector[i][j];
			current++;
		}
	}

	vector<vector<int>> graphlets;
	int k = 3;
	graphlets = enumerateSubgraphs(graphVector, delimiterInfo, k);
	vector<vector<unsigned int>> graphletAdjacencies = adjacencyNSFromGraphlets(graphVector, delimiterInfo, graphlets);
	vector<unsigned int *> canons = getCanons(k,graphletAdjacencies);
	unsigned int * canonical = ComputeLabel(k, &(graphletAdjacencies[0])[0]);

	return 0;
}



const vector<vector<int>> enumerateSubgraphs(vector<int>& graph, vector<int>& delimiter, const int k)
{

	const int nVertices = delimiter.size() - 1;		//number of vertices in the graph
	vector<vector<int>> answer;					//stores graphlets info

	for (int i = 0; i < nVertices; i++)
	{
		vector<int> extensionV = vertexNeigbors(graph, delimiter, i);	//stores neigbors for each vertex
		vector<int> subgraphV(1, i);
		extendSubgraph(subgraphV, extensionV, i, graph, delimiter, k, answer);
	}

	return answer;
}

void extendSubgraph(vector<int> subgraphVertices, vector<int> extensionVertices, int v, vector<int>& graph, vector<int>& delimiter, const int k, vector<vector<int>>& answer)
{
	if (subgraphVertices.size() == k)
	{
		answer.push_back(subgraphVertices);
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

		vector<int> neigbors = findNeighborhood(graph, delimiter, extVer);	//find the neighborhood of the remaining vertices


		//add neighboring vertices of w if they are bigger than v and in exclusive neighborhood of other vertices in extensionVertices
		for (int j = delimiter[w]; j < delimiter[w + 1]; j++)
		{
			if ((graph[j] > w) && (!ifContains(neigbors, graph[j])) && (!ifContains(extVer, graph[j])))
				extVer.push_back(graph[j]);
		}

		extendSubgraph(subVer, extVer, v, graph, delimiter, k, answer);
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
vector<vector<unsigned int>> adjacencyNSFromGraphlets(vector<int>& graph, vector<int>& delimiter, vector<vector<int>> & graphlets)
{
	unsigned int a = 1 << (sizeof(unsigned int)*8 -1);
	vector<vector<unsigned int>> graphletsAdjacency(graphlets.size());
	for (int i = 0; i < graphlets.size(); i++)
	{
		for (int j = 0; j < graphlets[0].size(); j++)
		{
			int bit = 0;
			int ifEdge;
			for (int k = 0; k < graphlets[0].size(); k++)
			{
				bit = bit >> 1;
				ifEdge = ifConnected(graphlets[i][j], graphlets[i][k], graph, delimiter);
				bit = bit | (ifEdge == 1 ? a : 0);
			}
			graphletsAdjacency[i].push_back(bit);
		}
		
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

//outputs canonical labling of an adjacency matrix of a graphlet
unsigned int *  ComputeLabel(unsigned int n, unsigned int *adjacencyMatrix)
{
	int m = 1;
	graph *canon = (graph *)malloc(n * sizeof(int) * 2);

	int lab[MAXN], ptn[MAXN], orbits[MAXN];
	DEFAULTOPTIONS(options);
	statsblk(stats);
	setword workspace[160 * MAXM];
	set *gv;

	options.writeautoms = FALSE;
	options.getcanon = TRUE;

	nauty(adjacencyMatrix, lab, ptn, NULL, orbits, &options, &stats, workspace, 160 * MAXM, m, n, canon);

	return canon;
}

//outputs canonical labling vector of an all canons of all graphlets
vector<unsigned int *> getCanons(unsigned int n, vector<vector<unsigned int>> graphletAdj){
	vector<unsigned int *> graphletCanon;
	for (int i = 0; i < graphletAdj.size(); i++)
		graphletCanon.push_back(ComputeLabel(n, &(graphletAdj[0])[0]));
	return graphletCanon;
}
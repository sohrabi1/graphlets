// graphlets.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "stdafx.h"

using namespace std;

const vector<vector<int>>& enumerateSubgraphs(vector<int>, vector<int> , int);
void extendSubgraph(vector<int>, vector<int>, int, vector<int>&, vector<int>&, const int, vector<vector<int>> &);

int main()
{
	vector<vector<int>> preVector(5, vector<int>(1));

	//read network from file
	ifstream readNetwork("network.txt", ios::in);
	int v1, v2, NumberEdges;	//v1 and v2 store the vertices for each edge entry
	NumberEdges = 0;			//NumberEdges counts the number of edges as input file streams

	while (readNetwork >> v1)
	{
		readNetwork >> v2;
		NumberEdges++;
		//if (preVector.size() < max(v1, v2))
		//{
		//	preVector.resize(max(v1, v2), 0);
		//}
		preVector[v1].push_back(v2);	//add corresponding info for an edge to the to vertices
		preVector[v2].push_back(v1);
	}
	//
	for (int i = 0; i < preVector.size(); i++)
	{
		cout << preVector[i][1];
	}

	vector<int> graphVector(NumberEdges * 2);		//graphVector is the finall vector used to represent the network
	vector<int> delimiterInfo(preVector.size()+1);	//the array used to store the information on where in graphVector each array starts
	delimiterInfo[0] = 0;
	int current = 0;								//to indicate where in the graphVector we are
	for (int i = 0; i < preVector.size(); i++)
	{
		delimiterInfo[i + 1] = delimiterInfo[i] + preVector[i].size()-1;
		for (int j = 1; j < preVector[i].size(); j++)
		{
			graphVector[current] = preVector[i][j];
			current++;
		}
	}

	return 0;
}


const vector<vector<int>>& enumerateSubgraphs(vector<int>& graph, vector<int>& delimiter, const int k)
{
	const int nEdges = graph.size() / 2;			//number of edges in the graph
	const int nVertices = delimiter.size() - 1;		//number of vertices in the graph
	vector<int> extensionV(1);						//stores neigbors for each vertex
	vector<vector<int>> answer(1);
	for (int i = 0; i < nVertices; i++)
	{
		for (int j = 0; j < delimiter[i + 1] - delimiter[i]; j++)
		{
			if (graph[i + j] > i)
				extensionV.push_back(graph[i + j]);
		}
		vector<int> subgraphV(1);
		subgraphV[0] = i;
		extendSubgraph(subgraphV, extensionV, i, graph, delimiter, k, answer);
	}

	return answer;
}

void extendSubgraph(vector<int> subgraphVertices, vector<int> extensionVertices, int v, vector<int>& graph, vector<int>& delimiter, const int k, vector<vector<int>>& answer)
{
	if (subgraphVertices.size() == k)
		answer.push_back( subgraphVertices);
	int w = 0;
	vector<int> neighborhood(1);
	neighborhood[1] = -1;

	while (extensionVertices.size() != 0)		//if not empty
	{
		for (int i = 0; i < extensionVertices.size(); i++)		//remove an arbitrary chosen vertex every time from extensionVertices called w
		{
			w = extensionVertices[i];
			extensionVertices.erase(extensionVertices.begin() + i);

			for (int j = 0; j < extensionVertices.size(); j++)	//find the neighborhood of the remaining vertices
				for (int k = delimiter[extensionVertices[j]]; k < delimiter[extensionVertices[j] + 1]; k++)
				{
					neighborhood.push_back(graph[k]);
				}

																//add neighboring vertices of w if they are bigger than v and in exclusive neighborhood of other vertices in extensionVertices
			for (int j = delimiter[w]; j < delimiter[w + 1]; j++)
			{
				if ((graph[w] > v) && (find(neighborhood.begin(), neighborhood.end(), graph[w]) == neighborhood.end()))
					extensionVertices.push_back(graph[w]);
			}
			subgraphVertices.push_back(w);
			extendSubgraph(subgraphVertices, extensionVertices, v, graph, delimiter, k, answer);
		}
	}
}
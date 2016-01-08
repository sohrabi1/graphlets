// graphlets.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "stdafx.h"

using namespace std;

const vector<vector<int>>& enumerateSubgraphs(vector<int> &, vector<int> & , int);
void extendSubgraph(vector<int>, vector<int>, int, vector<int>&, vector<int>&, const int, vector<vector<int>> &);

int main()
{
	vector<vector<int>> preVector(4, vector<int>(1));

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

	vector<vector<int>> graphlets;
	graphlets = enumerateSubgraphs(graphVector, delimiterInfo, 3);

	return 0;
}


const vector<vector<int>>& enumerateSubgraphs(vector<int>& graph, vector<int>& delimiter, const int k)
{
	const int nEdges = graph.size() / 2;			//number of edges in the graph
	const int nVertices = delimiter.size() - 1;		//number of vertices in the graph
	
	vector<vector<int>> answer;					//stores the graphlet info
	for (int i = 0; i < nVertices; i++)
	{
		vector<int> extensionV;						//stores neigbors for each vertex
		for (int j = delimiter[i]; j < delimiter[i + 1]; j++)
		{
			if (graph[j] > i)
				extensionV.push_back(graph[j]);
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
	{
		answer.push_back(subgraphVertices);
		return;
	}
	int w = 0;

	while (extensionVertices.size() != 0)		//if not empty
	{
		for (int i = 0; i < extensionVertices.size(); i++)		//remove an arbitrary chosen vertex every time from extensionVertices called w
		{
			w = extensionVertices[i];
			vector<int> extVer;
			extVer = extensionVertices;

			vector<int> subVer;
			subVer = subgraphVertices;

			extVer.erase(extVer.begin() + i);					//remove w from extVer

			for (int j = 0; j < extVer.size(); j++)				//remove elements in extVer if they are smaller than w (avoid graphlet repitition)
			{
				if (extVer[j]<w)
					extVer.erase(extVer.begin() + j);
			}
			vector<int> neighborhood;

			for (int j = 0; j < extVer.size(); j++)	//find the neighborhood of the remaining vertices
				for (int k = delimiter[extVer[j]]; k < delimiter[extVer[j] + 1]; k++)
				{
					neighborhood.push_back(graph[k]);
				}
			

																//add neighboring vertices of w if they are bigger than v and in exclusive neighborhood of other vertices in extensionVertices
			for (int j = delimiter[w]; j < delimiter[w + 1]; j++)
			{
				if ((graph[j] > w) && (find(neighborhood.begin(), neighborhood.end(), graph[j]) == neighborhood.end()) && (find(extVer.begin(), extVer.end(), graph[j]) == extVer.end()))
					extVer.push_back(graph[j]);
			}
			
			subVer.push_back(w);
			
			extendSubgraph(subVer, extVer, v, graph, delimiter, k, answer);
		}
		return;
	}
	return;
}
// graphlets.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "stdafx.h"

using namespace std;


int main()
{
	vector<vector<int>> preVector(4, vector<int>(1));

	//read network from file
	ifstream readNetwork("graphlets\\network.txt", ios::in);
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

	for (int i = 0; i < preVector.size(); i++)
	{
		cout << preVector[i][1];
	}

	vector<int> graphVector(NumberEdges * 2);		//graphVector is the finall vector used to represent the network
	vector<int> delimiterInfo(preVector.size()+1);	//the array used to store the information on where in graphVector each array ends
	delimiterInfo[0] = 0;
	int current = 0;
	for (int i = 0; i < preVector.size(); i++)
	{
		delimiterInfo[i + 1] = delimiterInfo[i] + preVector[i].size();
		for (int j = 0; j < preVector[i].size(); j++)
		{
			graphVector[current] = preVector[i][j];
			current++;
		}
	}
	return 0;
}


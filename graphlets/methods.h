// methods.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#include <vector>
#include <algorithm>

using namespace std;

const vector<vector<int>> enumerateSubgraphs(vector<int> &, vector<int> & , int);
void extendSubgraph(vector<int>, vector<int>, int, vector<int>&, vector<int>&, const int, vector<vector<int>> &);
void removeSmallerThanW(vector<int>&, int & );
vector<int> findNeighborhood(vector<int>&, vector<int>&, vector<int>&);
bool ifContains(vector<int>&, int);
vector<int> vertexNeigbors(vector<int>&, vector<int>&, int ); 


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

//remove elements in inputVector if they are smaller than w (avoid graphlet repitition)
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


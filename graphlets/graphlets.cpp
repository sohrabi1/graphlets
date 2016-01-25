// graphlets.cpp : Defines the entry point for the console application.
//
// graphlets.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <ctime>
#include "methods.h"




using namespace std;
//using web::json;

int main()
{
	clock_t startTime;
	clock_t endTime;
	startTime = clock();

	vector<vector<int>> dictionary = read_dictionary_from_file("dictionary.bin");
	network net = read_network_from_file("exmpl100.in");
	vector<vector<int>> network_adjacency = sparse_network(net);

	vector<vector<int>> gdd(net.Number_of_Vertices, vector<int>(73, 0));
	enumerateSubgraphs(net, network_adjacency, 5, gdd, dictionary);
	enumerateSubgraphs(net, network_adjacency, 4, gdd, dictionary);
	enumerateSubgraphs(net, network_adjacency, 3, gdd, dictionary);
	enumerateSubgraphs(net, network_adjacency, 2, gdd, dictionary);
	endTime = clock();
	double t = (endTime - startTime) / CLOCKS_PER_SEC;
	bool check = check_answers(gdd);
	cout << check << endl << t << endl;

	free(net.gVector);
	free(net.dInfo);

	return 0;
}


vector<vector<int>> sparse_network(network& net)
{
	int r = (net.Number_of_Vertices) % 32;
	int fraction = (net.Number_of_Vertices) / 32;
	vector<vector<int>> network_adjacency(net.Number_of_Vertices, vector<int>((fraction + 1), 0));
	int if_connected;
	for (int i = 0; i < net.Number_of_Vertices; i++)
		for (int j = 0; j < net.Number_of_Vertices; j++)
		{
			if_connected = ifConnected(i, j, net);
			network_adjacency[i][j / 32] = network_adjacency[i][j / 32] << 1;
			network_adjacency[i][j / 32] = network_adjacency[i][j / 32] | if_connected;
		}
	for (int i = 0; i < net.Number_of_Vertices; i++)
		network_adjacency[i][fraction] = network_adjacency[i][fraction] << (32 - r);
	return network_adjacency;
}

void enumerateSubgraphs(network& n, vector<vector<int>>& net_adjacency, const int k, vector<vector<int>>& GDD, vector<vector<int>>& orbit_dict)
{
	const int number_of_servers = 3;
	int itr = 0;
	for (int i = 0; i < number_of_servers; i++)
	{
		vector<int> initial_vertices((n.Number_of_Vertices % number_of_servers) ? n.Number_of_Vertices / number_of_servers + 1 : n.Number_of_Vertices / number_of_servers);
		for (int j = 0; j <initial_vertices.size(); j++)
		{
			if (itr<n.Number_of_Vertices)
			{
				initial_vertices[j]=itr;
				itr++;
			}
			else {
				initial_vertices.resize(j);
				break;
			}
		}
		vector<vector<int>> gdd_temp = call_extendSubgraph(initial_vertices, n, net_adjacency, k, orbit_dict);
		update_gdd_temp(gdd_temp, GDD);
	}

}

vector<vector<int>> call_extendSubgraph(vector<int> initial_vertices, network& n, vector<vector<int>>& net_adjacency, const int k, vector<vector<int>>& orbit_dict)
{
	vector<vector<int>> gdd(n.Number_of_Vertices, vector<int>(73, 0));
	for (int i = 0; i < initial_vertices.size(); i++)
	{
		vector<int> extensionV = vertexNeigbors(n.gVector, n.dInfo, initial_vertices[i]);	//stores neigbors for each vertex
		vector<int> subgraphV(1, initial_vertices[i]);
		extendSubgraph(subgraphV, extensionV, initial_vertices[i], n, net_adjacency, k, gdd, orbit_dict);
	}
	return gdd;
}


void update_gdd_temp(vector<vector<int>>& gdd_temp, vector<vector<int>>& GDD)
{
	for (int i = 0; i < GDD.size();i++)
		for (int j=0; j < 73; j++)
			GDD[i][j] += gdd_temp[i][j];
}


void extendSubgraph(vector<int>& subgraphVertices, vector<int>& extensionVertices, int v, network& n, vector<vector<int>>& net_adjacency, const int k, vector<vector<int>>& GDD, vector<vector<int>>& orbit_dict)
{
	if (subgraphVertices.size() == k)
	{
		update_GDD(subgraphVertices, n, net_adjacency, GDD, orbit_dict);
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

		vector<int> neigbors = findNeighborhood(n, subgraphVertices);	//find the neighborhood of the remaining vertices

		//add neighboring vertices of w if they are bigger than v and in exclusive neighborhood of other vertices in subgraphVertices

		for (int j = n.dInfo[w]; j < n.dInfo[w + 1]; j++)
		{
			if ((n.gVector[j] > v) && (!ifContains(neigbors, n.gVector[j])) && (!ifContains(extVer, n.gVector[j])))
				extVer.push_back(n.gVector[j]);
		}

		extendSubgraph(subVer, extVer, v, n, net_adjacency, k, GDD, orbit_dict);
	}
	return;
}

//stores neigbors for each vertex
vector<int> vertexNeigbors(int* graph, int* delimiter, int vertex)
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
	int* ptr;
	for (int i = 0; i < inputVector.size(); i++)
	{
		ptr = &inputVector[0];
		if (*(ptr+i) < w)
		{
			inputVector.erase(inputVector.begin() + i);
			i--;
		}
	}
}


//find the neighborhood of the remaining vertices in ExtensionVertices(EV)
vector<int> findNeighborhood(network& n, vector<int>& EV)
{
	vector<int> neighborhood;
	int* m = &EV[0];
	for (int i = 0; i < EV.size(); i++)
		for (int j = n.dInfo[m[i]]; j < n.dInfo[m[i] + 1]; j++)
		{
			neighborhood.push_back(n.gVector[j]);
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

/*outputs sparse adjacency matrix
vector<int> adjacencyFromGraphlets(vector<vector<int>>& net_adjacency, vector<vector<int>> & graphlets)
{
	vector<int> graphletsAdjacency;
	for (int i = 0; i < graphlets.size(); i++)
	{
		int bit = 0;
		int ifEdge;
		for (int j = 1; j < graphlets[0].size(); j++)
		{
			for (int k = 0; k < j; k++)
			{
				bit = bit << 1;
				ifEdge = ifConnected_fast(graphlets[i][j], graphlets[i][k], net_adjacency);
				bit = bit | ifEdge;
			}
		}
		graphletsAdjacency.push_back(bit);
	}
	return graphletsAdjacency;
}
*/
//outputs nonsparse adjacency matrix
vector<unsigned int> adjacencyNSFromGraphlets(vector<vector<int>>& net_adjacency, vector<int> & graphlet)
{
	int s = graphlet.size();
	unsigned int a = 1 << (sizeof(unsigned int) * 8 - 1);
	vector<unsigned int> graphletsAdjacency(s);
	for (int j = 0; j < s; j++)
	{
		unsigned int bit = 0;
		int ifEdge;
		for (int k = s - 1; k >= 0; k--)
		{
			bit = bit >> 1;
			ifEdge = ifConnected_fast(graphlet[j], graphlet[k], net_adjacency);
			bit = bit | (ifEdge == 1 ? a : 0);
		}
		graphletsAdjacency[j] = bit;
	}
	return graphletsAdjacency;
}

int ifConnected(int vertex1, int vertex2, network& n)
{
	for (int j = n.dInfo[vertex1]; j < n.dInfo[vertex1 + 1]; j++)
		if (n.gVector[j] == vertex2)
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

	nauty(NCAdjacency, output.label, output.ptn, nullptr, output.orbits, &options, &stats, workspace, 160 * MAXM, m, n, output.canonical);

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

	net.Number_of_Vertices = preVector.size();
	//net.graphVector.resize(net.Number_of_Edges * 2);		//graphVector is the finall vector used to represent the network
	//net.delimiterInfo.resize(preVector.size() + 1);			//the array used to store the information on where in graphVector each array starts
	net.dInfo[0] = 0;
	int current = 0;								//to indicate where in the graphVector we are
	for (int i = 0; i < preVector.size(); i++)
	{
		net.dInfo[i + 1] = net.dInfo[i] + preVector[i].size();
		for (int j = 0; j < preVector[i].size(); j++)
		{
			net.gVector[current] = preVector[i][j];
			current++;
		}
	}


	return net;
}

//reads the dictionary of orbits from a binary file
vector<vector<int>> read_dictionary_from_file(string dict_name){

	vector<vector<int>> dictionary(1024, vector<int>(5, -1));
	ifstream inputFileDict(dict_name, ios::in | ios::binary);
	for (int i = 0; i < 1024; i++)
		for (int j = 0; j < 5; j++)
			inputFileDict.read(reinterpret_cast<char*>(&dictionary[i][j]), sizeof(1));
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
	output = output << (32 - count);
	return output;
}



//compute Graphlet Degree Distribution
void update_GDD(vector<int>& graphlet, network& n, vector<vector<int>>& net_adjacency, vector<vector<int>>& GDD, vector<vector<int>>& orbit_dict){

	int k = graphlet.size();
	vector<unsigned int> graphletAdjacency = adjacencyNSFromGraphlets(net_adjacency, graphlet);
	nauty_output graphlet_canon = getCanons(k, graphletAdjacency);

	unsigned int canon_bitset = convert_Redundant_Adjacency_to_Bitset(graphlet_canon.canonical, k);
	canon_bitset = canon_bitset >> 22;
	int p, orb;
	for (int i = 0; i < k; i++){
		p = find_value(graphlet_canon.label, i, k);
		orb = orbit_dict[canon_bitset][p];
		GDD[graphlet[i]][orb]++;
	}

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

/*vector<bool> if_in_extVer (n.delimiterInfo.size(), FALSE);
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



extendSubgraph(subVer, extVer2, v, n, net_adjacency, k, GDD, orbit_dict);  */

vector<vector<int>> connect_server(network net)
{
	vector<vector<int>> gdd;


	return gdd;
}

//creates GDD object from jason values 
vector<vector<int>> GDD_generator(web::json::array jsonValue, network net)
{
	vector<vector<int>> igdd(73);
	int i = 0;
	for (auto itr = jsonValue.begin(); itr!=jsonValue.end(); ++itr)
		{	
			igdd[i/73][i%73] = jsonValue.at(i).as_integer();
			i++;
		}

	return igdd;
}

//creates jason object from the network
web::json::value network_to_jason(network net){
	//web::json::value json_network = web::json::value::array();

	std::vector<web::json::value> arrayNet(net.Number_of_Vertices + net.Number_of_Edges * 2 + 3);


	utility::stringstream_t ss1;
	ss1 << net.Number_of_Edges;
	web::json::value numberEdge = web::json::value::parse(ss1);	
	arrayNet.push_back(numberEdge);

	utility::stringstream_t ss2;
	ss1 << net.Number_of_Vertices;
	web::json::value numberVer = web::json::value::parse(ss2);
	arrayNet.push_back(numberVer);

	for (int i = 0; i < (net.Number_of_Edges * 2); i++)
	{
		utility::stringstream_t ssi;
		ssi << net.gVector[i];
		web::json::value inserted = web::json::value::parse(ssi);
		arrayNet.push_back(inserted);
	}

	for (int i = 0; i <= (net.Number_of_Vertices); i++)
	{
		utility::stringstream_t ssi;
		ssi << net.dInfo[i];
		web::json::value inserted = web::json::value::parse(ssi);
		arrayNet.push_back(inserted);
	}

	web::json::value njson = web::json::value::array(arrayNet);
	return njson;

}

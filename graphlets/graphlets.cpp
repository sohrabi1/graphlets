// graphlets.cpp : Defines the entry point for the console application.
//
// graphlets.cpp : Defines the entry point for the console application.
//

//#include <algorithm>

#include "methods.h"



int main()
{
	clock_t startTime;
	clock_t endTime;
	startTime = clock();

	

	try
	{
		http_listener _listener(L"http://localhost:8010");
		_listener.support(methods::POST, handle_post);
		_listener.open()
			.wait();
		while (true)
		{
			this_thread::sleep_for(chrono::milliseconds(2000));
		}
	}

	catch (exception const & e)
	{
		wcout << e.what() << endl;
	}


	endTime = clock();
	double t = (endTime - startTime) / CLOCKS_PER_SEC;
	cout << endl << t << endl;


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

		//call to server should be made here
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


json::value  convert_gdd_to_json(vector<vector<int>> gdd_vector)
{
	vector<json::value> arrayGDD(gdd_vector.size() * 73);

	for (int i = 0; i < (gdd_vector.size() * 73); i++)
	{
		utility::stringstream_t ssi;
		ssi << gdd_vector[i / 73][i % 73];
		arrayGDD[i] = web::json::value::parse(ssi);
	}

	web::json::value gdd_json = web::json::value::array(arrayGDD);


	return gdd_json;
}




void handle_post(http_request message)
	{
		json::value input = message.extract_json().get();
		json::value gdd = processRequest(input);
		respond(message, status_codes::OK, gdd);
	};

	//gets the json value sent by client and processes the gdd request. retuen value is the gdd response packed as json value
json::value  processRequest(json::value json_input)
	{
		network n;
		int graphlet_size = json_input.at(0).as_integer();
		int start_node = json_input.at(1).as_integer();
		int end_node = json_input.at(2).as_integer();
		n.Number_of_Edges = json_input.at(3).as_integer();
		n.Number_of_Vertices = json_input.at(4).as_integer();

		for (int i = 0; i < n.Number_of_Edges * 2; i++)
			n.gVector[i] = json_input.at(i + 5).as_integer();

		for (int i = 0; i < n.Number_of_Vertices + 1; i++)
			n.dInfo[i] = json_input.at(i + 5 + n.Number_of_Edges * 2).as_integer();

		vector<vector<int>> network_adjacency = sparse_network(n);
		vector<vector<int>> dictionary = read_dictionary_from_file("dictionary.bin");

		vector<int> initial_vertices(end_node - start_node + 1);
		for (int i = 0; i < initial_vertices.size(); i++)
			initial_vertices[i] = start_node + i;

		vector<vector<int>> gdd_temp = call_extendSubgraph(initial_vertices, n, network_adjacency, graphlet_size, dictionary);
		free(n.gVector);
		free(n.dInfo);
		json::value gdds=convert_gdd_to_json(gdd_temp) ;

		return gdds;
	}

	void respond(const http_request& request, const status_code& status, const json::value& response) {
		request.reply(status, response);
	}	


	//vector<vector<int>> connect_server(network net)
	//{
	//	vector<vector<int>> gdd;
	//	pplx::task<void> a = Get_GDD(net, gdd);

	//	return gdd;
	//}
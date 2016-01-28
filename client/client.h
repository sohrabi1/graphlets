#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <ctime>
#include <stdlib.h>
#include "http_client.h"
#include "ppltasks.h"
#include "json.h"



using namespace std;
using namespace web::json;
using namespace web::http;
using namespace web;


struct network {

	int* gVector = static_cast<int *>(malloc(350 * sizeof(int) * 2));;		//graphVector is the finall vector used to represent the network
	int* dInfo = static_cast<int *>(malloc(150 * sizeof(int)));	//the array used to store the information on where in graphVector each array starts

	int Number_of_Edges;
	int Number_of_Vertices;
};


vector<vector<int>> connect_server(network, int, int, int);
vector<vector<int>> GDD_generator(web::json::value, network);
web::json::value network_to_jason(network&, int, int, int);
pplx::task<void> Get_GDD(json::value, network, vector<vector<int>>&);
network read_network_from_file(string);

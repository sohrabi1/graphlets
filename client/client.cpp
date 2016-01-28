// client.cpp : Defines the entry point for the console application.
//

#include "client.h"
#include "ppltasks.h"

int main()
{
	//clock_t startTime;
	//clock_t endTime;
	//startTime = clock();

	network net = read_network_from_file("exmpl100.in");

	web::json::value js = network_to_jason(net,3,0,99);

	vector<vector<int>> Gdds = connect_server(net, 3, 0, 99);

	//vector<vector<int>> gdd(net.Number_of_Vertices, vector<int>(73, 0));
	//enumerateSubgraphs(net, network_adjacency, 5, gdd, dictionary);
	//enumerateSubgraphs(net, network_adjacency, 4, gdd, dictionary);
	//enumerateSubgraphs(net, network_adjacency, 3, gdd, dictionary);
	//enumerateSubgraphs(net, network_adjacency, 2, gdd, dictionary);
	//endTime = clock();
	//double t = (endTime - startTime) / CLOCKS_PER_SEC;
	//bool check = check_answers(gdd);
	//cout << check << endl << t << endl;

	//free(net.gVector);
	//free(net.dInfo);

	return 0;

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


//creates jason object from the network
web::json::value network_to_jason(network& net, int graphlet_size, int start_vertex, int end_vertex){
	//web::json::value json_network = web::json::value::array();

	std::vector<web::json::value> arrayNet(net.Number_of_Vertices+6+net.Number_of_Edges*2);

	utility::stringstream_t ss1;
	ss1 << graphlet_size;
	arrayNet[0] = web::json::value::parse(ss1);

	utility::stringstream_t ss2;
	ss2 << start_vertex;
	arrayNet[1] = web::json::value::parse(ss2);

	utility::stringstream_t ss3;
	ss3 << end_vertex;
	arrayNet[2] = web::json::value::parse(ss3);

	utility::stringstream_t ss4;
	ss4 << net.Number_of_Edges;
	arrayNet[3] = web::json::value::parse(ss4);

	utility::stringstream_t ss5;
	ss5 << net.Number_of_Vertices;
	arrayNet[4] = web::json::value::parse(ss5);

	for (int i = 0; i < (net.Number_of_Edges * 2); i++)
	{
		utility::stringstream_t ssi;
		ssi << net.gVector[i];
		arrayNet[5+i] = web::json::value::parse(ssi);
	}

	for (int i = 0; i <= (net.Number_of_Vertices); i++)
	{
		utility::stringstream_t ssi;
		ssi << net.dInfo[i];
		arrayNet[5 + net.Number_of_Edges * 2 + i] = web::json::value::parse(ssi);
	}

	web::json::value njson = web::json::value::array(arrayNet);


	return njson;

}


vector<vector<int>> connect_server(network net, int graphlet_size, int start_vertex, int end_vertex)
{
	web::json::value js = network_to_jason(net, graphlet_size, start_vertex, end_vertex);
	vector<vector<int>> gdd;
	pplx::task<void> a = Get_GDD(js, net, gdd);

	return gdd;
}

//creates GDD object from jason values 
vector<vector<int>> GDD_generator(web::json::value jsonValue, network net)
{
	vector<vector<int>> igdd(net.Number_of_Vertices, vector<int>(73));
	for (int i = 0; i<jsonValue.size(); ++i)
	{
		igdd[i / 73][i % 73] = jsonValue.at(i).as_integer();
		i++;
	}

	return igdd;
}


pplx::task<void> Get_GDD(json::value json_request, network net, vector<vector<int>>& gdds)
{

	return pplx::create_task([&]
	{
		web::http::client::http_client client(L"http://localhost:8010");
		return client.request(web::http::methods::GET, L"http://localhost:8010", json_request); //returns An asynchronous operation type task<http_response> that is completed once a response from the request is received.  
	})

		.then([&](web::http::http_response response)
	{
		if (response.status_code() == web::http::status_codes::OK)
		{
			return response.extract_json();		//returns>JSON value from the body of this message.  pplx::task<json::value>
		}

		return pplx::create_task([]{return web::json::value(); });
		//return web::json::value();
	})

		.then([&](web::json::value jsonValue)
	{
		if (jsonValue.is_null())
			return;
		try{
			gdds = GDD_generator(jsonValue, net);
			//writeGdds(gdds);
		}
		catch(http_exception const & e)
	{
		wcout << e.what() << endl;
	}

	});

}




// Build_Dictionary.cpp : Defines the entry point for the console application.
//

#include <vector>
#include<fstream>
#include "Build_Dictionary.h"


using namespace std;




int main()
{
	//make library of all possible adjacency matrices with size 5
	vector<unsigned int> allAdj = mkLibAdj();
	vector<vector<int>> Dict = mkDict();
	ofstream dictFile("dictionary.bin", ios::out | ios::binary);
	for (int i = 0; i < Dict.size(); i++)
		for (int j = 0; j < 5; j++)
			dictFile.write(reinterpret_cast<const char*>(&Dict[i][j]), sizeof(Dict[i][j]));
	dictFile.close();
	///

	return 0;
}

//make library of all possible adjacency matrices with size 5
vector<unsigned int> mkLibAdj(){
	vector<unsigned int> library;
	unsigned int bit = 0;
	int counter = 0;
	for (int j = 0; j < 2; j++)
		mkLibAdjHelper(library, counter, bit, j);

	return library;
}

//mkLibAdj helper function, recursive function
void mkLibAdjHelper(vector<unsigned int> & lib, int n, unsigned int bit, int zo){
	//there are ten nonredundant members in an adj matrix of size 5*5
	//unsigned int bit = 0;
	const unsigned int a = 1 << (sizeof(unsigned int) * 8 - 1);
	bit = bit >> 1;
	bit = bit | (zo == 1 ? a : 0);
	++n;
	if (n < 10)
		for (int i = 0; i < 2; i++)
			mkLibAdjHelper(lib, n, bit, i);
	else lib.push_back(bit);

	return;
}

//convert a nonreduntant adjacency bitset to a reduntand array of sets.
unsigned int* convert_NR_to_Set_array(unsigned int input){
	unsigned int* output = (unsigned int *)malloc(sizeof(unsigned int) * 5); //an array of bit sets presenting the adjacency matrix(the form acceptable to nauty)
	vector<vector<int>> temp(5, vector<int>(5));
	const unsigned int a = 1 << (sizeof(unsigned int) * 8 - 1);
	int b = 0;

	//fill temp with adjacency matrix information
	for (int i = 0; i < 10; i++){
		b = a & input;
		input = input << 1;

		if (b != 0)
			switch (i){
			case 0:
				temp[1][0] = temp[0][1] = 1;
				break;
			case 1:
				temp[2][0] = temp[0][2] = 1;
				break;
			case 2:
				temp[2][1] = temp[1][2] = 1;
				break;
			case 3:
				temp[3][0] = temp[0][3] = 1;
				break;
			case 4:
				temp[3][1] = temp[1][3] = 1;
				break;
			case 5:
				temp[3][2] = temp[2][3] = 1;
				break;
			case 6:
				temp[4][0] = temp[0][4] = 1;
				break;
			case 7:
				temp[4][1] = temp[1][4] = 1;
				break;
			case 8:
				temp[4][2] = temp[2][4] = 1;
				break;
			case 9:
				temp[4][3] = temp[3][4] = 1;
				break;
		}
	}

	//use temp to construct the output(bit set arrays representation of adj matrix)
	for (int i = 0; i < 5; i++){
		output[i] = 0;
		for (int j = 4; j >= 0; j--){
			output[i] = output[i] >> 1;
			output[i] = output[i] | (temp[i][j] == 1 ? a : 0);
		}
	}

	return output;

}

//outputs canonical adjacency matrix of a noncanonical adjacency matrix
nauty_output  ComputeCanonAdj(unsigned int n, unsigned int *NCAdjacencyMatrix)
{
	nauty_output output;
	int m = 1;
	//graph *canon = (graph *)malloc(n * sizeof(int) * 2);

	//int lab[MAXN], ptn[MAXN], orbits[MAXN];
	DEFAULTOPTIONS(options);
	statsblk(stats);
	setword workspace[160 * MAXM];
	set *gv;

	options.writeautoms = FALSE;
	options.getcanon = TRUE;

	nauty(NCAdjacencyMatrix, output.label, output.ptn, NULL, output.orbits, &options, &stats, workspace, 160 * MAXM, m, n, output.canonical);

	return output;
}

//convert a reduntand array of sets to a nonreduntant adjacency bitset.
unsigned int convert_R_Adjacency_to_Bitset(unsigned int* input){
	unsigned int output = 0;  //a bit set representing the adjacency matrix
	const unsigned int a = 1 << (sizeof(unsigned int) * 8 - 1);
	unsigned int b = 0;

	//fill adjacency matrix information into the bitset
	for (int i = 1; i < 5; i++)
		for (int j = 0; j < i; j++){
			output = output << 1;
			b = input[i] & a;
			output = output | (b != 0 ? 1 : 0);
			input[i] = input[i] << 1;
		}
	output = output << 22;
	return output;

}


vector<int> Graphlet_Lib(unsigned int adjacency){
	vector<int> orbits(5, -1);
	switch (adjacency){
	case 4194304:
		orbits = { -1, -1, -1, 0, 0 };
		return orbits;
	case 12582912:
		orbits = { -1, -1, 1, 1, 2 };
		return orbits;
	case 79691776:
		orbits = { -1, -1, 3, 3, 3 };
		return orbits;
	case 146800640:
		orbits = { -1, 4, 4, 5, 5 };
		return orbits;
	case 29360128:
		orbits = { -1, 6, 6, 6, 7 };
		return orbits;
	case 683671552:
		orbits = { -1, 8, 8, 8, 8 };
		return orbits;
	case 96468992:
		orbits = { -1, 9, 10, 10, 11 };
		return orbits;
	case 230686720:
		orbits = { -1, 12, 12, 13, 13 };
		return orbits;
	case 767557632:
		orbits = { -1, 14, 14, 14, 14 };
		return orbits;
	case 360710144:
		orbits = { 15, 15, 17, 16, 16 };
		return orbits;
	case 121634816:
		orbits = { 19, 19, 18, 20, 21 };
		return orbits;
	case 62914560:
		orbits = { 22, 22, 22, 22, 23 };
		return orbits;
	case 364904448:
		orbits = { 24, 24, 25, 26, 26 };
		return orbits;
	case 2243952640:
		orbits = { 27, 28, 29, 29, 30 };
		return orbits;
	case 130023424:
		orbits = { 31, 31, 32, 32, 33 };
		return orbits;
	case 3368026112:
		orbits = { 34, 34, 34, 34, 34 };
		return orbits;
	case 260046848:
		orbits = { 35, 37, 37, 36, 38 };
		return orbits;
	case 264241152:
		orbits = { 39, 40, 40, 41, 42 };
		return orbits;
	case 2277507072:
		orbits = { 43, 43, 43, 43, 44 };
		return orbits;
	case 784334848:
		orbits = { 45, 46, 48, 48, 47 };
		return orbits;
	case 528482304:
		orbits = { 49, 49, 49, 50, 50 };
		return orbits;
	case 2512388096:
		orbits = { 51, 51, 52, 53, 53 };
		return orbits;
	case 532676608:
		orbits = { 54, 54, 54, 55, 55 };
		return orbits;
	case 801112064:
		orbits = { 56, 57, 57, 57, 58 };
		return orbits;
	case 1337982976:
		orbits = { 59, 59, 60, 60, 61 };
		return orbits;
	case 3451912192:
		orbits = { 62, 63, 63, 64, 64 };
		return orbits;
	case 1069547520:
		orbits = { 65, 66, 66, 67, 67 };
		return orbits;
	case 3485466624:
		orbits = { 68, 68, 68, 68, 69 };
		return orbits;
	case 2143289344:
		orbits = { 70, 70, 71, 71, 71 };
		return orbits;
	case 4290772992:
		orbits = { 72, 72, 72, 72, 72 };
		return orbits;




	default:
		return orbits;
	}

}

//fined the assigned number for orbits for a given graph as bitset
vector<int> Find_Orbits(unsigned int adjacencyNR) {
	vector<int> orbits(5);
	unsigned int* outAdj = convert_NR_to_Set_array(adjacencyNR);
	nauty_output nautyRes = ComputeCanonAdj(5, outAdj);
	unsigned int bitset = convert_R_Adjacency_to_Bitset(nautyRes.canonical);
	free(nautyRes.canonical);
	vector<int> orb = Graphlet_Lib(bitset);

	for (int i = 0; i < 5; i++)
		orbits[nautyRes.label[i]] = orb[i];
	return orbits;
}

//global function that makes the Dictionary of orbits NOTE: the index of the output vector represents the adjacency right shifted by 22.
vector<vector<int>> mkDict() {
	vector<unsigned int> lib = mkLibAdj();
	vector<vector<int>> Dictionary(1024, vector<int>(5,-1));
	for (int i = 0; i < lib.size(); i++) {
		vector<int> orbits(5);
		orbits = Find_Orbits(lib[i]);
		unsigned int p = lib[i] >> 22;
		for (int j = 0; j < 5; j++){
			Dictionary[p][j] = orbits[j];
		}
	}
	return Dictionary;
}
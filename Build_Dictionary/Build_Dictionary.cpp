// Build_Dictionary.cpp : Defines the entry point for the console application.
//

#include <vector>
#include "Build_Dictionary.h"

using namespace std;

int main()
{
	//make library of all possible adjacency matrices with size 5
	vector<unsigned int> allAdj = mkLibAdj();
	unsigned int* out=convert_NR_to_Set_array(allAdj[10]);

	return 0;
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
		input=input << 1;

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
		for (int j = 0; j < 5; j++){
			output[i] = output[i] >> 1;
			output[i] = output[i] | (temp[i][j] == 1 ? a : 0);
		}
	}

	return output;
		
}


//make library of all possible adjacency matrices with size 5
vector<unsigned int> mkLibAdj(){
	vector<unsigned int> library;
	unsigned int bit = 0 ;
	int counter = 0;
	for (int j = 0; j < 2; j++)
		mkLibAdjHelper(library, counter, bit, j);

	return library;
}

//helps mkLibAdj, recursive function
void mkLibAdjHelper(vector<unsigned int> & lib , int n, unsigned int bit, int zo){
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
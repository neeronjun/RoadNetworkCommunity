//=======================================================================
// Copyright 2007 Aaron Windsor
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>

using namespace std;
using namespace boost;

int main(int argc, char** argv)
{

	// This program illustrates a simple use of boyer_myrvold_planar_embedding
	// as a simple yes/no test for planarity.
	

	typedef adjacency_list<vecS,
		vecS,
		undirectedS,
		property<vertex_index_t, int>
	> graph;


	graph K_4(4);
	fstream file;

	string line, temp, word; int traj[10];
	file.open("Cali_Edge_info.txt", ios::out | ios::in);
		
	while (file >> temp) {
		int count = 0;
		getline(file, line);
		stringstream s(line);
		while (getline(s, word, ' ')) {
			if (word == "")
				continue;
			traj[count] = atof(word.c_str());
			//cout << word << endl;
			count++;
		}
		add_edge(traj[0], traj[1], K_4);
		//printf("%d, %d \n", traj[0], traj[1]);
	}
	
	if (boyer_myrvold_planarity_test(K_4))
		std::cout << "K_4 is planar." << std::endl;
	else
		std::cout << "ERROR! K_4 should have been recognized as planar!"
		<< std::endl;

	return 0;
}
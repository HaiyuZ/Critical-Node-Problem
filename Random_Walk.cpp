#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <stack>
#include <algorithm>
#include <iomanip>
#include <math.h>

#define MAXV 30000
#define MAXE 200000
#define NIL -1

using namespace std;

int vertex_num = 2585568, edge_num = 33140018;
int steps = 0.0001 * vertex_num;
int D = 0;//设置为const 
int** adjlist;
int* vertex_degree = new int[vertex_num + 1];

double alpha = 1;
int* w = new int[vertex_num + 1];

struct pair{
	int first;
	int second;
	pair* next;
};

void read_adj(string file_name) {
	ifstream FIC;
	FIC.open(file_name);
	if (FIC.fail()) {
		cout << "### Erreur open, File_Name %s" << file_name;
		getchar();
		exit(0);
	}
	if (FIC.eof()) {
		cout << "### Erreur open, File_Name %s" << file_name;
		exit(0);
	}

	adjlist = new int* [vertex_num + 1];
	//adjlist[0] = NULL;
	for(int i = 0; i < vertex_num + 1; ++i)
		adjlist[i] = NULL;

	string line;
	int last_node = 1;
	int node_degree = 0;
	vector<int> neibors;
	
	pair* head = pair();
	pair* cur = head;
	for (int line_num = 0; line_num < edge_num; ++line_num) {
		getline(FIC, line);
		istringstream iss(line);
		int cur_node;
		iss >> cur_node;
		int neibor_node;
		iss >> neibor_node;
		if (line_num % 1000000 == 0) 
			cout << line_num << " " << cur_node << " " << neibor_node << endl;
		/*	
		if(cur_node != last_node) {//记得加上最后一行 
			
		}
		else {
			node_degree++;
			neibors.push_back(neibor_node);
		}
		*/ 
		pair* this_pair = 
		if (cur_node != last_node) {
			D += node_degree;
			if(adjlist[last_node] == NULL) {
				adjlist[last_node] = new int[node_degree];
				for (int i = 0; i < node_degree; ++i)
					adjlist[last_node][i] = neibors[i];
				vertex_degree[last_node] = node_degree;
			}
			else {
				int* old_array = adjlist[last_node];
				int old_degree = vertex_degree[last_node];
				int new_degree = old_degree + node_degree;
				adjlist[last_node] = new int[new_degree];
				for(int i = 0; i < old_degree; ++i) adjlist[last_node][i] = old_array[i];
				for(int i = 0; i < node_degree; ++i) adjlist[last_node][i + old_degree] = neibors[i];
				delete [] old_array;
			}
			neibors.clear();
			neibors.push_back(neibor_node);
			last_node = cur_node;
			node_degree = 1;
			/* 
			if (line_num == edge_num - 1) {
				adjlist[cur_node] = new int[1];
				adjlist[cur_node][0] = neibor_node;
				vertex_degree[cur_node] = 1;
				D++;
			}*/ 
		}
		else {
			node_degree++;
			neibors.push_back(neibor_node);
			/* 
			if (line_num == edge_num - 1) {
				adjlist[last_node] = new int[node_degree];
				for (int i = 0; i < node_degree; ++i)
					adjlist[last_node][i] = neibors[i];
				vertex_degree[last_node] = node_degree;
				D += node_degree;
			}*/ 
		}
	}
}

void show_adjlist() {
	for (int i = 1; i < vertex_num + 1; ++i) {
		for (int j = 0; j < vertex_degree[i]; ++j) {
			cout << adjlist[i][j] << " ";
		}
		cout << endl;
	}
}


int random_walk() {
	//for(int i = 1; i < vertex_num + 1; ++i)
		//if(adjlist[i] == NULL) cout << i << endl;
	
	//for (int i = 0; i < 10; ++i) cout << rand() % D << endl;
	//cout << "here" << endl;
	for(int i = 0; i < vertex_degree[514503]; i++) cout << adjlist[514503][i] << " ";
	cout << endl;
	for(int i = 0; i < vertex_degree[48082]; i++) cout << adjlist[48082][i] << " ";
	cout << endl;
	for(int i = 0; i < vertex_degree[132883]; i++) cout << adjlist[132883][i] << " ";
	cout << endl;
	int rand_int = rand() % D;
	int sum_degree = 0;
	int first_node = -1;
	for (int i = 1; i < vertex_num + 1; ++i) {
		if (rand_int >= sum_degree && rand_int < sum_degree + vertex_degree[i]) {
			first_node = i;
			break;
		}
		else sum_degree += vertex_degree[i];
	}
	int* R = new int[steps];
	R[0] = first_node;
	cout << "first_node : " << first_node << endl;
	int last_node = first_node;
	for (int i = 1; i < steps; ++i) {
		if(adjlist[last_node] == NULL) {
			cout << "error" << endl;
			return 0;
		}
		last_node = adjlist[last_node][rand() % vertex_degree[last_node]];
		R[i] = last_node;
		cout << i << "th" << last_node << endl;
	}
	int m = 0.025 * steps;
	int P = 0, Q = 0;
	int size = 0;
	cout << "here" << endl; 
	for (int i = 0; i < vertex_num; ++i) {
		w[i] = vertex_degree[i];
	}

	for (int i = 0; i < steps - m + 1; ++i) {
		for (int j = i + m; j < steps; ++j) {
			P += w[R[i]] * (R[i] == R[j] ? 1 : 0);
			Q += w[R[i]] * vertex_degree[R[i]] / vertex_degree[R[j]];
			size++;
		}
	}
	P /= size;
	Q /= size;
	return Q / P;
}

int main() {
	string file_name = "flickr-growth.txt/flickr-growth-sorted.txt";
	read_adj(file_name);
	//show_adjlist();
	cout << random_walk() - vertex_num << endl;
}

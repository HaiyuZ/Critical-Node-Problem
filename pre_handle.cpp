#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <vector>
#include <set>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <stack>
#include <algorithm>
#include <iomanip>
#include <math.h>
using namespace std;
/*
const int vertex_num = 3100000, edge_num = 250000000;

vector<int> tmpv;
vector<vector<int> > adjlist(vertex_num + 1, tmpv);
//set<int> tmps;
//vector<set<int> > s(vertex_num + 1, tmps);

void read_adj(string file_name) {
	ifstream FIC;
	FIC.open(file_name);
	string line;
	int max_ = -1;
	int last_cur = -1;
	int last_neibor = -1;
	//for(int i = 0; i < 100000000; ++i) getline(FIC, line);
	for (int line_num = 0; line_num < edge_num; ++line_num) {
		getline(FIC, line);
		istringstream iss(line);
		int cur_node, neibor_node;
		iss >> cur_node >> neibor_node;
		
		if(cur_node == neibor_node) continue;
		
		if(cur_node == last_cur && neibor_node == last_neibor) break;
		last_cur = cur_node;
		last_neibor = neibor_node;
		
		max_ = max(max_, neibor_node);
		max_ = max(max_, cur_node);
		//if(s[cur_node].find(neibor_node) == s[cur_node].end()) {
			//s[cur_node].insert(neibor_node);
			adjlist[cur_node].push_back(neibor_node);
		//}
		//if(s[neibor_node].find(cur_node) == s[neibor_node].end()) {
			//s[neibor_node].insert(cur_node);
			adjlist[neibor_node].push_back(cur_node);
		//}
		if (line_num % 1000000 == 0) 
			cout << line_num << " " << cur_node << " " << neibor_node << endl;
	}
	cout << "max vertex ID : " << max_ << endl;
}

void pre_handle() {
	ofstream out;
	out.open("orkut_Handled_total0.txt");
	int count = 0;
	for(int i = 0; i < vertex_num + 1; ++i) {
		if(adjlist[i].size() == 0) continue;
		count++;
		out << i << " ";
		for(int j = 0; j < adjlist[i].size(); ++j) {
			out << adjlist[i][j] << " ";
		}
		if(i % 100000 == 0) cout << i << endl;
		out << endl;
	}
	cout << "count vertex number : " << count << endl; 
}

void rehandle(string file_name) {
	ifstream FIC;
	FIC.open(file_name);
	ofstream out;
	out.open("orkut_Handled_total_rehandled0.txt");
	string line;
	set<int> s;
	int last = -1;
	for(int i = 0; i < vertex_num; ++i) {
		getline(FIC, line);
		istringstream iss(line);
		int cur, neibor;
		iss >> cur;
		if(cur == last) break;
		last = cur;
		while(iss >> neibor) {
			s.insert(neibor);
		}
		out << cur << " ";
		for(auto it = s.begin(); it != s.end(); it++) {
			out << *it << " ";
		}
		out << endl;
		s.clear();
		if(i % 100000 == 0) cout << i << " " << cur << endl;
	}
	FIC.close();
	/*
	for(int i = 0; i < 3100000; ++i) {
		if(i == curr.size()) break;
		out << curr[i] << " ";
		for(auto it = vs[i].begin(); it != vs[i].end(); it++) {
			out << *it << " ";
		}
		out << endl;
		if(i % 100000 == 0) cout << i << " " << curr[i] << endl;
	}
	out.close();
	//cout << cur << endl;
}

void handle(string file_name) {
	ifstream in;
	in.open(file_name);
	string line;
	int last = -1;
	for(int i = 0; i < vertex_num; ++i) {
		getline(in, line);
		istringstream iss(line);
		int cur;
		iss >> cur;
		if(cur == last) {
			cout << i << " " << cur << endl;
			//cout << line << endl;
			break;
		}
		last = cur;
	}
}
*/
ofstream out;

void find(string file_name) {
	ifstream in;
	in.open(file_name);
	int max_ = -1;
	int min_ = 100;
	int count = 0;
	int last = -1;
	vector<int> re;
	string line;
	for(int i = 0; i < 5000000; ++i) {
		int cur;
		getline(in, line);
		istringstream iss(line);
		iss >> cur;
		if(cur == last) break;
		last = cur;
		max_ = max(max_, cur);
		min_ = min(min_, cur);
		count++;
	}
	
	cout << file_name << endl;
	cout << "max_ : " << max_ << endl;
	cout << "min_ : " << min_ << endl;
	cout << " count : " << count << endl;
	//cout << "re : ";
	for(int i = 0; i < re.size(); ++i) cout << re[i] << " ";
	cout << endl << endl;
	
	out << file_name << endl;
	out << "max vertex ID : " << max_ << endl;
	out << "min vertex ID : " << min_ << endl;
	out << "total vertex number : " << count << endl;
	//out << "re : ";
	for(int i = 0; i < re.size(); ++i) out << re[i] << " ";
	out << endl;
	out << endl;
}

int main() {
	//string file_name = "orkut-links.txt/release-orkut-links.txt";
	//read_adj(file_name);
	//pre_handle();

	
	out.open("dataset_information.txt");
	string file_name = "orkut_Handled.txt";
	find(file_name);
	
	file_name = "Flickr_Handled.txt";
	find(file_name);
	
	file_name = "Live_Journal_Handled.txt";
	find(file_name);
}

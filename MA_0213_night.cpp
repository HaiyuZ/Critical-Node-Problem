#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <queue>
using namespace std;

#define maxn 30100
#define BETA 0.6
#define SP0 85
#define MaxIdleIters 1000//是否太大了
#define parent_num 20

int kbv;
int vertex_num, K, visited_vertex_num, run_times;
vector<vector<int>> adjlist;
int visited_vertex[maxn];
int weight[maxn], component_nid[maxn];
int XB[maxn];
bool deleted[maxn], visited[maxn];
vector<int> component[maxn];
vector<int> solutions[parent_num + 10];

vector<int> solution_operting;//当前操作的解
vector<int> best_solution_once;//单次最优解
int min_connectivity_once;//单次最优解的连通度
vector<int> best_solution;//所有次数里的最佳解
int f_best;//最优解的连通度
int f_avg;//所有单次最优解的平均连通度
int succeed_times;//达到最优解的次数
int exch;//交换次数
vector<int> tmpS;

stack<int> unused_nid;
set<int> used_nid;

double t_avg, init_cost_time, func_time;
long long iter;
clock_t start_time;
double time_out;

vector<int> connectivities(parent_num + 1, 0);//每个解的连通度
int* sizes = new int[maxn];
int* key_vertexs = new int[maxn];

vector<int> tmp_dis(parent_num + 1, 0);
vector<vector<int>> distances(parent_num + 1, tmp_dis);
vector<int> total_distances(parent_num + 1, 0);
vector<int> rank_connectivity(parent_num + 1, 0);

inline double time_cost() {
	return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

inline bool timeout() {
	return time_cost() > time_out;
}

void read_graph(string filename) {
	ifstream in;
	in.open(filename);
	string line;
	getline(in, line);
	istringstream iss(line);
	iss >> vertex_num;
	adjlist.resize(vertex_num);
	int tmp;
	char c;
	for (int i = 0; i < vertex_num; ++i) {
		getline(in, line);
		istringstream iss(line);
		iss >> tmp >> c;
		while(iss >> tmp) adjlist[i].push_back(tmp);
	}
}

void dfs(int cur) {
	visited_vertex[visited_vertex_num++] = cur;
	visited[cur] = true;
	for (auto &num:adjlist[cur])
		if (!deleted[num] && !visited[num]) dfs(num);//递归深度优先遍历改为循环遍历？
}

void dfs_split(int start) {//component信息记录
	visited_vertex_num=0;
	dfs(start);
	int cpn_nid = unused_nid.top();//component能用的id
	unused_nid.pop();
	used_nid.insert(cpn_nid);
	vector<int>& cur = component[cpn_nid];
	cur.resize(visited_vertex_num);
	key_vertexs[cpn_nid] = visited_vertex[0];
	sizes[cpn_nid] = visited_vertex_num;
	for (int j = 0; j < visited_vertex_num; ++j) {
		component_nid[visited_vertex[j]] = cpn_nid;//每个顶点所属的component的id
		cur[j] = visited_vertex[j];//每个component包含的顶点
		if (weight[visited_vertex[j]] > weight[key_vertexs[cpn_nid]])
			key_vertexs[cpn_nid] = visited_vertex[j];
		visited[visited_vertex[j]] = false;
	}
}

int Connect() {//可以采用蔡泽杰的方法求连通度
	while (!unused_nid.empty()) unused_nid.pop();
	for (int i=0; i<=vertex_num; ++i) unused_nid.push(i);  //add one more
	used_nid.clear();
	fill(visited, visited+vertex_num, false);
	int connectivity=0;
	for (int i=0; i<vertex_num; ++i)
		if (!deleted[i] && !visited[i]) {
			dfs_split(i);
			for (int j=0; j<visited_vertex_num; ++j)
				visited[visited_vertex[j]] = true;
			connectivity += (visited_vertex_num - 1) * visited_vertex_num / 2;
		}
	fill(visited, visited + vertex_num, false);
	return connectivity;
}

int find_large_component() {
	int L = vertex_num, R = 0;
	for (auto &num:used_nid){
		L = min(L, sizes[num]);
		R = max(R, sizes[num]);
	}
	int bound = (L + R) >> 1;//找阈值
	vector<int> large_component_id;
	for (auto &num:used_nid)
		if (sizes[num] >= bound)
			large_component_id.push_back(num);
	return large_component_id[rand() % large_component_id.size()];
}

void put_back() {                        //remove the optimal one 
	int increase=maxn*maxn, back_vertex=-1;
	int idx=-1;
	for (auto &num:solution_operting) {//遍历解内所有成员，找放回后使连通度上升最少的
		++idx;
		int delta=0, cnt=0, cur;
		int sub_vertex_num=0;
		for (auto &adj:adjlist[num])
			if (!deleted[adj])
				if (!visited[cur = component_nid[adj]]) {//这里的visited有没有被其他函数修改了，修改为参数传递过来而不是全局变量？
					visited[cur] = true;//直接每次要用visited都初始化为false，如果是同一个visited，就作为参数传递过去
					delta += component[cur].size() * sub_vertex_num;
					sub_vertex_num += sizes[cur];
				}
		delta += sub_vertex_num;
		for (auto &adj:adjlist[num])
			if (!deleted[adj]) visited[component_nid[adj]] = false;
		if (delta <= increase) increase = delta, back_vertex = idx;
	}
	swap(solution_operting[back_vertex], solution_operting.back());
	back_vertex = solution_operting.back();
	solution_operting.pop_back();

	deleted[back_vertex]=false;
	weight[back_vertex]=0;
	int cur;
	for (auto &adj:adjlist[back_vertex])
		if (!deleted[adj])
			if (!visited[cur=component_nid[adj]]) {//更新component的id使用情况
				visited[cur] = true;
				unused_nid.push(cur);
				used_nid.erase(cur);
			}
	for (auto &adj:adjlist[back_vertex])
		if (!deleted[adj]) visited[component_nid[adj]]=false;//更新visited数组
	dfs_split(back_vertex);//更新component
}

void remove(int remove_vertex) {         //remove the specific one
	int cur = component_nid[remove_vertex];
	deleted[remove_vertex] = true;
	for (auto &num:adjlist[remove_vertex])
		if (!deleted[num] && component_nid[num] == cur) {
			dfs_split(num);//重新分割原来的component
			for (int i = 0; i < visited_vertex_num; ++i)
				weight[visited_vertex[i]]++;//不用在每个循环都这样做，直接在那个component的list里面就可以找到所有的顶点
		}
	unused_nid.push(cur);
	used_nid.erase(cur);
}

void component_based_neghborhood_search() {
    clock_t func_start_time = clock();
    fill(weight, weight + vertex_num, 0);
	tmpS = solution_operting;
	int connectivity = Connect();
	int idle_iter = 0, tmp;
	while (idle_iter < MaxIdleIters) {
		int remove_vertex = key_vertexs[find_large_component()];
		solution_operting.push_back(remove_vertex);//增加顶点
		remove(remove_vertex);
		put_back();
        ++iter;
		if ((tmp=Connect()) < connectivity) {
			tmpS = solution_operting;
			connectivity = tmp;
			idle_iter = 0;
		}
		else idle_iter++;
	}
	for (auto &num:solution_operting) deleted[num] = false;
	//solution_operting = tmpS;
	for (auto &num:solution_operting) deleted[num] = true;
    func_time += (double)(clock() - func_start_time) / CLOCKS_PER_SEC;
}

double get_distance(vector<int> &S1, vector<int> &S2) {
	int common = 0;
	fill(visited, visited + vertex_num, false);
	for (auto &num:S1) visited[num] = true;
	for (auto &num:S2) common += visited[num];
	return K - common;
}

int deth_first_search(int start) {
	stack<int> s;
	s.push(start);
	int size = 0;
	while(!s.empty()) {
		int vertex = s.top();
		visited[vertex] = true;
		s.pop();
		size++;
		for(int i = 0; i < adjlist[vertex].size(); ++i) {
			int neibor = adjlist[vertex][i];
			if(!deleted[neibor] && !visited[neibor]) s.push(neibor);
		}
	}
	return size;
}

//只需要得到连通度，不需要将原图分割、记录子图信息时，使用该函数
int get_connectivity() {
	int connectivity = 0;
	fill(visited, visited + vertex_num, false);
	for(int i = 0; i < vertex_num; ++i) {
		if(!deleted[i] && !visited[i]) {
			int size = deth_first_search(i);
			connectivity += size * (size - 1) / 2;
		}
	}
	return connectivity;
}

void population_initialization() {
	for (int i = 0; i < parent_num; ++i) {
		fill(deleted, deleted + vertex_num, false);
		solution_operting.resize(K);
		for (int j = 0; j < K; ++j) {
			solution_operting[j] = rand()%vertex_num;
			while (deleted[solution_operting[j]]) solution_operting[j] = rand() % vertex_num;
			deleted[solution_operting[j]] = true;
		}
		component_based_neghborhood_search();
		while (true) {
			bool differ = true;
			for (int j = 0; j < i; ++j)
				if (get_distance(solution_operting, solutions[j]) == 0) {//可以设为小于某个阈值，则进行交换操作？
					differ = false;
					break;
				}
			if (differ) break;
			int idx = rand() % K;
			int add_vertex = rand() % vertex_num;
			while(deleted[add_vertex]) add_vertex = rand() % vertex_num;
			deleted[solution_operting[idx]] = false;
			deleted[add_vertex] = true;
			solution_operting[idx] = add_vertex;
		}
		solutions[i] = solution_operting;
		int cnty = (connectivities[i] = get_connectivity());
		cout << "initialization connectivity : " << cnty << endl;

		//记录单次最优解
		if(cnty < min_connectivity_once) {
			min_connectivity_once = cnty;
			best_solution_once = solution_operting;
		}
		
		//初始化距离矩阵与总距离
		for(int j = 0; j < i; ++j) {
            int common = 0;
            for(int m = 0; m < K; ++m) common += deleted[solutions[j][m]];
            int dis = K - common;
            distances[i][j] = (distances[j][i] = dis);
            total_distances[j] += dis;
            total_distances[i] += dis;
        }
	}

	//初始化各个解的连通度排名
	vector<int> v_connectivity(connectivities.begin(), connectivities.end() - 1);
	sort(v_connectivity.begin(), v_connectivity.end());//手写用二分插入方法而不用排序。
	for(int i = 0; i < parent_num; ++i) {
		auto it = find(v_connectivity.begin(), v_connectivity.end(), connectivities[i]);
		rank_connectivity[i] = distance(v_connectivity.begin(), it);
	}
}

void double_backbone_based_crossover(vector<int> &S1, vector<int> &S2) {
	int cnt = 0;
	fill(deleted, deleted + vertex_num, false);//重置del数组，表示所有点都未放进S
	solution_operting.clear();
	for (auto &num:S1) visited[num] = true;
	for (auto &num:S2)
		if (visited[num]) {
			solution_operting.push_back(num);//在s1也在s2的点
			visited[num] = false;
		}
		else XB[cnt++] = num;//在s2但是不在s1里的点
	for (auto &num:S1)
		if (visited[num]) XB[cnt++] = num;//在s1但是不在s2里的点
	for (int i=0; i<cnt; ++i)
		if (rand() % 100 <= SP0) solution_operting.push_back(XB[i]);

	for (auto &num:solution_operting) deleted[num] = true;
	if (solution_operting.size() != K) Connect();
	if (solution_operting.size() < K) {
		for (auto &num:S1) visited[num] = true;
		for (auto &num:S2) visited[num] = true;
		while (solution_operting.size() < K) {
			int large_idx = find_large_component();
			vector<int>& cur = component[large_idx];
			int idx = rand() % sizes[large_idx];
			if (!visited[cur[idx]]) {
				solution_operting.push_back(cur[idx]);
				remove(solution_operting.back());
			}
		}
	}
	while (solution_operting.size() > K) put_back();
}

bool compare(int a,int b) {return a > b;}

void rank_based_pool_updating(int connectivity) {//尝试不要直接加进去再算要不要加，而是根据原本的解，构建一个阈值，不用每次计算这么多
	solutions[parent_num] = solution_operting;

	//更新连通度以及连通度排名
	connectivities[parent_num] = connectivity;
	int count = 0;
	for(int i = 0; i < parent_num; ++i)
		if(connectivities[i] > connectivity) {
			rank_connectivity[i]++;
			count++;
		}
	rank_connectivity[parent_num] = parent_num - count;

	//求总距离
	for(int i = 0; i < parent_num; ++i) {
		int dis = 0;
		for(int j = 0; j < K; ++j) dis += deleted[solutions[i][j]];
		dis = K - dis;
		distances[parent_num][i] = dis;
		distances[i][parent_num] = dis;
		total_distances[parent_num] += dis;
		total_distances[i] += dis;
	}//get_distance已经优化为，总是求同一个解与其他不同的解时的距离时
	vector<int> v_distance = total_distances;
	sort(v_distance.begin(), v_distance.end(), compare);

	//找出加权排名最低的解
	int idx = -1;
	double max_score = -1;
	for(int i = 0; i < parent_num + 1; ++i) {
		auto it = find(v_distance.begin(), v_distance.end(), total_distances[i]);//返回该total_distance在排序后的数组中的迭代器
		double score = BETA * rank_connectivity[i] + (1 - BETA) * distance(v_distance.begin(), it);//求加权排名
		if(score > max_score) {
			idx = i;
			max_score = score;
		}
	}

	if(idx == parent_num) return;

	//更新距离矩阵distances，总距离total_distances
	//total_distance需要更新，distances上三角与下三角都要更新
	vector<int> delete_distances = distances[idx];
	for(int i = 0; i < idx; ++i) 
		distances[i][idx] = (distances[idx][i] = distances[parent_num][i]);
	for(int i = idx + 1; i < parent_num; ++i) 
		distances[idx][i] = (distances[i][idx] = distances[parent_num][i]);
	total_distances[idx] = total_distances[parent_num];
	total_distances[parent_num] = 0;
	for(int i = 0; i < parent_num; ++i) 
		if(i != idx)
			total_distances[i] -= delete_distances[i];
		else total_distances[i] -= delete_distances[parent_num];
	swap(solutions[idx], solutions[parent_num]);

	//交换connectivity，交换rank
	int delete_connectivtiy = connectivities[idx];
	connectivities[idx] = connectivities[parent_num];
	rank_connectivity[idx] = rank_connectivity[parent_num];
	for(int i = 0; i < parent_num + 1; ++i) 
		if(connectivities[i] >= delete_connectivtiy) rank_connectivity[i]--;
}

void Critical_Node_Problem() {
	srand(time(0));
	min_connectivity_once = vertex_num * vertex_num;//该次最小的连通度
	best_solution_once.clear();

	//初始化所有解
	start_time = clock();
	population_initialization();
	init_cost_time += time_cost();

	start_time = clock();
	double cost_time;

	fill(weight, weight + vertex_num, 0);
	while (!timeout()) {
		int Si = rand() % parent_num, Sj = rand() % parent_num;
		while (Si == Sj) Si = rand() % parent_num, Sj = rand() % parent_num;
		double_backbone_based_crossover(solutions[Si], solutions[Sj]);
		component_based_neghborhood_search();
		int offspring_connectvity = Connect();
		cout << "find solution objective value so far this time : " << min_connectivity_once << endl;
		if (offspring_connectvity <= kbv) {
			cout << "find best solution " << kbv << endl;
			cout << "time_cost : " << time_cost() << endl;
			break;
		}
		if (offspring_connectvity < min_connectivity_once) {
			best_solution_once = solution_operting;
			min_connectivity_once = offspring_connectvity;
			cost_time = time_cost();
		}
		rank_based_pool_updating(offspring_connectvity);
	}

	//更新所有次数最优解
	if (min_connectivity_once < f_best) {
        best_solution = best_solution_once;
        f_best = min_connectivity_once;
        t_avg = cost_time;
        succeed_times = 1;
    }
    else if (min_connectivity_once == f_best) {
		t_avg += cost_time;
		succeed_times++;
	}
    f_avg += min_connectivity_once;
}

void print(vector<int> &S){
	sort(S.begin(), S.end());
	for (auto &vertex:S) cout << vertex << " ";
	fill(deleted, deleted + vertex_num, false);
	for (auto &num:S) deleted[num]=true;
	cout << "\n" << Connect() << endl;
}

void show_adjlist() {
	for(int i = 0; i < vertex_num; ++i) {
		cout << i << " : ";
		for(int j = 0; j < adjlist[i].size(); ++j) {
			cout << adjlist[i][j] << " ";
		}
		cout << endl;
	}
}

int main() {
	K = 125;
	kbv = 2078;
	string filename = "BenchMarks/cnd/WattsStrogatz_n500.txt";
	read_graph(filename);
	time_out = 3600;
	run_times = 10;
    f_best = vertex_num * vertex_num, f_avg = 0, succeed_times = 0;
    t_avg = 0; init_cost_time = 0; func_time = 0; iter = 0;
    for (int i = 0; i < run_times; ++i) Critical_Node_Problem();
    cout << fixed;
    cout << "f_best " << f_best << endl;
    cout << "f_avg " << setprecision(1) << double(f_avg)/run_times << endl;
    cout << "succeed_times " << succeed_times << endl; 
    cout << "t_avg " << setprecision(1) << t_avg / succeed_times << endl;
    cout << "init_cost_time " << setprecision(1) << init_cost_time / run_times << endl;
    cout << "func_time " << setprecision(1) << func_time / run_times << endl;
    cout << "iter " << iter << endl;
    return 0;
}

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

int vertex_num, K, visited_vertex_num, run_times;
vector<int> *adjlist, *component;
int *vertex_degree, *visited_vertex, *weight, *component_nid, *nid;
bool *del, *visited;
vector<int> S, S0, S_best, tmpS;
stack<int> unused_nid;
set<int> used_nid;
int f_best, f_avg, succ;
double t_avg, init_cost_time, func_time;
long long iter;
clock_t start_time;
double time_out;
ifstream in;

int* connectivities = new int[parent_num + 1];//每个解的连通度
int* sizes = new int[maxn];
int* key_vertexs = new int[maxn];

vector<int> parent[parent_num + 1];
int** distances = new int*[parent_num + 1];
int* total_distance = new int[parent_num + 1];
vector<int> v_connectivity(parent_num + 1, 0);
vector<int> v_distance(parent_num + 1, 0);
int* rank_connectivity = new int[parent_num + 1];
int* rank_distance = new int[parent_num + 1];


inline double time_cost() {
	return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

inline bool timeout() {
	return time_cost() > time_out;
}

void read() {
	string line;
	int num;
	char c;
	for (int i = 0; i < vertex_num; ++i) {
		getline(in, line);
		istringstream iss(line);
		iss >> num >> c;
		while(iss >> num) adjlist[i].push_back(num);
	}
}

void dfs(int cur) {
	visited_vertex[visited_vertex_num++] = cur;
	visited[cur] = true;
	for (auto &num:adjlist[cur])
		if (!del[num] && !visited[num]) dfs(num);//递归深度优先遍历改为循环遍历？
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
	for (int i = 0; i < vertex_num; ++i) unused_nid.push(i);  //add one more
	used_nid.clear();
	fill(visited, visited+vertex_num, false);
	int connectivity=0;
	for (int i=0; i<vertex_num; ++i)
		if (!del[i] && !visited[i]) {
			dfs_split(i);
			for (int j=0; j<visited_vertex_num; ++j)
				visited[visited_vertex[j]] = true;
			connectivity += (visited_vertex_num - 1) * visited_vertex_num / 2;
		}
	fill(visited, visited + vertex_num, false);
	return connectivity;
}

int find_large_component() {
	int L=vertex_num, R=0;
	for (auto &num:used_nid)
	{
		L=min(L, sizes[num]);
		R=max(R, sizes[num]);
	}
	int bound = (L + R) >> 1;//找阈值
	int cnt = 0;
	for (auto &num:used_nid)
		if (sizes[num] >= bound)
			nid[cnt++]=num;
	return nid[rand()%cnt];
}

void put_back() {                        //remove the optimal one 
	int increase=maxn*maxn, back_vertex=-1;
	int idx=-1;
	for (auto &num:S) {//遍历解内所有成员，找放回后使连通度上升最少的
		++idx;
		int delta=0, cnt=0, cur;
		int sub_vertex_num=0;
		for (auto &adj:adjlist[num])
			if (!del[adj])
				if (!visited[cur = component_nid[adj]]) {//这里的visited有没有被其他函数修改了，修改为参数传递过来而不是全局变量？
					visited[cur] = true;//直接每次要用visited都初始化为false，如果是同一个visited，就作为参数传递过去
					delta -= sizes[cur] * sizes[cur];
					sub_vertex_num += sizes[cur];
				}
		delta += sub_vertex_num * sub_vertex_num - 1;
		delta >> 1;//HaiyuZ推导的公式
		for (auto &adj:adjlist[num])
			if (!del[adj]) visited[component_nid[adj]] = false;
		if (delta <= increase) increase = delta, back_vertex = idx;
	}
	swap(S[back_vertex], S.back());
	back_vertex=S.back();
	S.pop_back();

	del[back_vertex]=false;
	weight[back_vertex]=0;
	int cur;
	for (auto &adj:adjlist[back_vertex])
		if (!del[adj])
			if (!visited[cur=component_nid[adj]]) {//更新component的id使用情况
				visited[cur] = true;
				unused_nid.push(cur);
				used_nid.erase(cur);
			}
	for (auto &adj:adjlist[back_vertex])
		if (!del[adj]) visited[component_nid[adj]]=false;//更新visited数组
	dfs_split(back_vertex);//更新component
}

void remove(int remove_vertex) {         //remove the specific one
	int cur = component_nid[remove_vertex];
	del[remove_vertex] = true;
	for (auto &num:adjlist[remove_vertex])
		if (!del[num] && component_nid[num] == cur) {
			dfs_split(num);//重新分割原来的component
			for (int i = 0; i < visited_vertex_num; ++i)
				weight[visited_vertex[i]]++;//不用在每个循环都这样做，直接在那个component的list里面就可以找到所有的顶点
		}
	unused_nid.push(cur);
	used_nid.erase(cur);
}

void ComponentBasedNeighborhoodSearch() {
    clock_t func_start_time=clock();
    fill(weight, weight+vertex_num, 0);
	tmpS=S;
	int connectivity = Connect();
	int idle_iter = 0, tmp;
	while (idle_iter<MaxIdleIters) {
		int remove_vertex = key_vertexs[find_large_component()];
		S.push_back(remove_vertex);//增加顶点
		remove(remove_vertex);
		put_back();//放回一个顶点到图中
        ++iter;
		if ((tmp=Connect()) < connectivity) {
			tmpS = S;
			connectivity = tmp;
			idle_iter = 0;
		}
		else idle_iter++;
	}
	for (auto &num:S) del[num] = false;
	S=tmpS;
	for (auto &num:S) del[num] = true;
    func_time+=(double)(clock() - func_start_time)/CLOCKS_PER_SEC;
}

double get_distance(vector<int> &S1, vector<int> &S2) {
	int Sim=0;
	for (auto &num:S1) visited[num] = true;
	for (auto &num:S2) Sim += visited[num];
	for (auto &num:S1) visited[num] = false;
	return K - Sim;
}

void PoolInitialize() {
	fill(del, del + vertex_num, false);
	for (int i = 0; i < parent_num; ++i) {
		S.resize(K);
		for (int j=0; j<K; ++j) {
			S[j]=rand()%vertex_num;
			while (del[S[j]]) S[j]=rand()%vertex_num;
			del[S[j]]=true;
		}
		ComponentBasedNeighborhoodSearch();
		while (true)
		{
			bool flag=true;
			for (int j = 0; j < i; ++j)
				if (get_distance(S, parent[j])==0) {//可以设为小于某个阈值，则进行交换操作？
					flag=false;
					break;
				}
			if (flag) break;
			int idx=rand()%K;
			swap(S[idx], S.back());//不用把顶点放到最后再进行交换操作吧？
			del[S.back()]=false;
			idx=rand()%vertex_num;
			while (del[idx]) idx=rand()%vertex_num;//会不会浪费很多时间？
			del[S.back()=idx]=true;//交换顶点
		}
		parent[i] = S;
		v_connectivity[i] = (connectivities[i] = Connect());//初始化v_connectivity
		for(int j = 0; j < i; ++j)
			distances[i][j] = get_distance(parent[i], S);//初始化距离矩阵，get_distance也还可以优化为，总是求同一个解与其他不同的解时的距离时
		for (auto &num:S) del[num]=false;
	}
	sort(v_connectivity.begin(), v_connectivity.end());//手写用二分插入方法而不用排序。
	for(int i = 0; i < parent_num; ++i) {
		auto it = find(v_connectivity.begin(), v_connectivity.end(), connectivities[i]);
		rank_connectivity[i] = distance(it, v_connectivity.begin());
	}
	for(int i = 0; i < parent_num; ++i) {//初始化p_distance
		for(int j = i + 1; j < parent_num; ++j) 
			total_distance[i] += distances[j][i];
		for(int j = 0; j < i; ++j)
			total_distance[i] += distances[i][j];
		v_distance[i] = total_distance[i];
	}
}

void DoubleBackboneBasedCrossover(vector<int> &S1, vector<int> &S2) {
	int cnt=0;
	for (auto &num:S) del[num] = false;//重置del数组，表示所有点都未放进S
	S.clear();
	for (auto &num:S1) visited[num] = true;
	for (auto &num:S2)
		if (visited[num]) {
			S.push_back(num);
			visited[num] = false;
		}
		else nid[cnt++] = num;//与前面large component id用了同一个名字
	for (auto &num:S1)
		if (visited[num]) nid[cnt++] = num;
	for (int i=0; i<cnt; ++i)
		if (rand()%100<=SP0) S.push_back(nid[i]);

	for (auto &num:S) del[num] = true;
	if (S.size() != K) Connect();
	if (S.size() < K) {
		for (auto &num:S1) visited[num] = true;
		for (auto &num:S2) visited[num] = true;
		while (S.size() < K) {
			int large_idx = find_large_component();
			vector<int>& cur = component[large_idx];
			int idx = rand() % sizes[large_idx];
			if (!visited[cur[idx]]) {
				S.push_back(cur[idx]);
				remove(S.back());
			}
		}
	}
	while (S.size() > K) put_back();
}

void RankBasedPoolUpdating(int connectivity){
	parent[parent_num] = S;
	connectivities[parent_num] = connectivity;
	int cnt = 0;
	for(int i = 0; i < parent_num; ++i)
		if(connectivities[i] > connectivity) {
			rank_connectivity[i]++;
			cnt++;
		}
	rank_connectivity[parent_num] = parent_num - cnt;
	int idx = -1, max_score = 0;
	int distance_ = 0;
	for(int i = 0; i < parent_num; ++i) 
		distance_ += (distances[parent_num][i] = get_distance(S, parent[i]));//get_distance也还可以优化为，总是求同一个解与其他不同的解时的距离时
	total_distance[parent_num] = distance_;
	for(int i = 0; i < parent_num; ++i) {
		total_distance[i] += distances[parent_num][i];
		v_distance[i] = total_distance[i];
	}
	sort(v_distance.begin(), v_distance.end());
	for(int i = 0; i < parent_num + 1; ++i) {
		auto it = find(v_distance.begin(), v_distance.end(), total_distance[i]);
		int score = BETA * rank_connectivity[i] + (1 - BETA) * distance(it, v_distance.begin());//返回该total_distance在排序后的数组中的迭代器
		if(score > max_score) {
			idx = i;
			max_score = score;
		}
	}
	int delete_connectivtiy = connectivities[idx];
	if(idx != parent_num) {
		for(int i = 0; i < parent_num + 1; ++i) 
			if(connectivities[i] >= delete_connectivtiy) rank_connectivity[i]--;
		for(int i = 0; i < idx; ++i) 
			distances[idx][i] = distances[parent_num][i];
		for(int i = idx + 1; i < parent_num; ++i) 
			distances[i][idx] = distances[parent_num][i];
		swap(parent[idx], parent[parent_num]);
		rank_connectivity[idx] = rank_connectivity[parent_num];
		connectivities[idx] = connectivities[parent_num];
	}
}

void solve() {
	srand(time(0));
	int best_connectivity=maxn*maxn, tmp;
	S0.clear();
	start_time=clock();
	PoolInitialize();
	init_cost_time+=time_cost();
	for (int i=0; i<parent_num; ++i)
		if (connectivities[i]<best_connectivity)
		{
			best_connectivity = connectivities[i];
			S0 = parent[i];//创建解时动态更新，不需要这个时候再遍历
		}

	start_time=clock();
	double cost_time=0;

	fill(weight, weight+vertex_num, 0);
	while (!timeout()) {
		int Si=rand()%parent_num, Sj=rand()%parent_num;
		while (Si==Sj) Si=rand()%parent_num, Sj=rand()%parent_num;
		DoubleBackboneBasedCrossover(parent[Si], parent[Sj]);
		ComponentBasedNeighborhoodSearch();
		//cout << S[0] << " " << S[1] << endl;
		if ((tmp=Connect())<best_connectivity) {
			S0=S;
			best_connectivity=tmp;
			cost_time=time_cost();
		}
		RankBasedPoolUpdating(tmp);
	}
	if (best_connectivity<f_best)
    {
        S_best=S0;
        f_best=best_connectivity;
        t_avg=cost_time;
        succ=1;
    }
    else if (best_connectivity==f_best) t_avg+=cost_time, succ++;
    f_avg+=best_connectivity;
}

void print(vector<int> &S){
	sort(S.begin(), S.end());
	for (auto &vertex:S) cout << vertex << " ";
	fill(del, del+vertex_num, false);
	for (auto &num:S) del[num]=true;
	cout << "\n" << Connect() << endl;
}

void show_adjlist() {
	for(int i = 0; i < vertex_num; ++i) {
		cout << i << " : ";
		for(int j = 0; j < vertex_degree[i]; ++j) {
			cout << adjlist[i][j] << " ";
		}
		cout << endl;
	}
}

int main() {
	K = 50;
	string filename = "BenchMarks/cnd/BarabasiAlbert_n500m1.txt";
	in.open(filename);
	string line;
	getline(in, line);
	istringstream iss(line);
	iss >> vertex_num;
	adjlist = new vector<int>[vertex_num];
	read();
	//show_adjlist();
	for(int i = 0; i < parent_num + 1; ++i)
		distances[i] = new int[parent_num + 1];
	vertex_degree = new int[vertex_num];
	visited_vertex = new int[vertex_num];
	weight = new int[vertex_num];
	component_nid = new int[vertex_num];
	nid = new int[vertex_num];
	del = new bool[vertex_num];
	visited = new bool[vertex_num];
	component = new vector<int>[vertex_num];

	time_out=60;
	run_times=1;
    f_best=vertex_num*vertex_num, f_avg=0, succ=0;
    t_avg=0; init_cost_time=0; func_time=0; iter=0;
    for (int i=0; i<run_times; ++i) solve();
    cout << fixed;
    cout << "f_best " << f_best << "\n";
    cout << "f_avg " << setprecision(1) << double(f_avg)/run_times << "\n";
    cout << "succ " << succ << "\n"; 
    cout << "t_avg " << setprecision(1) << t_avg/succ << "\n";
    cout << "init_cost_time " << setprecision(1) << init_cost_time/run_times << "\n";
    cout << "func_time " << setprecision(1) << func_time/run_times << "\n";
    cout << "iter " << iter << "\n";
    //print(S_best);
    return 0;
}

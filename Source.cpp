#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <ctime>
#include <iomanip>

using namespace std;

// 1. Простые структуры
struct Point { double x, y; };
struct Edge { int to; double weight; };
struct Node {
    int id; double f;
    bool operator>(const Node& other) const { return f > other.f; }
};

// 2. Эвристика
double get_h(int a, int b, const vector<Point>& p) {
    return sqrt(pow(p[a].x - p[b].x, 2) + pow(p[a].y - p[b].y, 2));
}

// 3. A*
void a_star(int V, const vector<vector<Edge>>& adj, const vector<Point>& p) {
    vector<double> g(V, 1e15);
    priority_queue<Node, vector<Node>, greater<Node>> pq;
    g[0] = 0;
    pq.push({ 0, get_h(0, V - 1, p) });

    while (!pq.empty()) {
        int u = pq.top().id; pq.pop();
        if (u == V - 1) break;
        for (auto& e : adj[u]) {
            if (g[u] + e.weight < g[e.to]) {
                g[e.to] = g[u] + e.weight;
                pq.push({ e.to, g[e.to] + get_h(e.to, V - 1, p) });
            }
        }
    }
}

// 4. Флойд 
void floyd(int V, vector<vector<double>>& dist) {
    for (int k = 0; k < V; k++)
        for (int i = 0; i < V; i++)
            for (int j = 0; j < V; j++)
                if (dist[i][k] + dist[k][j] < dist[i][j])
                    dist[i][j] = dist[i][k] + dist[k][j];
}


void run(int V, double dens) {
    vector<vector<Edge>> adj(V);
    vector<vector<double>> mat(V, vector<double>(V, 1e15));
    vector<Point> p(V);

    for (int i = 0; i < V; i++) {
        p[i] = { (double)(rand() % 100), (double)(rand() % 100) };
        mat[i][i] = 0;
    }

    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (i != j && (double)rand() / RAND_MAX < dens) {
                double w = get_h(i, j, p) + 1;
                adj[i].push_back({ j, w });
                mat[i][j] = w;
            }
        }
    }

    // Замер A* 
    clock_t s1 = clock();
    for (int i = 0; i < 1000; i++) a_star(V, adj, p);
    clock_t e1 = clock();
    double t_a = ((double)(e1 - s1) * 1000000.0 / CLOCKS_PER_SEC) / 1000.0;

    // Замер Флойда 
    clock_t s2 = clock();
    floyd(V, mat);
    clock_t e2 = clock();
    double t_f = (double)(e2 - s2) * 1000.0 / CLOCKS_PER_SEC;

    cout << setw(5) << V << " | " << setw(5) << dens << " | "
        << setw(10) << fixed << setprecision(2) << t_a << " us | "
        << setw(10) << t_f << " ms" << endl;
}

int main() {
    srand(time(0));
    int counts[] = { 10, 20, 50, 100, 200, 500, 1000, 2000 };

    cout << "  V   |  Dens |   A* Time (us) | Floyd Time (ms)" << endl;
    cout << "------|-------|----------------|----------------" << endl;
    for (int v : counts) {
        run(v, 0.05);
        run(v, 0.5);
        cout << "------|-------|----------------|----------------" << endl;
    }
    return 0;
}
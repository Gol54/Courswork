#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <limits>
#include <cmath>
#include <random>
#include <chrono>
#include <numeric>
#include <functional> 
#include <iomanip>

using namespace std;

// --- 1. Œ¡Ÿ»≈ —“–” “”–€ ƒ¿ÕÕ€’ » “»œ€ ---

const double INF = numeric_limits<double>::max() / 2.0;

struct Edge {
    int target_node;
    double weight;
};

using Graph = vector<vector<Edge>>;
using DistanceMatrix = vector<vector<double>>;

// --- 2. ¿À√Œ–»“Ã A* (A-STAR) ---

struct NodeState {
    int node_id;
    double f_score;
    bool operator>(const NodeState& other) const {
        return f_score > other.f_score;
    }
};

vector<int> reconstruct_path(const vector<int>& came_from, int goal) {
    vector<int> path;
    int current = goal;
    while (current != -1) {
        path.push_back(current);
        current = came_from[current];
    }
    reverse(path.begin(), path.end());
    return path;
}

vector<int> a_star_search(const Graph& graph, int start, int goal,
    const vector<double>& h_scores, double& total_distance) {

    int V = graph.size();
    vector<double> g_score(V, INF);
    g_score[start] = 0.0;

    priority_queue<NodeState, vector<NodeState>, greater<NodeState>> open_set;
    open_set.push({ start, h_scores[start] });

    vector<int> came_from(V, -1);

    while (!open_set.empty()) {
        NodeState current_state = open_set.top();
        open_set.pop();
        int u = current_state.node_id;

        if (u == goal) {
            total_distance = g_score[u];
            return reconstruct_path(came_from, goal);
        }

        for (const auto& edge : graph[u]) {
            int v = edge.target_node;
            double edge_weight = edge.weight;


            if (edge_weight < 0) continue;

            double tentative_g_score = g_score[u] + edge_weight;

            if (tentative_g_score < g_score[v]) {
                came_from[v] = u;
                g_score[v] = tentative_g_score;

                double f_score = tentative_g_score + h_scores[v];
                open_set.push({ v, f_score });
            }
        }
    }

    total_distance = INF;
    return {};
}

// --- 3. ¿À√Œ–»“Ã ‘ÀŒ…ƒ¿-”Œ–ÿ≈ÀÀ¿ ---

void floyd_warshall(DistanceMatrix& dist) {
    int V = dist.size();

    for (int k = 0; k < V; ++k) {
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {

                if (dist[i][k] != INF && dist[k][j] != INF) {
                    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
                }
            }
        }
    }
}

bool check_negative_cycle(const DistanceMatrix& dist) {
    int V = dist.size();
    for (int i = 0; i < V; ++i) {
        if (dist[i][i] < 0) {
            return true;
        }
    }
    return false;
}

// --- 4. ¬—œŒÃŒ√¿“≈À‹Õ€≈ ‘”Õ ÷»» ƒÀﬂ “≈—“»–Œ¬¿Õ»ﬂ ---

pair<Graph, DistanceMatrix> generate_graphs(int V, double density, bool allow_negative) {
    Graph graph_astar(V);
    DistanceMatrix matrix_fw(V, vector<double>(V, INF));

    random_device rd;
    mt19937 gen(rd());

    uniform_real_distribution<> pos_weight_dist(1.0, 10.0);
    uniform_real_distribution<> fw_weight_dist(allow_negative ? -5.0 : 1.0, 10.0);

    uniform_real_distribution<> density_dist(0.0, 1.0);

    for (int i = 0; i < V; ++i) {
        matrix_fw[i][i] = 0.0;
        for (int j = 0; j < V; ++j) {
            if (i == j) continue;

            if (density_dist(gen) < density) {
                // √‡Ù ‰Îˇ A* (ÔÓÎÓÊËÚÂÎ¸Ì˚Â ‚ÂÒ‡)
                double pos_weight = pos_weight_dist(gen);
                graph_astar[i].push_back({ j, pos_weight });

                // Ã‡ÚËˆ‡ ‰Îˇ ‘ÎÓÈ‰‡-”Ó¯ÂÎÎ‡
                double fw_weight = fw_weight_dist(gen);
                matrix_fw[i][j] = fw_weight;
            }
        }
    }
    return { graph_astar, matrix_fw };
}



template<typename Func, typename Data>
double measure_time_single(Func algorithm_func, const Data& graph_data) {
    Data data_copy = graph_data;


    auto start = chrono::high_resolution_clock::now();

    algorithm_func(data_copy);

    auto end = chrono::high_resolution_clock::now();


    auto duration_nanos = chrono::duration_cast<chrono::nanoseconds>(end - start).count();

    return (double)duration_nanos / 1000000.0;
}


// --- 5. Œ—ÕŒ¬Õ¿ﬂ ‘”Õ ÷»ﬂ ƒÀﬂ “≈—“»–Œ¬¿Õ»ﬂ ---

void run_experiments() {

    cout << "=== —‡‚ÌËÚÂÎ¸Ì˚È ‡Ì‡ÎËÁ A* Ë ‘ÎÓÈ‰‡-”Ó¯ÂÎÎ‡ (C++) ===" << endl;
    cout << "V\tDensity\tA* (ms)\t\tFloyd-Warshall (ms)" << endl;
    cout << "--------------------------------------------------------" << endl;


    vector<int> vertex_counts = { 50, 100, 200, 500, 1000, 2000, 3000 };

    for (int V : vertex_counts) {


        double density_sparse = 0.05;
        auto [graph_astar_sparse, matrix_fw_sparse] = generate_graphs(V, density_sparse, true);

        vector<double> h_scores(V, 0.0);

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> node_dist(0, V - 1);
        int start_node = node_dist(gen);
        int end_node = node_dist(gen);
        double a_star_dist;



        auto measure_a_star = [&](const Graph& g) {
            a_star_search(g, start_node, end_node, h_scores, a_star_dist);
            };

        double time_a_star = measure_time_single(measure_a_star, graph_astar_sparse);

        auto measure_fw = [](DistanceMatrix& m) {
            floyd_warshall(m);
            };
        double time_fw = measure_time_single(measure_fw, matrix_fw_sparse);


        DistanceMatrix matrix_check_sparse = matrix_fw_sparse;
        floyd_warshall(matrix_check_sparse);
        bool negative_cycle_sparse = check_negative_cycle(matrix_check_sparse);

        cout << V << "\t" << density_sparse << "\t" << fixed << setprecision(3)
            << time_a_star << "\t\t" << time_fw << (negative_cycle_sparse ? " (NegCycle)" : "") << endl;


        if (V <= 2000) {
            double density_dense = 0.5;
            auto [graph_dense_astar, matrix_dense_fw] = generate_graphs(V, density_dense, false);
            double time_a_star_dense = measure_time_single(measure_a_star, graph_dense_astar);
            double time_fw_dense = measure_time_single(measure_fw, matrix_dense_fw);
            cout << V << "\t" << density_dense << "\t" << fixed << setprecision(3)
                << time_a_star_dense << "\t\t" << time_fw_dense << " (Dense)" << endl;
        }

    }
}

int main() {

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    cout.precision(6);
    run_experiments();

    return 0;
}
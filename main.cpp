#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>
#include <fstream>
#include "nlohmann/json.hpp"
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <chrono>
#include <cstdint>

#define START 0

using json = nlohmann::json;
using namespace lemon;
using namespace std;

struct Point {
    int x, y;
};

double distance(Point a, Point b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

struct Results
{
    vector<int> path;
    std::string name;
    double MST;
    double totalDistance;
    double firstEdge;
    double MDA;
    double time;
};

void to_json(nlohmann::json& j, const Results& r) {    
    j = nlohmann::json{{"name", r.name}, {"MDA", r.MDA}, {"Path", r.path}, {"Time", r.time}};
}

class Graph {
public:
    int V;
    vector<vector<double>> adjMatrix;
    double mstWeight;
    int targetNode;

    Graph(int V) : V(V), adjMatrix(V, vector<double>(V, 0)){
        mstWeight = 0;        
        targetNode = 0;
    }

    void addEdge(int u, int v, double w) {
        adjMatrix[u][v] = w;
        adjMatrix[v][u] = w;
    }

    vector<int> hamiltonianCircuit_Greedy() {
        vector<int> path;
        vector<bool> visited(V, false);        

        int currentNode = 0;
        visited[currentNode] = true;
        path.push_back(currentNode);
        for (int i = 1; i < V; i++) {
            int nextNode = -1;
            double minWeight = numeric_limits<double>::infinity();
            for (int j = 0; j < V; j++) {
                if (!visited[j] && adjMatrix[currentNode][j] < minWeight) {
                    minWeight = adjMatrix[currentNode][j];
                    nextNode = j;
                }
            }            
            visited[nextNode] = true;
            path.push_back(nextNode);
            currentNode = nextNode;
        }        
        path.push_back(START);        
        return path;
    }    

    vector<pair<int, int>> minimumSpanningTree(int forcedEdge=START) {
        vector<bool> inMST(V, false);
        vector<double> key(V, numeric_limits<double>::infinity());
        vector<int> parent(V, -1);
        key[0] = 0;
        parent[0] = -1;
        mstWeight = 0.0;

        // Use negative weight to force inclusion in MST
        if(forcedEdge != START){
            adjMatrix[START][forcedEdge] *= -1;
            adjMatrix[forcedEdge][START] *= -1;
            targetNode = forcedEdge;
        }

        for (int count = 0; count < V; ++count) {
            double minKey = numeric_limits<double>::infinity();
            int u;
            for (int v = 0; v < V; ++v) {
                if (!inMST[v] && key[v] < minKey) {
                    minKey = key[v];
                    u = v;
                }
            }

            inMST[u] = true;            
            mstWeight += abs(key[u]);

            for (int v = 0; v < V; ++v) {
                if (adjMatrix[u][v] && !inMST[v] && adjMatrix[u][v] < key[v]) {
                    key[v] = adjMatrix[u][v];
                    parent[v] = u;
                }
            }
        }

        // Cancel the forced edge
        if(forcedEdge != START){
            adjMatrix[START][forcedEdge] *= -1;
            adjMatrix[forcedEdge][START] *= -1;            
        }

        vector<pair<int, int>> mstEdges;
        for (int i = 1; i < V; ++i) {
            mstEdges.push_back({parent[i], i});
        }        
        return mstEdges;
    }

    vector<int> findOddDegreeVertices(const vector<pair<int, int>> &mstEdges) {
        vector<int> degree(V, 0);
        for (auto edge : mstEdges) {
            degree[edge.first]++;
            degree[edge.second]++;
        }
        vector<int> oddVertices;
        for (int i = 0; i < V; ++i) {
            if (degree[i] % 2 == 1) {
                oddVertices.push_back(i);
            }
        }
        return oddVertices;
    }

    vector<pair<int, int>> minimumWeightPerfectMatching (const vector<int> &oddVertices) {        
        int N = oddVertices.size();
        
        ListGraph graph;
        vector<ListGraph::Node> nodes(N);
        for (int i = 0; i < N; ++i) {
            nodes[i] = graph.addNode();
        }

        ListGraph::EdgeMap<double> weight(graph);
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                ListGraph::Edge e = graph.addEdge(nodes[i], nodes[j]);                
                weight[e] = adjMatrix[oddVertices[i]][oddVertices[j]]*(-1);
            }
        }
        
        MaxWeightedPerfectMatching<ListGraph, ListGraph::EdgeMap<double>> mwm(graph, weight);
        mwm.run();

        vector<pair<int, int>> result;
        for (ListGraph::EdgeIt e(graph); e != INVALID; ++e) {
            if (mwm.matching(e)) {
                int u = graph.id(graph.u(e));
                int v = graph.id(graph.v(e));
                result.push_back({oddVertices[u], oddVertices[v]});
            }
        }        

        return result;
    }
    vector<int> eulerianCircuit(const vector<pair<int, int>> &multigraphEdges, int firstNode=START) {
        struct CompareByValue {
            bool operator()(const std::pair<int, double>& lhs, const std::pair<int, double>& rhs) const {
                 // Sort in descending order of values
                return lhs.second > rhs.second;
            }
        };

        unordered_map<int, vector<pair<int, double>>> adjList;
        for (const auto& edge : multigraphEdges) {
            double weight = adjMatrix[edge.first][edge.second];
            adjList[edge.first].push_back(std::make_pair(edge.second, weight));                        
            adjList[edge.second].push_back(std::make_pair(edge.first, weight));
        }
        for (auto& node: adjList){
            std::sort(node.second.begin(), node.second.end(), CompareByValue());
        }

        // Vector to store the Eulerian circuit
        vector<int> circuit;

        // Stack to manage the current path in the graph        
        vector<int> currentPath;
        currentPath.push_back(START);
        if(firstNode != START){
            // manually add the first node to the path
            currentPath.push_back(firstNode);
            for(auto iter = adjList[START].begin(); iter != adjList[START].end(); ++iter){
                if(iter->first == firstNode){
                    adjList[START].erase(iter);                    
                    for(auto iter2 = adjList[firstNode].begin(); iter2 != adjList[firstNode].end(); ++iter2){
                        if(iter2->first == START){
                            adjList[firstNode].erase(iter2);
                            break;
                        }
                    }                    
                    break;
                }
            }
        }

        while (!currentPath.empty()) {
            int currentVertex = currentPath.back();

            // If the current vertex has unvisited edges
            if(!adjList[currentVertex].empty()){                
                int nextVertex = adjList[currentVertex].back().first;                
                currentPath.push_back(nextVertex);

                // Remove the edge from the graph
                adjList[currentVertex].pop_back();
                for(auto iter = adjList[nextVertex].begin(); iter != adjList[nextVertex].end(); ++iter){
                    if(iter->first == currentVertex){
                        adjList[nextVertex].erase(iter);
                        break;
                    }
                }                
            } else {
                // If the current vertex has no unvisited edges
                circuit.push_back(currentVertex);
                currentPath.pop_back();
            }
        }

        // Reverse the circuit to get the correct order
        reverse(circuit.begin(), circuit.end());
        return circuit;
    }   
    // vector<int> eulerianCircuit(const vector<pair<int, int>> &multigraphEdges, int firstNode) {
    //     unordered_map<int, list<int>> adjList;
    //     for (const auto& edge : multigraphEdges) {
    //         adjList[edge.first].push_back(edge.second);
    //         adjList[edge.second].push_back(edge.first);
    //     }

    //     // Vector to store the Eulerian circuit
    //     vector<int> circuit;

    //     // Stack to manage the current path in the graph        
    //     vector<int> currentPath;
    //     currentPath.push_back(START);
    //     if(firstNode != START){
    //         // manually add the first node to the path
    //         currentPath.push_back(firstNode);
    //         adjList[START].erase(find(adjList[START].begin(), adjList[START].end(), firstNode));
    //         adjList[firstNode].erase(find(adjList[firstNode].begin(), adjList[firstNode].end(), START));
    //     }

    //     while (!currentPath.empty()) {
    //         int currentVertex = currentPath.back();

    //         // If the current vertex has unvisited edges
    //         if(!adjList[currentVertex].empty()){
    //             int nextVertex = adjList[currentVertex].back();
    //             currentPath.push_back(nextVertex);

    //             // Remove the edge from the graph                
    //             adjList[currentVertex].pop_back();
    //             adjList[nextVertex].erase(find(adjList[nextVertex].begin(), adjList[nextVertex].end(), currentVertex));
    //         } else {
    //             // If the current vertex has no unvisited edges
    //             circuit.push_back(currentVertex);
    //             currentPath.pop_back();
    //         }
    //     }

    //     // Reverse the circuit to get the correct order
    //     reverse(circuit.begin(), circuit.end());
    //     return circuit;
    // }    

    vector<int> hamiltonianCircuit(const vector<int> &eulerianCircuit) {
        vector<bool> visited(V, false);
        vector<int> hamiltonianCircuit;
        for (int v : eulerianCircuit) {
            if (!visited[v]) {
                visited[v] = true;                
                hamiltonianCircuit.push_back(v);                
            }
        }        
        hamiltonianCircuit.push_back(hamiltonianCircuit.front());
        return hamiltonianCircuit;
    }

    Results analyzePath(vector<int> &path){
        Results results;
        results.MST = mstWeight;
        results.path = path;
        results.totalDistance = 0;
        results.firstEdge = 0;
        results.MDA = 0;
        
        for (int i = 0; i < path.size() - 1; ++i) {
            results.totalDistance += adjMatrix[path[i]][path[i + 1]];            
        }
        results.firstEdge = adjMatrix[path[0]][path[1]];
        results.MDA = 2*results.totalDistance - results.firstEdge;
        return results;
    }       
    int getNumNodes(){
        return V;    
    }    
};

void generatePathGreedy(Graph &g, nlohmann::json &output){

    auto start = std::chrono::high_resolution_clock::now();
    vector<int> greedyPath = g.hamiltonianCircuit_Greedy();
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate the duration in milliseconds
    std::chrono::duration<double, std::milli> duration = end - start;

    Results r = g.analyzePath(greedyPath);
    r.name = "Greedy";
    r.time = duration.count();    
    output["Greedy"] = {{"MDA", r.MDA}, {"Path", r.path}, {"Time", r.time}};
    //output.push_back(r);

    #ifdef DISPLAY
    cout << "Greedy Path: ";
    for (int v : greedyPath) {
        cout << v << " ";
    }
    cout << endl;
    cout << "Hamiltonian Weight: " << r.totalDistance << endl;
    cout << "Length of the first edge:" << r.firstEdge << endl;
    cout << "MDA: " << r.MDA << endl;
    cout << endl;
    #endif
}

void generatePathMST(Graph &g, nlohmann::json &output){
    auto start = std::chrono::high_resolution_clock::now();
    vector<pair<int, int>> mstEdges = g.minimumSpanningTree();
    vector<int> oddVertices = g.findOddDegreeVertices(mstEdges);    
    vector<pair<int, int>> matching = g.minimumWeightPerfectMatching(oddVertices);
    
    vector<pair<int, int>> multigraphEdges = mstEdges;
    multigraphEdges.insert(multigraphEdges.end(), matching.begin(), matching.end());
       
    vector<int> eulerianCircuit = g.eulerianCircuit(multigraphEdges);
    vector<int> hamiltonianCircuit = g.hamiltonianCircuit(eulerianCircuit);

    reverse(hamiltonianCircuit.begin(), hamiltonianCircuit.end());
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate the duration in milliseconds
    std::chrono::duration<double, std::milli> duration = end - start;

    Results r = g.analyzePath(hamiltonianCircuit);
    r.name = "MST";
    r.time = duration.count();
    output["MST"] = {{"MDA", r.MDA}, {"Path", r.path}, {"Time", r.time}};    
    //output.push_back(r);

    #ifdef DISPLAY
    cout << "Hamiltonian Cycle: ";
    for (int v : hamiltonianCircuit) {
        cout << v << " ";
    }
    cout << endl;
    cout << "MST Weight: " << r.MST << endl;
    cout << "Hamiltonian Weight: " << r.totalDistance << endl;
    cout << "Length of the first edge:" << r.firstEdge << endl;
    cout << "MDA: " << r.MDA << endl;
    cout << endl;
    #endif
}

void generatePathEnforce(Graph &g, nlohmann::json &output){

    double minMDA = numeric_limits<double>::infinity();
    int optimalTarget = 0;
    Results optimal;

    auto start = std::chrono::high_resolution_clock::now();
    for (int target = 1; target < g.getNumNodes(); ++target) {                
        #ifdef DISPLAY
        cout << "==Enforcing edge " << START << " -> " << target << "==" << endl;
        #endif
        
        vector<pair<int, int>> mstEdges = g.minimumSpanningTree(target);
        vector<int> oddVertices = g.findOddDegreeVertices(mstEdges);    
        vector<pair<int, int>> matching = g.minimumWeightPerfectMatching(oddVertices);
    

        vector<pair<int, int>> multigraphEdges = mstEdges;
        multigraphEdges.insert(multigraphEdges.end(), matching.begin(), matching.end());
       
        vector<int> eulerianCircuit = g.eulerianCircuit(multigraphEdges, target);
        vector<int> hamiltonianCircuit = g.hamiltonianCircuit(eulerianCircuit);

        Results r = g.analyzePath(hamiltonianCircuit);        

        #ifdef DISPLAY
        cout << "Hamiltonian Cycle: ";
        for (int v : hamiltonianCircuit) {
            cout << v << " ";
        }
        cout << endl;
        cout << "MST Weight: " << r.MST << endl;
        cout << "Hamiltonian Weight: " << r.totalDistance << endl;
        cout << "Length of the first edge:" << r.firstEdge << endl;
        cout << "MDA: " << r.MDA << endl;
        cout << endl;
        #endif
        
        if (r.MDA < minMDA) {
            minMDA = r.MDA;
            optimalTarget = target;
            optimal = r;
        }        
    }
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate the duration in milliseconds
    std::chrono::duration<double, std::milli> duration = end - start;

    optimal.name = "Enforced";
    optimal.time = duration.count();    
    output["Enforced"] = {{"MDA", optimal.MDA}, {"Path", optimal.path}, {"Time", optimal.time}};
    //output.push_back(optimal);
    #ifdef DISPLAY
    cout << "Optimal first node: " << optimalTarget;    
    cout << endl << endl;
    #endif
}

pair<double, vector<int>> heldKarp(vector<vector<double>> &dist, int n, int start, int end) {
    uint64_t subset_size = 1ULL << n; // 2^n subsets
    vector<vector<double>> dp(subset_size, vector<double>(n, numeric_limits<double>::infinity()));
    vector<vector<int>> parent(subset_size, vector<int>(n, -1));

    dp[1ULL << start][start] = 0;

    for (int mask = 0; mask < subset_size; ++mask) {
        for (int i = 0; i < n; ++i) {
            if (!(mask & (1 << i))) continue;
            for (int j = 0; j < n; ++j) {
                if (mask & (1 << j) || dist[i][j] == numeric_limits<double>::infinity()) continue;
                int new_mask = mask | (1 << j);
                if (dp[new_mask][j] > dp[mask][i] + dist[i][j]) {
                    dp[new_mask][j] = dp[mask][i] + dist[i][j];
                    parent[new_mask][j] = i;
                }
            }
        }
    }

    // for n >= 32: 
    // uint64_t subset_size = 1ULL << n; // 2^n subsets
    // unordered_map<uint64_t, vector<double>> dp;
    // unordered_map<uint64_t, vector<int>> parent;
    
    // dp.insert({1ULL << start, vector<double>(n, numeric_limits<double>::infinity())});
    // parent.insert({1ULL << start, vector<int>(n, -1)});
    // dp[1ULL << start][start] = 0;

    // for (uint64_t mask = 0; mask < subset_size; ++mask) {
    //     for (int i = 0; i < n; ++i) {
    //         if (!(mask & (1ULL << i))) continue;
    //         for (int j = 0; j < n; ++j) {
    //             if (mask & (1ULL << j) || dist[i][j] == numeric_limits<double>::infinity()) continue;
    //             if (dp.find(mask) == dp.end()) continue;
    //             uint64_t new_mask = mask | (1ULL << j);
    //             dp.insert({new_mask, vector<double>(n, numeric_limits<double>::infinity())});
    //             parent.insert({new_mask, vector<int>(n, -1)});
    //             if (dp[new_mask][j] > dp[mask][i] + dist[i][j]) {
    //                 dp[new_mask][j] = dp[mask][i] + dist[i][j];
    //                 parent[new_mask][j] = i;
    //             }
    //         }
    //     }
    // }


    // Retrieve the shortest path
    vector<int> path;
    uint64_t mask = (1ULL << n) - 1;
    int cur = end;
    while (cur != start) {
        path.push_back(cur);
        cur = parent[mask][cur];
        mask &= ~(1ULL << path.back());
    }
    path.push_back(start);
    path.push_back(0);
    reverse(path.begin(), path.end());

    return {dp[(1ULL << n) - 1][end], path};
}

void generateOptimal(Graph &g, vector<vector<double>> &dist, int n, nlohmann::json &output){
    double minCost = numeric_limits<double>::infinity();
    int bestFirstNode = 0;
    vector<int> optPath;
    auto start = std::chrono::high_resolution_clock::now();
    for(int target=1; target < n; target++){
        auto path = heldKarp(dist, n, target, 0); // alwasys end at 0
        double MDA = dist[0][target] + 2*path.first;
        if(minCost > MDA){
            minCost = MDA;
            bestFirstNode = target;
            optPath = path.second;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate the duration in milliseconds
    std::chrono::duration<double, std::milli> duration = end - start;

    // output
    Results optimal = g.analyzePath(optPath);
    optimal.name = "Optimal";
    optimal.time = duration.count();    
    output["Optimal"] = {{"MDA", optimal.MDA}, {"Path", optimal.path}, {"Time", optimal.time}};
    //output.push_back(optimal);

    #ifdef DISPLAY
    cout << "Optimal MDA: " << minCost << endl;
    cout << "Optimal Path: ";
    for (int node : optPath) {
        cout << node << " ";
    }
    cout << endl << endl;
    #endif
}

int main(int argc, char *argv[]) {
    // Read the JSON file
    //ifstream file("even-3-10.json");    
    ifstream file(argv[1]);
    json j;
    file >> j;

    std::ofstream output("output.json", std::ios::out);
    nlohmann::json results = nlohmann::json::array();

    int M = j["cases"];

    for (int scenarioN = 0; scenarioN < M; scenarioN++){        
        cout << "[Scenario " << scenarioN+1 << "]" << endl;                
        nlohmann::json result;

        int N = j["scenarios"][scenarioN]["N"];
        vector<Point> points(N);
        for (int i = 0; i < N; ++i) {
            points[i].x = j["scenarios"][scenarioN]["points"][i]["x"];
            points[i].y = j["scenarios"][scenarioN]["points"][i]["y"];
        }

        Graph g(N);
        vector<vector<double>> dist(N, vector<double>(N)); // Distance matrix for the Held-Karp algorithm
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                double w = distance(points[i], points[j]);
                g.addEdge(i, j, w);
                dist[i][j] = w;
                dist[j][i] = w;
            }
        }

        generatePathGreedy(g, result);
        generatePathMST(g, result);
        generatePathEnforce(g, result);
        generateOptimal(g, dist, N, result);
        results.push_back(result);
    }

    output << results.dump(4) << endl;
    output.close();

    return 0;
}
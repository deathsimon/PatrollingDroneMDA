#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>
#include <tuple>
#include <fstream>
#include "nlohmann/json.hpp"
#include <lemon/list_graph.h>
#include <lemon/matching.h>

using json = nlohmann::json;
using namespace lemon;
using namespace std;

struct Point {
    int x, y;
};

double distance(Point a, Point b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

class Graph {
public:
    int V;
    vector<vector<double>> adjMatrix;
    double mstWeight, HamiltonianWeight;

    Graph(int V) : V(V), adjMatrix(V, vector<double>(V, 0)){
        mstWeight = 0;
        HamiltonianWeight = 0;
    }

    void addEdge(int u, int v, double w) {
        adjMatrix[u][v] = w;
        adjMatrix[v][u] = w;
    }

    void addForcedEdge(int u, int v, double w) {
        adjMatrix[u][v] = -w; // Use negative weight to force inclusion in MST
        adjMatrix[v][u] = -w;
    }

    vector<pair<int, int>> minimumSpanningTree() {
        vector<bool> inMST(V, false);
        vector<double> key(V, numeric_limits<double>::infinity());
        vector<int> parent(V, -1);
        key[0] = 0;
        parent[0] = -1;
        mstWeight = 0.0;

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
                // assert(adjMatrix[oddVertices[i]][oddVertices[j]] >= 0);
                weight[e] = abs(adjMatrix[oddVertices[i]][oddVertices[j]])*(-1);
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

    vector<int> eulerianCircuit(const vector<pair<int, int>> &multigraphEdges) {
        vector<vector<int>> adjList(V);
        for (auto edge : multigraphEdges) {
            adjList[edge.first].push_back(edge.second);
            adjList[edge.second].push_back(edge.first);
        }

        vector<int> circuit;
        vector<int> currentPath;
        vector<bool> visitedEdges(multigraphEdges.size(), false);
        currentPath.push_back(0);
        while (!currentPath.empty()) {
            int u = currentPath.back();
            bool found = false;

            for (auto it = adjList[u].begin(); it != adjList[u].end(); ++it) {
                int v = *it;
                adjList[u].erase(it);
                currentPath.push_back(v);
                found = true;
                break;
            }

            if (!found) {
                circuit.push_back(u);
                currentPath.pop_back();
            }
        }
        //reverse(circuit.begin(), circuit.end());
        return circuit;
    }

    vector<int> hamiltonianCircuit(const vector<int> &eulerianCircuit) {
        vector<bool> visited(V, false);
        vector<int> hamiltonianCircuit;
        HamiltonianWeight = 0.0;
        for (int v : eulerianCircuit) {
            if (!visited[v]) {
                visited[v] = true;
                if(!hamiltonianCircuit.empty()){
                    HamiltonianWeight += abs(adjMatrix[hamiltonianCircuit.back()][v]);
                }
                hamiltonianCircuit.push_back(v);                
            }
        }
        HamiltonianWeight += abs(adjMatrix[hamiltonianCircuit.back()][hamiltonianCircuit.front()]);
        hamiltonianCircuit.push_back(hamiltonianCircuit.front());
        return hamiltonianCircuit;
    }

    double getMSTWeight() {
        return mstWeight;
    }
    double getHamiltonianWeight() {
        return HamiltonianWeight;
    }
    double getMDA() {
        double firstEdge = 0;
        for(int i = 0; i < V; i++){
            if(adjMatrix[0][i] < 0){
                firstEdge = abs(adjMatrix[0][i]);
                break;
            }            
        }        
        
        return 2*HamiltonianWeight - firstEdge;
    }
};

int main() {
    // Read the JSON file
    ifstream file("input.json");
    json j;
    file >> j;

    int N = j["N"];
    int start = j["start"];
    int end = j["end"];
    vector<Point> points(N);
    for (int i = 0; i < N; ++i) {
        points[i].x = j["points"][i]["x"];
        points[i].y = j["points"][i]["y"];
    }

    Graph g(N);
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double w = distance(points[i], points[j]);
            g.addEdge(i, j, w);
        }
    }

    // Add the forced edge with high priority (negative weight)
    g.addForcedEdge(start, end, distance(points[start], points[end]));

    vector<pair<int, int>> mstEdges = g.minimumSpanningTree();
    vector<int> oddVertices = g.findOddDegreeVertices(mstEdges);    
    vector<pair<int, int>> matching = g.minimumWeightPerfectMatching(oddVertices);
    

    vector<pair<int, int>> multigraphEdges = mstEdges;
    multigraphEdges.insert(multigraphEdges.end(), matching.begin(), matching.end());

    // Add the forced edge again in the multigraph to ensure it's part of the circuit
    multigraphEdges.push_back({start, end});

    vector<int> eulerianCircuit = g.eulerianCircuit(multigraphEdges);
    vector<int> hamiltonianCircuit = g.hamiltonianCircuit(eulerianCircuit);

    cout << "Hamiltonian Cycle: ";
    for (int v : hamiltonianCircuit) {
        cout << v << " ";
    }
    cout << endl;

    // assert(g.adjMatrix[hamiltonianCircuit[0]][hamiltonianCircuit[1]] == distance(points[start], points[end]));
    cout << "MST Weight: " << g.getMSTWeight() << endl;
    cout << "Hamiltonian Weight: " << g.getHamiltonianWeight() << endl;
    cout << "Length of the first edge:" << distance(points[start], points[end]) << endl;
    cout << "MDA: " << g.getMDA() << endl;

    return 0;
}

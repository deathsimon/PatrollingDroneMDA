#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>
#include <tuple>
#include <utility>
#include <fstream>
#include "nlohmann/json.hpp"
#include <lemon/list_graph.h>
#include <lemon/matching.h>

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

class Graph {
public:
    int V;
    vector<vector<double>> adjMatrix;
    double mstWeight, HamiltonianWeight, GreedyWeight;
    int targetNode;

    Graph(int V) : V(V), adjMatrix(V, vector<double>(V, 0)){
        mstWeight = 0;
        HamiltonianWeight = 0;
        targetNode = 0;
    }

    void addEdge(int u, int v, double w) {
        adjMatrix[u][v] = w;
        adjMatrix[v][u] = w;
    }

    vector<int> hamiltonianCircuit_Greedy() {
        vector<int> path;
        vector<bool> visited(V, false);
        GreedyWeight = 0;

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
            GreedyWeight += minWeight;            
            visited[nextNode] = true;
            path.push_back(nextNode);
            currentNode = nextNode;
        }
        GreedyWeight += adjMatrix[path.back()][START];
        path.push_back(START);
        cout << "Greedy Cycle: " << GreedyWeight << endl;
        return path;
    }    

    vector<pair<int, int>> minimumSpanningTree(int forcedEdge) {
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
    vector<int> eulerianCircuit(const vector<pair<int, int>> &multigraphEdges, int firstNode) {
        struct CompareByValue {
            bool operator()(const std::pair<int, double>& lhs, const std::pair<int, double>& rhs) const {
                 // Sort in descending order of values
                return lhs.second > rhs.second;
            }
        };

        unordered_map<int, vector<pair<int, double>>> adjList;
        for (const auto& edge : multigraphEdges) {
            double weight = abs(adjMatrix[edge.first][edge.second]);            
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
    double getGreedyWeight() {
        return GreedyWeight;
    }
    double getMDA() {
        return 2*HamiltonianWeight - adjMatrix[0][targetNode];
    }
    int getNumNodes(){
        return V;    
    }
    double getDistance(int u, int v){
        return adjMatrix[u][v];
    }
};

void generatePath_greedy(Graph &g){
    vector<int> greedyPath = g.hamiltonianCircuit_Greedy();
    cout << "Greedy Path: ";
    for (int v : greedyPath) {
        cout << v << " ";
    }
    cout << endl;
    cout << "Length of the first edge:" << g.getDistance(START,greedyPath[1]) << endl;            
    cout << "MDA: " << 2*g.getGreedyWeight() - g.getDistance(START,greedyPath[1]) << endl;
}

void generatePath(Graph &g){

    double minMDA = numeric_limits<double>::infinity();
    int optimalTarget = 0;

    for (int target = 0; target < g.getNumNodes(); ++target) {
        
        // Add the forced edge with high priority (negative weight)
        //g.addForcedEdge(START, target);
        if (target == 0) {
            cout << "Without enforcing any edge"  << endl;
        }
        else{
            cout << "==Enforcing edge " << START << " -> " << target << "==" << endl;
        }

        vector<pair<int, int>> mstEdges = g.minimumSpanningTree(target);
        vector<int> oddVertices = g.findOddDegreeVertices(mstEdges);    
        vector<pair<int, int>> matching = g.minimumWeightPerfectMatching(oddVertices);
    

        vector<pair<int, int>> multigraphEdges = mstEdges;
        multigraphEdges.insert(multigraphEdges.end(), matching.begin(), matching.end());
       
        vector<int> eulerianCircuit = g.eulerianCircuit(multigraphEdges, target);
        vector<int> hamiltonianCircuit = g.hamiltonianCircuit(eulerianCircuit);

        cout << "Hamiltonian Cycle: ";
        for (int v : hamiltonianCircuit) {
            cout << v << " ";
        }
        cout << endl;

        cout << "MST Weight: " << g.getMSTWeight() << endl;
        cout << "Hamiltonian Weight: " << g.getHamiltonianWeight() << endl;
        double MDA = g.getMDA();
        if(target != 0){
            cout << "Length of the first edge:" << g.getDistance(START,target) << endl;            
        }
        else{            
            MDA -= g.getDistance(START, hamiltonianCircuit[1]);
        }
        cout << "MDA: " << MDA << endl;

        if (MDA < minMDA) {
            minMDA = MDA;
            optimalTarget = target;
        }
        
    }

    if (optimalTarget == START)
    {
        cout << endl << "TSP leads to the optimal result" << endl;
    }
    else{
        cout << endl << "Optimal first node: " << optimalTarget << endl;
    }
}

int main() {
    // Read the JSON file
    ifstream file("input.json");
    json j;
    file >> j;

    int M = j["cases"];

    for (int scenarioN = 0; scenarioN < M; scenarioN++){         

        cout << "[Scenario " << scenarioN << "]" << endl;

        int N = j["scenarios"][scenarioN]["N"];
        vector<Point> points(N);
        for (int i = 0; i < N; ++i) {
            points[i].x = j["scenarios"][scenarioN]["points"][i]["x"];
            points[i].y = j["scenarios"][scenarioN]["points"][i]["y"];
        }

        Graph g(N);
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                double w = distance(points[i], points[j]);
                g.addEdge(i, j, w);
            }
        }
        generatePath_greedy(g);
        generatePath(g);
    }

    

    return 0;
}

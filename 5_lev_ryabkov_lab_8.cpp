#include <iostream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>
#include <limits>
#include <cmath>
#include <sstream>
#include <fstream>
#include <chrono>
#include <functional>
#include <cassert>

using namespace std;
using namespace chrono;

// Узел графа
struct Node {
    double lat, lon;
    vector<pair<Node*, double>> connections;
};

// Граф
class Graph {
public:
    unordered_map<string, Node> nodes;

    void loadFromFile(const string& filepath) {
        ifstream file(filepath);  // Открытие файла — O(1) по памяти
        if (!file.is_open()) {
            cerr << "Ошибка: не удалось открыть файл " << filepath << "\n";  // O(1) по сложности
            return;
        }

        string line;  // O(1) по памяти
        while (getline(file, line)) {  // O(n) по сложности, где n — количество строк в файле // Память: O(n * m), где m — длина строки
            size_t pos = line.find(':');  // O(m) по сложности для поиска символа ':' в строке
            if (pos == string::npos) continue;  // O(1) по сложности

            string mainNodeData = line.substr(0, pos);  // O(m) по сложности для извлечения подстроки // Память: O(m)
            double lon, lat;
            sscanf(mainNodeData.c_str(), "%lf,%lf", &lon, &lat);  // O(1) по сложности, память: 16 байт

            string key = toKey(lat, lon);  // O(1) по сложности // Память: O(1)
            Node& mainNode = nodes[key];  // O(1) по сложности для доступа и/или вставки в хеш-таблицу

            mainNode.lat = lat;  // O(1) по сложности // Память: 8 байт
            mainNode.lon = lon;  // O(1) по сложности // Память: 8 байт

            string neighborsData = line.substr(pos + 1);  // O(m) по сложности // Память: O(m)
            stringstream ss(neighborsData);  // O(m) по сложности // Память: O(m)
            string neighborEntry;

            while (getline(ss, neighborEntry, ';')) {  // O(k) по сложности для всех соседей // Память: O(k * l), где l — длина каждой записи
                double nLon, nLat, dist;
                sscanf(neighborEntry.c_str(), "%lf,%lf,%lf", &nLon, &nLat, &dist);  // O(1) по сложности // Память: 24 байта
                string nKey = toKey(nLat, nLon);  // O(1) по сложности // Память: O(1)
                Node& neighbor = nodes[nKey];  // O(1) по сложности // Память: O(1)
                neighbor.lat = nLat;  // O(1) по сложности // Память: 8 байт
                neighbor.lon = nLon;  // O(1) по сложности // Память: 8 байт

                mainNode.connections.push_back({&neighbor, dist});  // O(1) по сложности для добавления соединения // Память: 16 байт на каждое соединение
            }
        }
    }


    Node* findNearest(double lat, double lon) {
        Node* closest = nullptr;
        double smallestDist = numeric_limits<double>::max();  // O(1) по сложности // Память: 8 байт

        for (auto& [key, node] : nodes) {  // O(n) по сложности — перебор всех узлов в графе // Память: O(n * 16)
            double dist = distance(lat, lon, node.lat, node.lon);  // O(1) по сложности // Память: 8 байт
            if (dist < smallestDist) {  // O(1) по сложности // Память: 8 байт
                smallestDist = dist;  // O(1) по сложности // Память: 8 байт
                closest = &node;  // O(1) по сложности // Память: 8 байт
            }
        }
        return closest;  // O(1) по сложности // Память: 8 байт
    }


private:
    static string toKey(double lat, double lon) {
        stringstream ss;
        ss << lat << "," << lon;
        return ss.str();
    }

    static double distance(double lat1, double lon1, double lat2, double lon2) {
        return sqrt((lat1 - lat2) * (lat1 - lat2) + (lon1 - lon2) * (lon1 - lon2));
    }
};

// Поиск в ширину
vector<Node*> bfs(Graph& graph, Node* start, Node* goal) {
    if (!start || !goal) return {};  // O(1) по сложности

    unordered_map<Node*, Node*> parentMap;  // O(n) по памяти // Память: O(n * 16)
    set<Node*> visited;  // O(n) по памяти // Память: O(n * 16)
    queue<Node*> queue;  // O(n) по памяти // Память: O(n * 8)

    queue.push(start);  // O(1) по сложности // Память: 8 байт
    visited.insert(start);  // O(1) по сложности // Память: 8 байт
    parentMap[start] = nullptr;  // O(1) по сложности // Память: 16 байт

    while (!queue.empty()) {  // O(n) по сложности, т.к. перебираются все узлы
        Node* current = queue.front();  // O(1) по сложности // Память: 8 байт
        queue.pop();  // O(1) по сложности

        if (current == goal) break;  // O(1) по сложности

        for (auto& [neighbor, weight] : current->connections) {  // O(k) по сложности для всех соседей
            if (visited.find(neighbor) == visited.end()) {  // O(1) по сложности
                visited.insert(neighbor);  // O(1) по сложности // Память: 8 байт
                queue.push(neighbor);  // O(1) по сложности // Память: 8 байт
                parentMap[neighbor] = current;  // O(1) по сложности // Память: 16 байт
            }
        }
    }

    vector<Node*> path;  // O(k) по памяти // Память: O(k * 8)
    for (Node* node = goal; node != nullptr; node = parentMap[node]) {  // O(k) по сложности
        path.push_back(node);  // O(1) по сложности // Память: 8 байт
    }
    reverse(path.begin(), path.end());  // O(k) по сложности
    return path;  // O(k) по сложности // Память: O(k * 8)
}

// Поиск в глубину
vector<Node*> dfs(Graph& graph, Node* start, Node* goal) {
    if (!start || !goal) return {};  // O(1) по сложности // Память: O(1)

    unordered_map<Node*, Node*> parentMap;  // O(n) по памяти // Память: O(n * 16)
    set<Node*> visited;  // O(n) по памяти // Память: O(n * 16)
    stack<Node*> stack;  // O(n) по памяти // Память: O(n * 8)

    stack.push(start);  // O(1) по сложности // Память: O(8)
    parentMap[start] = nullptr;  // O(1) по сложности // Память: O(16)

    while (!stack.empty()) {  // O(n) по сложности // Память: O(1) для каждой итерации
        Node* current = stack.top();  // O(1) по сложности // Память: O(8)
        stack.pop();  // O(1) по сложности

        if (visited.find(current) != visited.end()) continue;  // O(1) по сложности
        visited.insert(current);  // O(1) по сложности // Память: O(8) для каждого элемента

        if (current == goal) break;  // O(1) по сложности

        for (auto& [neighbor, weight] : current->connections) {  // O(k) по сложности для всех соседей
            if (visited.find(neighbor) == visited.end()) {  // O(1) по сложности
                stack.push(neighbor);  // O(1) по сложности // Память: O(8)
                parentMap[neighbor] = current;  // O(1) по сложности // Память: O(16)
            }
        }
    }

    vector<Node*> path;  // O(k) по памяти // Память: O(k * 8)
    for (Node* node = goal; node != nullptr; node = parentMap[node]) {  // O(k) по сложности
        path.push_back(node);  // O(1) по сложности // Память: O(8) для каждого элемента
    }
    reverse(path.begin(), path.end());  // O(k) по сложности
    return path;  // O(k) по сложности // Память: O(k * 8)
}

// Алгоритм Дейкстры
vector<Node*> dijkstra(Graph& graph, Node* start, Node* goal) {
    if (!start || !goal) return {};  // O(1) по сложности // Память: O(1)

    unordered_map<Node*, double> costs;  // O(n) по памяти // Память: O(n * 8)
    unordered_map<Node*, Node*> parentMap;  // O(n) по памяти // Память: O(n * 16)
    priority_queue<pair<double, Node*>, vector<pair<double, Node*>>, greater<>> pq;  // O(n) по памяти // Память: O(n * 16)

    for (auto& [key, node] : graph.nodes) {  // O(n) по сложности // Память: O(n * 8)
        costs[&node] = numeric_limits<double>::infinity();
    }
    costs[start] = 0;
    pq.push({0, start});  // O(1) по сложности // Память: O(8)

    while (!pq.empty()) {  // O(n * log n) по сложности, т.к. очередь приоритетов
        auto [currentCost, current] = pq.top();  // O(1) по сложности // Память: O(8)
        pq.pop();  // O(log n) по сложности

        if (current == goal) break;  // O(1) по сложности

        for (auto& [neighbor, weight] : current->connections) {  // O(k) по сложности для каждого соседа
            double newCost = currentCost + weight;  // O(1) по сложности
            if (newCost < costs[neighbor]) {  // O(1) по сложности
                costs[neighbor] = newCost;  // O(1) по сложности
                parentMap[neighbor] = current;  // O(1) по сложности
                pq.push({newCost, neighbor});  // O(log n) по сложности
            }
        }
    }

    vector<Node*> path;  // O(k) по памяти // Память: O(k * 8)
    for (Node* node = goal; node != nullptr; node = parentMap[node]) {  // O(k) по сложности
        path.push_back(node);  // O(1) по сложности // Память: O(8)
    }
    reverse(path.begin(), path.end());  // O(k) по сложности
    return path;  // O(k) по сложности // Память: O(k * 8)
}


double measureTime(function<void()> func) {
    auto startTime = high_resolution_clock::now();
    func();
    auto endTime = high_resolution_clock::now();
    return duration<double>(endTime - startTime).count();
}

void bfsTime(Graph& graph, Node* start, Node* target) {
    auto bfsPath = bfs(graph, start, target);
}

void dfsTime(Graph& graph, Node* start, Node* target) {
    auto dfsPath = dfs(graph, start, target);
}

void dijkstraTime(Graph& graph, Node* start, Node* target) {
    auto dijkstraPath = dijkstra(graph, start, target);
}


void test_bfs() {

    Graph graph;
    graph.loadFromFile("test.txt");

    double startLon = 1.0, startLat = 1.0;
    double goalLon = 2.0, goalLat = 2.0;

    Node* start = graph.findNearest(startLat, startLon);
    Node* goal = graph.findNearest(goalLat, goalLon);

    vector<Node*> bfsPath = bfs(graph, start, goal);

    assert(!bfsPath.empty());

    double weight_bfs = 0.0;
    for (size_t i = 0; i < bfsPath.size() - 1; ++i) {
        for (auto& [neighbor, weight] : bfsPath[i]->connections) {
            if (neighbor == bfsPath[i + 1]) {
                weight_bfs += weight;
                break;
            }
        }
    }

    assert(weight_bfs == 22.0);
    cout << "BFS passed!" << endl;
}

void test_dfs() {
    Graph graph;
    graph.loadFromFile("test.txt");

    double startLon = 1.0, startLat = 1.0;
    double goalLon = 2.0, goalLat = 2.0;

    Node* start = graph.findNearest(startLat, startLon);
    Node* goal = graph.findNearest(goalLat, goalLon);

    vector<Node*> dfsPath = dfs(graph, start, goal);

    assert(!dfsPath.empty());

    double weight_dfs = 0.0;
    for (size_t i = 0; i < dfsPath.size() - 1; ++i) {
        for (auto& [neighbor, weight] : dfsPath[i]->connections) {
            if (neighbor == dfsPath[i + 1]) {
                weight_dfs += weight;
                break;
            }
        }
    }

    assert(weight_dfs == 74.0);
    cout << "DFS passed!" << endl;
}


void test_dijkstra() {

    Graph graph;
    graph.loadFromFile("test.txt");

    double startLon = 1.0, startLat = 1.0;
    double goalLon = 2.0, goalLat = 2.0;

    Node* start = graph.findNearest(startLat, startLon);
    Node* goal = graph.findNearest(goalLat, goalLon);

    vector<Node*> dijkstraPath = dijkstra(graph, start, goal);

    assert(!dijkstraPath.empty());

    double weight_dijkstra = 0.0;
    for (size_t i = 0; i < dijkstraPath.size() - 1; ++i) {
        for (auto& [neighbor, weight] : dijkstraPath[i]->connections) {
            if (neighbor == dijkstraPath[i + 1]) {
                weight_dijkstra += weight;
                break;
            }
        }
    }

    assert(weight_dijkstra == 18.0);
    cout << "Dijkstra passed!" << endl;
}


void testPerformance(Graph& graph, Node* start, Node* target) {
    cout << "Testing Performance...\n";
    double bfsTimeResult = measureTime([&]() { bfsTime(graph, start, target); });
    double dfsTimeResult = measureTime([&]() { dfsTime(graph, start, target); });
    double dijkstraTimeResult = measureTime([&]() { dijkstraTime(graph, start, target); });
}


void test_algorithms_time() {

    Graph graph;
    graph.loadFromFile("test.txt");

    double startLon = 1.0, startLat = 1.0;
    double goalLon = 2.0, goalLat = 2.0;

    Node* start = graph.findNearest(startLat, startLon);
    Node* goal = graph.findNearest(goalLat, goalLon);

    auto bfs_time = measureTime([&]() { bfs(graph, start, goal); });
    auto dfs_time = measureTime([&]() { dfs(graph, start, goal); });
    auto dijkstra_time = measureTime([&]() { dijkstra(graph, start, goal); });

    cout << "BFS Time: " << bfs_time << " seconds\n";
    cout << "DFS Time: " << dfs_time << " seconds\n";
    cout << "Dijkstra Time: " << dijkstra_time << " seconds\n";

    cout << "Algorithm time tests passed!" << endl;
}


int main() {

    test_algorithms_time();

    Graph graph;
    graph.loadFromFile("/Users/levr/CLionProjects/aboba/spb_graph.txt");

    Node* start = graph.findNearest(59.845388, 30.361925);
    Node* target = graph.findNearest(59.956566, 30.308424);

    cout << "BFS Time: " << measureTime([&]() { bfsTime(graph, start, target); }) << " seconds\n";
    cout << "DFS Time: " << measureTime([&]() { dfsTime(graph, start, target); }) << " seconds\n";
    cout << "Dijkstra Time: " << measureTime([&]() { dijkstraTime(graph, start, target); }) << " seconds\n";

    return 0;
}

// findNearest(): Время: O(n), Память: O(1)
// BFS/DFS: Время: O(n + k), Память: O(n)
// Dijkstra: Время: O((n + k) * log(n)), Память: O(n)
// main(): Время: Зависит от алгоритмов, Память: O(n)


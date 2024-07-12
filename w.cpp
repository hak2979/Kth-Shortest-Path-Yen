#include <iostream>
#include <climits>
#include <queue>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include<unordered_set>
#include<omp.h>
#include<random>
#include<ctime>
#include<mpi.h>

#include <chrono>

#include <thread>
using namespace std::chrono;

using namespace std;

class Graph {
public:
    unordered_map<string, unordered_map<string, int>> adjacency_list;
    int index;

    Graph() : index(0) {}

    bool find(string to_want) {
        return adjacency_list.find(to_want) != adjacency_list.end();
    }

    void addedge(string src, string dest, int cost) {
        if (!find(src)) {
            adjacency_list[src] = unordered_map<string, int>();
        }
        if (!find(dest)) {
            adjacency_list[dest] = unordered_map<string, int>();
        }
        adjacency_list[src][dest] = cost;
    }

    bool is_visited(const string& node, const string& to_want) {
        return adjacency_list[node].find(to_want) != adjacency_list[node].end();
    }

    pair<vector<string>, int> dijkstra(const unordered_map<string, unordered_map<string, int>>& adj_list, string src, string dest) {
        unordered_map<string, int> dist;
        unordered_map<string, string> path;
        unordered_set<string> visited; // New: to track visited nodes
        for (const auto& pair : adj_list) {
            dist[pair.first] = INT_MAX;
        }
        dist[src] = 0;
        priority_queue<pair<int, string>, vector<pair<int, string>>, greater<pair<int, string>>> queue;
        queue.push({ 0, src });
        while (!queue.empty()) {
            int u_d = queue.top().first;
            string u = queue.top().second;
            queue.pop();
            if (visited.count(u)) continue;
            visited.insert(u);
            if (u != dest) {
                for (const auto& neighbor : adj_list.at(u)) {
                    int v_dist = u_d + neighbor.second;
                    if (v_dist < dist[neighbor.first]) {
                        dist[neighbor.first] = v_dist;
                        path[neighbor.first] = u;
                        queue.push({ v_dist, neighbor.first });
                    }
                }
            }
            else {
                string x = dest;
                vector<string> new_path;
                while (x != "") {
                    new_path.push_back(x);
                    x = path[x];
                }
                reverse(new_path.begin(), new_path.end());
                return { new_path, u_d };
            }
        }
        return { {}, INT_MAX };
    }

    int getCost(string src, string dest) {
        if (is_visited(src, dest)) {
            return adjacency_list[src][dest];
        }
        return INT_MAX;
    }

    void yenKShortestPaths(string src, string dest, int k) {
        vector<pair<vector<string>, int>> paths;
        unordered_map<string, unordered_map<string, int>> adjacency_list_copys = adjacency_list;
        pair<vector<string>, int> p = dijkstra(adjacency_list_copys,src, dest);
        paths.push_back(p);

        int x = 0;
        bool spur_found = true;
        

        int index=0;
                while (spur_found && k>1) {
                    spur_found = false;
                    if(index==0){
                        #pragma omp parallel num_threads(10) shared(spur_found, paths,index)
                        {
                            #pragma omp for schedule(static,(1))
                                for (int i = 0; i < paths[index].first.size() - 1; i++) {
                                string cur = paths[index].first[i];
                                vector<pair<string, int>> red;
                                unordered_map<string, unordered_map<string, int>> adjacency_list_copy = adjacency_list;
                                // Remove all visited edges
                                for (int j = 0; j < paths.size(); j++) {
                                    const auto& pa = paths[j];
                                    for (int l = 0; l < pa.first.size() - 1; l++) {
                                        if (pa.first[l] == cur) {
                                            int cost = adjacency_list_copy[cur][pa.first[l + 1]];
                                            adjacency_list_copy[cur].erase(pa.first[l + 1]);
                                            red.push_back({ pa.first[l + 1], cost });
                                            adjacency_list_copy[pa.first[l + 1]].erase(cur);
                                        }
                                    }
                                }

                                pair<vector<string>, int> new_pa = dijkstra(adjacency_list_copy,cur, dest);
                                if (!new_pa.first.empty()) {
                                    // Modify existing path by removing visited edges and adding the spur path
                                    vector<string> updated_path;
                                    for (int j = 0; j <= i; j++) {
                                        updated_path.push_back(paths[index].first[j]);
                                    }
                                    for (int j = 1; j < new_pa.first.size(); j++) {
                                        updated_path.push_back(new_pa.first[j]);
                                    }

                                    int combined_cost = 0;
                                    int cost = 0;
                                    for (int i = 0; i < updated_path.size() - 1; i++) {
                                        const string& src = updated_path[i];
                                        const string& dest = updated_path[i + 1];
                                        // Check if the edge exists in the adjacency list
                                        if (adjacency_list.count(src) && adjacency_list[src].count(dest)) {
                                            cost += adjacency_list_copy[src][dest];
                                        }
                                        else {
                                        
                                        }
                                    }
                                    combined_cost=cost;

                                    // Check uniqueness and add new path with the updated cost
                                    #pragma omp critical
                                    {
                                        bool is_unique = true;
                                        #pragma omp parallel for shared(is_unique) num_threads(10)
                                        for (int j = 0; j < paths.size(); ++j) {
                                            if (is_unique) {
                                                const auto& pa = paths[j];
                                                if (pa.first == updated_path) {
                                                    is_unique = false;
                                                }
                                            }
                                        }
                                        if (is_unique) {
                                            paths.push_back({ updated_path, combined_cost });
                                            spur_found = true; // Set the flag to true if a spur node was found
                                        }
                                    }
                                }

                            }
                        }
                    }
                    index++;

                }
            sort(paths.begin(), paths.end(), [&](const auto& path1, const auto& path2) {
                        return comparePaths(path1, path2);
                        });

                    //Remove loops and excessive node revisits
                    while (true) {
                        bool removed = false;
                        for (size_t i = 1; i < paths.size(); ++i) {
                            auto& path = paths[i].first;
                            auto& prev_path = paths[i - 1].first;
                            int j = 0;
                            while (j < path.size() && j < prev_path.size() && path[j] == prev_path[j]) {
                                ++j;
                            }
                            if (j < path.size() && j < prev_path.size() && path[j] == prev_path[j]) {
                                paths.erase(paths.begin() + i);
                                removed = true;
                                break;
                            }
                        }
                        if (!removed || paths.size() <= k) break;
                    }

        // Print K shortest paths
        // int num_printed_paths = min(k, static_cast<int>(paths.size()));
        // for (int i = 0; i < num_printed_paths; i++) {
        //     const auto& path = paths[i];
        //     for (const auto& node : path.first) {
        //         cout << node << "->";
        //     }
        //     cout << " (cost: " << path.second << ")" << endl;
        // }

        if(k<paths.size()){
            const auto& path= paths[k];
            for (const auto& node : path.first) {
                cout << node << "->";
            }
            cout << " (cost: " << path.second << ")" << endl;
        }
        else{
            const auto& path= paths[paths.size()-1];
            for (const auto& node : path.first) {
                cout << node << "->";
            }
            cout << " (cost: " << path.second << ")" << endl;
        }
    }

    void yenKShortestPaths_p(string src, string dest, int k) {
        vector<pair<vector<string>, int>> paths;
        pair<vector<string>, int> p = dijkstra(adjacency_list,src, dest);
        if (p.first.empty()) {
            cout << "No path found between " << src << " and " << dest << endl;
            return;
        }
        paths.push_back(p);

        int x = 0;
        bool spur_found = true;
        

        int index=0;
                while (spur_found && k>1) {
                    spur_found = false;
                    if(index==0){
                                for (int i = 0; i < paths[index].first.size() - 1; i++) {
                                string cur = paths[index].first[i];
                                vector<pair<string, int>> red;
                                // Remove all visited edges
                                for (int j = 0; j < paths.size(); j++) {
                                    const auto& pa = paths[j];
                                    for (int l = 0; l < pa.first.size() - 1; l++) {
                                        if (pa.first[l] == cur) {
                                            int cost = adjacency_list[cur][pa.first[l + 1]];
                                            adjacency_list[cur].erase(pa.first[l + 1]);
                                            red.push_back({ pa.first[l + 1], cost });
                                            adjacency_list[pa.first[l + 1]].erase(cur);
                                        }
                                    }
                                }

                                pair<vector<string>, int> new_pa = dijkstra(adjacency_list,cur, dest);
                                if (!new_pa.first.empty()) {
                                    // Modify existing path by removing visited edges and adding the spur path
                                    vector<string> updated_path;
                                    for (int j = 0; j <= i; j++) {
                                        updated_path.push_back(paths[index].first[j]);
                                    }
                                    for (int j = 1; j < new_pa.first.size(); j++) {
                                        updated_path.push_back(new_pa.first[j]);
                                    }

                                    int combined_cost = 0;
                                    int cost = 0;
                                    for (int i = 0; i < updated_path.size() - 1; i++) {
                                        const string& src = updated_path[i];
                                        const string& dest = updated_path[i + 1];
                                        // Check if the edge exists in the adjacency list
                                        if (adjacency_list.count(src) && adjacency_list[src].count(dest)) {
                                            cost += adjacency_list[src][dest];
                                        }
                                        else {
                                        
                                        }
                                    }
                                    combined_cost=cost;

                                    // Check uniqueness and add new path with the updated cost
                                 
                                        bool is_unique = true;
                                        for (int j = 0; j < paths.size(); ++j) {
                                            if (is_unique) {
                                                const auto& pa = paths[j];
                                                if (pa.first == updated_path) {
                                                    is_unique = false;
                                                }
                                            }
                                        }
                                        if (is_unique) {
                                            paths.push_back({ updated_path, combined_cost });
                                            spur_found = true; // Set the flag to true if a spur node was found
                                        }
                                    
                                }

                            }
                        
                    }
                    index++;

                }
            sort(paths.begin(), paths.end(), [&](const auto& path1, const auto& path2) {
                        return comparePaths(path1, path2);
                        });

                    //Remove loops and excessive node revisits
                    while (true) {
                        bool removed = false;
                        for (size_t i = 1; i < paths.size(); ++i) {
                            auto& path = paths[i].first;
                            auto& prev_path = paths[i - 1].first;
                            int j = 0;
                            while (j < path.size() && j < prev_path.size() && path[j] == prev_path[j]) {
                                ++j;
                            }
                            if (j < path.size() && j < prev_path.size() && path[j] == prev_path[j]) {
                                paths.erase(paths.begin() + i);
                                removed = true;
                                break;
                            }
                        }
                        if (!removed || paths.size() <= k) break;
                    }

        // Print K shortest paths
        // int num_printed_paths = min(k, static_cast<int>(paths.size()));
        // for (int i = 0; i < num_printed_paths; i++) {
        //     const auto& path = paths[i];
        //     for (const auto& node : path.first) {
        //         cout << node << "->";
        //     }
        //     cout << " (cost: " << path.second << ")" << endl;
        // }

        if(k<paths.size()){
            const auto& path= paths[k];
            for (const auto& node : path.first) {
                cout << node << "->";
            }
            cout << " (cost: " << path.second << ")" << endl;
        }
        else{
            const auto& path= paths[paths.size()-1];
            for (const auto& node : path.first) {
                cout << node << "->";
            }
            cout << " (cost: " << path.second << ")" << endl;
        }
    }

    bool comparePaths(const pair<vector<string>, int>& path1, const pair<vector<string>, int>& path2) {
        return path1.second < path2.second;
    }

    void print() {
        for (const auto& pair : adjacency_list) {

                cout << pair.first << " -> ";
                for (const auto& neighbor : pair.second) {
                    cout << neighbor.first << "(" << neighbor.second << ")" << " -> ";
                }
                cout << "nullptr" << endl;
        }
    }
};

int main(int argc, char* argv[]) {

    int rank,size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        Graph g;

        /*g.addedge("C", "D", 3);
        g.addedge("C", "E", 2);
        g.addedge("D", "F", 4);
        g.addedge("E", "D", 1);
        g.addedge("E", "F", 2);
        g.addedge("E", "G", 3);
        g.addedge("F", "G", 2);
        g.addedge("F", "H", 1);
        g.addedge("G", "H", 2);*/
        ifstream file("classic-who.csv");
        if (!file.is_open()) {
            cerr << "Failed to open the file." << endl;
            return 1;
        }

        string line;
        getline(file, line);
        while (getline(file, line)) {
            stringstream ss(line);
            string source, dest;
            string weight;
            string gw;
            // Extracting data from CSV
            getline(ss, source, ',');
            getline(ss, dest, ',');
            getline(ss, weight, ',');
            getline(ss, gw, '\n');
            // Adding undirected edge to the graph
            g.addedge(source, dest, stoi(weight));
        }

        file.close();


        // ifstream file("Email-EuAll.txt");
        // if (!file.is_open()) {
        //    cerr << "Failed to open the file." << endl;
        //    return 0;
        // }

        // string line;

        // while (getline(file, line)) {
        //    stringstream ss(line);
        //    string fromNode, toNode;
        //    // Extracting data from the line
        //    ss >> fromNode >> toNode;
        //    // Add directed edge to the graph
        //    g.addedge(fromNode, toNode,1);
        // }

        // file.close();

        srand(time(NULL));
        int *s =new int[11];
        int *d=new int[11];
        for(int i=0;i<11;i++){
            s[i] = rand() % g.adjacency_list.size() - 1;
            d[i] = s[i];
            while (d[i] == s[i]) {
                d[i] = rand() % g.adjacency_list.size() - 1;
            }
        }
        auto start = high_resolution_clock::now();
        for(int i=0;i<size;i++){
            MPI_Send(s, 11, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(d, 11, MPI_INT, i, 0, MPI_COMM_WORLD);

        }

        for(int i=0;i<size;i++){
            char rdata[]="";
	        MPI_Recv(rdata, 100, MPI_CHAR, i, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        std::this_thread::sleep_for(std::chrono::seconds(1));
        cout << "Execution time: " << duration.count() << " milliseconds" << endl;

        // start = high_resolution_clock::now();
        // for (int i = 0; i < 10; i++) {
        //         string s_S = "";
        //         string d_S = "";
        //         int ind = 0;
        //         for (const auto& pair : g.adjacency_list) {

        //             if (ind == s[i]) {
        //                 s_S = pair.first;
        //             }

        //             ind++;
        //         }
        //         ind = 0;
        //         for (const auto& pair : g.adjacency_list) {

        //             if (ind == d[i]) {
        //                 d_S = pair.first;
        //             }

        //             ind++;
        //         }
        //         g.yenKShortestPaths_p(s_S, d_S, 10);
        //     }
        //     stop = high_resolution_clock::now();
        //     duration = duration_cast<milliseconds>(stop - start);
        //     std::this_thread::sleep_for(std::chrono::seconds(1));
        //     cout << "Execution time: " << duration.count() << " milliseconds" << endl;

    }
    else {
        int* recv_s = new int[11];
        int* recv_d = new int[11];
        MPI_Recv(recv_s, 11, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(recv_d, 11, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

        Graph g;

        /*g.addedge("C", "D", 3);
        g.addedge("C", "E", 2);
        g.addedge("D", "F", 4);
        g.addedge("E", "D", 1);
        g.addedge("E", "F", 2);
        g.addedge("E", "G", 3);
        g.addedge("F", "G", 2);
        g.addedge("F", "H", 1);
        g.addedge("G", "H", 2);*/
        ifstream file("classic-who.csv");
        if (!file.is_open()) {
            cerr << "Failed to open the file." << endl;
            return 1;
        }

        string line;
        getline(file, line);
        while (getline(file, line)) {
            stringstream ss(line);
            string source, dest;
            string weight;
            string gw;
            // Extracting data from CSV
            getline(ss, source, ',');
            getline(ss, dest, ',');
            getline(ss, weight, ',');
            getline(ss, gw, '\n');
            // Adding undirected edge to the graph
            g.addedge(source, dest, stoi(weight));
        }

        file.close();


        // ifstream file("Email-EuAll.txt");
        // if (!file.is_open()) {
        //    cerr << "Failed to open the file." << endl;
        //    return 0;
        // }

        // string line;

        // while (getline(file, line)) {
        //    stringstream ss(line);
        //    string fromNode, toNode;
        //    // Extracting data from the line
        //    ss >> fromNode >> toNode;
        //    // Add directed edge to the graph
        //    g.addedge(fromNode, toNode,1);
        // }

        // file.close();

        for (int i = 1; i <= 11; i++) {
            if (i % size == rank) {
                string s_S = "";
                string d_S = "";
                int ind = 0;
                for (const auto& pair : g.adjacency_list) {

                    if (ind == recv_s[i-1]) {
                        s_S = pair.first;
                    }

                    ind++;
                }
                ind = 0;
                for (const auto& pair : g.adjacency_list) {

                    if (ind == recv_d[i-1]) {
                        d_S = pair.first;
                    }

                    ind++;
                }
                g.yenKShortestPaths(s_S, d_S, 310);
                char sdata[] = "Hello PDC";
                MPI_Send(sdata, 9, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
            }
        }
        //MPI_Send
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
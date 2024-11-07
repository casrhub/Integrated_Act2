#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>
#include <cmath>
#include <map>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

/* Team members: 
César Alán Silva Ramos
José María Soto Valenzuela
German Avelino Del Rio Guzman 

*/

/*
 * Steps to compile and run this program on a Linux system:
 *
 * 1. Navigate to the project directory (where CMakeLists.txt is located):
 *    cd project_directory
 *
 * 2. Create and enter the build directory:
 *    mkdir build
 *    cd build
 *
 * 3. Run CMake to configure the project:
 *    cmake ..
 *
 * 4. Compile the project:
 *    make
 *
 * 5. Run the program with an input file:
 *    ./NeighborhoodFiberOptics < ../input.txt
 */




using namespace std;

// Part of CGAL
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef K::Point_2 Point;
typedef K::Vector_2 Vector;
typedef Delaunay::Face_handle Face_handle;
typedef Delaunay::Vertex_handle Vertex_handle;


// Structure to represent an edge in the graph
struct Edge
{
    int from;   // Starting node of the edge
    int to;     // Ending node of the edge
    int weight; // Weight of the edge (e.g., distance)
};


// Union-Find Disjoint Set class for Kruskal's Algorithm
class UnionFind
{
private:
    vector<int> parent; // Parent of each node in the disjoint set
    vector<int> rank;   // Rank (approximate depth) of each tree in the forest

public:
    // Constructor to initialize the Union-Find structure
    UnionFind(int n)
    {
        parent.resize(n);
        rank.resize(n, 0);
        // Initially, each node is its own parent (makeset operation)
        for (int i = 0; i < n; i++)
            parent[i] = i;
    }

    // Find operation with path compression
    int find(int x)
    {
        // If x is not its own parent, recursively find the root and compress the path
        if (parent[x] != x)
            parent[x] = find(parent[x]); // Path compression
        return parent[x];
    }

    // Union operation by rank
    void unite(int x, int y)
    {
        int rootX = find(x); // Find root of x
        int rootY = find(y); // Find root of y

        if (rootX == rootY)
            return; // x and y are already in the same set

        // Attach the smaller rank tree under the root of the higher rank tree
        if (rank[rootX] < rank[rootY])
            parent[rootX] = rootY; // RootY becomes parent of RootX
        else if (rank[rootX] > rank[rootY])
            parent[rootY] = rootX; // RootX becomes parent of RootY
        else
        {
            parent[rootY] = rootX; // RootX becomes parent of RootY
            rank[rootX]++;         // Increase rank of RootX
        }
    }
};

// Part 1: Kruskal's Algorithm for Minimum Spanning Tree
vector<Edge> kruskalMST(int numNodes, vector<vector<int>> &graph)
{
    vector<Edge> edges;

    // Build the list of edges from the adjacency matrix
    for (int i = 0; i < numNodes; i++)
    {
        for (int j = i + 1; j < numNodes; j++)
        {
            if (graph[i][j] > 0)
            {
                // Add edge to the list if there is a connection
                edges.push_back({i, j, graph[i][j]});
            }
        }
    }

    // Sort edges by weight in ascending order
    sort(edges.begin(), edges.end(), [](Edge a, Edge b)
         { return a.weight < b.weight; });

    UnionFind uf(numNodes); // Initialize Union-Find structure
    vector<Edge> mst;       // Vector to store edges of the MST

    // Iterate over sorted edges and build MST
    for (auto &edge : edges)
    {
        // If adding this edge doesn't form a cycle
        if (uf.find(edge.from) != uf.find(edge.to))
        {
            uf.unite(edge.from, edge.to); // Union the sets
            mst.push_back(edge);          // Add edge to MST
        }
    }
    return mst;
}

// Part 2: Traveling Salesman Problem using Held-Karp Algorithm (Dynamic Programming)
int tspHeldKarp(int numNodes, vector<vector<int>> &graph, vector<int> &optimalPath)
{
    const int INF = INT_MAX;
    int N = numNodes;
    int VISITED_ALL = (1 << N) - 1; // All nodes have been visited when the bitmask is all ones

    // dp[mask][i]: minimum cost to reach node i with visited nodes represented by mask
    vector<vector<int>> dp(1 << N, vector<int>(N, INF));

    // Base case: starting at node 0, cost is 0
    dp[1 << 0][0] = 0;

    // Iterate over all subsets of nodes (represented by masks)
    for (int mask = 0; mask < (1 << N); mask++)
    {
        for (int u = 0; u < N; u++)
        {
            // If u is included in the current subset represented by mask
            if (mask & (1 << u))
            {
                // Try to find the minimum cost to reach each node v from u
                for (int v = 0; v < N; v++)
                {
                    // Skip if v is already in the mask or u and v are the same
                    if ((mask & (1 << v)) || u == v)
                        continue;

                    // Update dp[mask | (1 << v)][v] if a better cost is found
                    if (dp[mask][u] != INF && dp[mask][u] + graph[u][v] < dp[mask | (1 << v)][v])
                    {
                        dp[mask | (1 << v)][v] = dp[mask][u] + graph[u][v];
                    }
                }
            }
        }
    }

    // Reconstruct the optimal path
    int mask = VISITED_ALL;
    int last = 0;
    int min_cost = INF;

    // Find the ending node that gives minimum cost
    for (int i = 1; i < N; i++)
    {
        if (dp[mask][i] != INF && dp[mask][i] + graph[i][0] < min_cost)
        {
            min_cost = dp[mask][i] + graph[i][0];
            last = i;
        }
    }

    // Backtracking to find the optimal path
    int current_mask = mask;
    int current_node = last;

    vector<int> path;
    path.push_back(last);

    while (current_mask != (1 << 0))
    {
        int prev_mask = current_mask ^ (1 << current_node);
        for (int i = 0; i < N; i++)
        {
            if ((prev_mask & (1 << i)) && dp[prev_mask][i] != INF && dp[prev_mask][i] + graph[i][current_node] == dp[current_mask][current_node])
            {
                path.push_back(i);
                current_node = i;
                current_mask = prev_mask;
                break;
            }
        }
    }

    // Reverse the path and add to optimalPath
    for (int i = path.size() - 1; i >= 0; i--)
    {
        optimalPath.push_back(path[i]);
    }

    optimalPath.push_back(0); // Return to starting node

    return min_cost;
}

// Part 3: Ford-Fulkerson Algorithm for Maximum Flow using DFS

// DFS function to find an augmenting path in the residual graph from 'source' to 'sink'
bool dfs(vector<vector<int>> &residualGraph, int current, int sink, vector<int> &parent)
{
    // If we've reached the sink node, a path has been found
    if (current == sink)
        return true; // Reached sink node

    // Explore each possible next node in the graph
    for (int next = 0; next < residualGraph.size(); ++next)
    {
        // Check if there's residual capacity (i.e., available capacity) and that 'next' is unvisited
        if (residualGraph[current][next] > 0 && parent[next] == -1)
        {
            parent[next] = current; // Set 'current' as the parent of 'next' to reconstruct the path

            // Recursively attempt to find a path to the sink node from 'next'
            if (dfs(residualGraph, next, sink, parent))
                return true; // Path to sink has been found
        }
    }

    // If no path to the sink node was found from the 'current' node
    return false;
}

bool arePointsColinear(const vector<Point>& points) {
    if (points.size() < 3) return true; // Less than 3 points are always colinear
    auto p0 = points[0];
    auto p1 = points[1];
    for (size_t i = 2; i < points.size(); ++i) {
        auto p = points[i];
        // Calculate the area of the triangle formed by p0, p1, and p
        double area = (p0.x() * (p1.y() - p.y()) + p1.x() * (p.y() - p0.y()) + p.x() * (p0.y() - p1.y())) / 2.0;
        if (abs(area) > 1e-6) return false; // Points are not colinear
    }
    return true; // All points are colinear
}

// Main function to calculate the maximum flow from 'source' to 'sink' in the graph
int maxFlow(vector<vector<int>> &graph, int source, int sink)
{
    int totalFlow = 0;                         // Initialize total flow to zero
    vector<vector<int>> residualGraph = graph; // Copy original graph to create the residual graph
    vector<int> parent(graph.size());          // Array to store the augmenting path (for each node, store its parent)

    // Loop until no more augmenting paths are found
    while (true)
    {
        fill(parent.begin(), parent.end(), -1); // Reset the parent array for a fresh search
        parent[source] = -2;                    // Mark the source node as visited

        // Use DFS to find an augmenting path from source to sink
        if (!dfs(residualGraph, source, sink, parent))
            break; // No more augmenting paths exist, so we terminate the loop

        // Find the bottleneck capacity (minimum residual capacity) along the augmenting path found
        int flow = INT_MAX;
        for (int v = sink; v != source; v = parent[v])
        {
            int u = parent[v];                     // Get the previous node in the path
            flow = min(flow, residualGraph[u][v]); // Update flow to the minimum residual capacity along the path
        }

        // Update residual capacities in the graph based on the bottleneck capacity
        for (int v = sink; v != source; v = parent[v])
        {
            int u = parent[v];
            residualGraph[u][v] -= flow; // Subtract bottleneck flow from the forward edge
            residualGraph[v][u] += flow; // Add bottleneck flow to the reverse edge (allowing for possible backflow)
        }

        // Add the bottleneck flow of the current path to the total maximum flow
        totalFlow += flow;
    }

    // Return the computed maximum flow from source to sink
    return totalFlow;
}

// Helper function to convert index to character (e.g., 0 -> 'A')
char indexToChar(int index)
{
    return 'A' + index;
}

int main()
{
    int numNodes;
    cin >> numNodes;

    vector<vector<int>> distanceMatrix(numNodes, vector<int>(numNodes)); // Distance between neighborhoods
    vector<vector<int>> capacityMatrix(numNodes, vector<int>(numNodes)); // Capacity between neighborhoods
    vector<pair<int, int>> coordinates(numNodes);                        // Coordinates of the exchanges

    // Read distance matrix
    for (int i = 0; i < numNodes; i++)
    {
        for (int j = 0; j < numNodes; j++)
        {
            cin >> distanceMatrix[i][j];
        }
    }

    // Read capacity matrix
    for (int i = 0; i < numNodes; i++)
    {
        for (int j = 0; j < numNodes; j++)
        {
            cin >> capacityMatrix[i][j];
        }
    }

    // Read coordinates of exchanges
    for (int i = 0; i < numNodes; i++)
    {
        char ch;
        cin >> ch; // Read '('
        int x;
        cin >> x;
        cin >> ch; // Read ','
        int y;
        cin >> y;
        cin >> ch; // Read ')'
        coordinates[i] = {x, y};
    }

    // Part 1: Optimal wiring (Minimum Spanning Tree)
    vector<Edge> mst = kruskalMST(numNodes, distanceMatrix);
    cout << "1. Way of wiring the neighborhoods with fiber (list of arcs):" << endl;
    for (auto &edge : mst)
    {
        cout << "(" << indexToChar(edge.from) << "," << indexToChar(edge.to) << ")" << endl;
    }

    // Part 2: Mail delivery route (Traveling Salesman Problem using Held-Karp Algorithm)
    vector<int> optimalRoute;
    int minTSPCost = tspHeldKarp(numNodes, distanceMatrix, optimalRoute);

    cout << "2. Route to be followed by the mail delivery personnel:" << endl;
    for (size_t i = 0; i < optimalRoute.size(); i++)
    {
        cout << indexToChar(optimalRoute[i]);
        if (i != optimalRoute.size() - 1)
            cout << " -> ";
    }
    cout << endl;

    // Part 3: Maximum information flow value from the initial node to the final node
    int sourceNode = 0;          // Initial node ('A')
    int sinkNode = numNodes - 1; // Final node (last node)
    int max_flow = maxFlow(capacityMatrix, sourceNode, sinkNode);
    cout << "3. Maximum information flow value from the initial node to the final node: " << max_flow << endl;


    
      // Part 4: Compute and Output Voronoi Polygons for Each Exchange
    
    vector<Point> points;
    for (const auto& coord : coordinates)
        points.emplace_back(coord.first, coord.second);


    if (arePointsColinear(points)) {
    cout << "4. Voronoi cells cannot be computed because all points are colinear." << endl;
} else {
    Delaunay delaunay;
    delaunay.insert(points.begin(), points.end());


    map<Point, vector<pair<double, double>>> voronoi_polygons;
    for (auto vertex = delaunay.finite_vertices_begin(); vertex != delaunay.finite_vertices_end(); ++vertex) {
        vector<pair<double, double>> polygon;
        Delaunay::Face_circulator circulator = delaunay.incident_faces(vertex), done = circulator;
        do {
            Face_handle face = circulator;
            if (!delaunay.is_infinite(face)) {
                Point voronoi_vertex = delaunay.dual(face);
                polygon.emplace_back(voronoi_vertex.x(), voronoi_vertex.y());
            } else {
                int i = face->index(vertex);
                Face_handle neighbor = face->neighbor(i);
                if (!delaunay.is_infinite(neighbor)) {
                    Point p1 = delaunay.dual(neighbor);
                    Vertex_handle v_target = face->vertex((i + 1) % 3);
                    Vertex_handle v_source = face->vertex((i + 2) % 3);
                    Vector vec = v_target->point() - v_source->point();
                    Vector dir = Vector(-vec.y(), vec.x()); // Perpendicular to the edge
dir = dir / std::sqrt(dir.squared_length()); // Normalize the direction vector
                    Point extended_point = delaunay.dual(face) + dir * 1000;  // Extend by arbitrary distance
                    polygon.emplace_back(extended_point.x(), extended_point.y());
                }
            }
            ++circulator;
        } while (circulator != done);
        voronoi_polygons[vertex->point()] = polygon;
    }

    cout << "4. Voronoi cells for exchanges:" << endl;
    for (const auto& [exchange, polygon] : voronoi_polygons) {
        cout << "Exchange at (" << exchange.x() << ", " << exchange.y() << "): ";
        for (const auto& vertex : polygon)
            cout << "(" << vertex.first << ", " << vertex.second << ") ";
        cout << endl;
    }
}


    return 0;
}
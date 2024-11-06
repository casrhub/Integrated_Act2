
# Neighborhood Fiber Optics Project

This project implements algorithms to optimize the planning of a fiber optics network in a neighborhood. The program performs the following tasks:

1. **Minimum Spanning Tree (MST)**: Calculates the optimal way to wire the neighborhoods with fiber, minimizing the total length of the cables using Kruskal's algorithm.

2. **Traveling Salesman Problem (TSP)**: Determines the most efficient route for mail delivery personnel to follow, ensuring the shortest possible route that visits all neighborhoods using the Held-Karp algorithm.

3. **Maximum Flow**: Computes the maximum information flow value from the initial node to the final node using the Ford-Fulkerson algorithm.

4. **Voronoi Diagram**: Generates Voronoi cells for exchanges based on the given coordinates, which can be used for further geographical analysis.

For detailed explanations of the algorithms and their implementations, please refer to the **[Project Report](int_act_2.pdf)** included in the repository.

---

## Steps to Compile and Run the Program on a Linux System

1. **Navigate to the project directory** (where `CMakeLists.txt` is located):
   ```bash
   cd project_directory
   ```

2. **Create and enter the build directory**:
   ```bash
   mkdir build
   cd build
   ```

3. **Run CMake to configure the project**:
   ```bash
   cmake ..
   ```

4. **Compile the project**:
   ```bash
   make
   ```

5. **Run the program with an input file**:
   ```bash
   ./NeighborhoodFiberOptics < ../input.txt
   ```

---

## Input Format

The program requires an input file with the following structure:

1. An integer, `N`, representing the number of neighborhoods.
2. An `N x N` matrix representing the distances between neighborhoods.
3. An `N x N` matrix representing the maximum data flow capacities between neighborhoods.
4. A list of `N` coordinates in the form `(x,y)` representing the geographic location of exchanges.

**Example Input**:
```
4

0 10 15 20
10 0 35 25
15 35 0 30
20 25 30 0

0 16 13 10
16 0 12 14
13 12 0 18
10 14 18 0

(100,200)
(200,100)
(300,200)
(400,300)
```

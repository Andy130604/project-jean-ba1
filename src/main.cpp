#include <iostream>
#include <iomanip>
#include <exception>
#include "pathsearch.hpp"
#include "utils.hpp"
#include "unit_test.hpp"

/*---------------- Helper test functions  -----------------------*/
void test_part1();
void test_part2();
void test_part3();
void test_larger_networks();

/*---------------- Main  -----------------------*/

int main(int argc, char *argv[])
{

    std::cout << "========= TESTING PART 1 ================" << std::endl;
    test_part1(); // UNCOMMENT WHEN READY TO TEST

    std::cout << "========= TESTING PART 2 ================" << std::endl;
    test_part2();

    std::cout << "========= TESTING PART 3 ================" << std::endl;
    test_part3();

    std::cout << "========= TESTING WITH LARGER NETWORKS  ================" << std::endl;
    // test_larger_networks();

    std::cout << "========= UNIT TESTS ================" << std::endl;
    run_unit_tests(1); // UNCOMMENT WHEN READY TO TEST PART 1
    run_unit_tests(2); // UNCOMMENT WHEN READY TO TEST PART 2
    run_unit_tests(3); // UNCOMMENT WHEN READY TO TEST PART 3

    return 0;
}

void test_part1()
{
    std::cout << " ======= Testing find_compoundID ======= " << std::endl;
    Network basic_network(read_network("data/basic.txt"));
    std::cout << find_compoundID(basic_network, "C4") << std::endl; // should display 4

    std::cout << " ======= Testing build_adjacency_graph ======= " << std::endl;

    /* Should display:
     *
    0:C0: {{1,0}}
    1:C1: {{0,0}{2,1}{4,2}}
    2:C2: {{1,1}{3,4}{5,5}}
    3:C3: {{2,4}{5,6}}
    4:C4: {{1,2}{5,3}}
    5:C5: {{2,5}{3,6}{4,3}}
   */
    std::cout << to_string(basic_network, build_adjacency_graph(basic_network));

    std::cout << " ======= Testing find_reactionID ======= " << std::endl;
    AdjacencyGraph graph{{{1, 0}},
                         {{0, 0}, {2, 1}, {4, 2}},
                         {{1, 1}, {3, 4}, {5, 5}},
                         {{2, 4}, {5, 6}},
                         {{1, 2}, {5, 3}},
                         {{2, 5}, {3, 6}, {4, 3}}};
    std::cout << find_reactionID(graph, 3, 5) << std::endl; // should display 6

    std::cout << " ======= Testing bfs ======= " << std::endl;
    /*
     * should diplay :
    start:0
    Node0: {-1 }, distance : 0
    Node1: {0 }, distance : 1
    Node2: {1 }, distance : 2
    Node3: {2 }, distance : 3
    Node4: {1 }, distance : 2
    Node5: {2 4 }, distance : 3
    */

    BFS bfs_result(bfs(graph, 0));
    std::cout << to_string(bfs_result) << std::endl;
}
void test_part2()
{

    AdjacencyGraph graph{{{1, 0}},
                         {{0, 0}, {2, 1}, {4, 2}},
                         {{1, 1}, {3, 4}, {5, 5}},
                         {{2, 4}, {5, 6}},
                         {{1, 2}, {5, 3}},
                         {{2, 5}, {3, 6}, {4, 3}}};

    std::cout << " ======= Testing find_shortest_path ======= " << std::endl;
    /* should display :
     Found path: Reactions: 0 1 5
     or
     Found path: Reactions: 0 2 3
      */

    auto path(find_shortest_path(graph, 0, 5));
    std::cout << "Found path: " << to_string(path) << std::endl;

    std::cout << " ======= Testing find_all_shortest_paths ======= " << std::endl;
    /* should display :
    Number of shortest paths: 2
    Path Number 1 - Reactions: 0 1 5
    Path Number 2 - Reactions: 0 2 3
    */

    Paths paths(find_all_shortest_paths(graph, 0, 5));
    std::cout << "Number of shortest paths: " << paths.size() << std::endl;
    std::cout << to_string(paths) << std::endl;
}
void test_part3()
{
    Network network = read_network("data/basic.txt");
    Concentrations initial = read_initial_concentrations(network, "data/basic_concentrations.txt");

    std::cout << "Testing with network basic.txt " << std::endl;

    std::cout << " ======= Testing test_michaelis_reversible_rate ======= " << std::endl;

    std::cout << michaelis_reversible_rate(network.reactions[0], 0.5, 0.1) << std::endl;  // 1.1348
    std::cout << michaelis_reversible_rate(network.reactions[0], 0.05, 0.8) << std::endl; // -3.97324
    std::cout << " ======= test_compute_ss_concentration ======= " << std::endl;

    std::cout << "concentrations (all compounds) before Euler simulation : " << std::endl;
    /* should  display:
    0 = 0.312
    1 = 0.762
    2 = 0.435
    3 = 0.963
    4 = 0.045
    5 = 0.187
    */
    std::cout << to_string(initial);
    std::cout << "concentrations after Euler simulation on shortest path {0,1,5} : " << std::endl;
    /* should display:
    0 = 0.937932
    1 = 0.403755
    2 = 0.467368
    5 = 0.310334
     */
    Concentrations ss_path1 = compute_ss_concentration(network, {0, 1, 5}, initial, 1e-3);
    std::cout << to_string(ss_path1) << std::endl;
    std::cout << "concentrations (all compounds) must remain unchanged: " << std::endl;
    /* should  display:
    0 = 0.312
    1 = 0.762
    2 = 0.435
    3 = 0.963
    4 = 0.045
    5 = 0.187
    */
    std::cout << to_string(initial);

    std::cout << "concentrations after Euler simulation on shortest path {0,2,3} : " << std::endl;
    /*should print:
    0 = 0.881971
    1 = 0.325478
    4 = 0.354201
    5 = 0.590138
    */
    Concentrations ss_path2 = compute_ss_concentration(network, {0, 2, 3}, initial, 1e-3);
    std::cout << to_string(ss_path2) << std::endl;

    std::cout << " ======= test_compute_path_rate ======= " << std::endl;

    std::cout << "path rate on shortest path {0,1,5} : " << std::endl;
    std::cout << compute_path_rate(network, {0, 1, 5}, ss_path1) << std::endl; // should display : 0.310335
    std::cout << "path rate on shortest path {0,2,3} : " << std::endl;
    std::cout << compute_path_rate(network, {0, 2, 3}, ss_path2) << std::endl; // should display : 0.590141

    std::cout << " ======= test_fastest_path ======= " << std::endl;

    Path fastest_path = find_fastest_path(network, {{0, 1, 5}, {0, 2, 3}}, initial, 1e-2);
    std::cout << to_string(fastest_path) << std::endl; // should display:  0 2 3
}
void test_larger_networks()
{
    /* should display:
    Path Number 1 - Reactions: 38 5 12
    Path Number 2 - Reactions: 18 22 12
    Path Number 3 - Reactions: 38 50 47

    Reactions: 38 5 12
      */
    Network network = read_network("data/C00025-C00148.txt");
    Concentrations initial = read_initial_concentrations(network, "data/C00025-C00148_concentrations.txt");
    Paths shortest_paths(find_all_shortest_paths(build_adjacency_graph(network), 0, 32));
    std::cout << to_string(shortest_paths) << std::endl;
    std::cout << to_string(find_fastest_path(network, shortest_paths, initial, 1e-2)) << std::endl;
}

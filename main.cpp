#include <iostream>
#include "drt.h"

int main() {
    DRT rt_model;

    int v1 = rt_model.add_vertex("v1", 1, 4);
    int v2 = rt_model.add_vertex("v2", 1, 5);
    int v3 = rt_model.add_vertex("v3", 1, 6);

    rt_model.add_edge(v1, v2, 10);
    rt_model.add_edge(v1, v3, 10);

    auto vertices = rt_model.get_vertices();

    for(auto vertex : vertices) {
        std::cout << vertex.name << " : <" << vertex.wcet << ", " << vertex.d << "> " << std::endl;
        for(auto edge : vertex.edges) {
            std::cout << "\t--- " << edge.p  << " ---> "  << std::endl;
        }
    }
    
    return 0;
}

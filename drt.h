#ifndef DRT_HPP
#define DRT_HPP

#include <vector>
#include <string>

struct edge {
    int p;
    std::size_t from;
    std::size_t to;
};

struct vertex {
    std::size_t         idx;
    std::string         name;
    int                 wcet;
    int                 d;
    std::vector<edge>   edges;
};

class DRT {
    private:
        std::vector<vertex> vertices;
    public:
        std::size_t add_vertex(std::string name, int e, int d) {
            vertices.push_back({0, std::move(name), e, d, {}});
            vertices.back().idx = vertices.size() - 1;
            return vertices.back().idx;
        }

        void add_edge(std::size_t from, std::size_t to, int p){
            vertices[from].edges.push_back({p, from, to});
        }
        
        const std::vector<vertex>& get_vertices() const { return vertices; }
};

#endif

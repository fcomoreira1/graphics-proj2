#include "mesh.h"
#include "random.h"
#include <cassert>
#include <cinttypes>
#include <limits>
#include <list>
#include <map>
#include <utility>

void TriangleMesh::transform(double scale_factor, Vector translation,
                             double rotation_angle) {
    // Vectors for rotation matrix (rotation on y-axis)
    Vector row1(cos(rotation_angle), 0, sin(rotation_angle));
    Vector row3(-sin(rotation_angle), 0, cos(rotation_angle));
    for (auto &v : vertices) {
        v = scale_factor * v + translation;
        v = Vector(dot(row1, v), v.data[1], dot(row3, v));
    }
}

bool operator<(const Edge &current, const Edge &other) {
    if (current.vtxi < other.vtxi) {
        return true;
    } else if (current.vtxi > other.vtxi) {
        return false;
    }
    return current.vtxj < other.vtxj;
}

TriangleMesh TriangleMesh::tutte(int N_iter) {
    std::map<Edge, std::vector<int>> edge_to_triangles;
    std::map<int, std::vector<Edge>> vtx_to_triangles;
    for (int i = 0; i < indices.size(); i++) {
        int vtxi = indices[i].vtxi;
        int vtxj = indices[i].vtxj;
        int vtxk = indices[i].vtxk;
        edge_to_triangles[Edge(vtxi, vtxj)].push_back(i);
        edge_to_triangles[Edge(vtxi, vtxk)].push_back(i);
        edge_to_triangles[Edge(vtxj, vtxk)].push_back(i);
        vtx_to_triangles[vtxi].push_back(Edge(vtxj, vtxk));
        vtx_to_triangles[vtxj].push_back(Edge(vtxi, vtxk));
        vtx_to_triangles[vtxk].push_back(Edge(vtxi, vtxj));
    }
    // get boundary vertices
    std::vector<Edge> boundary_edges;
    std::map<int, std::vector<int>> unordered_boundary_edges;
    for (auto it : edge_to_triangles) {
        if (it.second.size() == 1) {
            unordered_boundary_edges[it.first.vtxi].push_back(it.first.vtxj);
            unordered_boundary_edges[it.first.vtxj].push_back(it.first.vtxi);
        }
    }
    int first = unordered_boundary_edges.begin()->first;
    boundary_edges.push_back(Edge(first, unordered_boundary_edges[first][0]));
    int next = unordered_boundary_edges[first][0];
    int prev = first;
    while (next != first) {
        auto cur = unordered_boundary_edges[next];
        if (cur[0] != prev) {
            prev = next;
            next = cur[0];
        } else {
            prev = next;
            next = cur[1];
        }
        boundary_edges.push_back(Edge(prev, next));
    }

    std::vector<bool> is_on_boundary(vertices.size(), false);
    TriangleMesh result = *this;
    int cur_vtx = boundary_edges[0].vtxi;
    for (int i = 0; i < boundary_edges.size(); i++) {
        Edge current_edge = boundary_edges[i];
        cur_vtx = (current_edge.vtxj != cur_vtx) ? current_edge.vtxj
                                                 : current_edge.vtxi;
        double theta = i / (double)boundary_edges.size() * 2.0 * M_PI;
        result.vertices[cur_vtx] = Vector(cos(theta), sin(theta), 0);
        is_on_boundary[cur_vtx] = true;
    }
    std::vector<Vector> new_positions(vertices.size());
    for (int i = 0; i < new_positions.size(); i++) {
        result.vertices[i] =
            is_on_boundary[i]
                ? result.vertices[i]
                : Vector(uniform_distribution(), uniform_distribution(), 0);
    }
    for (int iter = 0; iter < N_iter; iter++) {
        for (int i = 0; i < vertices.size(); i++) {
            if (is_on_boundary[i]) {
                new_positions[i] = result.vertices[i];
                continue;
            }
            Vector average_vtx(0, 0, 0);
            for (int e = 0; e < vtx_to_triangles[i].size(); e++) {
                const Edge &opposite_e = vtx_to_triangles[i][e];
                average_vtx = average_vtx + result.vertices[opposite_e.vtxi] / 2.0;
                average_vtx = average_vtx + result.vertices[opposite_e.vtxj] / 2.0;
            }
            new_positions[i] = average_vtx / vtx_to_triangles[i].size();
        }
        result.vertices = new_positions;
    }
    // for (int i = 0; i < vertices.size(); i++) {
    //     std::cout << i << ": " << result.vertices[i][0] << ", "
    //               << result.vertices[i][1] << ", " << result.vertices[i][2]
    //               << std::endl;
    // }
    return result;
}

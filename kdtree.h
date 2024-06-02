#pragma once
#include "nanoflann.hpp"
#include "vector.h"
#include <iostream>
#include <vector>

const int init_k = 100;

struct PointCloud {

    std::vector<Vector> pts;

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0)
            return pts[idx].data[0];
        else if (dim == 1)
            return pts[idx].data[1];
        else
            return pts[idx].data[2];
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX> bool kdtree_get_bbox(BBOX & /* bb */) const {
        return false;
    }
};

class KDTree {
    using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud,
        3 /* dim */
        >;

  public:
    KDTree(const std::vector<Vector> &points) {
        const PointCloud cloud{points};
        index = new my_kd_tree_t(3 /*dim*/, cloud, {10 /* max leaf */});
    }
    ~KDTree() { delete index; }
    void findNeighbors(const double query_pt[3], const int k,
                       std::vector<u_long> &ret_indexes) {
        ret_indexes.resize(k);
        std::vector<double> out_dist_sqr(k);
        nanoflann::KNNResultSet<double> resultSet(k);
        resultSet.init(&ret_indexes[0], &out_dist_sqr[0]);
        index->findNeighbors(resultSet, query_pt);
        // u_long num_results =
        //     index->knnSearch(query_pt, k, &ret_indexes[0], &out_dist_sqr[0]);
        ret_indexes.resize(resultSet.size());
        std::cout << "Dumping the distances for debugging" << std::endl;
        for (int i = 0; i < 4; i++) {
            std::cout << out_dist_sqr[i] << ", ";
        }
        std::cout << std::endl;
    }

  private:
    my_kd_tree_t *index;
};

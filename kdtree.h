#pragma once
#include "nanoflann.hpp"
#include "vector.h"
#include <iostream>
#include <vector>

const int init_k = 10;

class PointCloud {
  public:
    std::vector<Vector> *pts_ptr;

    PointCloud(std::vector<Vector> *points) : pts_ptr(points) {}

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return (*pts_ptr).size(); }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0)
            return (*pts_ptr)[idx].data[0];
        else if (dim == 1)
            return (*pts_ptr)[idx].data[1];
        else
            return (*pts_ptr)[idx].data[2];
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
    KDTree(std::vector<Vector> &points) {
        const PointCloud cloud{&points};
        index = new my_kd_tree_t(3 /*dim*/, cloud, {10 /* max leaf */});
    }
    ~KDTree() { delete index; }
    std::pair<double, int> findNeighbors(const double query_pt[3], const int k,
                       std::vector<u_long> &ret_indexes) {
        // ret_indexes.resize(k);
        std::vector<double> out_dist_sqr(k);
        nanoflann::KNNResultSet<double> result_set(k);
        result_set.init(&ret_indexes[0], &out_dist_sqr[0]);
        index->findNeighbors(result_set, query_pt);
        // u_long num_results =
        //     index->knnSearch(query_pt, k, &ret_indexes[0], &out_dist_sqr[0]);
        // ret_indexes.resize(result_set.size());
        return {out_dist_sqr[result_set.size() - 1], result_set.size()};
    }

  private:
    my_kd_tree_t *index;
};

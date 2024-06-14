#include "random.h"
#include "vector.h"
#include <algorithm>
#include <vector>

typedef std::vector<Vector> Image;

void color_match(Image &I, Image &M, int num_iter = 100) {
    std::vector<std::pair<double, int>> projI(I.size());
    std::vector<std::pair<double, int>> projM(M.size());
    while (num_iter--) {
        auto v = random_direction();
        for (int i = 0; I.size(); i++) {
            projI[i] = std::make_pair(dot(v, I[i]), i);
            projM[i] = std::make_pair(dot(v, M[i]), i);
        }
        std::sort(projI.begin(), projI.end());
        std::sort(projM.begin(), projM.end());
        for (int i = 0; i < I.size(); i++) {
            I[projI[i].second] += (projM[i].first - projI[i].first) * v;
        }
    }
}

int main() { return 0; }

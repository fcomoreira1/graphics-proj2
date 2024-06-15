#include "random.h"
#include "vector.h"
#include <algorithm>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include "../stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"

typedef std::vector<Vector> Image;

void color_match(Image &I, Image &M, int num_iter = 100) {
    printf("Color matching\n");
    std::vector<std::pair<double, int>> projI(I.size());
    std::vector<std::pair<double, int>> projM(M.size());
    while (num_iter--) {
        printf("Iteration %d\n", num_iter);
        auto v = random_direction();
        printf("Random direction: %f, %f, %f\n", v[0], v[1], v[2]);
        for (int i = 0; i < I.size(); i++) {
            projI[i] = std::make_pair(dot(I[i], v), i);
            projM[i] = std::make_pair(dot(M[i], v), i);
        }
        std::sort(projI.begin(), projI.end());
        std::sort(projM.begin(), projM.end());
        for (int i = 0; i < I.size(); i++) {
            I[projI[i].second] += (projM[i].first - projI[i].first) * v;
        }
    }
    printf("Color matching done\n");
}

Image process_image(const char *filename, int &w, int &h) {
    int c;
    unsigned char *data = stbi_load(filename, &w, &h, &c, 3);
    if (data == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }
    if (c != 3) {
        printf("Warning: only RGB images are supported, not %d channels\n", c);
    }
    Image img(w * h);
    for (int i = 0; i < w * h; i++) {
        img[i] = Vector(data[i * 3] / 255.0, data[i * 3 + 1] / 255.0,
                        data[i * 3 + 2] / 255.0);
    }
    stbi_image_free(data);
    return img;
}

std::vector<unsigned char> export_image(const Image &img) {
    std::vector<unsigned char> data(img.size() * 3);
    for (int i = 0; i < img.size(); i++) {
        for (int j = 0; j < 3; j++) {
            data[i * 3 + j] =
                (unsigned char)(255 * std::max(0.0, std::min(1.0, img[i][j])));
        }
    }
    return data;
}

int main() {
    int w_input, h_input, w_model, h_model;
    Image I = process_image("davy.png", w_input, h_input);
    Image M = process_image("brad.png", w_model, h_model);
    if (w_input != w_model || h_input != h_model || I.size() != M.size()) {
        printf("Error: input and model images should have the same size\n");
        exit(1);
    }
    color_match(I, M);
    // CHECK THIS
    stbi_write_png("out.png", w_input, h_input, 3, &export_image(I)[0], w_input * 3);
    return 0;
}

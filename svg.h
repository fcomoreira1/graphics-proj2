#include "polygon.h"
#define _CRT_SECURE_NO_WARNINGS 1

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
// if the Polygon class name conflicts with a class in wingdi.h on Windows, use
// a namespace or change the name

// saves a static svg file. The polygon vertices are supposed to be in the range
// [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename,
              std::string fillcol = "none") {
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" "
               "height = \"1000\">\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i][j][0] * 1000),
                    (1000 - polygons[i][j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

void save_svg(const std::vector<Polygon> &polygons,
              const std::vector<Vector> &vertices, std::string filename,
              std::string fillcol = "none") {
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" "
               "height = \"1000\">\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i][j][0] * 1000),
                    (1000 - polygons[i][j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    for (int i = 0; i < vertices.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<circle cx= \"");
        fprintf(f, "%3.3f\" cy=\"%3.3f", vertices[i][0] * 1000,
                1000 - vertices[i][1] * 1000);
        fprintf(f, "\" r=\"3\"\nfill = \"red\" stroke = \"red\"/>\n",
                fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

// Adds one frame of an animated svg file. frameid is the frame number (between
// 0 and nbframes-1). polygons is a list of polygons, describing the current
// frame. The polygon vertices are supposed to be in the range [0..1], and a
// canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons,
                       std::string filename, int frameid, int nbframes) {
    FILE *f;
    if (frameid == 0) {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = "
                   "\"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    } else {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i][j][0] * 1000),
                    (1000 - polygons[i][j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"blue\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "    id = \"frame%u\"\n", frameid);
    fprintf(f, "    attributeName = \"display\"\n");
    fprintf(f, "    values = \"");
    for (int j = 0; j < nbframes; j++) {
        if (frameid == j) {
            fprintf(f, "inline");
        } else {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n    keyTimes = \"");
    for (int j = 0; j < nbframes; j++) {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n   dur = \"5s\"\n");
    fprintf(f, "    begin = \"0s\"\n");
    fprintf(f, "    repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1) {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}

static int sgn(double x) { return x < 0 ? -1 : x > 0; }

void save_frame(const std::vector<Polygon> &cells, std::string filename,
                int frameid = 0) {
    int W = 1000, H = 1000;
    std::cerr << "Saving frame " << frameid << std::endl;
    std::vector<unsigned char> image(W * H * 3, 255);
    std::cerr << "Created image" << std::endl;
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {
        std::cerr << "Processing cell " << i << std::endl;

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].size(); j++) {
            bminx = std::min(bminx, cells[i][j][0]);
            bminy = std::min(bminy, cells[i][j][1]);
            bmaxx = std::max(bmaxx, cells[i][j][0]);
            bmaxy = std::max(bmaxy, cells[i][j][1]);
        }
        bminx = std::min(W - 1., std::max(0., W * bminx));
        bminy = std::min(H - 1., std::max(0., H * bminy));
        bmaxx = std::max(W - 1., std::max(0., W * bmaxx));
        bmaxy = std::max(H - 1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].size(); j++) {
                    double x0 = cells[i][j][0] * W;
                    double y0 = cells[i][j][1] * H;
                    double x1 = cells[i][(j + 1) % cells[i].size()][0] * W;
                    double y1 = cells[i][(j + 1) % cells[i].size()][1] * H;
                    double det = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0);
                    int sign = sgn(det);
                    if (prevSign == 0)
                        prevSign = sign;
                    else if (sign == 0)
                        sign = prevSign;
                    else if (sign != prevSign) {
                        isInside = false;
                        break;
                    }
                    prevSign = sign;
                    double edgeLen = std::sqrt((x1 - x0) * (x1 - x0) +
                                               (y1 - y0) * (y1 - y0));
                    double distEdge = std::abs(det) / edgeLen;
                    double dotp = (x - x0) * (x1 - x0) + (y - y0) * (y1 - y0);
                    if (dotp < 0 || dotp > edgeLen * edgeLen)
                        distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    if (i < cells.size() - 1) {
                        // the N first particles may represent
                        // fluid, displayed in blue
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 255;
                    } else if (mindistEdge <= 2) {
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 0;
                    }
                }
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    std::cerr << "Finished Processing" << std::endl;
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

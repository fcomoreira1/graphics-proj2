#include "mesh.h"

int main() {
    TriangleMesh init;
    std::cout << "Reading Obj" << std::endl;
    init.readOBJ("goethe.obj");
    std::cout << "Running Tutte" << std::endl;
    TriangleMesh res = init.tutte(1000);
    std::cout << "Writting obj" << std::endl;
    res.writeOBJ("out.obj");
    return 0;
}

#include <stdlib.h>
#include <stdio.h>

#include "Model.h"
extern int g_sharp;
int main(int argc, char** argv)
{
    Model obj;
    int resolution = 20000;
    float min_hole_size = 0.02;
    if (argc < 3)
    {
        cout << "./manifold input.obj output.obj [resolution=20000] [-s]\n";
        return 0;
    }
    obj.Load(argv[1]);
    if (obj.vertices.empty() || obj.faces.empty()) {
        std::cout << "Empty mesh, nothing to do." << std::endl;
        exit(0);
    }

    if (argc > 3)
    {
        if (strcmp(argv[3], "-s") == 0) {
            g_sharp = 1;
        } else {
            sscanf(argv[3], "%d", &resolution);
            if (argc > 4) {
                if ( strcmp(argv[4], "-s") == 0) {
                    g_sharp = 1;
                } else {
                    sscanf(argv[4], "%f", &min_hole_size);
                    if (argc > 5 && strcmp(argv[5], "-s") == 0) {
                        g_sharp = 1;
                    }
                }
            }
        }
    }
    printf("manifold %s %s %d %f\n", argv[1], argv[2], resolution, min_hole_size);
    std::cout << "faces: " << obj.faces.size() << std::endl << std::flush;
    std::cout << "vertices: " << obj.vertices.size() << std::endl << std::flush;
    obj.Process_Manifold(resolution, min_hole_size);
    obj.SaveOBJ(argv[2]);
    return 0; 
}

#include "tiny_obj_loader.h"
#include "obj_writer.h"
#include <vector>
#include <iostream>


// g++ tiny_obj_loader_example.cpp tiny_obj_loader.cc obj_writer.cc

int main()
{
    std::string inputfile = "cube.obj";
    // load file and display to meshlab
    tinyobj::attrib_t attribs;
    std::vector<tinyobj::material_t> materials;
    std::vector<tinyobj::shape_t> shapes;

    std::string warn;
    std::string err;
    bool ret = tinyobj::LoadObj(&attribs, &shapes, &materials, &warn, &err, inputfile.c_str());
    
    const std::string outputfile = "output_cube.obj";

    bool write_ret = WriteObj(outputfile, attribs, shapes, materials, false);
    std::cout << write_ret;
    
    return 0;
}
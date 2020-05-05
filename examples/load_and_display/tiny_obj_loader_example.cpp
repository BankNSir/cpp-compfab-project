#include "../../myincludes/tiny_obj_loader.h"
#include "../../myincludes/obj_writer.h"
#include <vector>
#include <iostream>

void displayVertices(tinyobj::attrib_t &t_attribs) {

    // for accessing vertex
    std::vector<tinyobj::real_t> vertices_holder = t_attribs.vertices;
    int index = 0;
    for (auto &val : vertices_holder) {

        std::cout << val << " ";

        if ( index % 3 == 2 ) {
            std::cout << std::endl;
        }

        index++;
    }
    // std::cout << t_attribs << std::endl;
    int nVertices = index / 3;
    std::cout << "Number of Points: " << nVertices << std::endl;
    std::cout << "This is after debug" << std::endl;
}

void displayNormals(tinyobj::attrib_t &t_attribs) {

    // for accessing vertex
    std::vector<tinyobj::real_t> vertices_holder = t_attribs.normals;
    int index = 0;
    for (auto &val : vertices_holder) {

        std::cout << val << " ";

        if ( index % 3 == 2 ) {
            std::cout << std::endl;
        }

        index++;
    }
    // std::cout << t_attribs << std::endl;
    int nNormals = index / 3;
    std::cout << "Number of point: " << nNormals << std::endl;
    std::cout << "This is after debug" << std::endl;
}

// g++ tiny_obj_loader_example.cpp tiny_obj_loader.cc obj_writer.cc
// this example download cube.obj file in current directory and write an obj file named output_cube.obj
int main()
{
    // load obj file and write output obj file
    std::string inputfile = "cube.obj";

    tinyobj::attrib_t attribs;
    std::vector<tinyobj::material_t> materials;
    std::vector<tinyobj::shape_t> shapes;
    std::string warn;
    std::string err;

    bool load_ret = tinyobj::LoadObj(&attribs, &shapes, &materials, &warn, &err, inputfile.c_str());

    // Debug
    displayVertices(attribs);
    displayNormals(attribs);
    std::cout << "endDebug" << std::endl;

    const std::string outputfile = "output_cube.obj";

    bool write_ret = WriteObj(outputfile, attribs, shapes, materials, false);
    std::string write_ret_text;

    if(write_ret == true){
        write_ret_text = "Successfully write obj file";
    }
    else{
        write_ret_text = "Failed to write obj file";
    }

    std::cout << "exit status: " << write_ret_text << std::endl;
    
    return 0;
}
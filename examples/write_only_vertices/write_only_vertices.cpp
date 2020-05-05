#include "../../myincludes/tiny_obj_loader.h"
#include "../../myincludes/obj_writer.h"
#include <vector>
#include <iostream>

// This example shows that obj file can be read and write using only vertex, without other info

// compile : g++ tiny_obj_loader_example.cpp tiny_obj_loader.cc obj_writer.cc

int main()
{
    // load obj file and write output obj file
    std::string inputfile = "cube_only_vert.obj";

    tinyobj::attrib_t attribs;
    std::vector<tinyobj::material_t> materials;
    std::vector<tinyobj::shape_t> shapes;
    std::string warn;
    std::string err;

    bool load_ret = tinyobj::LoadObj(&attribs, &shapes, &materials, &warn, &err, inputfile.c_str());

    // This is for sending vertices only
    // std::vector<tinyobj::real_t> vertices_holder = attribs.vertices;
    // std::cout << "endDebug" << std::endl;
    
    // tinyobj::attrib_t onlyVertAttribs;
    // onlyVertAttribs.vertices = vertices_holder;

    const std::string outputfile = "output_cube_only_vert.obj";

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

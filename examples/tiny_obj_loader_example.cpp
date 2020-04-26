#include "tiny_obj_loader.h"
#include "obj_writer.h"
#include <vector>
#include <iostream>


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
#include <iostream>
#include <vector>
#include <string>

#include "Eigen/Dense"       // for linear algebra operation
// #include "tiny_obj_loader.h"    // load and save obj files

// using namespace Eigen;
// using Vector2d = Vec2;
using Vec2 = Eigen::Vector2d;
using Mat2 = Eigen::Matrix2d;
// using Eigen::Matrix2d = Mat2;

struct Particle
{
    // velocity waiting for eigen matrix vector 2x2
    Vec2 v;
    // position vector 2x2
    Vec2 x;
    // deformation gradient matrix 2x2
    Mat2 F;
    // APIC coefficient matrix 2x2
    Mat2 C;
    float mass;
};

// need to create grid struct here


int main()
{

    // Eigen::MatrixXd test_matrix(2, 2);
    // test_matrix(0, 0) = 3;
    // test_matrix(0, 1) = 2;
    // test_matrix(1, 0) = 2;
    // test_matrix(1, 1) = 3;
    // std::cout << test_matrix << std::endl;

    Eigen::Vector2d test_vector(2, 2);
    std::cout << test_vector << std::endl;

    // main program

    // I need Particle and Grid struct first
    // So i need to use eigen for vector and matrix
    
    // initialize() 

    // sendinfo from particle to grid
    // particleToGrid()

    // update grid info
    // gridUpdate()

    // sendinfo from grid to particle
    // gridToParticle()
    
    // write an objfile by particle vertices
    // visualize() 
}
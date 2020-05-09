// rignt now finished particle constructor
// Next test particle and P2G

#include <cmath>            // for floor()
#include <cstdlib>          // for rand()
#include <iostream>         // cout, cin
#include <string>           // string
#include <vector>           // vector<T>

#include "Eigen/Dense"       // for linear algebra operation
#include "myincludes/tiny_obj_loader.h"    // load and save obj files

// using namespace Eigen;
// using Vector2d = Vec2;
using Vec2 = Eigen::Vector2d;
using Mat2 = Eigen::Matrix2d;
// using Eigen::Matrix2d = Mat2;

struct Particle
{
    // velocity waiting for eigen matrix vector 2x2
    Vec2 velocity;
    // position vector 2x2
    Vec2 position;
    // deformation gradient matrix 2x2
    Mat2 F;
    // APIC coefficient matrix 2x2
    Mat2 C;
    // need precision (not float)
    double mass;

    Particle() = default;
    Particle(Vec2 t_position, Vec2 t_velocity = Vec2(0, 0))
    : position(t_position), velocity(t_velocity) // if this bug initialize vector by Eigen format v << 1, 2
    {
        // Eigen initialize format
        F << 1, 1, 1, 1;
        C << 0, 0, 0, 0;
    }
};

// need to create grid struct here
// with grid struct I need to have mass and velocity

struct Grid
{
    Vec2 momentum;
    Vec2 velocity;
    double mass;
};

// set particle data pass by reference -> I want to change particles without returning it
void initialize(std::vector<Particle> &particles, const std::vector<Vec2> &positions)
{
    int i = 0;
    for (auto &p: particles)
    {
        p.velocity = positions[i];
        i++;
    }
}

// This function create particle of Box bounding by bottom_left and top_right
// with each particle spacing is particle_space
std::vector<Vec2> createBox(Vec2 bottom_left, Vec2 top_right, double particle_space)
{   
    // 1) calc v2-v1 --> look Eigen
    // 2) nParticle_create = floor((v2-v1)/particle_space) --> Vec2/constant Eigen
    // backup: Vec2 nParticle_in_box = Vec2((top_right - bottom_left)/particle_space)
    
    Vec2 nParticle_in_box = Vec2((top_right - bottom_left)/particle_space);
    int xParticle_in_box = floor(nParticle_in_box(0));
    int yParticle_in_box = floor(nParticle_in_box(1));

    std::vector<Vec2> box_Coords;

    // 3) v1 + i*interval from i = 0 to i = nParticle_create use normal for loop
    // <= is correct
    for (int i = 0; i <= xParticle_in_box; i++)
    {
        for (int j = 0; j <= yParticle_in_box; j++)
        {
            box_Coords.push_back(bottom_left + Vec2(i, j)*particle_space);        
        }
    }
    return box_Coords;
}
// create and return a series of 2D-vector random position
// std::vector<Vec2> createRandomPosition(int min_range, int max_range)
// {
//     // 
//     srand()
//     std::vector<Vec2> rand_positions;
    
//     rand_positions.push_back(Vec2 )
//     for (auto &x: positions)
//     {
//         x = Vec2(0, 0);
//         std::cout << x;
//     }

//     return rand_positions;
// }


int main()
{

    Particle test_particle(Vec2(1, 1));


    std::vector<Particle> particles;
    std::vector<Vec2> box = createBox(Vec2(-1, -1), Vec2(1, 1), 0.1);

    particles.assign(box.begin(), box.end());


    for (auto it = particles.begin(); it != particles.end(); it++)
    {
        std::cout << (*it).position << std::endl;
    }
    
    // finished assign a box to particle position


    // main program
    
    // I need Particle and Grid struct first: CLEAR
    // So i need to use eigen for vector and matrix: CLEAR
    
    // createBox 
    // initialize(particles,) 

    // sendinfo from particle to grid
    // particleToGrid()

    // update grid info
    // gridUpdate()

    // sendinfo from grid to particle
    // gridToParticle()
    
    // write an objfile by particle vertices
    // visualize() 

    std::cout << "End of program" << std::endl;
}


// Eigen::Vector2d test_vector(2, 2);
// std::cout << test_vector << std::endl;

// Eigen::MatrixXd test_matrix(2, 2);
    // test_matrix(0, 0) = 3;
    // test_matrix(0, 1) = 2;
    // test_matrix(1, 0) = 2;
    // test_matrix(1, 1) = 3;
    // std::cout << test_matrix << std::endl;
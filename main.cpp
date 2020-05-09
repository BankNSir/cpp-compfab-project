// ParticleToGridCoord still doesn't handle boundary point -1 and 33

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

// number of grid boxes
const int N_GRID_BOX = 32;
// width of each grid boxes
const int dx = 1/N_GRID_BOX;
// number of Grid lines = number of box + 1
Grid GRID_LINE[N_GRID_BOX + 1][N_GRID_BOX + 1];

struct Particle
{
    // velocity waiting for eigen matrix vector 2x1
    Vec2 velocity;
    // position vector 2x1
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
    for (int index = 0; index < positions.size(); index++)
    {
        // Particle p = Particle(positions[index]);
        particles.push_back(Particle(positions[index]));
    }
    // for (auto &p: particles)
    // {
    //     p.position = positions[i];
    //     i++;
    // }
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





// pos coord in [0, 1)
// find which grid are in this particle kernel 
Vec2* ParticleToGridCoord(Particle particle)
{
    Vec2 *index = new Vec2[9];

    // pos [0, 1)
    // which box this particle in
    Vec2 grid_box_coord = particle.position * N_GRID_BOX;     //Vec20.5 * 32
        // if grid box in [0.5, 1.5) return 0 ans so onnn
    // DEBUG: std::cout << "grid_box_coord: " << grid_box_coord << std::endl;

    // pos coord to bottomleft grid line coords
    Vec2 tmp_grid_box_coord = grid_box_coord - Vec2(0.5, 0.5);
    // base_coord(x, y)
    Vec2 base_coord = Vec2(floor(tmp_grid_box_coord(0)), floor(tmp_grid_box_coord(1)));

    // set index3x3 from base_coord
    int count = 0;
    for (int i = 0; i <= 2; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            index[count] = base_coord + Vec2(i, j); 
            count++;
        }
    }

    return index;
}


int main()
{

    Particle test_particle(Vec2(1, 1));


    std::vector<Particle> particles;
    std::vector<Vec2> boxCoords = createBox(Vec2(-1, -1), Vec2(1, 1), 0.1);

    // particles.assign(box.begin(), box.end());


    // for (auto it = particles.begin(); it != particles.end(); it++)
    // {
    //     std::cout << (*it).position << std::endl;
    // }

    

    // main program
    
    // I need Particle and Grid struct first: CLEAR
    // So i need to use eigen for vector and matrix: CLEAR
    
    // createBox 
    // initialize(particles,) : CLEAR
    initialize(particles, boxCoords);
    // DEBUG
    // std::cout << "B4 loop";
    // for (int i = 0; i < particles.size(); i++)
    // {
    //     std::cout << particles[i].position << std::endl;
    // }


    // sendinfo from particle to grid
    // particleToGrid()

    // correct caveat: 0, 0 return -1, -1
    for (int i = 0; i < 9; i++)
    {
        Vec2 *gridCoords = ParticleToGridCoord(Particle(Vec2(0/N_GRID_BOX, 0/N_GRID_BOX)));
        std::cout << *(gridCoords + i) << std::endl;
    }


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
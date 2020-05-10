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
    // Vec2 momentum;
    Vec2 velocity;
    double mass;
};

// Grid Coord(v, m)
Grid GRID_LINE[N_GRID_BOX + 1][N_GRID_BOX + 1];

// set particle data pass by reference -> I want to change particles without returning it
void Initialize(std::vector<Particle> &particles, const std::vector<Vec2> &positions)
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
std::vector<Vec2> CreateBox(Vec2 bottom_left, Vec2 top_right, double particle_space)
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

// pos coord in [0, 1)
// find which grid are in this particle kernel 
std::vector<Vec2> GetGridCoord(Particle const &particle)
{
    // calculate pos x and pos y separately
    // transform to grid coords
    double dPosX = particle.position[0] * N_GRID_BOX;
    std::cout << "particle pos x: " << dPosX << std::endl;
    double dPosY = particle.position[1] * N_GRID_BOX;

    // minus 0.5
    // floor to get result
    int dBaseCoordX = floor(dPosX - 0.5);
    int dBaseCoordY = floor(dPosY - 0.5);

    std::vector<Vec2> surroundGridCoords;

    // create 3x3 by +1 +2 and push vec2(x, y) into vector<Vec2>
    for (int indexX = 0; indexX < 3; indexX++)
    {
        for (int indexY = 0; indexY < 3; indexY++)
        {
            int coordX = dBaseCoordX + indexX;
            int coordY = dBaseCoordY + indexY;

            // hack for boundary -1, N_GRID_BOX + 1
            if (coordX == -1)
            {
                coordX = 0; 
            }
            else if (coordX == N_GRID_BOX + 1)
            {
                coordX = N_GRID_BOX;
            }
            if (coordY == -1)
            {
                coordY = 0;
            }
            else if (coordY == N_GRID_BOX + 1)
            {
                coordY = N_GRID_BOX;
            }

            surroundGridCoords.push_back(Vec2(coordX, coordY));
        }
    }
    
    return surroundGridCoords;
}

// This function send info from particle to grid


int main()
{

    Particle test_particle(Vec2(1, 1));


    std::vector<Particle> particles;
    std::vector<Vec2> boxCoords = CreateBox(Vec2(-1, -1), Vec2(1, 1), 0.1);

    // main program
    
    // I need Particle and Grid struct first: CLEAR
    // So i need to use eigen for vector and matrix: CLEAR
    
    // createBox 
    // 1) Initialize(particles,) : CLEAR
    Initialize(particles, boxCoords);

    // sendinfo from particle to grid
    // 2) particleToGrid() : DO NOW

    // correct caveat:  -1, -1 and 33, 33
    // test GetGridCoord
    std::vector<Vec2> gridCoords;
    gridCoords = GetGridCoord(Particle(Vec2(0.0, 0.0)));
    // std::cout << gridCoords[0] << std::endl;
    for (auto &gridCoord: gridCoords)
    {
        std::cout << gridCoord << std::endl;
    }

    // create next function to call GetGridCoord and update momentum w.r.t gridCoord
    // ParticleToGrid function input: vector<Particle> function output: void just update grid data
    


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
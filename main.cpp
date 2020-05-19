// ParticleToGridCoord still doesn't handle boundary point -1 and 33

#include <cmath>            // for floor()
#include <cstdlib>          // for rand()
#include <iostream>         // cout, cin
#include <string>           // string
#include <vector>           // vector<T>
#include <typeinfo>

#include "Eigen/Dense"       // for linear algebra operation  compile: -I .
#include "myincludes/tiny_obj_loader.h"    // load and save obj files

// using namespace Eigen;
// using Vector2d = Vec2;
using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Mat2 = Eigen::Matrix2d;
// using Eigen::Matrix2d = Mat2;

// number of grid boxes
const int N_GRID_BOX = 32;
// width of each grid boxes
const double dx = 1/N_GRID_BOX;
// deltaTime for each simulation step
const double dt = 1e-4;
const double frame_dt = 1e-3;

// snow hardening
const auto particle_mass = 1.0;
const auto vol = 1.0;        // Particle Volume
const auto hardening = 10.0; 
const auto E = 1e4;          // Young's Modulus
const auto nu = 0.2;         // Poisson ratio

// Initial Lamé parameters
const double mu_0 = E / (2 * (1 + nu));
const double lambda_0 = E * nu / ((1+nu) * (1 - 2 * nu));

struct Particle
{
    // velocity waiting for eigen matrix vector 2x1
    Vec2 velocity;
    // position vector 2x1
    Vec2 position;
    // deformation gradient matrix 2x2
    Mat2 F;
    // J = det(F), current volume
    double J;
    // APIC coefficient matrix 2x2
    Mat2 C;
    // need precision (not float)
    double mass;

    Particle() = default;
    Particle(Vec2 t_position, Vec2 t_velocity = Vec2(0, 0), double J = 1)
    : position(t_position), velocity(t_velocity) // if this bug initialize vector by Eigen format v << 1, 2
    {
        // Eigen initialize format
        F << 1, 1, 1, 1;
        C << 0, 0, 0, 0;
    }
};

// need to create grid struct here
// with grid struct I need to have mass and velocity

// Grid Coord(v, m)
Vec3 GRID[N_GRID_BOX + 1][N_GRID_BOX + 1];

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

// element-wise power
Vec2 PowVec2(Vec2 vec, int exponent) 
{
    return Vec2(pow(vec[0], exponent), pow(vec[1], exponent));
}

double DetMat2(Mat2 mat) 
{
    return mat(0, 0)*mat(1, 1) - mat(1, 0)*mat(0, 1);
}

Eigen::Vector2i FloorVec2(Vec2 vec) 
{
    return Eigen::Vector2i(static_cast<int>(vec[0]), static_cast<int>(vec[1]));
}

// This function send info from particle to grid
void ParticleToGrid(const std::vector<Particle> &particles)
{
    for (auto &p: particles) {
        
        // compute fx
        // std::cout << p.position << std::endl;
        Eigen::Vector2i baseCoords = FloorVec2(p.position * N_GRID_BOX - Vec2(0.5, 0.5));
        // std::cout << baseCoords << std::endl;
        Vec2 fx = p.position * N_GRID_BOX - baseCoords.cast<double>();
        // weight
        Vec2 weight[3] = 
        { 
            0.5 * PowVec2(Vec2(1.5, 1.5) - fx, 2),
            Vec2(0.75, 0.75) - PowVec2(fx - Vec2(1.0, 1.0), 2),
            0.5 * PowVec2(fx - Vec2(0.5, 0.5), 2)
        };

        double e = std::exp(hardening * (1.0f - p.J));
        double mu = mu_0 * e;
        double lambda = lambda_0 * e;
        
        double J = DetMat2(p.F);
        double Dinv = 4*N_GRID_BOX*N_GRID_BOX;

        // polar decomposition of F
        Mat2 r = p.F.colPivHouseholderQr().matrixR();
       
        // usage m.transpose() ;

        // DEBUG HERE
        // std::cout << "Debug: " << typeid(2*mu*(p.F - r)).name() << std::endl;
        
        //BUG
        Mat2 rhMat;
        double rhTerm = lambda*(J-1)*J;
        rhMat << rhTerm, rhTerm, rhTerm, rhTerm;
        Mat2 PF = (2*mu*(p.F-r))*p.F.transpose() + rhMat;

        Mat2 stress = - (dt * vol) * (Dinv * PF);

        Mat2 affine = stress + particle_mass * p.C;

        // BUG HERE
        for (int i = 0; i < 3; i++) { 
            for (int j = 0; j < 3; j++) { 
                Vec2 dParticleToBaseCoords = (Vec2(i, j) - fx) * dx;
                Vec2 momentum = particle_mass * p.velocity;
                Vec3 momentumMass(momentum.x(), momentum.y(), particle_mass);
                Vec2 affMuldParticle = affine * dParticleToBaseCoords;
                GRID[baseCoords.x() + i][baseCoords.y() + j] += weight[i].x()*weight[j].y() 
                    * (momentumMass + Vec3(affMuldParticle.x(), affMuldParticle.y(), 0));
            }    
        }
    }
}

void UpdateGrid() {

    for (int i = 0; i < N_GRID_BOX+1; i++) {
        for (int j = 0; j < N_GRID_BOX+1; j++) {
            // update
            // // by mass then add by gravity minus at y coord * dt
            Vec3 &grid = GRID[i][j];
            if (grid[2] > 0) {

                grid /= grid[2]; // divide by mass
            
                grid += dt * Vec3(0, -200, 0); // gravity

                double boundary = 0.05;
                double nodeCoordsX = static_cast<double>(i) / N_GRID_BOX; 
                double nodeCoordsY = static_cast<double>(j) / N_GRID_BOX;
                printf("%.2f %.2f  ", nodeCoordsX, nodeCoordsY);
                // sticky
                if (nodeCoordsX < boundary || nodeCoordsX > 1-boundary || nodeCoordsY > 1-boundary) {
                    grid = Vec3(0, 0, 0);
                    std::cout << "STICKY" << std::endl;
                }            

                // seperate
                if (nodeCoordsY < boundary) {
                    grid[1] = std::max(0.0, grid[1]);
                }
            }
        }
    }
} 

// void GridToParticle(std::vector<Particle> &particles) {

//     for (auto &p : particles) {
        
//     }
// }

void ShowGrid(const Vec3 grid[][N_GRID_BOX+1]) {
    for (int i = 0; i < N_GRID_BOX + 1; i++) {
        std::cout << "Row: " << i << " ";
        for (int j = 0; j < N_GRID_BOX + 1; j++) {
            Vec3 temp(grid[i][j]);
            // std::cout << temp.x() << " " << temp.y() << " " << temp.z() << "| ";
            std::printf("%.3f %.3f %.3f | ", temp.x(), temp.y(), temp.z());
        }
        std::cout << std::endl << "________________________________________________________________" 
            << std::endl << std::endl;
    }
}

// Current volume
    // real J = determinant(p.F);

    // // Polar decomposition for fixed corotated model
    // Mat r, s;
    // polar_decomp(p.F, r, s);

    // // [http://mpm.graphics Paragraph after Eqn. 176]
    // real Dinv = 4 * inv_dx * inv_dx;
    // // [http://mpm.graphics Eqn. 52]
    // auto PF = (2 * mu * (p.F-r) * transposed(p.F) + lambda * (J-1) * J);

    // // Cauchy stress times dt and inv_dx
    // auto stress = - (dt * vol) * (Dinv * PF);

    // // Fused APIC momentum + MLS-MPM stress contribution
    // // See http://taichi.graphics/wp-content/uploads/2019/03/mls-mpm-cpic.pdf
    // // Eqn 29
    // auto affine = stress + particle_mass * p.C;

    // // P2G
    // for (int i = 0; i < 3; i++) {
    //   for (int j = 0; j < 3; j++) {
    //     auto dpos = (Vec(i, j) - fx) * dx;
    //     // Translational momentum
    //     Vector3 mass_x_velocity(p.v * particle_mass, particle_mass);
    //     grid[base_coord.x + i][base_coord.y + j] += (
    //       w[i].x*w[j].y * (mass_x_velocity + Vector3(affine * dpos, 0))
    //     );
    //   }
    // }

int main()
{

    Particle test_particle(Vec2(1, 1));


    std::vector<Particle> particles;
    // boxCoord, positoin must be in range [0, 1)
    std::vector<Vec2> boxCoords = CreateBox(Vec2(0.25, 0.25), Vec2(0.75, 0.75), 0.05);

    // main program
    
    // I need Particle and Grid struct first: CLEAR
    // So i need to use eigen for vector and matrix: CLEAR
    
    // createBox 
    // 1) Initialize(particles,) : CLEAR
    Initialize(particles, boxCoords);

    // sendinfo from particle to grid
    // 2) particleToGrid() : DO NOW
    // ParticleToGrid(particles); // I"M HERE
    ParticleToGrid(particles);

    ShowGrid(GRID);
    UpdateGrid();
    ShowGrid(GRID);
    // create next function to call GetGridCoord and update momentum w.r.t gridCoord
    // ParticleToGrid function input: vector<Particle> function output: void just update grid data
    
    // update grid info
    // gridUpdate()

    // sendinfo from grid to particle
    // gridToParticle()
    
    // write an objfile by particle vertices
    // visualize() 
    //DEBUG
    // Vec2 a(1.0);
    // std::cout << temp * 10 << std::endl;
    // std::cout << typeid(temp(0)).name() << std::endl;
    // std::cout << temp[0] << " " << temp[1] << std::endl << temp(0) << " " << temp(1) << " " << temp(0, 0) << std::endl;
    // Mat2 a(2, 2);
    // a << 1, 1, 1, 1;
    // Vec2 a(1.0, 2.0);
    // std::cout << a << std::endl;
    // std::cout << a.x() << " " << a.y() << std::endl;

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
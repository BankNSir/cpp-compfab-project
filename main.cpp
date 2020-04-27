#include <iostream>
#include <vector>
#include <string>

#include "Eigen/Dense"       // for linear algebra operation
// #include "tiny_obj_loader.h"    // load and save obj files


int main()
{

    Eigen::MatrixXd test_matrix(2, 2);
    test_matrix(0, 0) = 3;
    test_matrix(0, 1) = 2;
    test_matrix(1, 0) = 2;
    test_matrix(1, 1) = 3;
    std::cout << test_matrix << std::endl;

}
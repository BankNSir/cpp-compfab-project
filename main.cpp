#include <iostream>
#include <vector>
#include <string>

#include "Eigen/Dense"
#include "tiny_obj_loader.h"

using namespace std;

int main()
{
    vector<string> msg {"Hello", "C++", "World", "from", "VS Code", "and the C++ extension!"};

    for (const string& word : msg)
    {
        cout << word << " ";
    }
    cout << endl;
    int x = 10;
    cout << x << endl;
    cout << x;
    cout << " ";
    cout << endl;

    Eigen::MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(0, 1) = 2;
    m(1, 0) = 2;
    m(1, 1) = 3;
    cout << m << endl;

    // tinyobj::ObjReader()
}
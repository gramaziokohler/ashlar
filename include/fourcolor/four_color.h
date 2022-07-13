#ifndef FOUR_COLOR_H
#define FOUR_COLOR_H

#include "ashlar/eigen_stl_vector_specialization.h"
#include <Eigen/Core>
#include <vector>

using namespace std;

/*
An eigen based implementation of the four color theorem.  Modified from okaydemir's implementation: https://github.com/okaydemir/4-color-theorem

Inputs:  
F is an NxN integer matrix of 1s and 0s representing the adjacency of our graph.  It is a symmetric matrix, where N is the number of elements
and the diagonal (i.e. self adjacency) components are all zero.  Otherwise, an adjacency is represented as a 1

Outputs:
C is an Nx1 vector with an integer value for the color of that node.  Usually 0-3, but in the case that four colors doesn't work, 0-4 is used as a loose 5 color solution
*/

namespace fourcolor {
struct Node {//using nodes for graph representation
  int number = -1;
  int color = -1;
  vector<Node *> *links = new vector<Node *>;//addresses of all vertices that this vertice is linked to
};

void four_color_theorem(const Eigen::MatrixXi &F, Eigen::VectorXi &C);
void input(const Eigen::MatrixXi &F, vector<Node *> &vertices, int &N); // Reads input and updates vertices and links
vector<Node *> sortbyValence(vector<Node *> *arr); //Sorts the vector with quicksort according to their adjacency list size;
}
#endif
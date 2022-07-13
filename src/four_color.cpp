#include "fourcolor/four_color.h"

using namespace fourcolor;

void fourcolor::four_color_theorem(const Eigen::MatrixXi &F, Eigen::VectorXi &C) {
  vector < Node * > vertices;
  int N = 0;
  input(F, vertices, N);
  //beginning of Welsh-Powell algorithm
  vector < Node * > vertices_sorted = sortbyValence(&vertices);
  for (int color = 1; color < 5; color++) {
    for (int i = 0; i < N; i++) {
      if (vertices_sorted[i]->color == -1) {
        for (size_t j = 0; j < static_cast<size_t>(N); j++) {
          if (j < vertices_sorted[i]->links->size()) {
            if (vertices_sorted[i]->links->at(j)->color == color) {
              goto skip;//if one of the links has the same color, dont color
            }
          }
        }
        vertices_sorted[i]->color = color;
        skip:;
      }
    }
  }
  //end of Welsh-Powell algorithm, fill in our colors into our color vector
  //colors are 1-4 above, instead of correcting (because I am lazy), we shift them down by one and swap -1 with 4 for a fifth color if necessary.

  C.setZero(N);
  for (int i = 0; i < N; i++) {
    if (vertices[i]->color < 5 && vertices[i]->color > 0) {
      C(i) = vertices[i]->color - 1;
    } else {
      //graph can't be four colored...give it a fifth!
      C(i) = 4;
    }
  }
}

void fourcolor::input(const Eigen::MatrixXi &F, vector<Node *> &vertices, int &N) {
  //input our adjacency matrix values into our node vector
  N = F.cols(); //how many items do we have?
  for (int i = 0; i < N; i++) {
    Node *x = new Node;
    x->number = i;
    vertices.push_back(x);
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (F(i, j) == 1) {
        vertices[i]->links->push_back(vertices[j]);
      }
    }
  }
}

vector<Node *> fourcolor::sortbyValence(vector<Node *> *arr) {
  vector < Node * > ls;
  vector < Node * > eq;
  vector < Node * > gr;
  if (arr->size() <= 1) {
    return *arr;
  }
  size_t pivot = arr->at(0)->links->size();
  for (size_t i = 0; i < arr->size(); i++) {
    if (arr->at(i)->links->size() > pivot) {
      ls.push_back(arr->at(i));
    }
    if (arr->at(i)->links->size() == pivot) {
      eq.push_back(arr->at(i));
    }

    if (arr->at(i)->links->size() < pivot) {
      gr.push_back(arr->at(i));
    }
  }

  vector < Node * > lsf = sortbyValence(&ls);
  vector < Node * > grf = sortbyValence(&gr);
  lsf.insert(lsf.end(), eq.begin(), eq.end());
  lsf.insert(lsf.end(), grf.begin(), grf.end());
  return lsf;
}
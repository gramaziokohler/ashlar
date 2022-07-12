#include <iostream>

#include "igl/opengl/glfw/Viewer.h"


#include "ashlar/stone.h" //simple class for a stone object

ashlar::Stone stone;

using Viewer = igl::opengl::glfw::Viewer;


int main(int argc, char *argv[])
{
  Viewer viewer;

  std::string filename;
  if (argc == 2) {
    filename = std::string(argv[1]);
  }
  else {
    std::cout<<"No input mesh.  Please provide filename argument.\n";
  }

  //create new instance of a stone by loading the vertex and face data from the provided file
  stone = ashlar::Stone(filename);
  //set the density if you want to get estimated mass properties
  stone.SetDensity(2800.0);
  //reorient the mesh such that it is centered at the origin
  stone.Reorient();
  stone.ComputeProperties();
  stone.Downsample();
  stone.PrintProperties();

  // Plot the mesh
  viewer.data().clear();
  viewer.data().set_mesh(stone.mesh_.V,stone.mesh_.F);
  viewer.data().compute_normals();
  viewer.data().set_face_based(true);
  viewer.launch();
}

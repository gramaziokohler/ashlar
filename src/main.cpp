#include <iostream>
#include <fstream>
#include "igl/opengl/glfw/Viewer.h"
#include "ashlar/stone.h" //simple class for a stone object

ashlar::Stone stone;
std::ofstream logger;  //log our stone properties for later use

using Viewer = igl::opengl::glfw::Viewer;

int main(int argc, char *argv[])
{
  // Process command line arguments
  std::vector<std::string> filenames;
  if (argc < 2) {
    std::cout<<"No input mesh.  Please provide at least one filename argument.\n";
  }
  else {
    for (int i=1; i<argc; i++) {
      filenames.push_back(std::string(argv[i]));
    }
  }

  //initialize the log csv with headers
  logger.open("../logger.csv", std::ios::trunc);
  //For each process stone, we will write these properties to a logfile in the main directory
  //short, intermediate, and long bounding box dimensions (DS,DI,DL) are otherwise called width,length,height
  logger << "FILE,VOLUME,DS,DI,DL,FLATNESS,ELONGATION,RECTANGULARITY,SPHERICITY,METAFACES,ASHLARNESS\n";

  //make some variables to store averages
  double average_ashlarness = 0.0;
  double average_volume = 0.0;
  double average_ds = 0.0;
  double average_di = 0.0;
  double average_dl = 0.0;
  double average_sphericity = 0.0;
  double average_rectangularity = 0.0;

  // For all given filenames, compute, print, and save properties
  for(int i=0; i<filenames.size(); i++) {
    std::cout << "PROCESSING MESH " << i+1 << " OF " << filenames.size() << "\n";
    //create new instance of a stone by loading the vertex and face data from the provided file
    stone = ashlar::Stone(filenames.at(i));
    //set the density if you want to get estimated mass properties
    stone.SetDensity(2800.0);
    //reorient the mesh such that it is centered at the origin
    stone.Reorient();
    //compute intrinsics such as the bounding box and mass properties
    stone.ComputeProperties();
    //downsample the mesh to speed up VSA computation
    //(for very large meshes, this can be moved to before reorientation for further speedup)
    stone.Downsample();
    //run iterative variational shape approximation to determine the face properties
    stone.ComputeVSA(25, 10, 0.26);
    //compute the ashlarness metric for this object
    stone.ComputeAshlarness();
    //print the relevant properties of this stone to the console
    stone.PrintProperties();
    //log the relevant properties to the log file
    stone.LogProperties(logger);

    //add to our running average
    average_ashlarness+=stone.ashlarness_;
    average_volume+=stone.mass_props_.volume;
    average_ds+=stone.aabb_.width;
    average_di+=stone.aabb_.length;
    average_dl+=stone.aabb_.height;
    average_sphericity+=stone.mass_props_.sphericity;
    average_rectangularity+=stone.aabb_.rectangularity;
  }

  //close the logger
  logger.close();

  //Print our average properties if relevant
  if(filenames.size() > 1) {
    average_ashlarness /= (double) filenames.size();
    average_volume /= (double) filenames.size();
    average_ds /= (double) filenames.size();
    average_di /= (double) filenames.size();
    average_dl /= (double) filenames.size();
    average_sphericity /= (double) filenames.size();
    average_rectangularity /= (double) filenames.size();
    std::cout << "AVERAGE SET ASHLARNESS: " << average_ashlarness << "\n";
    std::cout << "AVERAGE SET VOLUME: " << average_volume << "\n";
    std::cout << "AVERAGE SET Ds (width): " << average_ds << "\n";
    std::cout << "AVERAGE SET Di (length): " << average_di << "\n";
    std::cout << "AVERAGE SET Dl (height): " << average_dl << "\n";
    std::cout << "AVERAGE SET SPHERICITY: " << average_sphericity << "\n";
    std::cout << "AVERAGE SET 3D RECTANGULARITY: " << average_rectangularity << "\n";
    std::cout<< "=========================================\n";
  }

  // Plot the most recently processed mesh
  Viewer viewer;
  viewer.core(viewer.core_list[0].id).background_color.setConstant(0.7);
  viewer.data().clear();
  viewer.data().set_mesh(stone.mesh_.V,stone.mesh_.F);
  viewer.data().set_colors(stone.vsa_.per_face_colors);

  // Plot the edges of the bounding box
  for (unsigned i=0;i<stone.aabb_.edges.rows(); ++i)
    viewer.data().add_edges
        (
            stone.aabb_.corners.row(stone.aabb_.edges(i,0)),
            stone.aabb_.corners.row(stone.aabb_.edges(i,1)),
            Eigen::RowVector3d(0.5,0.5,0.5)
        );

  // Plot the edge midpoints
  viewer.data().add_points(stone.vsa_.flattened_edge_centers,Eigen::RowVector3d(1,0,0));

  //Launch Plotter
  viewer.data().compute_normals();
  viewer.data().set_face_based(true);
  viewer.launch_init(true, false, "ashlar");
  viewer.launch_rendering();
  viewer.launch_shut();
}

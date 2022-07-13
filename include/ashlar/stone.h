#ifndef EXAMPLE_SRC_STONE_H_
#define EXAMPLE_SRC_STONE_H_

#include <fstream>
#include <Eigen/Core>
#include "eigen_stl_vector_specialization.h"
#include "ashlar/properties.h"

namespace ashlar {
class Stone {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW



  Stone();
  Stone(const std::string filename);

  void SetDensity(double density);
  void Reorient();   //reorient the stone such that its eigenvectors are aligned to the world reference frame, and it is centered on its centroid
  void ComputeProperties(); //compute key properties such as bounding box, surface area, and mass properties
  void ComputeVSA(int max_clusters, int subiterations, double max_error); //compute variational shape approximation for this rock
  void ComputeAshlarness();  //compute the ashlarness metric for this stone
  void PrintProperties(); //print the key stone properties to the console
  void LogProperties(std::ofstream& logger);
  void Downsample();  //downsample the mesh with consistent face sizes

  /*!
   * Compute the moment of intertia and mass from a triangle mesh
   * @param V vertices of the input mesh
   * @param F faces of the input mesh
   * @param density specified density
   * @param mass_props [output] mass properties to set
   */
  void MassProps(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double& density, ashlar::MassProperties& mass_props);

  /*!
   * Compute the volume of a triangle mesh
   * @param verts vertices of the input mesh
   * @param faces faces of the input mesh
   * @return the volume of the input mesh (i.e. m3 if mesh units are meters)
   */
  double ComputeVolume(const Eigen::MatrixXd& verts, const Eigen::MatrixXi& faces);

  /*!
   * Get a rotation matrix given a #Vx3 matrix of points, such that the geometry is aligned with the principal axes
   * @param verts vertices of the input mesh
   * @param rotation [output] rotation matrix to align the mesh
   */
  void ComputeRotation(const Eigen::MatrixXd& verts, Eigen::Matrix3d& rotation);  //compute rotation for reorienting stone using PCA



  //public stone properties

  //store the mesh properties
  ashlar::MeshProperties mesh_;
  //store the dimension properties of our axis aligned bounding box
  ashlar::BoxProperties aabb_;
  //store our mass properties
  ashlar::MassProperties mass_props_;
  //store our VSA properties
  ashlar::VsaProperties vsa_;
  double ashlarness_ = 0.0;  //store ashlarness for this stone

 private:
  double density_ = 2800.0;  //density in kg/m3 (asuming mesh units are meters)
  std::string filename_;
};
}

#endif //EXAMPLE_SRC_STONE_H_

#include "ashlar/stone.h"

#include "igl/read_triangle_mesh.h"

#include <igl/edges.h>
#include <igl/decimate.h>
#include "igl/volume.h"
#include "igl/barycenter.h"
#include "igl/adjacency_list.h"
#include "igl/triangle_triangle_adjacency.h"
#include "igl/loop.h"
#include "igl/centroid.h"
#include "igl/per_corner_normals.h"
#include "igl/per_vertex_normals.h"
#include "igl/per_face_normals.h"
#include "igl/random_points_on_mesh.h"

ashlar::Stone::Stone(){

}

ashlar::Stone::Stone(std::string filename){
  igl::read_triangle_mesh(filename,mesh_.V,mesh_.F);
}

void ashlar::Stone::SetDensity(double density){
  density_ = density;
}
double ashlar::Stone::ComputeVolume(const Eigen::MatrixXd& verts, const Eigen::MatrixXi& faces){
  // from https://github.com/libigl/libigl/issues/694
  Eigen::MatrixXd V2(verts.rows() + 1, verts.cols());
  V2.topRows(verts.rows()) = verts;
  V2.bottomRows(1).setZero();
  Eigen::MatrixXi T(faces.rows(), 4);
  T.leftCols(3) = faces;
  T.rightCols(1).setConstant(verts.rows());
  Eigen::VectorXd vol;
  igl::volume(V2, T, vol);
  return std::abs(vol.sum());
}

void ashlar::Stone::ComputeProperties(){
  //get mesh volume
  mass_props_.volume = ComputeVolume(mesh_.V,mesh_.F);

  // get surface area of mesh
  Eigen::VectorXd face_areas;
  igl::doublearea(mesh_.V,mesh_.F, face_areas);  // 2* the area of each face #F x 1
  mass_props_.surface_area = face_areas.sum() / 2.0;

  //get mass properties
  MassProps(mesh_.V, mesh_.F, density_, mass_props_);

  //get bounding box properties
  aabb_.min = mesh_.V.colwise().minCoeff();
  aabb_.max = mesh_.V.colwise().maxCoeff();
  aabb_.center = (aabb_.min + aabb_.max) / 2;
  Eigen::RowVector3d dim = aabb_.max - aabb_.min;
  aabb_.width = dim[0];
  aabb_.length = dim[1];
  aabb_.height = dim[2];
  aabb_.diagonal = dim.norm();  // get the length of our aabb_diagonal to give us another metric of general rock size

  //compute form factors (3D rectangularity, elongation, flatness)
  double aabb_volume = aabb_.width * aabb_.length * aabb_.height;
  aabb_.rectangularity = mass_props_.volume / aabb_volume;
  aabb_.elongation = aabb_.length / aabb_.height;
  aabb_.flatness = aabb_.width / aabb_.length;
}

void ashlar::Stone::ComputeVSA() {

}

void ashlar::Stone::PrintProperties() {
  std::cout<< "Surface Area: " << mass_props_.surface_area << "\n";
  std::cout<< "Volume: " << mass_props_.volume << "\n";
  std::cout<< "Mass: " << mass_props_.mass << " (from estimated density)\n";
  std::cout<< "Width: " << aabb_.width << "\n";
  std::cout<< "Length: " << aabb_.length << "\n";
  std::cout<< "Height: " << aabb_.height << "\n";
  std::cout<< "3D Rectangularity: " << aabb_.rectangularity << "\n";
  std::cout<< "Flatness: " << aabb_.flatness << "\n";
  std::cout<< "Elongation: " << aabb_.elongation << "\n";
}

void ashlar::Stone::ComputeRotation(const Eigen::MatrixXd& verts, Eigen::Matrix3d& rotation) {
  // use PCA to reorient our rock
  Eigen::MatrixXd Y = verts.rowwise() - verts.colwise().mean();  // each point minus the centroid m
  Eigen::MatrixXd S = Y.adjoint() * Y;                           // scatterMatrix S = Y * Y.transpose();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs(S,
                                                      Eigen::ComputeEigenvectors);  // get the eigen vectors and values

  Eigen::Matrix3d xform_a = Eigen::Matrix3d::Identity() * eigs.eigenvectors();
  Eigen::Vector3d x_axis = xform_a.col(0).normalized();
  Eigen::Vector3d z_axis = x_axis.cross(xform_a.col(1).normalized());
  Eigen::Vector3d y_axis = (z_axis.cross(x_axis)).normalized();

  rotation.setZero();
  rotation.col(0) = x_axis;
  rotation.col(1) = y_axis;
  rotation.col(2) = z_axis;
}

void ashlar::Stone::Reorient(){
  //usually we reorient based on mesh vertices using PCA,
  //but in this example we can't be guaranteed an even distribution.
  //Instead, we get the mesh orientation by populating the mesh with points
  int num_samples = 2000;
  Eigen::VectorXi FI;
  Eigen::SparseMatrix<double> BC; //the barycentric coords of each point on the mesh
  igl::random_points_on_mesh(num_samples, mesh_.V,mesh_.F, BC, FI);
  Eigen::MatrixXd P; //#num_samples x 3 matrix of points populated on the mesh surface
  P = BC * mesh_.V; //convert bc coords to global with libigl magic https://github.com/libigl/libigl/issues/808

  //compute rotation using PCA on populated points
  Eigen::Matrix3d rotation;
  ComputeRotation(P,rotation);

  //compute translation using mesh centroid
  double volume = 0.0;
  Eigen::Vector3d centroid;
  igl::centroid(mesh_.V,mesh_.F, centroid, volume);

  // create a combined transformation
  Eigen::Transform<double, 3, Eigen::Affine> transform(rotation);
  transform = Eigen::Translation3d(centroid) * transform;

  //apply transform to mesh vertices
  std::cout << "Centering mesh by applying translation: (" << centroid.transpose() << ")\n";
  std::cout << "and rotation:\n" << rotation << "\n";
  for(int i=0; i<mesh_.V.rows(); i++){
    Eigen::Vector3d pt = mesh_.V.row(i).transpose();
    mesh_.V.row(i) = (transform * pt).transpose();
  }
}

void ashlar::Stone::Downsample(){
    //downsampling to somewhat evenly sized faces speeds up the VSA calculation
    //and improves the consistency of the results with the exisiting dataset
    //for the purpose of this example, downsample to a given number of faces
    //determined by scaling the volume to our typical volume of 1/3 m3
    //and then getting a number of faces based on a fixed density of 500 faces/m2
    double target_volume = 0.33;
    double face_density = 500.0;
    double scale_factor = cbrt(target_volume / mass_props_.volume);
    double scaled_area = scale_factor*mass_props_.surface_area;
    int num_samples = ceil(face_density*scaled_area);
    Eigen::VectorXi J, I; //references to original birth faces and birth vertices
    Eigen::MatrixXd v_temp = mesh_.V;
    Eigen::MatrixXi f_temp = mesh_.F;
    igl::decimate(v_temp, f_temp, num_samples, mesh_.V, mesh_.F, J, I);
    std::cout << "Mesh decimated from " << f_temp.rows() << " to " << mesh_.F.rows() << " faces for VSA approximation.\n";
}

// compute mass properties from mesh
void ashlar::Stone::MassProps(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double& density, ashlar::MassProperties& mass_props) {
  // Method from supplemental material for "Spin-it: Optimizing moment of inertia for spinnable objects"
  // http://www.baecher.info/publications/spin_it_sup_mat_sig14.pdf
  // s10 implementation by Simon Huber <simon.huber@inf.ethz.ch>, M_I from RLJ <rjohns@ethz.ch>
  // phi is density
  // triangle vertices are counter clockwise: a, b, c

  double phi = density;  // set our phi to our rock density
  Eigen::VectorXd s10;   // used to store mass properties
  s10.resize(10);
  s10.setZero();

  for (int i = 0; i < F.rows(); i++) {                      // for each face
    Eigen::Vector3d a = V.row(F.row(i)[0]).transpose();  // get triangle verts of this face
    Eigen::Vector3d b = V.row(F.row(i)[1]).transpose();
    Eigen::Vector3d c = V.row(F.row(i)[2]).transpose();

    Eigen::Vector3d u = (b - a), v = (c - a);  // get edge vectors
    Eigen::Vector3d n = u.cross(v);      // get normal

    Eigen::Vector3d h1 = a + b + c;
    Eigen::Vector3d h2 = a.array() * a.array() + b.array() * (a.array() + b.array());  // .array();
    Eigen::Vector3d h3 = h2.array() + c.array() * h1.array();
    Eigen::Vector3d h4 = a.array() * a.array() * a.array() + b.array() * h2.array() + c.array() * h3.array();
    Eigen::Vector3d h5 = h3.array() + a.array() * (h1.array() + a.array());
    Eigen::Vector3d h6 = h3.array() + b.array() * (h1.array() + b.array());
    Eigen::Vector3d h7 = h3.array() + c.array() * (h1.array() + c.array());
    Eigen::Vector3d a_bar;
    a_bar << a[1], a[2], a[0];  // see notation reference from http://baecher.info/publications/spin_it_sup_mat_sig14.pdf
    Eigen::Vector3d b_bar;
    b_bar << b[1], b[2], b[0];
    Eigen::Vector3d c_bar;
    c_bar << c[1], c[2], c[0];

    Eigen::Vector3d h8 = a_bar.array() * h5.array() + b_bar.array() * h6.array() + c_bar.array() * h7.array();
    s10[0] += (n.array() * h1.array())[0];

    Eigen::Vector3d stemp = n.array() * h3.array();
    s10[1] += stemp[0];
    s10[2] += stemp[1];
    s10[3] += stemp[2];

    stemp = n.array() * h8.array();
    s10[4] += stemp[0];
    s10[5] += stemp[1];
    s10[6] += stemp[2];

    stemp = n.array() * h4.array();
    s10[7] += stemp[0];
    s10[8] += stemp[1];
    s10[9] += stemp[2];
  }

  s10[0] *= 1. / 6.;
  s10[1] *= 1. / 24.;
  s10[2] *= 1. / 24.;
  s10[3] *= 1. / 24.;
  s10[4] *= 1. / 120.;
  s10[5] *= 1. / 120.;
  s10[6] *= 1. / 120.;
  s10[7] *= 1. / 60.;
  s10[8] *= 1. / 60.;
  s10[9] *= 1. / 60.;

  s10 *= phi;

  if (s10[0] < 0) {  //fix randomly inverted masses...?
    s10 *= -1.0;
  }

  // get the moment of inertia from our s vector
  /*
  for quick reference:
  0 = S_1
  1 = S_x
  2 = S_y
  3 = S_z
  4 = S_xy
  5 = S_yz
  6 = S_xz
  7 = S_x^2
  8 = S_y^2
  9 = S_z^2
  */

  Eigen::Matrix3d M_I;
  M_I.setZero();
  M_I(0, 0) = s10[8] + s10[9];
  M_I(0, 1) = -1.0 * s10[4];
  M_I(0, 2) = -1.0 * s10[6];

  M_I(1, 0) = -1.0 * s10[4];
  M_I(1, 1) = s10[7] + s10[9];
  M_I(1, 2) = -1.0 * s10[5];

  M_I(2, 0) = -1.0 * s10[6];
  M_I(2, 1) = -1.0 * s10[5];
  M_I(2, 2) = s10[7] + s10[8];

  mass_props.moment_of_inertia = M_I;  // save our moment of inertia
  mass_props.mass = s10[0];   // get our computed mass
  mass_props.centroid << s10[1], s10[2], s10[3];   //computed a different way from igl centroid redundantly, but should be close to zero
}

#include "ashlar/stone.h"

#include <math.h>

#include "igl/volume.h"
#include "igl/decimate.h"
#include "igl/centroid.h"
#include "igl/doublearea.h"
#include "igl/barycenter.h"
#include "igl/adjacency_list.h"
#include "igl/per_face_normals.h"
#include "igl/read_triangle_mesh.h"
#include "igl/random_points_on_mesh.h"
#include "igl/triangle_triangle_adjacency.h"

#include "ashlar/variational_shape_approximation.h"

ashlar::Stone::Stone(){
}

ashlar::Stone::Stone(std::string filename){
  filename_ = filename;
  std::cout<< "=========================================\n";
  std::cout<< "File: " << filename_ << "\n";
  igl::read_triangle_mesh(filename_,mesh_.V,mesh_.F);
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
  std::cout << "Computing mesh properties...\n";
  //get mesh volume
  mass_props_.volume = ComputeVolume(mesh_.V,mesh_.F);

  // get surface area of mesh
  Eigen::VectorXd face_areas;
  igl::doublearea(mesh_.V,mesh_.F, face_areas);  // 2* the area of each face #F x 1
  mass_props_.surface_area = face_areas.sum() / 2.0;
  mass_props_.sphericity = (pow(M_PI,(1.0/3.0))*pow(6.0*mass_props_.volume,(2.0/3.0)))/mass_props_.surface_area;

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

  //store corners of our box for plotting
  aabb_.corners.setZero(8,3);
  aabb_.corners <<
      aabb_.min(0), aabb_.min(1), aabb_.min(2),
      aabb_.max(0), aabb_.min(1), aabb_.min(2),
      aabb_.max(0), aabb_.max(1), aabb_.min(2),
      aabb_.min(0), aabb_.max(1), aabb_.min(2),
      aabb_.min(0), aabb_.min(1), aabb_.max(2),
      aabb_.max(0), aabb_.min(1), aabb_.max(2),
      aabb_.max(0), aabb_.max(1), aabb_.max(2),
      aabb_.min(0), aabb_.max(1), aabb_.max(2);

  //edges always reference the same corners
  aabb_.edges.setZero(12,2);
  aabb_.edges<<
      0, 1,
      1, 2,
      2, 3,
      3, 0,
      4, 5,
      5, 6,
      6, 7,
      7, 4,
      0, 4,
      1, 5,
      2, 6,
      7 ,3;

  //compute form factors (3D rectangularity, elongation, flatness)
  double aabb_volume = aabb_.width * aabb_.length * aabb_.height;
  aabb_.rectangularity = mass_props_.volume / aabb_volume;
  aabb_.elongation = aabb_.length / aabb_.height;
  aabb_.flatness = aabb_.width / aabb_.length;

}

void ashlar::Stone::ComputeVSA(int max_clusters, int subiterations, double max_error) {
  std::cout << "Computing iterative shape approximation...\n";
  ashlar::VsaParams params;
  params.max_clusters = max_clusters; //if our resolution is set super high, stop before we get here
  params.subiterations = subiterations; //how many subiterations per loop
  params.max_error = max_error; //target weighted normal error per region
  ashlar::VariationalShapeApproximation(params,mesh_,vsa_);
}

void ashlar::Stone::ComputeAshlarness() {
  double min_size = mass_props_.surface_area * .02;  //arbitrary smallest allowable surface area for a primary meta-face in ashlarness metric
  int min_count = 2;  //fewest triangles allowed to consider a meta-face as primary ashlarness face (in case there are small islands)
  int skip_count = 0; //count how many faces we skipped due to the above restrictions

  std::vector<Eigen::RowVector3d> normal_vectors;  //store the normal vectors of all non-skipped regions
  std::vector<double> pseudoface_areas;  //store the total area of all non-skipped pseudoface regions
  int proxy_pseudoface = 0;  //store the index of the face that is most closely aligned to the broad axis of the bounding box (either side), weighted by area
  double closest_alignment = -1.0;  //get the area weighted alignment of the proxy pseudoface with the "broad normal" of the bounding box
  Eigen::RowVector3d broad_normal = Eigen::RowVector3d::UnitX();  //it is assumed that our objects are aligned such that the largest bounding box face is perpendicular to the x axis
  Eigen::RowVector3d aligned_normal = Eigen::RowVector3d::UnitX();  //the normal of the face that is most closely aligned to the x axis considering area weighting

  //precompute properties of vsa meta-faces, and find the region that is most closely aligned with the broad face of the bounding box
  for (int i = 0; i < vsa_.region_faces.size(); i++) {
    //Extract the mesh of each meta-face region
    Eigen::VectorXi IM, inRows;
    Eigen::MatrixXd VC;
    Eigen::MatrixXi FC;
    igl::remove_unreferenced(mesh_.V, vsa_.region_faces.at(i), VC, FC, IM, inRows); //remove all unused verts/faces

    //get the area of this region
    Eigen::VectorXd ps_area;
    igl::doublearea(VC,FC, ps_area); //2* the area of each face #F x 1
    double total_area = ps_area.sum() / 2.0;

    //recompute the normals for this region and get the weighted average normal
    Eigen::MatrixXd recomputed_norms;
    igl::per_face_normals(VC, FC, recomputed_norms);
    recomputed_norms = recomputed_norms.array().colwise() * ps_area.array();
    Eigen::RowVector3d average_norm = recomputed_norms.colwise().mean();
    average_norm.normalize();

    if(total_area > min_size && VC.rows() > min_count){
      //store the area and normal of this region for later use
      normal_vectors.push_back(vsa_.region_normals.row(i));
      pseudoface_areas.push_back(total_area);

      //get how closely this face approximates our broad face (area weighted)
      double alignment_score = total_area * abs(broad_normal.dot(vsa_.region_normals.row(i)));
      if(alignment_score > closest_alignment){
        //this is our best-yet meta-face in terms of being aligned with the bounding box broad-face
        proxy_pseudoface = normal_vectors.size()-1;  //our current object is the best
        closest_alignment = alignment_score;
        aligned_normal = vsa_.region_normals.row(i);
      }
    }
    else{
      skip_count++;
    }
  }

  //get area of whole stone again (in case this was done before downsampling)
  Eigen::VectorXd rock_face_areas;
  igl::doublearea(mesh_.V,mesh_.F,rock_face_areas); //2* the area of each face #F x 1
  double rock_area = rock_face_areas.sum() / 2.0;

  //get area of most aligned area weighted pseudoface
  double area_1 = pseudoface_areas.at(proxy_pseudoface);

  //find best ashlarness value
  double best_ashlarness = -1.0;

  for(int i=0; i<normal_vectors.size(); i++){
    if(i != proxy_pseudoface){  //don't compare with self
      double alignment = aligned_normal.dot(normal_vectors.at(i));  //get alignment between this face and the proxy meta face
      if(alignment < 0){  //they need to be opposing faces for ashlarness metric
        double curr_area = pseudoface_areas.at(i);  //area of this meta face
        double curr_area_ratio = area_1 > curr_area ? curr_area / area_1 : area_1 / curr_area;  //area ratio between these two meta faces
        double curr_flat_ratio = (area_1 + curr_area) / rock_area;  //what percentage of the total rock area are these flat faces?
        double curr_alignment = abs(alignment);  //how aligned are the faces? (absolute value of dot product)
        double curr_ashlarness = (curr_area_ratio + curr_flat_ratio + curr_alignment*curr_alignment)/3;  //ashlarness
        if(curr_ashlarness > best_ashlarness){
          best_ashlarness = curr_ashlarness;
        }
      }
    }
  }

  //save our computed value
  ashlarness_ = best_ashlarness;
}

void ashlar::Stone::PrintProperties() {
  std::cout<< "Surface Area: " << mass_props_.surface_area << "\n";
  std::cout<< "Volume: " << mass_props_.volume << "\n";
  std::cout<< "Mass: " << mass_props_.mass << " (from estimated density)\n";
  std::cout<< "Width: " << aabb_.width << "\n";
  std::cout<< "Length: " << aabb_.length << "\n";
  std::cout<< "Height: " << aabb_.height << "\n";
  std::cout<< "3D Rectangularity: " << aabb_.rectangularity << "\n";
  std::cout<< "Sphericity: " << mass_props_.sphericity << "\n";
  std::cout<< "Flatness: " << aabb_.flatness << "\n";
  std::cout<< "Elongation: " << aabb_.elongation << "\n";
  std::cout<< "Ashlarness: " << ashlarness_ << "\n";
  std::cout<< "=========================================\n";
}

void ashlar::Stone::LogProperties(std::ofstream& logger){
  //FILE,VOLUME,WIDTH,LENGTH,HEIGHT,FLATNESS,ELONGATION,RECTANGULARITY,SPHERICITY,METAFACES,ASHLARNESS
  logger<<filename_<<",";
  logger<<mass_props_.volume<<",";
  logger<<aabb_.width<<",";
  logger<<aabb_.length<<",";
  logger<<aabb_.height<<",";
  logger<<aabb_.flatness<<",";
  logger<<aabb_.elongation<<",";
  logger<<aabb_.rectangularity<<",";
  logger<<mass_props_.sphericity<<",";
  logger<<vsa_.region_count<<",";
  logger<<ashlarness_<<"\n";
}

void ashlar::Stone::ComputeRotation(const Eigen::MatrixXd& verts, Eigen::Matrix3d& rotation) {
  // use PCA to reorient our rock
  Eigen::MatrixXd Y = verts.rowwise() - verts.colwise().mean();
  Eigen::MatrixXd S = Y.adjoint() * Y;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs(S,Eigen::ComputeEigenvectors);

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
    std::cout << "Downsampling...\n";
    //downsampling to somewhat evenly sized faces speeds up the VSA calculation
    //and improves the consistency of the results with the exisiting dataset
    //for the purpose of this example, I downsample to a given number of faces
    //determined by scaling the volume to our typical volume of about 1/3 m3
    //and then getting a number of faces that object would have given surface area
    //and based on a fixed density of 500 faces/m2
    double target_volume = 0.33;
    double face_density = 500.0;
    double scale_factor = cbrt(target_volume / mass_props_.volume);
    double scaled_area = scale_factor*scale_factor*mass_props_.surface_area;
    int num_samples = ceil(face_density*scaled_area);
    Eigen::VectorXi J, I; //references to original birth faces and birth vertices
    Eigen::MatrixXd v_temp = mesh_.V;
    Eigen::MatrixXi f_temp = mesh_.F;
    igl::decimate(v_temp, f_temp, num_samples, mesh_.V, mesh_.F, J, I);
    std::cout << "Mesh decimated from " << f_temp.rows() << " to " << mesh_.F.rows() << " faces.\n";
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
    a_bar << a[1], a[2], a[0];
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

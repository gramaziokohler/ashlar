#ifndef ASHLAR_INCLUDE_ASHLAR_PROPERTIES_H_
#define ASHLAR_INCLUDE_ASHLAR_PROPERTIES_H_

#include <vector>
#include <Eigen/Core>
#include "eigen_stl_vector_specialization.h"

namespace ashlar {

struct MeshProperties{
  // Vertex array, #V x3
  Eigen::MatrixXd V;
  // Face array, #F x3
  Eigen::MatrixXi F;
  // Per-face normal array, #F x3
  Eigen::MatrixXd FN;
  // Per-vertex normal array, #V x3
  Eigen::MatrixXd VN;
  // Per-corner normal array, (3#F) x3
  Eigen::MatrixXd CN;
  // Vectors of indices for adjacency relations
  std::vector<std::vector<int>> VF, VFi, VV;
  // Integer vector of component IDs per face, #F x1
  Eigen::VectorXi cid;
  // Per-face color array, #F x3
  Eigen::MatrixXd component_colors_per_face;
};

struct MassProperties{
  double volume = 0.0;               //volume, computed from mesh
  double surface_area = 0.0;         //surface area of mesh
  Eigen::Vector3d centroid;          //centroid
  Eigen::Matrix3d moment_of_inertia; //moment of intertia, computed from mesh
  double mass = 0.0;                 //computed mass based on specified density
  double sphericity = 0.0;           //compute sphericity from area and volume
};

struct BoxProperties{
  Eigen::Vector3d min;    //min corner of the axis aligned bounding box
  Eigen::Vector3d max;    //min corner of the axis aligned bounding box
  Eigen::Vector3d center; //center of the axis aligned bounding box
  Eigen::MatrixXd corners;  //corners of the bounding box
  Eigen::MatrixXi edges;    //12x2 matrix of edges referencing corners of box (used for plotting
  double width;   //S1 of the axis aligned bounding box
  double length;  //S2 of the axis aligned bounding box
  double height;  //S3 of the axis aligned bounding box
  double diagonal;//diagonal length of the axis aligned bounding box
  double elongation = 0.0;           //d_i / d_l
  double flatness = 0.0;             //d_s / d_i
  double rectangularity = 0.0;       //volume object / volume bounding box
};

//store computed variational shape approximation properties
struct VsaProperties {
  // http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  //Each region is assigned a color index (0-4) stored here
  Eigen::VectorXi region_color_key;
  //For each region, a list of faces referencing the original mesh faces (F)
  std::vector<Eigen::MatrixXi> region_faces;
  std::vector<std::vector<int>> region_face_indices; //region_face_indices_; //list of lists where rfi.at(i).at(j) indexes faces in F for each region i
  Eigen::VectorXi region_key; //region_key_; //#F x 1 vector where face_region_indices_(i) provides the region that face i belongs to
  Eigen::MatrixXd region_centers; //regionCenters; //#_region x 3, 3d location of each region center (on the plane)
  Eigen::MatrixXd region_normals; //regionNormals; //#_region x 3, normal direction of face approximation
  Eigen::MatrixXi adjacencies; //flattened list of adjacency regions with valid edges/midpoints, where MPA(0,0) and MPA(0,1) are the indices of the adjacent regions (useful for looping through edges)
  std::map<int, std::map<int, Eigen::RowVector3d>> edge_centers; //map of 3d midpoints of each region edge, where edgeCenters.at(4).at(1) is a row vector of the x,y,z location of that shared region midpoint
  std::map<int, eigen_map<int, Eigen::Quaterniond>> edge_frames; //edgeFrames;
  Eigen::MatrixXd flattened_edge_centers; //store a #midpoints x 3 matrix of 3d coordinates that can be referenced in the map with vsa_.adjacencies

  std::map<int, std::map<int, int>> edge_center_proxies; //map of index in V of nearest vert to the region region edge center
  Eigen::VectorXi proxy_faces; //#_region x 1 index of proxy face for each region (the triangle index on the mesh that best approximates the clustered face)
  //map of edges per region combination, where edgePaths.at(4).at(1).row(3).col(0) is the index of the first vert in V of the 4th triangle edge on the edge between region 4 and 1
  std::map<int, std::map<int, Eigen::VectorXi>>edge_paths;

  Eigen::MatrixXi stable_adjacencies;  //#_regions x #_regions adjacency matrix where 1 indicates an edge that the centroid projects onto (Legacy, was used for flipping stones)
  std::map<int, std::map<int, double>> edge_lengths; //length from the start to end point of a given edge
  std::map<int, std::map<int, Eigen::RowVector3d>> edge_dirs;  //direction of a given edge between two regions

  std::map<int, std::map<int, double>> inter_angles; //map of angles between regions, where interAngles.at(2).at(4) is the angle between region 2 and 4
  int region_count = 0;  //how many regions do we have

  Eigen::MatrixXd per_face_colors;  //per face colors for the original mesh, derived from the four color region map
};

}


#endif //ASHLAR_INCLUDE_ASHLAR_PROPERTIES_H_

#ifndef SHAPE_APPROXIMATION_H
#define SHAPE_APPROXIMATION_H

#include "eigen_stl_vector_specialization.h"
#include "igl/igl_inline.h"
#include "Eigen/Core"
#include <Eigen/Eigen>
#include <vector>
#include <map>

namespace ashlar{
// Divide a mesh into subsurface clusters using Variational Shape Approximation.
// David Cohen-Steiner, Pierre Alliez, and Mathieu Desbrun. Variational shape approximation. ACM Transactions on Graphics, 23(3):905-914, August 2004.
// With additional example information from the python script kmeans.py by Jesus Galvez.
//
struct VsaParams {
  int max_clusters = 25; //if our resolution is set super high, stop before we get here
  int subiterations = 10; //how many subiterations per loop
  double max_error = 0.26; //target weighted normal error per region
  bool verbose = false;
};

// Inputs:
//	 V	#V by 3 list of vertices
//   F  #F by 3 list of triangle indices
//	 k_clusters	Integer indicating number of clusters to find
//	 n_iterations	Integer indicating number of iterations to complete
// Outputs:
//   C  #k list of cluster component ids
void variational_shape_approximation(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const int k_clusters,
    const int n_iterations,
    Eigen::VectorXi &C,
    Eigen::MatrixXi &adjacency_mx);

void variational_shape_approximation(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const VsaParams &params,
    Eigen::VectorXi &C,
    Eigen::MatrixXi &adjacency_mx,
    Eigen::VectorXi &proxyFaces,
    Eigen::MatrixXd &regionCenters,
    Eigen::MatrixXd &regionNormals,
    std::vector<Eigen::MatrixXi> &regionFaces,
    std::vector<std::vector<int>> &regionFaceIndices
); //region faces is a list of all faces in a give region as a #region faces x 3 referencing V

struct face_error {
  double _error;
  int region_id;
  int face_id;
};

struct prox_region {
  std::vector<int>
      faces; //create a vector to hold our region faces, as a list of the face indices (of F) that live in this region
  Eigen::RowVector3d center; //hold the center point of each region, #K x 3
  Eigen::RowVector3d normal; //hold the normalised normal of each region, #K x 3
  Eigen::RowVector3d weightedNormal; //hold the area weighted normalof each region, #K x 3
  int proxyFace; //hold the index of the face that is the best seed of this region
  int id;
};

bool evaluate_mean_error(std::vector<prox_region> &proxies,
                         const Eigen::MatrixXd &WFN,
                         const Eigen::MatrixXd &FN,
                         const Eigen::VectorXd &FA,
                         const Eigen::MatrixXd &FC,
                         const double max_error);
void update_proxies(std::vector<prox_region> &proxies,
                    const Eigen::MatrixXd &WFN,
                    const Eigen::VectorXd &FA,
                    const Eigen::MatrixXd &FC);
void update_proxy(prox_region &proxy, const Eigen::MatrixXd &WFN, const Eigen::VectorXd &FA, const Eigen::MatrixXd &FC);
void get_proxy_seed(std::vector<prox_region> &proxies, const Eigen::MatrixXd &FN, const Eigen::VectorXd &FA);
void metric_error(const prox_region &PR,
                  std::vector<int> faces,
                  const Eigen::MatrixXd &FN,
                  const Eigen::VectorXd &FA,
                  Eigen::VectorXd &errors,
                  std::vector<face_error> &faceErrors);
void build_queue(const std::vector<prox_region> &regions,
                 const Eigen::MatrixXd &FN,
                 const Eigen::VectorXd &FA,
                 const Eigen::MatrixXi &TT,
                 std::vector<face_error> &queue,
                 std::vector<int> &assignedIndices);
void update_queue(const prox_region &region,
                  const Eigen::MatrixXd &FN,
                  const Eigen::VectorXd &FA,
                  std::vector<int> &seed_neighbors,
                  std::vector<face_error> &queue);
void update_queue_new(const prox_region &region,
                      const Eigen::MatrixXd &FN,
                      const Eigen::VectorXd &FA,
                      std::vector<int> &seed_neighbors,
                      std::vector<face_error> &queue);
void assign_regions(const Eigen::MatrixXd &FN,
                    const Eigen::VectorXd &FA,
                    const Eigen::MatrixXi &TT,
                    const std::vector<face_error> &queue,
                    std::vector<int> &assignedIndices,
                    std::vector<prox_region> &regions,
                    face_error &worst);
void split_regions(const Eigen::MatrixXd &FC,
                   const Eigen::MatrixXd &WFN,
                   const Eigen::MatrixXd &FN,
                   const Eigen::VectorXd &FA,
                   const Eigen::MatrixXi &TT,
                   std::vector<prox_region> &regions,
                   face_error &worst);
void assign_to_wost_regions(
    const Eigen::MatrixXd &FN,
    const Eigen::VectorXd &FA,
    const Eigen::MatrixXi &TT,
    const face_error &worst,
    const std::vector<face_error> &queue,
    std::vector<prox_region> &regions,
    std::vector<prox_region> &splitRegions,
    std::vector<int> &assignedIndices
);

void find_adjacent_regions(
    const Eigen::MatrixXi &TT,
    std::vector<prox_region> &regions,
    Eigen::MatrixXi &adjacency_matrix
);

void combine_regions(
    const Eigen::MatrixXd &FC,
    const Eigen::MatrixXd &WFN,
    const Eigen::MatrixXd &FN,
    const Eigen::VectorXd &FA,
    const Eigen::MatrixXi &TT,
    std::vector<prox_region> &regions,
    Eigen::MatrixXi &adjacency_matrix
);

//returns angles in degrees
void adjacency_angles(
    const Eigen::MatrixXi &adjacency_matrix,
    const Eigen::MatrixXd &regionNormals,
    std::map<int, std::map<int, double>> &interAngles
);

//given a point P, a plane origin O, and a plane normal N, return the point projected onto the plane, PP
void project_point_to_plane(
    const Eigen::RowVector3d &P,
    const Eigen::RowVector3d &O,
    const Eigen::RowVector3d &N,
    Eigen::RowVector3d &PP
);

//stable adjacency matrix is a num_regions x num_regions matrix where a given index is 1 if the centroid projects onto the edge
void find_anchors_and_edges(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &E,
    const Eigen::VectorXi &C,
    const Eigen::MatrixXd &regionCenters,
    const Eigen::MatrixXd &regionNormals,
    const Eigen::MatrixXi &adjacency_matrix,
    Eigen::MatrixXd &MP,
    Eigen::MatrixXi &MPA,
    std::map<int, std::map<int, Eigen::VectorXi>> &edgePaths,
    std::map<int, std::map<int, int>> &edgeCenterProxies,
    std::map<int, std::map<int, Eigen::RowVector3d>> &edgeCenters,
    std::map<int, std::map<int, Eigen::RowVector3d>> &frameXAxes,
    std::map<int, std::map<int, Eigen::RowVector3d>> &frameYAxes,
    std::map<int, eigen_map<int, Eigen::Quaterniond>> & edgeFrames,
Eigen::MatrixXi &stable_adjacency_matrix,
    std::map<int, std::map < int, double>> & pseudo_edge_lengths,
std::map<int, std::map<int, Eigen::RowVector3d>> &edge_dirs
);
}
#endif

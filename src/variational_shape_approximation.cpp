#include "ashlar/variational_shape_approximation.h"

#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/per_face_normals.h>
#include <igl/edges_to_path.h>
#include <igl/edge_lengths.h>
#include <igl/doublearea.h>
#include <igl/barycenter.h>
#include <igl/edges.h>

#include <algorithm>
#include <iostream>
//#include <vector>
#include <float.h>
#include <math.h>       /* acos */
#include <set>

//this version solves for a fixed number of k_clusters with a given number of n_iterations
void ashlar::variational_shape_approximation(const Eigen::MatrixXd &V,
                                          const Eigen::MatrixXi &F,
                                          const int k_clusters,
                                          const int n_iterations,
                                          Eigen::VectorXi &C,
                                          Eigen::MatrixXi &adjacency_mx) {

  //Get mesh face normals, areas, centers, and adjacencies
  Eigen::MatrixXd FN;
  igl::per_face_normals(V, F, FN); //compute per face normals, #F x 3
  Eigen::VectorXd FA;
  igl::doublearea(V, F, FA);
  FA /= 2.;  //compute our face areas #Fx1
  Eigen::MatrixXd FC;
  igl::barycenter(V, F, FC); //compute the global centers of each face #Fx3
  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);//hold our adjacent triangles and adjacent edges of these triangles
  Eigen::MatrixXd WFN;
  WFN = FN.array().colwise() * FA.array();  //weight our face normals by our face areas

  Eigen::VectorXi seedFaces;
  seedFaces = Eigen::VectorXi::LinSpaced(k_clusters, 0, F.rows() - 1); //set our seed faces as evenly spaced through our face indices (or try random)
  C.setZero(F.rows()); //empty our c vector

  //a proxy is a plane that best approximates a cluster of faces, we start with a proxy beign randomly selected faces from the mesh
  //a region is a grouping of faces/properties that best approximate a proxy
  std::vector<prox_region> proxies, regions;
  proxies.resize(k_clusters);
  regions.resize(k_clusters);

  //std::cout << "V: " << V.rows() << " F " << F.rows() << " FN " << FN.rows() << " F0 " << FN.row(0) << "\n";
  face_error worst = face_error();

  for (int i = 0; i < k_clusters; i++) {
    proxies.at(i).id = i;
    proxies.at(i).proxyFace = seedFaces(i);
    proxies.at(i).faces.push_back(seedFaces(i));
  }
  //for each iteration
  for (int n = 0; n < n_iterations; n++) {
    //updates the fit normals, center, weights for a proxy plane for a given region
    update_proxies(proxies, WFN, FA, FC);
    //update our region properties to match our proxy properties
    regions = proxies; //copy our proxies into our regions
    for (int i = 0; i < regions.size(); i++) {
      regions.at(i).faces.clear(); //clear our faces from this region
      regions.at(i).faces.push_back(proxies.at(i).proxyFace);
    }

    //build our queue
    std::vector<face_error> queue; //our error list for all faces in regions
    std::vector<int> assignedIndices; //indices that have already been assigned to a region
    build_queue(regions, FN, FA, TT, queue, assignedIndices);

    //grow regions from seeds according to the queue
    worst = face_error(); //store the error, region and face ID of our worst fit
    assign_regions(FN, FA, TT, queue, assignedIndices, regions, worst);
    //split worst regions
    split_regions(FC, WFN, FN, FA, TT, regions, worst);


    //find adjacent regions
    Eigen::MatrixXi adjacency_matrix;
    adjacency_matrix.setZero(regions.size(), regions.size());
    find_adjacent_regions(TT, regions, adjacency_matrix);
    //std::cout << "ADJACENCY MATRIX: \n" << adjacency_matrix << "\n";

    //find and combine the most similar
    combine_regions(FC, WFN, FN, FA, TT, regions, adjacency_matrix);
    //updates the fit normals, center, weights for a proxy plane for a given region
    update_proxies(proxies, WFN, FA, FC);
    get_proxy_seed(regions, FN, FA);  //update the proxy seed face for our regions
    //save for the next time around
    proxies = regions;
  }

  find_adjacent_regions(TT, regions, adjacency_mx); //get our final adjacency matrix
  Eigen::MatrixXi amt;
  amt.setZero(adjacency_mx.rows(), adjacency_mx.cols());
  amt = adjacency_mx.transpose();
  adjacency_mx = adjacency_mx + amt; //fill in the lower left corner to make a symmetric adjacency matrix

  for (int i = 0; i < regions.size(); i++) {
    //for each region, put our region index into the corresponding row of C
    for (int j = 0; j < regions.at(i).faces.size(); j++) {
      //for each face in this region, add the region index to C
      int f_idx = regions.at(i).faces.at(j);
      C(f_idx) = regions.at(i).id; //set this face to belong to this region
      //std::cout << "region: " << regions.at(i).id << "face : " << f_idx << "\n";
    }
  }
}

//this version takes an error metric and iteratively increases the size of K until that value is met, or the max number of faces/clusters/regions is reached
void ashlar::variational_shape_approximation(const Eigen::MatrixXd &V,
                                          const Eigen::MatrixXi &F,
                                          const VsaParams &params,
                                          Eigen::VectorXi &C,
                                          Eigen::MatrixXi &adjacency_mx,
                                          Eigen::VectorXi &proxyFaces,
                                          Eigen::MatrixXd &regionCenters,
                                          Eigen::MatrixXd &regionNormals,
                                          std::vector<Eigen::MatrixXi> &regionFaces,
                                          std::vector<std::vector<int>> &regionFaceIndices) {

  //Get mesh face normals, areas, centers, and adjacencies
  Eigen::MatrixXd FN;
  igl::per_face_normals(V, F, FN); //compute per face normals, #F x 3
  Eigen::VectorXd FA;
  igl::doublearea(V, F, FA);
  FA /= 2.;  //compute our face areas #Fx1
  Eigen::MatrixXd FC;
  igl::barycenter(V, F, FC); //compute the global centers of each face #Fx3
  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);//hold our adjacent triangles and adjacent edges of these triangles
  Eigen::MatrixXd WFN;
  WFN = FN.array().colwise() * FA.array();  //weight our face normals by our face areas

  int k_clusters = 4; //we start with this many clusters
  Eigen::VectorXi seedFaces;
  //set our seed faces as evenly spaced through our face indices (or try random)
  seedFaces = Eigen::VectorXi::LinSpaced(k_clusters, 0, F.rows() - 1);
  C.setZero(F.rows()); //empty our c vector

  //a proxy is a plane that best approximates a cluster of faces, we start with a proxy being randomly selected faces from the mesh
  //a region is a grouping of faces/properties that best approximate a proxy
  std::vector<prox_region> proxies, regions;

  face_error worst = face_error();

  //std::cout << "V: " << V.rows() << " F " << F.rows() << " FN " << FN.rows() << " F0 " << FN.row(0) << "\n";
  //std::cout << TT;

  bool withinErrorThreshold = false; //are we within the specified error threshold?

  while (!withinErrorThreshold && k_clusters < params.max_clusters) {  //while we haven't hit our error threshold, and we are below our max number of clusters

    proxies.clear();
    regions.clear(); //reset our proxies and our regions
    proxies.resize(k_clusters);
    regions.resize(k_clusters);  //resize our proxies and our regions
    for (int i = 0; i < k_clusters; i++) { //set our proxies to contain 1 seed face per proxy
      proxies.at(i).id = i;
      proxies.at(i).proxyFace = seedFaces(i);
      proxies.at(i).faces.push_back(seedFaces(i));
    }

    for (int n = 0; n < params.subiterations; n++) {
      //for each iteration
      update_proxies(proxies, WFN, FA, FC);  //updates the fit normals, center, weights for a proxy plane for a given region
      //update our region properties to match our proxy properties
      regions = proxies; //copy our proxies into our regions
      for (int i = 0; i < regions.size(); i++) {
        regions.at(i).faces.clear(); //clear our faces from this region
        regions.at(i).faces.push_back(proxies.at(i).proxyFace);
      }

      //build our queue
      std::vector<face_error> queue; //our error list for all faces in regions
      std::vector<int> assignedIndices; //indices that have already been assigned to a region
      build_queue(regions, FN, FA, TT, queue, assignedIndices);

      //grow regions from seeds according to the queue
      worst = face_error(); //store the error, region and face ID of our worst fit
      assign_regions(FN, FA, TT, queue, assignedIndices, regions, worst);
      //split worst regions
      split_regions(FC, WFN, FN, FA, TT, regions, worst);

      //find adjacent regions
      Eigen::MatrixXi adjacency_matrix;
      adjacency_matrix.setZero(regions.size(), regions.size());
      find_adjacent_regions(TT, regions, adjacency_matrix);
      //std::cout << "ADJACENCY MATRIX: \n" << adjacency_matrix << "\n";

      //find and combine the most similar
      combine_regions(FC, WFN, FN, FA, TT, regions, adjacency_matrix);
      update_proxies(proxies, WFN, FA, FC);  //updates the fit normals, center, weights for a proxy plane for a given region
      get_proxy_seed(regions, FN, FA);  //update the proxy seed face for our regions
      //save for the next time around
      proxies = regions;
    }

    withinErrorThreshold = evaluate_mean_error(proxies, WFN, FN, FA, FC, params.max_error); //returns true if within threshold, otherwise proxies is modified with an extra seed face from the worst region

    if (!withinErrorThreshold) {
      k_clusters++;  //increment the number of clusters we have
      seedFaces.setZero(k_clusters);
      for (int i = 0; i < seedFaces.rows(); i++) {
        seedFaces(i) = proxies.at(i).proxyFace;
      }
    }
  }

  //zero our region properties to the correct size
  proxyFaces.setZero(regions.size());
  regionCenters.setZero(regions.size(), 3);
  regionNormals.setZero(regions.size(), 3);

  //at the very end we find our adjacency matrix
  find_adjacent_regions(TT, regions, adjacency_mx); //get our final adjacency matrix
  Eigen::MatrixXi amt;
  amt.setZero(adjacency_mx.rows(), adjacency_mx.cols());
  amt = adjacency_mx.transpose();
  adjacency_mx = adjacency_mx + amt; //fill in the lower left corner to make a symmetric adjacency matrix
  regionFaces.clear(); //empty our matrix of region faces in case it was full
  regionFaceIndices.clear();
  for (int i = 0; i < regions.size(); i++) {
    //save our proxy properties as an output
    proxyFaces(i) = regions.at(i).proxyFace;
    regionCenters.row(i) = regions.at(i).center;
    regionNormals.row(i) = regions.at(i).normal;
    Eigen::MatrixXi rgnFaces;
    rgnFaces.setZero(regions.at(i).faces.size(), 3); //store the faces in a per region face list
    //for each region, put our region index into the corresponding row of C
    for (int j = 0; j < regions.at(i).faces.size(); j++) {
      //for each face in this region, add the region index to C
      int f_idx = regions.at(i).faces.at(j);
      rgnFaces.row(j) = F.row(f_idx); //save the values of this face in F to our F list for this region
      C(f_idx) = regions.at(i).id; //set this face to belong to this region
    }
    regionFaces.push_back(rgnFaces); //save this to our master list of region faces
    regionFaceIndices.push_back(regions.at(i).faces); //add our indexed face list as yet another thing to store
  }
}

//this is used in the case of recursive reprocessing where we need to calculate the error of each region to decide if we want to keep going with more faces
bool ashlar::evaluate_mean_error(std::vector<prox_region> &proxies,
                              const Eigen::MatrixXd &WFN,
                              const Eigen::MatrixXd &FN,
                              const Eigen::VectorXd &FA,
                              const Eigen::MatrixXd &FC,
                              const double max_error) {
  face_error worstError;
  worstError._error = 0.0;

  for (int i = 0; i < proxies.size(); i++) { //for each proxy
    Eigen::RowVector3d wfn;
    wfn.setZero(); //hold our weighted normal
    Eigen::RowVector3d pbc;
    pbc.setZero(); //hold our approximate centroid calculated by weighted face centers
    double proxyArea = 0;  //store the total area of this region
    double totalError = 0.0; //store the region error

    face_error worstRegionError; //keep track of our worst region for placing seeds next time
    worstRegionError._error = 0.0;

    for (int j = 0; j < proxies.at(i).faces.size(); j++) { //for each face in each proxy
      Eigen::RowVector3d normError;
      normError = FN.row(proxies.at(i).faces[j]) - proxies.at(i).normal;
      double errorLength = normError.squaredNorm();
      errorLength *= FA(proxies.at(i).faces[j]);
      wfn += WFN.row(proxies.at(i).faces[j]);
      pbc += FC.row(proxies.at(i).faces[j]) * FA(proxies.at(i).faces[j]);

      if (errorLength > worstRegionError._error) {  //save our record for the worst face of this region
        worstRegionError._error = errorLength;
        worstRegionError.face_id = proxies.at(i).faces[j];
        worstRegionError.region_id = i;
      }

      proxyArea += FA(proxies.at(i).faces[j]);
      totalError += errorLength;
    }

    //update our proxy properties while we're at it
    proxies.at(i).weightedNormal = wfn;
    proxies.at(i).normal = wfn.normalized();
    proxies.at(i).center = pbc / proxyArea; //normalize our weighted average;

    double proxyError = totalError / proxyArea; //get the normalized mean error

    if (proxyError > worstError._error) { //if this region has the worst proxy error so far, save it
      worstError._error = proxyError;
      worstError.face_id = worstRegionError.face_id;
      worstError.region_id = worstRegionError.region_id;
    }
  }

  if (worstError._error > max_error) {
    prox_region newSeedRegion; //create a new proxy region so we start our next iteration with one extra seed face that is better than randomly placed
    newSeedRegion.id = proxies.size();
    newSeedRegion.proxyFace = worstError.face_id;
    proxies.push_back(newSeedRegion); //add our split regions to our region
    return false; //we didn't meet our min requirements, so we need to add one to K and try again with our new seeds
  } else {
    return true;
  }
}

//Update the normal properties of proxy regions based on the weighted sum of all faces in each region (does not update seed face)
void ashlar::update_proxies(std::vector<prox_region> &proxies, const Eigen::MatrixXd &WFN, const Eigen::VectorXd &FA, const Eigen::MatrixXd &FC) {
  //update the center, normals, and weighted normals of our proxy given our index of faces, their weighted normals, areas and face centers
  //proxyNormals.setZero(proxyFaces.size(), 3); //empty our proxy normals
  //proxyWeightedNormals.setZero(proxyFaces.size(), 3);
  for (int i = 0; i < proxies.size(); i++) { //for each proxy
    Eigen::RowVector3d wfn;
    wfn.setZero(); //hold our weighted normal
    Eigen::RowVector3d pbc;
    pbc.setZero(); //hold our approximate centroid calculated by weighted face centers
    double proxyArea = 0;
    for (int j = 0; j < proxies.at(i).faces.size(); j++) { //for each face in each proxy
      wfn += WFN.row(proxies.at(i).faces[j]);
      pbc += FC.row(proxies.at(i).faces[j]) * FA(proxies.at(i).faces[j]);
      proxyArea += FA(proxies.at(i).faces[j]);
    }
    proxies.at(i).weightedNormal = wfn;
    proxies.at(i).normal = wfn.normalized();
    proxies.at(i).center = pbc / proxyArea; //normalize our weighted average;
  }
}

//Update the normal properties of a proxy region based on the weighted sum of all faces in the region (does not update seed face)
void ashlar::update_proxy(prox_region &proxy, const Eigen::MatrixXd &WFN, const Eigen::VectorXd &FA, const Eigen::MatrixXd &FC) {
  //update the center, normals, and weighted normals of our proxy given our index of faces, their weighted normals, areas and face centers
  //proxyNormals.setZero(proxyFaces.size(), 3); //empty our proxy normals
  //proxyWeightedNormals.setZero(proxyFaces.size(), 3);
  Eigen::RowVector3d wfn;
  wfn.setZero(); //hold our weighted normal
  Eigen::RowVector3d pbc;
  pbc.setZero(); //hold our approximate centroid calculated by weighted face centers
  double proxyArea = 0;
  for (int j = 0; j < proxy.faces.size(); j++) { //for each face in each proxy
    wfn += WFN.row(proxy.faces[j]);
    pbc += FC.row(proxy.faces[j]) * FA(proxy.faces[j]);
    proxyArea += FA(proxy.faces[j]);
  }
  proxy.weightedNormal = wfn;
  proxy.normal = wfn.normalized();
  proxy.center = pbc / proxyArea; //normalize our weighted average;
}

//update the seed face in each proxy region as the face with the least error when compared to the proxy normal
void ashlar::get_proxy_seed(std::vector<prox_region> &proxies, const Eigen::MatrixXd &FN, const Eigen::VectorXd &FA) {
  //iterate through each face listed in a proxy region to determine which one has the least metric error
  //this becomes our seed face from which the proxy will grow
  //regionProxyFaces.setZero(proxyFaces.size());
  //std::cout << "\nGET PROXY SEED: \n";
  for (int i = 0; i < proxies.size(); i++) { //for each proxy
    Eigen::VectorXd errors;
    std::vector<face_error> faceErrors;

    metric_error(proxies.at(i), proxies.at(i).faces, FN, FA, errors, faceErrors); //returns the errors of a proxy region as a list when compared to the proxy normal
    Eigen::VectorXd::Index seedProxyFaceIndex;
    double leastError = errors.minCoeff(&seedProxyFaceIndex); //get the least error and the index of our proxyfaces with the least error
    int seedFaceIndex = proxies.at(i).faces.at(seedProxyFaceIndex); //get the actual face index of our whole mesh from the proxy face with the least error
    proxies.at(i).proxyFace = seedFaceIndex;
    //std::cout << seedFaceIndex << "\n";
  }
}

//given a proxy region and a list of faces, return a list of all face errors between the faces and the normal of the proxy
void ashlar::metric_error(const prox_region &PR,
                       std::vector<int> faces,
                       const Eigen::MatrixXd &FN,
                       const Eigen::VectorXd &FA,
                       Eigen::VectorXd &errors,
                       std::vector<face_error> &faceErrors) {
  //return the metric error of all faces in a proxy
  errors.setZero(faces.size()); //empty our proxy normals
  faceErrors.clear(); //empty our face errors
  for (int i = 0; i < faces.size(); i++) { //for each proxy
    Eigen::RowVector3d normError;
    normError = FN.row(faces[i]) - PR.normal;
    double errorLength = normError.squaredNorm();
    errorLength *= FA(faces[i]);
    errors(i) = errorLength; //set our error row vector

    face_error structError; //hold our error as a struct;
    structError._error = errorLength;
    structError.face_id = faces[i];
    structError.region_id = PR.id;
    faceErrors.push_back(structError); //add our face error to this list
  }
}

//builds a new queue from scratch that contains all the neighbors of seed faces for all regions in one error list
void ashlar::build_queue(const std::vector<prox_region> &regions,
                      const Eigen::MatrixXd &FN,
                      const Eigen::VectorXd &FA,
                      const Eigen::MatrixXi &TT,
                      std::vector<face_error> &queue,
                      std::vector<int> &assignedIndices) {
  queue.clear();
  assignedIndices.clear();
  for (int i = 0; i < regions.size(); i++) {
    //for each region
    //int seedIndex =
    int seedIndex = regions.at(i).faces.at(0);
    //std::cout << "SEED INDEX: " << seedIndex << "\n";
    assignedIndices.push_back(seedIndex); //the seed faces are already accounted for
    std::vector<int> seed_neighbors;
    for (int j = 0; j < TT.row(seedIndex).size(); j++) { //convert eigen row to vector (maybe find an eigen-y way of doing this?)
      if (TT(seedIndex, j) != -1) {
        seed_neighbors.push_back(TT(seedIndex, j));//TODO: TT returns -1 if no neighbor for a given face, not sure if there will be consequences for disjointed pieces
      }
    }
    update_queue(regions.at(i), FN, FA, seed_neighbors, queue);
  }
}

//for a list of faces "seed neighbors", calculate the error from each face to the region proxy and insert these errors into a sorted queue.
void ashlar::update_queue(const prox_region &region, const Eigen::MatrixXd &FN, const Eigen::VectorXd &FA, std::vector<int> &seed_neighbors, std::vector<face_error> &queue) {
  //int regionIndex = region.id;
  Eigen::RowVector3d proxyNormal;
  proxyNormal = region.normal;
  Eigen::VectorXd newErrors;
  std::vector<face_error> newFaceErrors;
  metric_error(region, seed_neighbors, FN, FA, newErrors, newFaceErrors);

  queue.insert(queue.end(), newFaceErrors.begin(), newFaceErrors.end()); //add our new errors into our old ones
  sort(queue.begin(), queue.end(), [](const face_error &a, const face_error &b) { return (a._error < b._error); });
  //std::cout << "\nUPDATED QUEUE.  Min Error is: " << queue.front()._error << " max error is: " << queue.back()._error << "\n";
}

//adds the errors of each face in the seed_neighbors face list to the queue, as measured against the normal of the region.  Does not sort the queue.
void ashlar::update_queue_new(const prox_region &region, const Eigen::MatrixXd &FN, const Eigen::VectorXd &FA, std::vector<int> &seed_neighbors, std::vector<face_error> &queue) {
  int regionIndex = region.id;
  Eigen::RowVector3d proxyNormal;
  proxyNormal = region.normal;
  for (int i = 0; i < seed_neighbors.size(); i++) {
    int index = seed_neighbors.at(i);
    Eigen::RowVector3d normError;
    normError = FN.row(index) - proxyNormal;
    double errorLength = normError.squaredNorm();
    errorLength *= FA(index);

    face_error structError; //hold our error as a struct;
    structError._error = errorLength;
    structError.face_id = index;
    structError.region_id = regionIndex;
    queue.push_back(structError); //add our face error to this list
    //std::push_heap(queue.begin(), queue.end()); //not sure if necessary, but should readjust heap
  }
}

//assign all of our faces to proxy regions, and return the region that has the face with the worst error
void ashlar::assign_regions(const Eigen::MatrixXd &FN,
                         const Eigen::VectorXd &FA,
                         const Eigen::MatrixXi &TT,
                         const std::vector<face_error> &queue,
                         std::vector<int> &assignedIndices,
                         std::vector<prox_region> &regions,
                         face_error &worst) {
//V, F, FN, FA, TT, queue, assignedIndices, regions, worst) {
  std::vector<face_error> heapq = queue;
  std::vector<face_error> globalQueue;
  std::make_heap(heapq.begin(), heapq.end(), [](const face_error &a, const face_error &b) { return (a._error > b._error); });
  while (heapq.size() > 0) {
    face_error mostPriority = heapq.front(); //get the smallest val
    std::pop_heap(heapq.begin(), heapq.end(), [](const face_error &a, const face_error &b) { return (a._error > b._error); });
    heapq.pop_back(); //remove the smallest val
    int faceIndex = mostPriority.face_id;

    if (std::find(assignedIndices.begin(), assignedIndices.end(), faceIndex) == assignedIndices.end()) { //if faceIndex is NOT found in our list of assigned indices
      globalQueue.push_back(mostPriority); //add this face to our global queue at the end.
      int regionIndex = mostPriority.region_id;
      regions.at(regionIndex).faces.push_back(faceIndex); //add this face to this region
      assignedIndices.push_back(faceIndex); //make sure we remember we've been to this one
      //get the adjacent faces to the popped face and add them to the queue
      std::vector<int> face_neighbors;
      for (int j = 0; j < TT.row(faceIndex).size(); j++) { //convert eigen row to vector (maybe find an eigen-y way of doing this?)
        if (TT(faceIndex, j) != -1) {
          face_neighbors.push_back(TT(faceIndex, j));//TODO: check if -1 check is adequate
        }
      }
      //if any of these neighbors have already been assigned, remove them
      std::vector<int> unique_neighbors;
      std::sort(face_neighbors.begin(), face_neighbors.end());
      std::sort(assignedIndices.begin(), assignedIndices.end());
      std::set_difference(face_neighbors.begin(), face_neighbors.end(), assignedIndices.begin(), assignedIndices.end(), back_inserter(unique_neighbors));
      //update our queue
      update_queue_new(regions.at(regionIndex), FN, FA, unique_neighbors, heapq);
      std::push_heap(heapq.begin(), heapq.end(), [](const face_error &a, const face_error &b) { return (a._error > b._error); });  //was a<b
    }
  }
  sort(globalQueue.begin(), globalQueue.end(), [](const face_error &a, const face_error &b) {
    return (a._error < b._error);
  }); //reverse comparator to get the worst value
  worst = globalQueue.back();
}

void ashlar::split_regions(const Eigen::MatrixXd &FC,
                        const Eigen::MatrixXd &WFN,
                        const Eigen::MatrixXd &FN,
                        const Eigen::VectorXd &FA,
                        const Eigen::MatrixXi &TT,
                        std::vector<prox_region> &regions,
                        face_error &worst) {
  //here we take the worst region and we split it into two, taking the worst face of a region and using it as a new seed
  prox_region *worstRegion = &regions.at(worst.region_id);

  //make our split regions and set their properties
  prox_region splitRegionA, splitRegionB;
  splitRegionA.id = regions.size();
  splitRegionB.id = regions.size() + 1;
  splitRegionA.proxyFace = worstRegion->proxyFace;
  splitRegionA.faces.push_back(worstRegion->proxyFace);
  splitRegionB.proxyFace = worst.face_id;
  splitRegionB.faces.push_back(worst.face_id);
  std::vector<prox_region> splitRegions = {splitRegionA, splitRegionB};
  update_proxies(splitRegions, WFN, FA, FC);

  //build our queue
  std::vector<face_error> queue; //our error list for all faces in regions
  std::vector<int> assignedIndices; //indices that have already been assigned to a region
  ashlar::build_queue(splitRegions,
                   FN,
                   FA,
                   TT,
                   queue,
                   assignedIndices); //create a sorted list of the errors associated with the seed faces for each region (does not guarantee seed faces are in region)

  //assign to regions
  assign_to_wost_regions(FN, FA, TT, worst, queue, regions, splitRegions, assignedIndices);

  //we also need to update the region IDs of our regions, because we did some funny adding and deleting
  for (int i = 0; i < regions.size(); i++) {
    regions.at(i).id = i; //reset our ID to be sequential!
  }

}

//take the worst region and replace it with two regions that are subcomponents of it
void ashlar::assign_to_wost_regions(const Eigen::MatrixXd &FN,
                                 const Eigen::VectorXd &FA,
                                 const Eigen::MatrixXi &TT,
                                 const face_error &worst,
                                 const std::vector<face_error> &queue,
                                 std::vector<prox_region> &regions,
                                 std::vector<prox_region> &splitRegions,
                                 std::vector<int> &assignedIndices) {

  //remove any elements in queue that are not also in regionDomain
  std::vector<face_error> heapq = queue;
  std::vector<int> regionDomain = regions.at(worst.region_id).faces;
  for (int i = heapq.size() - 1; i >= 0;
       i--) { //iterate back through queue, removing elements that are not in regionDomain (we might have some initial neighbors to our seed that are not in this region)
    if (std::find(regionDomain.begin(), regionDomain.end(), heapq.at(i).face_id) == regionDomain.end()) { //if the element is not in the assigned indices
      heapq.erase(heapq.begin() + i); //erase this element in the queue
    }
  }

  //heapify
  std::sort(regionDomain.begin(), regionDomain.end()); //all the indices of all faces in our worst region
  std::make_heap(heapq.begin(), heapq.end(), [](const face_error &a, const face_error &b) { return (a._error > b._error); });//used to be a<b

  while (heapq.size() > 0) {
    face_error mostPriority = heapq.front(); //get the smallest val
    std::pop_heap(heapq.begin(), heapq.end(), [](const face_error &a, const face_error &b) {
      return (a._error > b._error);
    }); //need to pop heap in order to pop back (remove) smallest element
    heapq.pop_back(); //remove the smallest val
    int faceIndex = mostPriority.face_id;
    if (std::find(assignedIndices.begin(), assignedIndices.end(), faceIndex) == assignedIndices.end()) { //if faceIndex is NOT found in our list of assigned indices
      int regionIndex = mostPriority.region_id; //this should be regions.size or regions.size+1
      for (int i = 0; i < splitRegions.size(); i++) { //our split regions are indexed as 0 and 1, but we need to get the priority regionID which is in our global region ID numbers
        if (regionIndex == splitRegions.at(i).id) {
          regionIndex = i;
          break;
        }
      }
      if (regionIndex > 1) {
        std::cout << "================SPLIT REGION INDEX TOO LARGE!  DEBUG! DEBUG!=================\n"; //just here in case there is a big hole in my logic...
      }

      splitRegions.at(regionIndex).faces.push_back(faceIndex); //add our face to this region
      assignedIndices.push_back(faceIndex); //make sure we remember we've been to this one
      std::vector<int> face_neighbors; //get the face neighbors
      for (int j = 0; j < TT.row(faceIndex).size(); j++) { //convert eigen row to vector (maybe find an eigen-y way of doing this?)
        if (TT(faceIndex, j) != -1) {
          face_neighbors.push_back(TT(faceIndex, j));//TODO:  check if -1 check is adequate for naked edges
        }
      }
      //convert face neighbors, assignedIndices, and  to a set
      std::sort(face_neighbors.begin(), face_neighbors.end());
      std::vector<int> intersected;
      //keep only the faces that are also in our master region that we've split
      std::set_intersection(face_neighbors.begin(), face_neighbors.end(), regionDomain.begin(), regionDomain.end(), back_inserter(intersected));
      //remove faces that have already been assigned
      std::sort(intersected.begin(), intersected.end());
      std::sort(assignedIndices.begin(), assignedIndices.end());
      std::vector<int> unique_neighbors;
      std::set_difference(intersected.begin(), intersected.end(), assignedIndices.begin(), assignedIndices.end(), back_inserter(unique_neighbors));
      if (unique_neighbors.size() > 0) {//if there are any left
        //update our queue with these neighbors
        update_queue_new(splitRegions.at(regionIndex), FN, FA, unique_neighbors, heapq);
        //update_queue(splitRegions.at(regionIndex), FN, FA, unique_neighbors, heapq);
        std::push_heap(heapq.begin(), heapq.end(), [](const face_error &a, const face_error &b) { return (a._error > b._error); });
      }
    }
  }
  regions.push_back(splitRegions.at(0)); //add our split regions to our region
  regions.push_back(splitRegions.at(1));
  regions.erase(regions.begin() + worst.region_id); //erase our old bad region
}

//construct an RxR adjacency matrix for region adjacency, where only the top right of the matrix is filled in (self-self is 0, but otherwise adjacencies in top right corner are represented as 1)
void ashlar::find_adjacent_regions(const Eigen::MatrixXi &TT, std::vector<prox_region> &regions, Eigen::MatrixXi &adjacency_matrix) {
  //construct a vector<int> for each region of all adjacent faces
  std::vector<std::vector<int>> region_face_neighbors; //store a list of lists where each sublist is a 3*nFaces long list of ints representing face indices of adjacent values
  for (int i = 0; i < regions.size(); i++) {
    std::vector<int> face_neighbors; //store all adjacencies for this region as a big list
    for (int j = 0; j < regions.at(i).faces.size(); j++) {
      //for each face in this region, get the adjacent faces
      int f_idx = regions.at(i).faces.at(j);
      for (int k = 0; k < TT.row(f_idx).size(); k++) { //for each adjacency
        face_neighbors.push_back(TT(f_idx, k));//add the value of this neighbor
      }
    }
    std::sort(face_neighbors.begin(), face_neighbors.end()); //sort this list
    region_face_neighbors.push_back(face_neighbors); //add the face neighbors of this region
  }

  //now we need to loop through the adjacency list for each region and see if it fits any face indices of any other regions
  //we store this adjacency list as a numRegions * numRegions square matrix where all matches are stored in the upper/right half of the diagonal (which is all zeros because we don't check for self matching)
  adjacency_matrix.setZero(regions.size(), regions.size());
  for (int i = 0; i < regions.size() - 1; i++) { //for each region except the last one (because we would have already checked for all matches with it)
    //we need to check against all other regions EXCEPT the current one, and any ones we have checked before
    std::vector<int> regionF = regions.at(i).faces; //this is the faces of our region
    std::sort(regionF.begin(), regionF.end());  //sort all of our faces (necessary for set intersection)
    for (int j = i + 1; j < regions.size(); j++) {
      std::vector<int> intersected;
      //keep only the faces that are also in our master region that we've split
      std::set_intersection(regionF.begin(), regionF.end(), region_face_neighbors.at(j).begin(), region_face_neighbors.at(j).end(), back_inserter(intersected));
      if (intersected.size() > 0) {
        //if we had any intersections, put a 1 in our intersection matrix
        adjacency_matrix(i, j) = 1;
      }
    }
  }
}

//combine the two adjacent regions that have the least mean error when combined
void ashlar::combine_regions(const Eigen::MatrixXd &FC,
                          const Eigen::MatrixXd &WFN,
                          const Eigen::MatrixXd &FN,
                          const Eigen::VectorXd &FA,
                          const Eigen::MatrixXi &TT,
                          std::vector<prox_region> &regions,
                          Eigen::MatrixXi &adjacency_matrix) {
  int bestI = 0; //store the best index of our corresponding regions to merge
  int bestJ = 0;
  double minMaxError = 0; //set this to the largest possible number
  int mergedCount = 0; //store how many regions we have tried to merge
  prox_region minMerged;
  for (int i = 0; i < regions.size(); i++) {
    for (int j = 0; j < regions.size(); j++) {
      int isAdjacency = adjacency_matrix(i, j); //if the adjacency matrix is 1, there is an intersection
      if (isAdjacency == 1) {
        prox_region *regionA = &regions.at(i);
        prox_region *regionB = &regions.at(j);
        prox_region merged;
        merged.faces.insert(merged.faces.begin(), begin(regionA->faces), end(regionA->faces)); //insert the faces from region A into our merged face region
        merged.faces.insert(end(merged.faces), begin(regionB->faces), end(regionB->faces)); //insert the faces from region B into our merged face region
        update_proxy(merged, WFN, FA, FC); //update our proxy properties for this combined region

        double regionError = 0; //store the total error of all faces in this region when compared to the region proxy
        for (int k = 0; k < merged.faces.size(); k++) {//for each face in this combined region, get the error
          Eigen::RowVector3d normError;
          normError = FN.row(merged.faces.at(k)) - merged.normal;
          double errorLength = normError.squaredNorm();
          errorLength *= FA(merged.faces.at(k)); //weight for face area
          regionError += errorLength;

          if (mergedCount > 0 && regionError > minMaxError) { //if we've alread gone beyond our max error before even finishing the loop, we can jump out of it
            break;
          }
        }

        if (mergedCount == 0 || regionError <= minMaxError) {
          minMaxError = regionError; //store our best result!
          bestI = i;
          bestJ = j;
          minMerged = merged;
        }
        mergedCount++; //add one to our counter to let us know we've gone beyond the first one

      } //end if adjacency
    } //for adjacency cols
  } //for adjacency rows

  //now that we have a best result, we need to edit our proxies accordingly by removing two and adding one.
  regions.erase(regions.begin() + bestJ); //because of the way our adjacency matrix is structured, J is always bigger than I, so we can safely delete J first
  regions.erase(regions.begin() + bestI); //because of the way our adjacency matrix is structured, J is always bigger than I, so we can safely delete J first
  regions.push_back(minMerged); //add our best region

  //we also need to update the region IDs of our regions, because we did some funny adding and deleting
  for (int i = 0; i < regions.size(); i++) {
    regions.at(i).id = i; //reset our ID to be sequential!
  }
}

//compute the angle between any two adjacent regions in our matrix, and return the angles as a map.  returns dot product result, to convert to degrees need to use acos / arctan
void ashlar::adjacency_angles(const Eigen::MatrixXi &adjacency_matrix, const Eigen::MatrixXd &regionNormals, std::map<int, std::map<int, double>> &interAngles) {

  for (int i = 0; i < adjacency_matrix.rows(); i++) { //for each region
    for (int j = 0; j < adjacency_matrix.cols(); j++) { //for each region
      if (adjacency_matrix(i, j) == 1 && i < j) { //enforce that this is an adjacent edge, and that we only consider the upper half of our adjacency mx to not count twice
        //double interregionalAngle = acos(regionNormals.row(i).dot(regionNormals.row(j)))*(180.0*M_PI);
        Eigen::RowVector3d v_a = regionNormals.row(i).normalized();
        Eigen::RowVector3d v_b = regionNormals.row(j).normalized();
        double interregionalAngle = acos(v_a.dot(v_b)) * (180.0 / M_PI);
        //ToDo:  Check if angle is acute or obtuse
        interregionalAngle = (90.0 - interregionalAngle / 2.0) * 2; //convert to INNER ANGLE
        interAngles[i][j] = interregionalAngle; //add this angle to our map for this face relationship
      }
    }
  }
}

void ashlar::project_point_to_plane(const Eigen::RowVector3d &P, const Eigen::RowVector3d &O, const Eigen::RowVector3d &N, Eigen::RowVector3d &PP) {
  //https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d?noredirect=1&lq=1
  //https://www.sciencedirect.com/topics/computer-science/orthographic-projection
  PP = P - (P - O).dot(N) * N;
}

//find all anchor vertices (verts that have adjacent faces connected to three or more proxy regions) and "edges" between proxy anchors
void ashlar::find_anchors_and_edges(const Eigen::MatrixXd &V,
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
                                 std::map<int, eigen_map<int, Eigen::Quaterniond>> &edgeFrames,
Eigen::MatrixXi &stable_adjacency_matrix,
    std::map<int, std::map<int, double>> &pseudo_edge_lengths,
std::map<int, std::map<int, Eigen::RowVector3d>> &edge_dirs) {

// set up a matrix that is the same size as V, but we will only fill in where our anchors are..actually, for safety we just use V
Eigen::MatrixXd VAdjustedAnchors(V);
std::vector<std::vector<int>> VF; //#V list of lists of incident faces(adjacency list)
std::vector<std::vector<int>> VFi; //#V list of lists of index of incidence within incident faces listed
//calculate the adjacent faces to each vertex
igl::vertex_triangle_adjacency(V, F, VF, VFi);
std::vector<std::vector<int>> VR; //#V list of lists of incident regions to each vert
//0,1,2 if the vertex is within a region, on the edge of two, or an anchor
Eigen::VectorXi V_type;
V_type.setZero(V.rows()); //set our type to zero
for (int i = 0; i < VF.size(); i++) { //for each vertex
std::vector<int> adjacentRegionsToVert;
for (int j = 0; j < VF.at(i).size(); j++) { //for each adjacent face
int adjacentFaceRegionIndex = C(VF.at(i).at(j));
adjacentRegionsToVert.push_back(adjacentFaceRegionIndex);
}
VR.push_back(adjacentRegionsToVert); //add this list of adjacent regions to vert to the region

//get the number of unique adjacent regions to this vert
//https://stackoverflow.com/questions/29321127/number-of-unique-elements-in-a-vector
std::set<double> arvSet(adjacentRegionsToVert.begin(), adjacentRegionsToVert.end()); // get a set version of our adjacent regions to remove duplicates
int numAdjRegions = arvSet.size(); //use the size of this vector converted to a set to give us the number of unique elements
if (numAdjRegions > 2) {
V_type(i) = 2; //this vertex is an anchor (3 or more neighboring regions)

//if it's an anchor, get the proxy approximation of this vert as the average of this anchor projected onto each connected proxy plane
Eigen::RowVector3d anchorSum = Eigen::RowVector3d::Zero();
for (int k = 0; k < arvSet.size(); k++) {
Eigen::RowVector3d projectedAnchor = Eigen::RowVector3d::Zero();
//https://thispointer.com/how-to-access-element-by-index-in-a-set-c/
int currRegionIdx = *(std::next(arvSet.begin(), k)); //need to get the current region index using K from our arvset
project_point_to_plane(V.row(i), regionCenters.row(currRegionIdx), regionNormals.row(currRegionIdx), projectedAnchor);
anchorSum += projectedAnchor; //add our current point to the sum of all projected points
}
//divide sum of projected anchors by number of anchors to get the mean projected anchor
anchorSum /= (double) arvSet.size();
VAdjustedAnchors.row(i) = anchorSum; //save this new calculated version of V into of normal verts and projected anchors.
} else if (numAdjRegions > 1) {
V_type(i) = 1; //this vertex is an edge (2 adjacent regions)
} else {
V_type(i) = 0; //this vertex is entirely within a region
}
}

//now we need to determine edges that lie along intersection lines
//each edge of the mesh can either divide two regions, or it lies entirely within a region.
//we need to take our list of (completely unsorted edges) and put the index of each edge into separate lists corresponding to our region adjacency matrix
//i.e. if an edge is adjacent to region 1 and region 2, it goes in the edge list at(1).at(2)


std::map < int, std::map < int, std::vector < int>>>
    edgeMap; //map of edges per region combination
//should hold in the format (region1,2Edge0 = region_edges.at(1).at(2).at(0)

//for each edge, we need to find the adjacent faces and their region types
//to my knowledge, there isn't a libigl function for this, so we can find the adjacent faces to the verts at each end of the line
//and then get the two common faces (or possibly 1 in edge/hole cases)
for (int i = 0; i < E.rows(); i++) {
//for each edge
if (V_type(E(i, 0)) > 0 && V_type(E(i, 1)) > 0) { //if both verts are either on an edge or are anchors, we have an edge
//we need to determine which face indices are part of this edge, which we will determine as the two face indices that are common to both edge verts

std::vector<int> VFE0(VF.at(E(i, 0))); //adjacent faces to the first vert of this edge, copied from our VF list element
std::vector<int> VFE1(VF.at(E(i, 1))); //adjacent faces to the second vert of this edge, copied from our VF list element

//sort adjacent faces to each vert
std::sort(VFE0.begin(), VFE0.end());
std::sort(VFE1.begin(), VFE1.end());

//make a new vector to store faces common to both
std::vector<int> commonFaces;
//keep only the faces that are also in our master region that we've split

std::set_intersection(VFE0.begin(), VFE0.end(), VFE1.begin(), VFE1.end(), back_inserter(commonFaces));
if (commonFaces.size() == 2) {
//if we have two common faces, we know what to do
//if we only have one common face, this is a half-edge or outer edge of an open mesh, a case that we do not deal with now
//if we have 0 or more than 2 common faces, something is wrong ;)
std::vector<int> commonRegions; //store the regions common to these faces
for (int j = 0; j < commonFaces.size(); j++) {
commonRegions.push_back(C(commonFaces.at(j)));
}
std::sort(commonRegions.begin(), commonRegions.end()); //sort our adjacent regions in order
//we now need to push this edge into the region edge map indicated by commonRegions(0),commonRegions(1)
edgeMap[commonRegions.at(0)][commonRegions.at(1)].push_back(i); //put this edge in the correct bin

}

}
}

//our edge map should now have correct edges stored in correct bins, but unordered.  Now we need to go through those bins and chain the edges together in order
//this can be tricky if there are multiple edges per region adjacency, but if we don't consider this case, there is a libigl function for this!
//should modify in the future to check if there are more than 2 unique values in the list, in which case we would need to make two or more edge lists per region adjaceny
std::vector<Eigen::RowVector3d> mps;
std::vector<Eigen::RowVector2i> mpAdjacencies; //this stores the actual adjacencies that we've found midpoints for, which might be a subset of adjacencies if we have a crenellated face, etc.
stable_adjacency_matrix.setZero(adjacency_matrix.rows(), adjacency_matrix.cols());  //set our stable matrix to zero (should be similar to adjacency matrix but with some missing 1s

//now we need to loop through the adjacency list for each region and see if it fits any face indices of any other regions
for (int i = 0; i < adjacency_matrix.rows(); i++) { //for each region except the last one (because we would have already checked for all matches with it)
for (int j = 0; j < adjacency_matrix.cols(); j++) {
if (adjacency_matrix(i, j) == 1 && i < j) { //enforce that this is an adjacent edge, and that we only consider the upper half of our adjacency mx to not count twice
//if key doesn't exist
if (edgeMap.find(i) == edgeMap.end() || edgeMap.at(i).find(j) == edgeMap.at(i).end()) {
std::cout << "Key doesn't exist at " << i << " , " << j << "\n";
} else if (edgeMap.at(i).at(j).size() > 0) {
//if there are edges in this adjacent region
//copy our data into an eigen matrix
Eigen::MatrixX2i unsortedTempEdges;
unsortedTempEdges.setZero(edgeMap.at(i).at(j).size(), 2); //create an empty 2xnumEdgesInRegion matrix, storing the index of V as copies of rows of E
for (int k = 0; k < unsortedTempEdges.rows(); k++) { //get the endpoint indices from each edge and put them into our temp edge matrix
unsortedTempEdges(k, 0) = E(edgeMap.at(i).at(j).at(k), 0);
unsortedTempEdges(k, 1) = E(edgeMap.at(i).at(j).at(k), 1);
}
//get the length of all of our triangle edges for this master edge
Eigen::VectorXd UEL; // unsorted edge lengths
igl::edge_lengths(V, unsortedTempEdges, UEL); //uel stores the length of each edge in unsortedTempEdges
double HEL = UEL.sum() / 2.0; //half the length of this total edge
//chain and sort our edges.  Warning, will break if there are disjointed edges!
Eigen::VectorXi vpath, epath, eend;
igl::edges_to_path(unsortedTempEdges, vpath, epath, eend); //todo deal with disjoint edges case!
edgePaths[i][j] = vpath; // add our polyline points to our map as indices in V
double LUM = 0; //length until midpoint
for (int k = 0; k < epath.rows(); k++) { //for each reordered index in unsortedTempEdges edge in our path
LUM += UEL(epath(k)); //add the length of this edge (the length that corresponds to the same index K as an index in unsortedTempEdges
if (LUM >= HEL) { //if we have passed our midlength
//get how far we have gone past it
double overshoot = LUM - HEL;
//get the percentage of the current edge that is overshot
double osFrac = overshoot / UEL(epath(k));
int firstIndex = eend(k); //determine which way this edge is pointing
int secondIndex = (firstIndex + 1) % 2;

Eigen::RowVector3d firstVert = V.row(unsortedTempEdges(epath(k), firstIndex));
Eigen::RowVector3d secondVert = V.row(unsortedTempEdges(epath(k), secondIndex));
Eigen::RowVector3d backTrack = (firstVert - secondVert) * osFrac;
Eigen::RowVector3d midpoint = secondVert + backTrack; // get our midpoint location as a 3D coord
mps.push_back(midpoint);

edgeCenters[i][j] = midpoint; //save our midpoint to this map
edgeCenters[j][i] = midpoint; //save our midpoint to this map

if (osFrac <= 0.5) {
//if we are less than halfway overshot, the nearest vert is the end of this triangle edge
edgeCenterProxies[i][j] = unsortedTempEdges(epath(k), secondIndex); //save the closest actual vert index to this map
} else {
edgeCenterProxies[i][j] = unsortedTempEdges(epath(k), firstIndex);
}

//if we get to here, we have found the midpoint of this region region edge.  Let's now also get its mean direction
//as the vector between its approximated anchor points
Eigen::RowVector3d projectedEdgeDirection, end_pt_a, end_pt_b;
end_pt_a = VAdjustedAnchors.row(vpath(0));
if (vpath(vpath.rows() - 1) != vpath(0)) { // if our edge is not a closed curve, we treat its direction as the vector between end points
projectedEdgeDirection = VAdjustedAnchors.row(vpath(vpath.rows() - 1)) - VAdjustedAnchors.row(vpath(0)); //t
end_pt_b = VAdjustedAnchors.row(vpath(vpath.rows() - 1));
} else {
//if this is a closed edge, we'll just use the first triangle edge direction as our direction vector
projectedEdgeDirection = VAdjustedAnchors.row(vpath(1)) - VAdjustedAnchors.row(vpath(0)); //t
end_pt_b = VAdjustedAnchors.row(vpath(1));
}

projectedEdgeDirection.normalize();

//calculate the properties (x and y axis, Z is existing region normal) of the reference frame for each region region edge
//this is done twice for each midpoint, so we have a frame where region A points out, and another where region B points out (and A is down)

Eigen::RowVector3d iCenter = regionCenters.row(i); // get the centerpoints of the regions that are adjacent to this edge
Eigen::RowVector3d jCenter = regionCenters.row(j);

Eigen::RowVector3d iNormal = regionNormals.row(i); // get the normals of the regions that are adjacent to this edge
Eigen::RowVector3d jNormal = regionNormals.row(j);

Eigen::RowVector3d iCenterDir = (iCenter - midpoint).normalized(); // get the vector from the edge midpoint to the face center
Eigen::RowVector3d jCenterDir = (jCenter - midpoint).normalized();

Eigen::RowVector3d Yi = iNormal.cross(projectedEdgeDirection); //get our local frame y axis as the cross product of our edge and face normal
Eigen::RowVector3d Yj = jNormal.cross(projectedEdgeDirection);

if (Yi.dot(iCenterDir) < 0)
Yi *= -1; // if our y axis is pointing the "wrong" way, reverse it
if (Yj.dot(jCenterDir) < 0)
Yj *= -1;

Eigen::RowVector3d Xi = Yi.cross(iNormal);
Eigen::RowVector3d Xj = Yj.cross(jNormal);

frameXAxes[i][j] = Xi;  //save our local reference frames to our map
frameXAxes[j][i] = Xj;

frameYAxes[i][j] = Yi;
frameYAxes[j][i] = Yj;

//convert transformation frame to quaternion
Eigen::Matrix3d frameMxI, frameMxJ;
frameMxI.col(0) = Xi.transpose();
frameMxI.col(1) = Yi.transpose();
frameMxI.col(2) = iNormal.transpose();

frameMxJ.col(0) = Xj.transpose();
frameMxJ.col(1) = Yj.transpose();
frameMxJ.col(2) = jNormal.transpose();

edgeFrames[i][j] = Eigen::Quaternion<double>(frameMxI);
edgeFrames[j][i] = Eigen::Quaternion<double>(frameMxJ);

mpAdjacencies.push_back(Eigen::RowVector2i(i, j)); //add this adjacency to our list so we know it has a valid midpoint for future work

//for purposes of stone flipping, we want to see if the centroid of the stone falls within the edge if projected onto it
//0<= (c - a).dot(b-a) <= (b-a).dot(b-a)  //in our case our centroid (c) is at zero, so this becomes:
Eigen::RowVector3d b_a = end_pt_b - end_pt_a;

edge_dirs[i][j] = projectedEdgeDirection;  //store our unitized edge direction
//edge_dirs[j][i] = projectedEdgeDirection;  //store our unitized edge direction

pseudo_edge_lengths[i][j] = b_a.norm(); //store the length of the pseudo edge (not the mesh edge length)
//pseudo_edge_lengths[j][i] = pseudo_edge_lengths[i][j]; //store the length of the pseudo edge (not the mesh edge length)

double ca_ba = (-1.0 * end_pt_a).dot(b_a);
if (ca_ba > 0 && ca_ba < b_a.dot(b_a)) {
//this is a "stable edge"
stable_adjacency_matrix(i, j) = 1;
stable_adjacency_matrix(j, i) = 1;
}

break; // leave this edge and go to the next adjacency
}
}
}
}
}
}
std::cout << "Number of Edge Midpoints Found: " << mps.size() << " \n";
//std::cout << "Face Adjacency Matrix:\n";
//std::cout << adjacency_matrix << "\n";

MP.setZero(mps.size(), 3); //convert our list of midpoints to a flattened eigen matrix
for (int i = 0; i < MP.rows(); i++) {
MP.row(i) = mps.at(i);
}

MPA.setZero(mpAdjacencies.size(), 2);
//convert our list of midpoint adjacencies to a flattened eigen matrix (same length as nonduplicate adjacency matrix or smaller)
for (int i = 0; i < MPA.rows(); i++) {
MPA.row(i) = mpAdjacencies.at(i);
}

}

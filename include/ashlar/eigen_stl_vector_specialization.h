#pragma once

/**
 * Define specializations for std::vector<...> for fixed-size vectorizable Eigen-type.s
 * Needed to avoid unnecessary SEGFAULTS.
 * This header needs to precede
 * @see https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html
 * Special thanks to Tom Lankhorst for this fix
 */

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Geometry>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Quaterniond)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::RowVector2i)

// alias for maps that use Eigen types
// e.g.: `eigen_map<Eigen::Quaterniond>`
#include <map>
template<typename Key, typename Value>
using eigen_map = std::map<Key, Value, std::less<Key>,
                           Eigen::aligned_allocator<std::pair<const Key, Value>>>;


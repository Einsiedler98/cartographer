/*
 * Copyright 2016 The Cartographer Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CARTOGRAPHER_MAPPING_3D_TSDFS_H_
#define CARTOGRAPHER_MAPPING_3D_TSDFS_H_

#include <memory>
#include <string>
#include <vector>

#include <open_chisel/Chisel.h>

#include "Eigen/Geometry"
#include "cartographer/common/port.h"
#include "cartographer/mapping/proto/submap_visualization.pb.h"
#include "cartographer/mapping/submaps.h"
#include "cartographer/mapping_2d/probability_grid.h"
#include "cartographer/mapping_2d/range_data_inserter.h"
#include "cartographer/mapping_3d/hybrid_grid.h"
#include "cartographer/mapping_3d/proto/submaps_options.pb.h"
#include "cartographer/mapping_3d/range_data_inserter.h"
#include "cartographer/sensor/range_data.h"
#include "cartographer/transform/transform.h"

namespace cartographer {
namespace mapping_3d {

struct TSDF : public mapping::Submap {
  TSDF(float high_resolution, float low_resolution,
         const Eigen::Vector3f& origin, int begin_range_data_index);

  chisel::ChiselPtr tsdf;
  bool finished = false;
  std::vector<int> trajectory_node_indices;
};


// A container for Truncated Signed Distance Fields similar to the Submaps container.
class TSDFs : public mapping::Submaps {
 public:
  TSDFs();
  TSDFs(const proto::SubmapsOptions& options);

  TSDFs(const TSDFs&) = delete;
  TSDFs& operator=(const TSDFs&) = delete;

  const TSDF* Get(int index) const override;
  const chisel::ProjectionIntegrator* GetIntegrator(int index) const;
  int size() const override;

    // Returns the indices of the Submap into which point clouds will
  // be inserted.
  std::vector<int> insertion_indices() const;
  //Inserts 'range_data' into the Submap collection.
   void InsertRangeData(const sensor::RangeData& range_data_in_tracking,
                        const transform::Rigid3d& pose_observation);

  /*void SubmapToProto(
      int index, const std::vector<mapping::TrajectoryNode>& trajectory_nodes,
      const transform::Rigid3d& global_submap_pose,
      mapping::proto::SubmapQuery::Response* response) const override;

  //

  // Returns the 'high_resolution' HybridGrid to be used for matching.
  const HybridGrid& high_resolution_matching_grid() const;

  // Returns the 'low_resolution' HybridGrid to be used for matching.
  const HybridGrid& low_resolution_matching_grid() const;*/

  // Adds a node to be used when visualizing the submap.
  void AddTrajectoryNodeIndex(int trajectory_node_index);

 private:

  void AddTSDF(const Eigen::Vector3f& origin);

  const proto::SubmapsOptions options_;

  std::vector<std::unique_ptr<TSDF>> submaps_;
  std::vector<chisel::ProjectionIntegrator> projection_integrators_;

  // Number of RangeData inserted.
  int num_range_data_ = 0;

  // Number of RangeData inserted since the last Submap was added.
  int num_range_data_in_last_submap_ = 0;
};

}  // namespace mapping_3d
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_3D_TSDFS_H_

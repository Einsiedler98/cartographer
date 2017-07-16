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

#ifndef CARTOGRAPHER_MAPPING_3D_KALMAN_TSDF_LOCAL_TRAJECTORY_BUILDER_H_
#define CARTOGRAPHER_MAPPING_3D_KALMAN_TSDF_LOCAL_TRAJECTORY_BUILDER_H_

#include <memory>

#include "cartographer/common/time.h"
#include "cartographer/kalman_filter/pose_tracker.h"
#include "cartographer/mapping/global_trajectory_builder_interface.h"
#include "cartographer/mapping_3d/local_tsdf_trajectory_builder_interface.h"
#include "cartographer/mapping_3d/motion_filter.h"
#include "cartographer/mapping_3d/proto/local_trajectory_builder_options.pb.h"
#include "cartographer/mapping_3d/scan_matching/ceres_tsdf_scan_matcher.h"
#include "cartographer/mapping_3d/tsdfs.h"
#include "cartographer/sensor/range_data.h"
#include "cartographer/sensor/voxel_filter.h"

namespace cartographer {
namespace mapping_3d {


// Wires up the local SLAM stack (i.e. UKF, scan matching, etc.) without loop
// closure.
class KalmanTSDFLocalTrajectoryBuilder : public LocalTSDFTrajectoryBuilderInterface {
 public:
  explicit KalmanTSDFLocalTrajectoryBuilder(
      const proto::LocalTrajectoryBuilderOptions& options);
  ~KalmanTSDFLocalTrajectoryBuilder() override;

  KalmanTSDFLocalTrajectoryBuilder(const KalmanTSDFLocalTrajectoryBuilder&) = delete;
  KalmanTSDFLocalTrajectoryBuilder& operator=(const KalmanTSDFLocalTrajectoryBuilder&) =
      delete;

  void AddImuData(common::Time time, const Eigen::Vector3d& linear_acceleration,
                  const Eigen::Vector3d& angular_velocity) override;
  std::unique_ptr<InsertionResult> AddRangefinderData(
      common::Time time, const Eigen::Vector3f& origin,
      const sensor::PointCloud& ranges) override;
  void AddOdometerData(common::Time time,
                       const transform::Rigid3d& pose) override;
  const mapping_3d::TSDFs* submaps() const override;
  const PoseEstimate& pose_estimate() const override;

 private:
  std::unique_ptr<InsertionResult> AddAccumulatedRangeData(
       common::Time time, const sensor::RangeData& range_data_in_tracking,
       const Eigen::Vector3f& sensor_origin);

  std::unique_ptr<InsertionResult> InsertIntoSubmap(
      const common::Time time, const sensor::RangeData& range_data_in_tracking,
      const transform::Rigid3d& pose_observation,
      const kalman_filter::PoseCovariance& covariance_estimate,
      const Eigen::Vector3f& sensor_origin);

  const proto::LocalTrajectoryBuilderOptions options_;
  std::unique_ptr<mapping_3d::TSDFs> submaps_;

  PoseEstimate last_pose_estimate_;

  // Pose of the last computed scan match.
  transform::Rigid3d scan_matcher_pose_estimate_;

  MotionFilter motion_filter_;
  std::unique_ptr<scan_matching::CeresTSDFScanMatcher> ceres_scan_matcher_;

  std::unique_ptr<kalman_filter::PoseTracker> pose_tracker_;  
  int num_accumulated_;
  int num_optimizations_;
  int num_optimization_iterations_;
  int num_map_updates_;
  float summed_duration_optimization;
  float summed_duration_map_update;
  double summed_residuals;
  int summed_solver_iterations;
  transform::Rigid3f first_pose_prediction_;
  sensor::RangeData accumulated_range_data_;
  std::vector<CombinedRangeData> accumulated_range_data_with_pose_;

};

}  // namespace mapping_3d
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_3D_KALMAN_TSDF_LOCAL_TRAJECTORY_BUILDER_H_

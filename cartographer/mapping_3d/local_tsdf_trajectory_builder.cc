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

#include "cartographer/mapping_3d/local_tsdf_trajectory_builder.h"

#include "cartographer/common/make_unique.h"
#include "cartographer/mapping_3d/kalman_tsdf_local_trajectory_builder.h"
#include "cartographer/mapping_3d/optimizing_tsdf_local_trajectory_builder.h"
#include "cartographer/mapping_3d/continuously_optimizing_tsdf_local_trajectory_builder.h"

namespace cartographer {
namespace mapping_3d {

std::unique_ptr<LocalTSDFTrajectoryBuilderInterface> CreateLocalTSDFTrajectoryBuilder(
    const proto::LocalTrajectoryBuilderOptions&
        local_trajectory_builder_options) {
  switch (local_trajectory_builder_options.use()) {
    case proto::LocalTrajectoryBuilderOptions::KALMAN:
      LOG(INFO) << "Initializing KalmanTSDFLocalTrajectoryBuilder";
      return common::make_unique<KalmanTSDFLocalTrajectoryBuilder>(
          local_trajectory_builder_options);
    case proto::LocalTrajectoryBuilderOptions::OPTIMIZING:
      LOG(INFO) << "Initializing OptimizingTSDFLocalTrajectoryBuilder";
      return common::make_unique<OptimizingTSDFLocalTrajectoryBuilder>(
          local_trajectory_builder_options);
    case proto::LocalTrajectoryBuilderOptions::CONTINUOUSLY_OPTIMIZING:
      LOG(INFO) << "Initializing OptimizingTSDFLocalTrajectoryBuilder";
      return common::make_unique<ContinuouslyOptimizingTSDFLocalTrajectoryBuilder>(
          local_trajectory_builder_options);
  }
  LOG(FATAL);
}

}  // namespace mapping_3d
}  // namespace cartographer

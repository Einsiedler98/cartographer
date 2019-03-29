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

#ifndef CARTOGRAPHER_SENSOR_INTERNAL_RANDOM_FILTER_H_
#define CARTOGRAPHER_SENSOR_INTERNAL_RANDOM_FILTER_H_

#include <bitset>
#include <unordered_set>

#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/sensor/point_cloud.h"
#include "cartographer/sensor/proto/adaptive_voxel_filter_options.pb.h"
#include "cartographer/sensor/timed_point_cloud_data.h"
#include "cartographer/sensor/internal/scan_matching_filter.h"
#include <iostream>

namespace cartographer {
namespace sensor {

// Voxel filter for point clouds. For each voxel, the assembled point cloud
// contains the first point that fell into it from any of the inserted point
// clouds.
class RandomFilter : public ScanMatchingFilter{
 public:
  // 'size' is the length of a voxel edge.
  explicit RandomFilter(float min_num_points) {

    number_of_points_ = min_num_points;
  }

  RandomFilter(const RandomFilter&) = delete;
  RandomFilter& operator=(const RandomFilter&) = delete;

  // Returns a voxel filtered copy of 'point_cloud'.
  PointCloud Filter(const PointCloud& point_cloud);

  // Same for TimedPointCloud.
  TimedPointCloud Filter(const TimedPointCloud& timed_point_cloud);

  // Same for RangeMeasurement.
  std::vector<TimedPointCloudOriginData::RangeMeasurement> Filter(
      const std::vector<TimedPointCloudOriginData::RangeMeasurement>&
          range_measurements);

 private:
  float number_of_points_;
};

}  // namespace sensor
}  // namespace cartographer

#endif  // CARTOGRAPHER_SENSOR_INTERNAL_RANDOM_FILTER_H_

/*
 * Copyright 2018 The Cartographer Authors
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

#ifndef CARTOGRAPHER_MAPPING_2D_EDF_2D_H_
#define CARTOGRAPHER_MAPPING_2D_EDF_2D_H_

#include <vector>

#include "cartographer/common/port.h"
#include "cartographer/mapping/2d/grid_2d.h"
#include "cartographer/mapping/2d/tsdf_2d.h"
#include "cartographer/mapping/2d/map_limits.h"
#include "cartographer/mapping/2d/xy_index.h"
#include "cartographer/mapping/internal/tsd_value_converter.h"

#include <queue>
#include "absl/memory/memory.h"
namespace cartographer {
namespace mapping {

// Represents a 2D grid of truncated signed distances and weights.
class EDF2D {
 public:
  EDF2D(const TSDF2D& tsdf, float max_distance)
      : max_distance_(max_distance),
        resolution_(tsdf.limits().resolution()),
        limits_(tsdf.limits().cell_limits()),
        cells_(limits_.num_x_cells * limits_.num_y_cells, 255) {
    std::queue<Eigen::Array2i> update_queue;

    // Find seeds
    const cartographer::mapping::MapLimits& limits = tsdf.limits();
    int num_x_cells = limits.cell_limits().num_y_cells;
    int num_y_cells = limits.cell_limits().num_x_cells;
    for (int ix = 0; ix < num_x_cells; ++ix) {
      for (int iy = 0; iy < num_y_cells; ++iy) {
        if (std::abs(tsdf.GetTSD({iy, ix})) <= tsdf.limits().resolution()) {
          update_queue.push({iy, ix});
          SetCell({iy, ix}, ComputeCellValue(std::abs(tsdf.GetTSD({iy, ix}))));
          //          LOG(INFO)<<GetValue({iy, ix});
          //          LOG(INFO)<<GetED({iy,
          //          ix})<<"\t"<<std::abs(tsdf.GetTSD({iy,
          //          ix}))<<"\t"<<max_distance/255.f<<"\t"<<max_distance;
        }
      }
    }

    while (!update_queue.empty()) {
      Eigen::Array2i cell_idx = update_queue.front();
      update_queue.pop();
      float center_tsd = GetED(cell_idx);
      for (int ix = -1; ix < 2; ix++) {
        for (int iy = -1; iy < 2; iy++) {
          if (ix == 0 && iy == 0) continue;
          Eigen::Array2i candidate = {cell_idx[0] + iy, cell_idx[1] + ix};
          if (!Contains(candidate)) continue;
          float candidate_tsd = GetED(candidate);
          float d = std::sqrt(std::abs(ix) + std::abs(iy)) * resolution_;
          if (std::abs(candidate_tsd) > center_tsd + d) {
            if (SetCell(candidate, ComputeCellValue(center_tsd + d)))
              update_queue.push(candidate);
          }
        }
      }
    }
  };

  float GetED(const Eigen::Array2i& cell_index) const {
    return ComputeDistance(GetValue(cell_index));
  };

  int GetValue(const Eigen::Array2i& xy_index) const {
    const Eigen::Array2i local_xy_index = xy_index;
    // The static_cast<unsigned> is for performance to check with 2 comparisons
    // xy_index.x() < offset_.x() || xy_index.y() < offset_.y() ||
    // local_xy_index.x() >= wide_limits_.num_x_cells ||
    // local_xy_index.y() >= wide_limits_.num_y_cells
    // instead of using 4 comparisons.
    if (static_cast<unsigned>(local_xy_index.x()) >=
            static_cast<unsigned>(limits_.num_x_cells) ||
        static_cast<unsigned>(local_xy_index.y()) >=
            static_cast<unsigned>(limits_.num_y_cells)) {
      return 255;
    }
    const int stride = limits_.num_x_cells;
    return cells_[local_xy_index.x() + local_xy_index.y() * stride];
  }

  bool SetCell(const Eigen::Array2i& xy_index, uint8 value) {
    const Eigen::Array2i local_xy_index = xy_index;
    // The static_cast<unsigned> is for performance to check with 2 comparisons
    // xy_index.x() < offset_.x() || xy_index.y() < offset_.y() ||
    // local_xy_index.x() >= wide_limits_.num_x_cells ||
    // local_xy_index.y() >= wide_limits_.num_y_cells
    // instead of using 4 comparisons.
    if (static_cast<unsigned>(local_xy_index.x()) >=
            static_cast<unsigned>(limits_.num_x_cells) ||
        static_cast<unsigned>(local_xy_index.y()) >=
            static_cast<unsigned>(limits_.num_y_cells)) {
      return false;
    }
    const int stride = limits_.num_x_cells;
    cells_[local_xy_index.x() + local_xy_index.y() * stride] = value;
    return true;
  }

  uint8 ComputeCellValue(const float distance) const {
    const int cell_value = common::RoundToInt(distance * 255.f / max_distance_);
    CHECK_GE(cell_value, 0);
    CHECK_LE(cell_value, 255);
    return cell_value;
  }

  float ComputeDistance(const uint8 value) const {
    const float distance = float(value) * max_distance_ / 255.f;
    CHECK_GE(distance, 0.);
    CHECK_LE(distance, max_distance_);
    return distance;
  }

  bool Contains(const Eigen::Array2i& xy_index) const {
    return !(static_cast<unsigned>(xy_index.x()) >=
                 static_cast<unsigned>(limits_.num_x_cells) ||
             static_cast<unsigned>(xy_index.y()) >=
                 static_cast<unsigned>(limits_.num_y_cells));
  };
  float max_distance_;
  float resolution_;

 private:
  const CellLimits limits_;
  std::vector<uint8> cells_;
};

EDF2D CreateEDFFromTSDF(float truncation_distance,
                          const TSDF2D& tsdf);

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_2D_EDF_2D_H_

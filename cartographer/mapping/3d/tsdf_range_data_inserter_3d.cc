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

#include "cartographer/mapping/3d/tsdf_range_data_inserter_3d.h"
#include "cairo/cairo.h"
#include "cartographer/mapping/internal/3d/scan_matching/interpolated_tsdf.h"

#include "Eigen/Core"
#include "cartographer/mapping/probability_values.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {
namespace {
}  // namespace

proto::TSDFRangeDataInserterOptions3D CreateTSDFRangeDataInserterOptions3D(
    common::LuaParameterDictionary* parameter_dictionary) {
  proto::TSDFRangeDataInserterOptions3D options;
  options.set_truncation_distance(
      parameter_dictionary->GetDouble("truncation_distance"));
  options.set_maximum_weight(parameter_dictionary->GetDouble("maximum_weight"));
  options.set_project_sdf_distance_to_scan_normal(
      parameter_dictionary->GetBool("project_sdf_distance_to_scan_normal"));
  options.set_num_free_space_voxels(
      parameter_dictionary->GetInt("num_free_space_voxels"));
  return options;
}

TSDFRangeDataInserter3D::TSDFRangeDataInserter3D(
    const proto::RangeDataInserterOptions3D& options)
    : options_(options) {}

void TSDFRangeDataInserter3D::InsertHit(const Eigen::Vector3f& hit,
                                        const Eigen::Vector3f& origin,
                                        HybridGridTSDF* tsdf) const {
  const Eigen::Vector3f ray = hit - origin;
  const float range = ray.norm();
  const float truncation_distance =
      options_.tsdf_range_data_inserter_options_3d().truncation_distance() *
      tsdf->resolution() /
      0.05f;  // todo(kdaun) clean up truncation distance scaling
  //      static_cast<float>(options_.truncation_distance());
  //  LOG(INFO)<<"tau "<<truncation_distance;
  if (range < truncation_distance) return;
  const float truncation_ratio = truncation_distance / range;
  bool update_free_space =
      options_.tsdf_range_data_inserter_options_3d().num_free_space_voxels() >
      0;  // todo(kdaun) use num free space
          // cells value instead of bool
  const Eigen::Vector3f ray_begin =
      update_free_space ? origin : origin + (1.0f - truncation_ratio) * ray;
  const Eigen::Vector3f ray_end = origin + (1.0f + truncation_ratio) * ray;

  bool use_default_raycast = false;

  if (use_default_raycast) {
    const Eigen::Array3i begin_cell = tsdf->GetCellIndex(ray_begin);
    const Eigen::Array3i end_cell = tsdf->GetCellIndex(ray_end);
    const Eigen::Array3i delta = end_cell - begin_cell;
    const int num_samples = delta.cwiseAbs().maxCoeff();
    CHECK_LT(num_samples, 1 << 15);
    // 'num_samples' is the number of samples we equi-distantly place on the
    // line between 'origin' and 'hit'. (including a fractional part for sub-
    // voxels) It is chosen so that between two samples we change from one voxel
    // to the next on the fastest changing dimension.

    for (int position = 0; position < num_samples; ++position) {
      const Eigen::Array3i update_cell_index =
          begin_cell + delta * position / num_samples;
      //    if (tsdf->CellIsUpdated(update_cell)) continue;
      Eigen::Vector3f cell_center = tsdf->GetCenterOfCell(update_cell_index);
      float distance_cell_to_origin = (cell_center - origin).norm();
      float update_tsd = range - distance_cell_to_origin;
      update_tsd =
          common::Clamp(update_tsd, -truncation_distance, truncation_distance);
      float update_weight = 1.0;
      UpdateCell(update_cell_index, update_tsd, update_weight, tsdf);
    }
  } else {
    // Based on Amanatides, John, and Andrew Woo. "A fast voxel traversal
    // algorithm for ray tracing." Eurographics. Vol. 87. No. 3. 1987.
    Eigen::Array3i update_cell = tsdf->GetCellIndex(ray_begin);
    const Eigen::Array3i end_cell = tsdf->GetCellIndex(ray_end);
    const Eigen::Vector3f ray_begin_to_end = ray_end - ray_begin;
    const Eigen::Vector3f begin_cell_center =
        tsdf->GetCenterOfCell(update_cell);
    int step_x = update_cell[0] <= end_cell[0] ? 1 : -1;
    int step_y = update_cell[1] <= end_cell[1] ? 1 : -1;
    int step_z = update_cell[2] <= end_cell[2] ? 1 : -1;

    const Eigen::Vector3f step_direction({step_x, step_y, step_z});
    Eigen::Vector3f t_max = begin_cell_center - ray_begin +
                            0.5f * float(tsdf->resolution()) *
                                step_direction.cwiseQuotient(ray_begin_to_end);
    const Eigen::Vector3f t_delta =
        float(tsdf->resolution()) *
        (ray_begin_to_end.cwiseInverse()).cwiseAbs();
    while (t_max[0] < 1.0 || t_max[1] < 1.0 || t_max[2] < 1.0) {
      const Eigen::Array3i update_cell_index(
          {update_cell[0], update_cell[1], update_cell[2]});
      Eigen::Vector3f cell_center = tsdf->GetCenterOfCell(update_cell_index);
      float distance_cell_to_origin = (cell_center - origin).norm();
      float update_tsd = range - distance_cell_to_origin;
      update_tsd =
          common::Clamp(update_tsd, -truncation_distance, truncation_distance);
      float update_weight = 1.0;
      UpdateCell(update_cell_index, update_tsd, update_weight, tsdf);

      if (t_max[0] < t_max[1]) {
        if (t_max[0] < t_max[2]) {
          update_cell[0] = update_cell[0] + step_x;
          t_max[0] = t_max[0] + t_delta[0];
        } else {
          update_cell[2] = update_cell[2] + step_z;
          t_max[2] = t_max[2] + t_delta[2];
        }
      } else {
        if (t_max[1] < t_max[2]) {
          update_cell[1] = update_cell[1] + step_y;
          t_max[1] = t_max[1] + t_delta[1];
        } else {
          update_cell[2] = update_cell[2] + step_z;
          t_max[2] = t_max[2] + t_delta[2];
        }
      }
    }
  }
}

void TSDFRangeDataInserter3D::Insert(const sensor::RangeData& range_data,
                                              GridInterface* grid) const {
  CHECK(grid != nullptr);
  CHECK(grid->GetGridType() == GridType::TSDF);
  HybridGridTSDF* tsdf = static_cast<HybridGridTSDF*>(grid);

  const Eigen::Vector3f origin = range_data.origin.head<3>();
  for (const sensor::RangefinderPoint& hit_point : range_data.returns) {
    const Eigen::Vector3f hit = hit_point.position.head<3>();
    InsertHit(hit, origin, tsdf);
  }
  tsdf->FinishUpdate();
}

void TSDFRangeDataInserter3D::UpdateCell(const Eigen::Array3i& cell,
                                         float update_sdf, float update_weight,
                                         HybridGridTSDF* tsdf) const {
  if (update_weight == 0.f) return;

  const float old_weight = tsdf->GetWeight(cell);
  const float old_sdf = tsdf->GetTSD(cell);
  float updated_weight = old_weight + update_weight;
  float updated_sdf =
      (old_sdf * old_weight + update_sdf * update_weight) / updated_weight;
  float maximum_weight = static_cast<float>(
      options_.tsdf_range_data_inserter_options_3d().maximum_weight());
  updated_weight = std::min(updated_weight, maximum_weight);
  tsdf->SetCell(cell, updated_sdf, updated_weight);
}

}  // namespace mapping
}  // namespace cartographer

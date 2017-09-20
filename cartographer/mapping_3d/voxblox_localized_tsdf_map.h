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

#ifndef CARTOGRAPHER_MAPPING_3D_VOXBLOX_LOCALIZED_TSDF_MAP_H_
#define CARTOGRAPHER_MAPPING_3D_VOXBLOX_LOCALIZED_TSDF_MAP_H_

#include "cartographer/transform/rigid_transform.h"

#include <voxblox/core/tsdf_map.h>

namespace cartographer {
namespace mapping_3d {

class LocalizedTsdfMap : public voxblox::TsdfMap {
public:
  explicit LocalizedTsdfMap(const voxblox::TsdfMap::Config& config, const transform::Rigid3d& origin = transform::Rigid3d::Identity())
    : voxblox::TsdfMap(config),
      origin_(origin){
    block_size_ = config.tsdf_voxel_size * config.tsdf_voxels_per_side;
  }
  transform::Rigid3d getOrigin() {
    return origin_;
  }
protected:
  const transform::Rigid3d origin_;
};


}  // namespace mapping_3d
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_3D_VOXBLOX_LOCALIZED_TSDF_MAP_H_

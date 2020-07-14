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

#include "cartographer/io/ply_writing_points_processor.h"

#include <iomanip>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include "absl/memory/memory.h"
#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/io/points_batch.h"
#include "glog/logging.h"
#include "Open3D/Open3D.h"

namespace cartographer {
  namespace io {

    namespace {

// Writes the PLY header claiming 'num_points' will follow it into
// 'output_file'.
      void WriteBinaryPlyHeader(const bool has_color, const bool has_intensities,
                                const std::vector<std::string>& comments,
                                const int64 num_points,
                                FileWriter* const file_writer) {
        const std::string color_header = !has_color ? ""
                                                    : "property uchar red\n"
                                                      "property uchar green\n"
                                                      "property uchar blue\n";
        const std::string intensity_header =
                !has_intensities ? "" : "property float intensity\n";
        std::ostringstream stream;
        stream << "ply\n"
               << "format binary_little_endian 1.0\n"
               << "comment generated by Cartographer\n";
        for (const std::string& comment : comments) {
          stream << "comment " << comment << "\n";
        }
        stream << "element vertex " << std::setw(15) << std::setfill('0')
               << num_points << "\n"
               << "property float x\n"
               << "property float y\n"
               << "property float z\n"
               << color_header << intensity_header << "end_header\n";
        const std::string out = stream.str();
        CHECK(file_writer->WriteHeader(out.data(), out.size()));
      }

      void WriteBinaryPlyPointCoordinate(const Eigen::Vector3f& point,
                                         std::shared_ptr<open3d::geometry::PointCloud> cloud) {
        // TODO(sirver): This ignores endianness.
        cloud->points_.push_back({point[0], point[1], point[2]});
      }

      void WriteBinaryIntensity(const float intensity,
                                FileWriter* const file_writer) {
        // TODO(sirver): This ignores endianness.
//  CHECK(file_writer->Write(reinterpret_cast<const char*>(&intensity),
//                           sizeof(float)));
      }

      void WriteBinaryPlyPointColor(const Uint8Color& color,
                                    FileWriter* const file_writer) {
//  CHECK(file_writer->Write(reinterpret_cast<const char*>(color.data()),
//                           color.size()));
      }

    }  // namespace

    std::unique_ptr<PlyWritingPointsProcessor>
    PlyWritingPointsProcessor::FromDictionary(
            const FileWriterFactory& file_writer_factory,
            common::LuaParameterDictionary* const dictionary,
            PointsProcessor* const next) {
      return absl::make_unique<PlyWritingPointsProcessor>(
              file_writer_factory(dictionary->GetString("filename")),
              dictionary->GetInt("aggregate"),
              dictionary->HasKey("poisson_depth") ? dictionary->GetInt("poisson_depth") : 0,
              dictionary->HasKey("trim_surface") ? dictionary->GetDouble("trim_surface") : 0,
              dictionary->HasKey("statistical_outlier_neighbours") ? dictionary->GetInt("statistical_outlier_neighbours") : 0,
              dictionary->HasKey("statistical_outlier_radius") ? dictionary->GetDouble("statistical_outlier_radius") : 0,
              std::vector<std::string>(), next);
    }

    PlyWritingPointsProcessor::PlyWritingPointsProcessor(std::unique_ptr<FileWriter> file_writer,
                                                         const size_t &aggregate,
                                                         const int64 &poisson_depth,
                                                         const double &trim_surface,
                                                         const int64 &statistical_outlier_neighbours,
                                                         const double &statistical_outlier_radius,
                                                         const std::vector<std::string> &comments,
                                                         PointsProcessor *const next)
            : next_(next),
              aggregate_(aggregate),
              poisson_depth_(poisson_depth),
              trim_surface_(trim_surface),
              statistical_outlier_neighbours_(statistical_outlier_neighbours),
              statistical_outlier_radius_(statistical_outlier_radius),
              comments_(comments),
              num_points_(0),
              currentTime_(common::FromUniversal(0)),
              has_colors_(false),
              file_(std::move(file_writer)) {
      name_ = file_->GetFilename();
      pc_ = std::make_shared<open3d::geometry::PointCloud>();
      resultpc_ = std::make_shared<open3d::geometry::PointCloud>();
    }

    PointsProcessor::FlushResult PlyWritingPointsProcessor::Flush() {
      if(statistical_outlier_neighbours_ != 0 && statistical_outlier_radius_ != 0) {
        LOG(INFO) << "Removing statistical outliers using: " << std::to_string(statistical_outlier_neighbours_) << " " << std::to_string(statistical_outlier_radius_);
        std::vector<size_t> outliers;
        std::tie(resultpc_, outliers) = resultpc_->RemoveStatisticalOutliers(statistical_outlier_neighbours_, statistical_outlier_radius_);
      }

      if(poisson_depth_ == 0) {
        LOG(INFO) << "Writing point cloud to file: " + name_;
        open3d::io::WritePointCloudToPLY(name_, *resultpc_);
      } else {
        LOG(INFO) << "Calculating mesh using poisson reconstruction with depth: " + std::to_string(poisson_depth_);
        std::shared_ptr<open3d::geometry::TriangleMesh> mesh_es;
        std::vector<double> densities_es;
        std::tie(mesh_es, densities_es) = open3d::geometry::TriangleMesh::CreateFromPointCloudPoisson(*resultpc_, poisson_depth_);

        std::vector<bool> density_mask(densities_es.size(), false);
        double max = 0;
        for (int i = 0 ; i < densities_es.size(); i++) {
          if(densities_es[i] > max) max = densities_es[i];
          if(densities_es[i] < trim_surface_) density_mask[i] = true;
        }
        std::vector<Eigen::Vector3d> densities;
        for (std::vector<double>::iterator it = densities_es.begin() ; it != densities_es.end(); ++it) {
          double val = *it/max;
          densities.push_back({ val, val, val });
        }
        mesh_es->vertex_colors_ = densities;
        if(trim_surface_ > 0) {
          LOG(INFO) << "Trimming Mesh below density: " + std::to_string(trim_surface_);
          mesh_es->RemoveVerticesByMask(density_mask);
        }
        mesh_es->RemoveDegenerateTriangles();
        mesh_es->RemoveDuplicatedTriangles();
        mesh_es->RemoveDuplicatedVertices();
        mesh_es->RemoveNonManifoldEdges();

        open3d::io::WriteTriangleMesh(name_, *mesh_es);
      }

      switch (next_->Flush()) {
        case FlushResult::kFinished:
          return FlushResult::kFinished;

        case FlushResult::kRestartStream:
          LOG(FATAL) << "PLY generation must be configured to occur after any "
                        "stages that require multiple passes.";
      }
      LOG(FATAL);
    }

    void PlyWritingPointsProcessor::Process(std::unique_ptr<PointsBatch> batch) {
      if (batch->points.empty()) {
        next_->Process(std::move(batch));
        return;
      }

//  std::cout << std::to_string(aggregation_counter_) << std::endl;
      if(aggregation_counter_ == 0) {
        num_points_ = 0;
//    std::string path = name_.substr(0, name_.find(".ply")) + std::to_string(common::ToUniversal(batch->start_time));
//    std::string plyPath = path + ".ply";
//
//    file_->UpdateFileName(path + ".pose");
//    WriteBinaryPlyPointCoordinate(batch->origin, file_.get());
//    std::cout << plyPath << std::endl;
//    file_->UpdateFileName(plyPath);
//    WriteBinaryPlyHeader(has_colors_, has_intensities_, comments_, 0,
//                         file_.get());
      }

      if (has_colors_) {
        CHECK_EQ(batch->points.size(), batch->colors.size())
                << "First PointsBatch had colors, but encountered one without. "
                   "frame_id: "
                << batch->frame_id;
      }
      if (has_intensities_) {
        CHECK_EQ(batch->points.size(), batch->intensities.size())
                << "First PointsBatch had intensities, but encountered one without. "
                   "frame_id: "
                << batch->frame_id;
      }

      for (size_t i = 0; i < batch->points.size(); ++i) {
        pc_->points_.push_back({
                                       batch->points[i].position[0],
                                       batch->points[i].position[1],
                                       batch->points[i].position[2]});
        ++num_points_;
      }
      ++aggregation_counter_;
      if(aggregation_counter_ >= aggregate_) {
        aggregation_counter_ = 0;
        std::string path = name_.substr(0, name_.find(".ply")) + std::to_string(common::ToUniversal(batch->start_time));

        size_t position = path.find_last_of("/");
        std::string folderPath = path.substr(0, position) + "/results";
        std::string filePath = path.substr(position, path.length());
        struct stat buffer;



        std::string plyPath = folderPath + filePath + ".ply";
        pc_->EstimateNormals(open3d::geometry::KDTreeSearchParamHybrid(0.5, 30));
        pc_->OrientNormalsTowardsCameraLocation(batch->origin.cast<double>());
        resultpc_->operator+=(*pc_);

//        if(stat(folderPath.c_str(), &buffer) != 0) {
//          mkdir(folderPath.c_str(), 0755);
//        }
//        open3d::io::WritePointCloudToPLY(plyPath, *pc_);
        pc_ = std::make_shared<open3d::geometry::PointCloud>();
      }
      next_->Process(std::move(batch));
    }

  }  // namespace io
}  // namespace cartographer

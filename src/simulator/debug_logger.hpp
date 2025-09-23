#pragma once

#include <fstream>
#include <string>
#include <glm/glm.hpp>

class DebugLogger {
public:
    static DebugLogger& instance() {
        static DebugLogger logger;
        return logger;
    }

    void initialize(const std::string& filename) {
#ifdef _DEBUG
        if (file_.is_open()) {
            file_.close();
        }
        file_.open(filename);
        if (file_.is_open()) {
            file_ << "PhotonID,Event,PosX,PosY,PosZ,DirX,DirY,DirZ,Weight,VoxelX,VoxelY,VoxelZ,MediumID,Energy,Description\n";
        }
#endif
    }

    void log_photon_event(int photon_id, const std::string& event, 
                         const glm::dvec3& position, const glm::dvec3& direction, 
                         double weight, const glm::ivec3& voxel_coords = glm::ivec3(-1),
                         int medium_id = -1, double energy = 0.0,
                         const std::string& description = "") {
#ifdef _DEBUG
        if (file_.is_open()) {
            file_ << photon_id << "," << event << ","
                  << position.x << "," << position.y << "," << position.z << ","
                  << direction.x << "," << direction.y << "," << direction.z << ","
                  << weight << ","
                  << voxel_coords.x << "," << voxel_coords.y << "," << voxel_coords.z << ","
                  << medium_id << "," << energy << "," << description << "\n";
            file_.flush();
        }
#endif
    }

    void log_voxel_emittance(int photon_id, const glm::dvec3& position, const glm::dvec3& direction,
                           double weight, const glm::ivec3& voxel_coords, double emittance, 
                           const std::string& surface_type = "") {
#ifdef _DEBUG
        if (file_.is_open()) {
            file_ << photon_id << ",VOXEL_EMITTANCE,"
                  << position.x << "," << position.y << "," << position.z << ","
                  << direction.x << "," << direction.y << "," << direction.z << ","
                  << weight << ","
                  << voxel_coords.x << "," << voxel_coords.y << "," << voxel_coords.z << ","
                  << "-1," << emittance << "," << surface_type << "\n";
            file_.flush();
        }
#endif
    }

    ~DebugLogger() {
#ifdef _DEBUG
        if (file_.is_open()) {
            file_.close();
        }
#endif
    }

private:
    std::ofstream file_;
    DebugLogger() = default;
};
#pragma once

enum class CameraMode
{
	Orbit,
	Free
};

enum class VoxelMode
{
	Layers,
	Absorption,
	Emittance
};

enum class MenuItem
{
	DisplayMode,
	CameraMode,
	TogglePaths,
	ToggleBounds
};

struct Settings
{
	bool draw_paths {true};
	bool draw_voxels {true}; // Enable/disable voxel rendering
	bool draw_volume {true}; // Geometry-aware bounds
	bool draw_labels {true}; // Energy labels

	CameraMode camera_mode {CameraMode::Orbit};
	VoxelMode voxel_mode {VoxelMode::Absorption};

	// Default constructor uses default member initialization
	Settings() = default;
};

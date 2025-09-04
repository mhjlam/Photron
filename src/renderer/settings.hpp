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
	bool draw_paths;
	bool draw_voxels; // Enable/disable voxel rendering
	bool draw_volume; // Geometry-aware bounds
	bool draw_labels; // Energy labels

	CameraMode camera_mode;
	VoxelMode voxel_mode;

	Settings() {
		draw_paths = true;
		draw_voxels = true;
		draw_volume = true;
		draw_labels = true;

		camera_mode = CameraMode::Orbit;
		voxel_mode = VoxelMode::Absorption;
	}
};

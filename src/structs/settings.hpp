#pragma once

enum class CameraMode
{
	Orbit,
	Free
};
enum class VoxelMode
{
	None,
	Absorption,
	Emittance
};
enum class TextMode
{
	None,
	HeadsUpDisplay,
	All
};

enum class MenuItem
{
	DisplayMode,
	CameraMode,
	TextMode,
	TogglePaths,
	ToggleBounds
};

struct Settings
{
	bool draw_paths;
	bool draw_bounds;
	bool draw_path_labels; // New setting for energy labels

	CameraMode camera_mode;
	VoxelMode voxel_mode;
	TextMode text_mode;

	Settings() {
		draw_paths = true;
		draw_bounds = true;
		draw_path_labels = false;

		camera_mode = CameraMode::Orbit;
		voxel_mode = VoxelMode::Absorption;
		text_mode = TextMode::All;
	}
};

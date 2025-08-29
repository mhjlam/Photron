#pragma once

enum class CameraMode { Arc, Free };
enum class VoxelMode { None, Absorption, Emittance };
enum class TextMode { None, HeadsUpDisplay, All };
enum class GridMode { None, Partial, All };
enum class GeometryMode { None, White, Color };

enum class MenuItem {
	DisplayMode,
	CameraMode,
	GridMode,
	TextMode,
	GeometryMode,
	ToggleFrame,
	TogglePaths,
	ToggleBounds
};

struct Settings {
	bool draw_paths;
	bool draw_frame;
	bool draw_bounds;

	CameraMode camera_mode;
	VoxelMode voxel_mode;
	TextMode text_mode;
	GridMode grid_mode;
	GeometryMode geometry_mode;

	Settings() {
		draw_paths = true;
		draw_frame = false;
		draw_bounds = true;

		camera_mode = CameraMode::Arc;
		voxel_mode = VoxelMode::None;
		text_mode = TextMode::All;
		grid_mode = GridMode::None;
		geometry_mode = GeometryMode::White;
	}
};

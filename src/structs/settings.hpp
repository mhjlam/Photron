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
		draw_frame = false; // User wants this removed
		draw_bounds = true; // Enable wireframe quads by default (user request 3)

		camera_mode = CameraMode::Arc;
		voxel_mode = VoxelMode::Absorption; // Restore voxel display (user wants to see them)
		text_mode = TextMode::All;
		grid_mode = GridMode::None;
		geometry_mode = GeometryMode::Color; // Use colored geometry like original
	}
};

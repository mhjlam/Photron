/**
 * @file settings.hpp
 * @brief Rendering and interaction settings for the 3D visualization system
 *
 * Defines configuration options for camera modes, rendering styles, and
 * visualization features. These settings are typically controlled through
 * the GUI interface and affect real-time rendering behavior.
 */

#pragma once

/**
 * @enum CameraMode
 * @brief Camera interaction paradigms for 3D navigation
 */
enum class CameraMode
{
	Orbit, ///< Orbital/arc-ball mode - camera orbits around target point
	Free   ///< Free-flight mode - first-person style movement
};

/**
 * @enum VoxelMode
 * @brief Voxel visualization and coloring modes
 */
enum class VoxelMode
{
	Layers,     ///< Color voxels by material layer assignment
	Absorption, ///< Color voxels by accumulated energy absorption
	Emittance   ///< Color voxels by energy emission/transmission
};

/**
 * @struct Settings
 * @brief Comprehensive rendering and visualization settings
 *
 * Contains all user-configurable options that affect the visual
 * presentation of simulation results. Settings are typically modified
 * through the GUI overlay and immediately affect rendering.
 *
 * All boolean flags default to enabled for maximum information display,
 * allowing users to selectively disable features for performance or clarity.
 */
struct Settings
{
	bool draw_paths {true};                       ///< Render photon paths and scattering points
	bool draw_voxels {true};                      ///< Render voxelized geometry with energy coloring
	bool draw_volume {true};                      ///< Display geometry bounds and wireframes
	bool draw_labels {true};                      ///< Show energy value labels and annotations

	CameraMode camera_mode {CameraMode::Orbit};   ///< Current camera interaction mode
	VoxelMode voxel_mode {VoxelMode::Absorption}; ///< Current voxel coloring scheme

	/**
	 * @brief Default constructor initializes all settings to sensible defaults
	 */
	Settings() = default;
};

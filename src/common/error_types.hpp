/**
 * @file error_types.hpp
 * @brief Centralized error type definitions and message formatting utilities
 *
 * Provides structured error types for all major system components, replacing
 * inconsistent error handling patterns with type-safe error reporting. Includes
 * utility classes for consistent error message formatting across the codebase.
 */

#pragma once

#include <sstream>
#include <string>

/**
 * @brief Structured error types for different system components
 *
 * Replaces inconsistent error handling patterns with structured,
 * type-safe error reporting across the entire codebase.
 */

/**
 * @enum ConfigError
 * @brief Configuration system error types
 *
 * Covers TOML parsing, validation, and configuration file issues.
 */
enum class ConfigError
{
	FileNotFound,         ///< Configuration file not found at specified path
	ParseError,           ///< TOML syntax or structure parsing error
	ValidationError,      ///< Configuration values failed validation rules
	MissingRequiredField, ///< Required configuration field is missing
	InvalidValue,         ///< Configuration value is outside valid range
	GeometryError         ///< Geometric configuration parameters are inconsistent
};

/**
 * @enum SimulationError
 * @brief Monte Carlo simulation engine error types
 *
 * Covers photon transport, medium initialization, and energy conservation issues.
 */
enum class SimulationError
{
	InitializationFailed,        ///< Simulator failed to initialize properly
	ConfigNotInitialized,        ///< Configuration not loaded before simulation start
	InvalidConfiguration,        ///< Configuration is structurally invalid
	NoMediumFound,               ///< No medium defined for simulation
	NoVoxelFound,                ///< Voxel grid is empty or invalid
	InvalidPhotonState,          ///< Photon position/direction is invalid
	ExceededMaxIterations,       ///< Photon exceeded maximum scattering iterations
	NoLightSources,              ///< No light sources defined for simulation
	NoSources,                   ///< Legacy: alias for NoLightSources
	NoPhotons,                   ///< Zero photons specified for simulation
	SourceIntersectionFailure,   ///< Light source geometry intersection failed
	MediumInitializationFailure, ///< Medium layer initialization failed
	EnergyConservationViolation  ///< Total energy does not sum to 100%
};

/**
 * @enum IOError
 * @brief File input/output operation error types
 *
 * Covers all file system operations including reading, writing, and permissions.
 */
enum class IOError
{
	FileNotFound,     ///< Requested file does not exist
	PermissionDenied, ///< Insufficient permissions for file operation
	WriteError,       ///< Failed to write data to file
	ReadError,        ///< Failed to read data from file
	InvalidFormat,    ///< File format is not supported or corrupted
	DiskSpaceError    ///< Insufficient disk space for write operation
};

/**
 * @enum RenderError
 * @brief OpenGL rendering system error types
 *
 * Covers shader compilation, texture loading, and OpenGL state management.
 */
enum class RenderError
{
	InitializationFailed, ///< OpenGL context or renderer initialization failed
	ShaderCompileError,   ///< GLSL shader compilation or linking failed
	TextureLoadError,     ///< Texture file loading or GPU upload failed
	BufferError,          ///< Vertex/index buffer creation or update failed
	OpenGLError           ///< General OpenGL API error occurred
};

/**
 * @class ErrorMessage
 * @brief Utility class for consistent error message formatting
 *
 * Provides standardized error message formatting across all system components.
 * Ensures consistent format, context information, and user-friendly messages
 * for debugging and error reporting.
 */
class ErrorMessage
{
public:
	static std::string format(ConfigError err, const std::string& context = "") {
		std::ostringstream oss;
		oss << "Config Error: ";

		switch (err) {
			case ConfigError::FileNotFound: oss << "Configuration file not found."; break;
			case ConfigError::ParseError: oss << "Failed to parse configuration file."; break;
			case ConfigError::ValidationError: oss << "Configuration validation failed."; break;
			case ConfigError::MissingRequiredField: oss << "Missing required configuration field."; break;
			case ConfigError::InvalidValue: oss << "Invalid configuration value."; break;
			case ConfigError::GeometryError: oss << "Geometry configuration error."; break;
		}

		if (!context.empty()) {
			oss << " - " << context;
		}

		return oss.str();
	}

	static std::string format(SimulationError err, const std::string& context = "") {
		std::ostringstream oss;
		oss << "Simulation Error: ";

		switch (err) {
			case SimulationError::InitializationFailed: oss << "Simulation initialization failed"; break;
			case SimulationError::ConfigNotInitialized: oss << "Configuration not initialized"; break;
			case SimulationError::InvalidConfiguration: oss << "Invalid simulation configuration"; break;
			case SimulationError::NoMediumFound: oss << "No medium found at photon position"; break;
			case SimulationError::NoVoxelFound: oss << "No voxel found at photon position"; break;
			case SimulationError::InvalidPhotonState: oss << "Invalid photon state detected"; break;
			case SimulationError::ExceededMaxIterations: oss << "Exceeded maximum iterations"; break;
			case SimulationError::NoLightSources: oss << "No light sources available"; break;
			case SimulationError::NoSources: oss << "No sources found in configuration"; break;
			case SimulationError::NoPhotons: oss << "No photons available for simulation"; break;
			case SimulationError::SourceIntersectionFailure:
				oss << "Light source does not intersect with medium geometry";
				break;
			case SimulationError::MediumInitializationFailure: oss << "Failed to initialize simulation medium"; break;
			case SimulationError::EnergyConservationViolation: oss << "Energy conservation violation detected"; break;
		}

		if (!context.empty()) {
			oss << " - " << context;
		}

		return oss.str();
	}

	static std::string format(IOError err, const std::string& context = "") {
		std::ostringstream oss;
		oss << "I/O Error: ";

		switch (err) {
			case IOError::FileNotFound: oss << "File not found"; break;
			case IOError::PermissionDenied: oss << "Permission denied"; break;
			case IOError::WriteError: oss << "Write operation failed"; break;
			case IOError::ReadError: oss << "Read operation failed"; break;
			case IOError::InvalidFormat: oss << "Invalid file format"; break;
			case IOError::DiskSpaceError: oss << "Insufficient disk space"; break;
		}

		if (!context.empty()) {
			oss << " - " << context;
		}

		return oss.str();
	}

	static std::string format(RenderError err, const std::string& context = "") {
		std::ostringstream oss;
		oss << "Render Error: ";

		switch (err) {
			case RenderError::InitializationFailed: oss << "Renderer initialization failed"; break;
			case RenderError::ShaderCompileError: oss << "Shader compilation failed"; break;
			case RenderError::TextureLoadError: oss << "Texture loading failed"; break;
			case RenderError::BufferError: oss << "Buffer operation failed"; break;
			case RenderError::OpenGLError: oss << "OpenGL error occurred"; break;
		}

		if (!context.empty()) {
			oss << " - " << context;
		}

		return oss.str();
	}
};

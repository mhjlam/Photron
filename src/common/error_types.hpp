#pragma once

#include <string>
#include <sstream>

/**
 * @brief Structured error types for different system components
 * 
 * Replaces inconsistent error handling patterns with structured, 
 * type-safe error reporting across the entire codebase.
 */

// Configuration-related errors
enum class ConfigError {
    FileNotFound,
    ParseError, 
    ValidationError,
    MissingRequiredField,
    InvalidValue,
    GeometryError
};

// Simulation-related errors  
enum class SimulationError {
    InitializationFailed,
    ConfigNotInitialized,
    InvalidConfiguration,
    NoMediumFound,
    NoVoxelFound,
    InvalidPhotonState,
    ExceededMaxIterations,
    NoLightSources,
    NoSources,
    NoPhotons,
    SourceIntersectionFailure,
    MediumInitializationFailure,
    EnergyConservationViolation
};

// File I/O errors
enum class IOError {
    FileNotFound,
    PermissionDenied,
    WriteError,
    ReadError,
    InvalidFormat,
    DiskSpaceError
};

// Rendering/Graphics errors
enum class RenderError {
    InitializationFailed,
    ShaderCompileError,
    TextureLoadError,
    BufferError,
    OpenGLError
};

// Error message formatting utilities
class ErrorMessage {
public:
    static std::string format(ConfigError err, const std::string& context = "") {
        std::ostringstream oss;
        oss << "Config Error: ";
        
        switch (err) {
            case ConfigError::FileNotFound:
                oss << "Configuration file not found.";
                break;
            case ConfigError::ParseError:
                oss << "Failed to parse configuration file.";
                break;
            case ConfigError::ValidationError:
                oss << "Configuration validation failed.";
                break;
            case ConfigError::MissingRequiredField:
                oss << "Missing required configuration field.";
                break;
            case ConfigError::InvalidValue:
                oss << "Invalid configuration value.";
                break;
            case ConfigError::GeometryError:
                oss << "Geometry configuration error.";
                break;
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
            case SimulationError::InitializationFailed:
                oss << "Simulation initialization failed";
                break;
            case SimulationError::ConfigNotInitialized:
                oss << "Configuration not initialized";
                break;
            case SimulationError::InvalidConfiguration:
                oss << "Invalid simulation configuration";
                break;
            case SimulationError::NoMediumFound:
                oss << "No medium found at photon position";
                break;
            case SimulationError::NoVoxelFound:
                oss << "No voxel found at photon position";
                break;
            case SimulationError::InvalidPhotonState:
                oss << "Invalid photon state detected";
                break;
            case SimulationError::ExceededMaxIterations:
                oss << "Exceeded maximum iterations";
                break;
            case SimulationError::NoLightSources:
                oss << "No light sources available";
                break;
            case SimulationError::NoSources:
                oss << "No sources found in configuration";
                break;
            case SimulationError::NoPhotons:
                oss << "No photons available for simulation";
                break;
            case SimulationError::SourceIntersectionFailure:
                oss << "Light source does not intersect with medium geometry";
                break;
            case SimulationError::MediumInitializationFailure:
                oss << "Failed to initialize simulation medium";
                break;
            case SimulationError::EnergyConservationViolation:
                oss << "Energy conservation violation detected";
                break;
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
            case IOError::FileNotFound:
                oss << "File not found";
                break;
            case IOError::PermissionDenied:
                oss << "Permission denied";
                break;
            case IOError::WriteError:
                oss << "Write operation failed";
                break;
            case IOError::ReadError:
                oss << "Read operation failed";
                break;
            case IOError::InvalidFormat:
                oss << "Invalid file format";
                break;
            case IOError::DiskSpaceError:
                oss << "Insufficient disk space";
                break;
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
            case RenderError::InitializationFailed:
                oss << "Renderer initialization failed";
                break;
            case RenderError::ShaderCompileError:
                oss << "Shader compilation failed";
                break;
            case RenderError::TextureLoadError:
                oss << "Texture loading failed";
                break;
            case RenderError::BufferError:
                oss << "Buffer operation failed";
                break;
            case RenderError::OpenGLError:
                oss << "OpenGL error occurred";
                break;
        }
        
        if (!context.empty()) {
            oss << " - " << context;
        }
        
        return oss.str();
    }
};
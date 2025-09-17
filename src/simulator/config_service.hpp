#pragma once

#include "config.hpp"
#include <memory>
#include <stdexcept>

/**
 * @brief ConfigService provides global access to the current configuration
 * 
 * This service pattern avoids global variables while making the config accessible
 * from any part of the codebase. The config is registered once during application
 * initialization and can be safely accessed throughout the application lifetime.
 * 
 * Thread-safe for read access after initialization.
 */
class ConfigService {
public:
    /**
     * @brief Register the config instance for global access
     * @param config Reference to the config instance (must remain valid)
     */
    static void initialize(const Config& config) {
        instance_ = &config;
    }
    
    /**
     * @brief Get the current config instance
     * @return Reference to the current config
     * @throws std::runtime_error if config has not been initialized
     */
    static const Config& get() {
        if (!instance_) {
            throw std::runtime_error("ConfigService not initialized. Call initialize() first.");
        }
        return *instance_;
    }
    
    /**
     * @brief Check if config service has been initialized
     * @return true if initialized, false otherwise
     */
    static bool is_initialized() {
        return instance_ != nullptr;
    }
    
    /**
     * @brief Reset the config service (for testing or shutdown)
     */
    static void reset() {
        instance_ = nullptr;
    }

private:
    static const Config* instance_;
    
    // Prevent instantiation
    ConfigService() = delete;
    ConfigService(const ConfigService&) = delete;
    ConfigService& operator=(const ConfigService&) = delete;
};
#pragma once

#include <string>
#include <iostream>
#include <sstream>

/**
 * @brief Centralized error handling and reporting system
 * 
 * Provides consistent error handling across the entire codebase with
 * appropriate routing to console, logs, or both based on error severity.
 */
class ErrorHandler {
public:
    enum class Level {
        Info,      // Informational messages (log only)
        Warning,   // Warnings (console + log)  
        Error,     // Errors (console + log)
        Critical   // Critical errors (console + log + immediate)
    };

    // Singleton access
    static ErrorHandler& instance();

    // Main error reporting methods
    void report_info(const std::string& message);
    void report_warning(const std::string& message);
    void report_error(const std::string& message);
    void report_critical(const std::string& message);

    // Formatted reporting with context
    void report_error(const std::string& component, const std::string& operation, const std::string& details);
    void report_warning(const std::string& component, const std::string& operation, const std::string& details);

    // Config and initialization errors (always to console)
    void report_config_error(const std::string& message);
    void report_initialization_error(const std::string& component, const std::string& message);
    
    // Configuration methods
    void set_logging_enabled(bool enabled) { logging_enabled_ = enabled; }
    bool is_logging_enabled() const { return logging_enabled_; }

private:
    ErrorHandler() = default;
    
    void write_to_console(Level level, const std::string& message);
    void write_to_log(Level level, const std::string& message);
    std::string format_message(Level level, const std::string& message) const;
    std::string get_level_string(Level level) const;
    
    bool logging_enabled_ = false; // Set by App during initialization
};

// Convenience macros for common patterns
#define REPORT_INFO(msg) ErrorHandler::instance().report_info(msg)
#define REPORT_WARNING(msg) ErrorHandler::instance().report_warning(msg)
#define REPORT_ERROR(msg) ErrorHandler::instance().report_error(msg)
#define REPORT_CRITICAL(msg) ErrorHandler::instance().report_critical(msg)

// Formatted reporting macros
#define REPORT_COMPONENT_ERROR(comp, op, details) ErrorHandler::instance().report_error(comp, op, details)
#define REPORT_COMPONENT_WARNING(comp, op, details) ErrorHandler::instance().report_warning(comp, op, details)

// Config/init specific macros (always console)
#define REPORT_CONFIG_ERROR(msg) ErrorHandler::instance().report_config_error(msg)
#define REPORT_INIT_ERROR(comp, msg) ErrorHandler::instance().report_initialization_error(comp, msg)
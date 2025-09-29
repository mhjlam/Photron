#include "error_handler.hpp"

// Add GLM include for complete type definitions
#include <glm/glm.hpp>
#include "../simulator/logger.hpp"
#include <chrono>
#include <iomanip>

ErrorHandler& ErrorHandler::instance() {
    static ErrorHandler handler;
    return handler;
}

void ErrorHandler::report_info(const std::string& message) {
    // Info only goes to log files (if logging enabled)
    if (logging_enabled_) {
        write_to_log(Level::Info, message);
    }
}

void ErrorHandler::report_warning(const std::string& message) {
    // Warnings go to both console and log
    write_to_console(Level::Warning, message);
    if (logging_enabled_) {
        write_to_log(Level::Warning, message);
    }
}

void ErrorHandler::report_error(const std::string& message) {
    // Errors go to both console and log
    write_to_console(Level::Error, message);
    if (logging_enabled_) {
        write_to_log(Level::Error, message);
    }
}

void ErrorHandler::report_critical(const std::string& message) {
    // Critical errors always go everywhere
    write_to_console(Level::Critical, message);
    if (logging_enabled_) {
        write_to_log(Level::Critical, message);
    }
}

void ErrorHandler::report_error(const std::string& component, const std::string& operation, const std::string& details) {
    std::ostringstream oss;
    oss << component << ": " << operation << " - " << details;
    report_error(oss.str());
}

void ErrorHandler::report_warning(const std::string& component, const std::string& operation, const std::string& details) {
    std::ostringstream oss;
    oss << component << ": " << operation << " - " << details;
    report_warning(oss.str());
}

void ErrorHandler::report_config_error(const std::string& message) {
    // Config errors always go to console (config might not be ready)
    write_to_console(Level::Error, message);
}

void ErrorHandler::report_initialization_error(const std::string& component, const std::string& message) {
    // Initialization errors always go to console
    std::ostringstream oss;
    oss << component << ": " << message;
    write_to_console(Level::Error, oss.str());
}

void ErrorHandler::write_to_console(Level level, const std::string& message) {
    std::string formatted_msg = format_message(level, message);
    
    // Route to appropriate stream
    if (level == Level::Info) {
        std::cout << formatted_msg << std::endl;
    } else {
        std::cerr << formatted_msg << std::endl;
    }
}

void ErrorHandler::write_to_log(Level level, const std::string& message) {
    // Delegate to existing Logger system
    switch (level) {
        case Level::Info:
            Logger::instance().log_info(message);
            break;
        case Level::Warning:
            Logger::instance().log_warning(message);
            break;
        case Level::Error:
        case Level::Critical:
            Logger::instance().log_error(message);
            break;
    }
}

std::string ErrorHandler::format_message(Level level, const std::string& message) const {
    // Simple formatting for console output
    return get_level_string(level) + ": " + message;
}

std::string ErrorHandler::get_level_string(Level level) const {
    switch (level) {
        case Level::Info:    return "Info";
        case Level::Warning: return "Warning";
        case Level::Error:   return "Error";
        case Level::Critical: return "CRITICAL";
        default:             return "Unknown";
    }
}
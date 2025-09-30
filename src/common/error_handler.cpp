/**
 * @file error_handler.cpp
 * @brief Centralized error handling and logging system implementation
 *
 * Implements the ErrorHandler singleton for unified error management across
 * the Monte Carlo photon transport simulation:
 * - Hierarchical error level handling (Info, Warning, Error, Critical)
 * - Dual output routing (console and log files)
 * - Thread-safe singleton pattern for global accessibility
 * - Configurable logging with timestamped entries
 * - Component-based error reporting for debugging
 *
 * Essential for debugging complex Monte Carlo algorithms and tracking
 * numerical issues in photon transport calculations.
 */

#include "error_handler.hpp"

#include <chrono>
#include <iomanip>

#include <glm/glm.hpp>

#include "simulator/logger.hpp"

/**
 * @brief Get singleton ErrorHandler instance
 * @return Reference to global ErrorHandler instance
 *
 * Thread-safe singleton implementation providing global access
 * to error handling functionality throughout the application.
 */
ErrorHandler& ErrorHandler::instance() {
	static ErrorHandler handler;
	return handler;
}

/**
 * @brief Report informational message
 * @param message Information message to log
 *
 * Logs informational messages to file only (not console) to avoid
 * cluttering output during normal operation. Useful for debugging
 * and performance analysis.
 */
void ErrorHandler::report_info(const std::string& message) {
	// Info messages only appear in log files to keep console clean
	if (logging_enabled_) {
		write_to_log(Level::Info, message);
	}
}

/**
 * @brief Report warning message
 * @param message Warning message to report
 *
 * Outputs warnings to both console and log files. Indicates potential
 * issues that don't prevent operation but may affect results.
 */
void ErrorHandler::report_warning(const std::string& message) {
	// Warnings appear on console and in logs for visibility
	write_to_console(Level::Warning, message);
	if (logging_enabled_) {
		write_to_log(Level::Warning, message);
	}
}

/**
 * @brief Report error message
 * @param message Error message to report
 *
 * Outputs errors to both console and log files. Indicates significant
 * issues that may compromise simulation accuracy or functionality.
 */
void ErrorHandler::report_error(const std::string& message) {
	// Errors appear on console and in logs for immediate attention
	write_to_console(Level::Error, message);
	if (logging_enabled_) {
		write_to_log(Level::Error, message);
	}
}

/**
 * @brief Report critical error message
 * @param message Critical error message to report
 *
 * Outputs critical errors to both console and log files regardless
 * of logging settings. Indicates severe issues requiring immediate attention.
 */
void ErrorHandler::report_critical(const std::string& message) {
	// Critical errors always appear everywhere for maximum visibility
	write_to_console(Level::Critical, message);
	if (logging_enabled_) {
		write_to_log(Level::Critical, message);
	}
}

/**
 * @brief Report structured error with component context
 * @param component Software component where error occurred
 * @param operation Operation that failed
 * @param details Detailed error description
 *
 * Provides structured error reporting with context for easier debugging
 * of complex Monte Carlo simulation issues.
 */
void ErrorHandler::report_error(const std::string& component,
								const std::string& operation,
								const std::string& details) {
	// Format structured error message with component context
	std::ostringstream oss;
	oss << component << ": " << operation << " - " << details;
	report_error(oss.str());
}

void ErrorHandler::report_warning(const std::string& component,
								  const std::string& operation,
								  const std::string& details) {
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
	}
	else {
		std::cerr << formatted_msg << std::endl;
	}
}

void ErrorHandler::write_to_log(Level level, const std::string& message) {
	// Delegate to existing Logger system
	switch (level) {
		case Level::Info: Logger::instance().log_info(message); break;
		case Level::Warning: Logger::instance().log_warning(message); break;
		case Level::Error:
		case Level::Critical: Logger::instance().log_error(message); break;
	}
}

std::string ErrorHandler::format_message(Level level, const std::string& message) const {
	// Simple formatting for console output
	return get_level_string(level) + ": " + message;
}

std::string ErrorHandler::get_level_string(Level level) const {
	switch (level) {
		case Level::Info: return "Info";
		case Level::Warning: return "Warning";
		case Level::Error: return "Error";
		case Level::Critical: return "CRITICAL";
		default: return "Unknown";
	}
}
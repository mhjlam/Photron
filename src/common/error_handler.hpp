/**
 * @file error_handler.hpp
 * @brief Centralized error handling and reporting system
 *
 * Provides consistent error handling across the entire application with
 * intelligent routing based on error severity and context. Supports both
 * console output for immediate user feedback and structured logging for
 * debugging and analysis.
 */

#pragma once

#include <iostream>
#include <sstream>
#include <string>

/**
 * @class ErrorHandler
 * @brief Singleton error reporting and logging system
 *
 * The ErrorHandler class provides centralized error management with:
 *
 * - **Severity Levels**: Info, Warning, Error, and Critical classifications
 * - **Smart Routing**: Automatic decision of console vs. log output based on context
 * - **Structured Messages**: Consistent formatting with timestamps and context
 * - **Specialized Handling**: Special treatment for configuration and initialization errors
 * - **Runtime Configuration**: Enable/disable logging as needed
 *
 * Design principles:
 * - Single point of truth for all error reporting
 * - Non-blocking operation (no exceptions thrown)
 * - Thread-safe singleton implementation
 * - Graceful degradation if logging fails
 *
 * Usage examples:
 * @code
 * ErrorHandler::instance().report_error("Failed to load configuration");
 * ErrorHandler::instance().report_warning("Renderer", "initialize", "OpenGL 4.5 not available, falling back to 4.1");
 * @endcode
 */
class ErrorHandler
{
public:
	/**
	 * @enum Level
	 * @brief Error severity levels for intelligent message routing
	 */
	enum class Level
	{
		Info,    ///< Informational messages (log only, no console output)
		Warning, ///< Warnings (both console and log output)
		Error,   ///< Errors (both console and log output)
		Critical ///< Critical errors (immediate console + log output)
	};

	/**
	 * @brief Get singleton instance of ErrorHandler
	 *
	 * Thread-safe singleton access using Meyer's singleton pattern.
	 *
	 * @return ErrorHandler& Reference to the global ErrorHandler instance
	 */
	static ErrorHandler& instance();

	// Basic error reporting methods

	/**
	 * @brief Report informational message (log only)
	 *
	 * @param message Informational message text
	 */
	void report_info(const std::string& message);

	/**
	 * @brief Report warning message (console + log)
	 *
	 * @param message Warning message text
	 */
	void report_warning(const std::string& message);

	/**
	 * @brief Report error message (console + log)
	 *
	 * @param message Error message text
	 */
	void report_error(const std::string& message);

	/**
	 * @brief Report critical error message (immediate console + log)
	 *
	 * @param message Critical error message text
	 */
	void report_critical(const std::string& message);

	// Formatted reporting with context

	/**
	 * @brief Report error with structured context information
	 *
	 * @param component Component or subsystem name where error occurred
	 * @param operation Operation or function name that failed
	 * @param details Detailed error description
	 */
	void report_error(const std::string& component, const std::string& operation, const std::string& details);

	/**
	 * @brief Report warning with structured context information
	 *
	 * @param component Component or subsystem name where warning occurred
	 * @param operation Operation or function name that generated warning
	 * @param details Detailed warning description
	 */
	void report_warning(const std::string& component, const std::string& operation, const std::string& details);

	// Specialized error types

	/**
	 * @brief Report configuration-related errors (always to console)
	 *
	 * Configuration errors are always displayed to console since they
	 * typically prevent the application from running properly.
	 *
	 * @param message Configuration error description
	 */
	void report_config_error(const std::string& message);

	/**
	 * @brief Report initialization errors (always to console)
	 *
	 * @param component Component that failed to initialize
	 * @param message Specific initialization error details
	 */
	void report_initialization_error(const std::string& component, const std::string& message);

	// Configuration methods

	/**
	 * @brief Enable or disable file logging
	 *
	 * @param enabled True to enable logging, false to disable
	 */
	void set_logging_enabled(bool enabled) { logging_enabled_ = enabled; }

	/**
	 * @brief Check if file logging is currently enabled
	 *
	 * @return true if logging is enabled, false otherwise
	 */
	bool is_logging_enabled() const { return logging_enabled_; }

private:
	ErrorHandler() = default;

	void write_to_console(Level level, const std::string& message);
	void write_to_log(Level level, const std::string& message);
	std::string format_message(Level level, const std::string& message) const;
	std::string get_level_string(Level level) const;

	bool logging_enabled_ {false}; // Set by App during initialization
};

// Convenience macros for common patterns
#define REPORT_INFO(msg)     ErrorHandler::instance().report_info(msg)
#define REPORT_WARNING(msg)  ErrorHandler::instance().report_warning(msg)
#define REPORT_ERROR(msg)    ErrorHandler::instance().report_error(msg)
#define REPORT_CRITICAL(msg) ErrorHandler::instance().report_critical(msg)

// Formatted reporting macros
#define REPORT_COMPONENT_ERROR(comp, op, details)   ErrorHandler::instance().report_error(comp, op, details)
#define REPORT_COMPONENT_WARNING(comp, op, details) ErrorHandler::instance().report_warning(comp, op, details)

// Config/init specific macros (always console)
#define REPORT_CONFIG_ERROR(msg)     ErrorHandler::instance().report_config_error(msg)
#define REPORT_INIT_ERROR(comp, msg) ErrorHandler::instance().report_initialization_error(comp, msg)
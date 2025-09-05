#pragma once

#include <numbers>

/**
 * @file numerical_constants.hpp
 * @brief Centralized numerical constants for consistent precision handling across the codebase
 * 
 * This file defines standard epsilon values and other numerical constants used throughout
 * the Monte Carlo simulation for robust floating-point comparisons and geometric calculations.
 */

namespace NumericalConstants {
    
    // Primary epsilon for high-precision geometric calculations
    constexpr double GEOMETRIC_EPSILON = 1e-12;
    
    // Epsilon for ray-triangle intersection self-intersection avoidance
    constexpr double RAY_EPSILON = 1e-12;
    
    // Epsilon for AABB ray intersection (slab method)
    constexpr double AABB_EPSILON = 1e-12;
    
    // Epsilon for Fresnel reflection calculations
    constexpr double FRESNEL_EPSILON = 1e-12;
    
    // Epsilon for Monte Carlo step size calculations  
    constexpr double MCMC_EPSILON = 1e-12;
    
    // Epsilon for photon position/direction comparisons
    constexpr double PHOTON_EPSILON = 1e-10;
    
    // Epsilon for voxel boundary calculations
    constexpr double VOXEL_EPSILON = 1e-6;
    
    // Mathematical constants
    namespace Math {
        constexpr double PI = std::numbers::pi;
        constexpr double TWO_PI = 2.0 * PI;
        constexpr double HALF_PI = PI / 2.0;
        constexpr double INV_PI = 1.0 / PI;
        constexpr double INV_TWO_PI = 1.0 / TWO_PI;
    }
    
    // Physical constants for Monte Carlo simulation
    namespace Physics {
        constexpr double MIN_PHOTON_WEIGHT = 1e-4;
        constexpr double DEFAULT_SURVIVAL_PROBABILITY = 0.1;
        constexpr double MAX_SURVIVAL_PROBABILITY = 0.5;
    }
}

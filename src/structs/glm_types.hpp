#pragma once

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/norm.hpp>

// Use double precision for mathematical computations (dvec) and float for OpenGL (vec)

// Conversion helpers between float and double precision
inline glm::vec3 toFloat(const glm::dvec3& v) {
    return glm::vec3(static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z));
}

inline glm::vec4 toFloat(const glm::dvec4& v) {
    return glm::vec4(static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z), static_cast<float>(v.w));
}

inline glm::dvec3 toDouble(const glm::vec3& v) {
    return glm::dvec3(static_cast<double>(v.x), static_cast<double>(v.y), static_cast<double>(v.z));
}

inline glm::dvec2 toDouble(const glm::vec2& v) {
    return glm::dvec2(static_cast<double>(v.x), static_cast<double>(v.y));
}

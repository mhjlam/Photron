#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Type aliases for gradual migration from custom types to GLM
// Use double precision for mathematical computations
using Vec2 = glm::dvec2;
using Vec3 = glm::dvec3;
using Vec4 = glm::dvec4;

// Matrix types
using Mat3 = glm::dmat3;
using Mat4 = glm::dmat4;

// Float precision for OpenGL (GPU-friendly)
using Vec2f = glm::vec2;
using Vec3f = glm::vec3;
using Vec4f = glm::vec4;
using Mat3f = glm::mat3;
using Mat4f = glm::mat4;

// Conversion helpers
inline Vec3f toFloat(const Vec3& v) {
    return Vec3f(static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z));
}

inline Vec4f toFloat(const Vec4& v) {
    return Vec4f(static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z), static_cast<float>(v.w));
}

inline Vec3 toDouble(const Vec3f& v) {
    return Vec3(static_cast<double>(v.x), static_cast<double>(v.y), static_cast<double>(v.z));
}

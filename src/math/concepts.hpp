#pragma once

#include <algorithm>
#include <concepts>
#include <ranges>

// C++20 Concepts for geometric types
template<typename T>
concept Point3D = requires(T t) {
	{ t.x } -> std::convertible_to<double>;
	{ t.y } -> std::convertible_to<double>;
	{ t.z } -> std::convertible_to<double>;
};

template<typename T>
concept Vector3D = Point3D<T>;

// Container concepts
template<typename Container>
concept TriangleContainer = std::ranges::range<Container>;

template<typename Container>
concept IndexContainer =
	std::ranges::range<Container> && std::convertible_to<std::ranges::range_value_t<Container>, int>;

template<typename Container>
concept ConfigContainer = std::ranges::range<Container>;

// Numeric concepts
template<typename T>
concept Numeric = std::integral<T> || std::floating_point<T>;

// Additional concepts for volume calculations
template<typename Vec>
concept Vector3Like = requires(Vec v) {
	{ v.x } -> std::convertible_to<double>;
	{ v.y } -> std::convertible_to<double>;
	{ v.z } -> std::convertible_to<double>;
};

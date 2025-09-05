#pragma once

#include <memory>

#include <glm/glm.hpp>

struct PhotonNode
{
	glm::dvec3 position;
	double value;

	std::shared_ptr<PhotonNode> prev = nullptr; // previous internal vertex
	std::shared_ptr<PhotonNode> next = nullptr; // next internal vertex
	std::shared_ptr<PhotonNode> emit = nullptr; // external vertex

	PhotonNode(double xx, double yy, double zz, double v) noexcept : position(xx, yy, zz), value(v) {}

	PhotonNode(const glm::dvec3& pos, double v) noexcept : position(pos), value(v) {}
};

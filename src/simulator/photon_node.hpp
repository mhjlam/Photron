#pragma once

#include <memory>

#include <glm/glm.hpp>

// Forward declarations to avoid circular dependency
struct Photon;
struct Emitter;

struct PhotonNode
{
	glm::dvec3 position;
	double value;

	// Exit classification for energy labeling
	enum class ExitType { NONE, REFLECTED, TRANSMITTED };
	ExitType exit_type {ExitType::NONE};

	std::shared_ptr<PhotonNode> prev = nullptr; // previous internal vertex
	std::shared_ptr<PhotonNode> next = nullptr; // next internal vertex
	std::shared_ptr<PhotonNode> emit = nullptr; // external vertex
	
	// Optional connection to associated emitter (when this node represents an exit point)
	std::shared_ptr<Emitter> emitter = nullptr;

	PhotonNode(double xx, double yy, double zz, double v) noexcept : position(xx, yy, zz), value(v) {}

	PhotonNode(const glm::dvec3& pos, double v) noexcept : position(pos), value(v) {}

	// Constructor with exit type for external vertices
	PhotonNode(const glm::dvec3& pos, double v, ExitType exit_classification) noexcept : 
		position(pos), value(v), exit_type(exit_classification) {}
};

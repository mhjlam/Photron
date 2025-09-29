#pragma once

// Standard library includes
#include <memory>
#include <vector>
#include <set>
#include <string>

// Third-party includes for interface types
#include <glm/glm.hpp>

// Project includes for member types
#include "math/triangle.hpp"
#include "simulator/voxel.hpp"
#include "math/concepts.hpp"

// Source properties (moved from Photon class to avoid circular dependencies)
struct Source
{
	uint64_t id {0};                     // identifier
	glm::dvec3 origin {0.0};             // origin
	glm::dvec3 direction {0.0};          // direction of incidence
	glm::dvec3 specular_direction {0.0}; // direction of specular reflectance
	glm::dvec3 intersect {0.0};          // intersection point
	Triangle triangle;                   // triangle at intersection point

	Source() = default;
	explicit Source(uint64_t i, const glm::dvec3& p, const glm::dvec3& v) noexcept : 
		id(i), origin(p), direction(v) {}
};

// Emitted photon
struct Emitter
{
	uint64_t id;
	glm::dvec3 position;
	glm::dvec3 direction;
	double weight;
	
	// Exit classification for accurate energy labeling
	enum class ExitType { NONE, REFLECTED, TRANSMITTED };
	ExitType exit_type {ExitType::NONE};

	template<Point3D P, Vector3D V>
	explicit constexpr Emitter(uint64_t i, const P& p, const V& d, double w, ExitType exit_classification = ExitType::NONE) noexcept :
		id(i), position(static_cast<glm::dvec3>(p)), direction(static_cast<glm::dvec3>(d)), weight(w), exit_type(exit_classification) {}
};

class PhotonNode
{
public:
	glm::dvec3 position;
	double value;

	// Exit classification for energy labeling
	enum class ExitType { NONE, REFLECTED, TRANSMITTED };
	ExitType exit_type {ExitType::NONE};

	std::shared_ptr<PhotonNode> prev{nullptr}; // previous internal vertex
	std::shared_ptr<PhotonNode> next{nullptr}; // next internal vertex
	std::shared_ptr<PhotonNode> emit{nullptr}; // external vertex
	
	// Optional connection to associated emitter (when this node represents an exit point)
	std::shared_ptr<Emitter> emitter{nullptr};

	PhotonNode(double xx, double yy, double zz, double v) noexcept : position(xx, yy, zz), value(v) {}
	PhotonNode(const glm::dvec3& pos, double v) noexcept : position(pos), value(v) {}
	PhotonNode(const glm::dvec3& pos, double v, ExitType exit_classification) noexcept : 
		position(pos), value(v), exit_type(exit_classification) {}
};

class Photon
{
public:
	// Core photon properties
	uint64_t id {0};               		// identifier
	bool alive {true};             		// true if photon still propagates
	bool cross {false};            		// true if photon crosses voxel boundary in substep
	bool exits {false};		   			// true if photon crosses medium boundary in substep

	// ENERGY CONSERVATION TRACKING
	double total_energy_budget {0.0};	// total energy available to this photon (includes Russian Roulette amplifications)
	double total_energy_radiated {0.0};	// cumulative energy radiated through all radiate() calls
	double total_energy_absorbed {0.0};	// cumulative energy absorbed through all absorption events
	int radiate_call_count {0};		// number of times radiate() was called for this photon

	double step {0.0};             		// step distance between scattering points
	double sub_step {0.0};         		// step distance inside a voxel
	double weight {0.0};           		// remaining weight of the photon packet

	glm::dvec3 position {0.0};     		// position at the start of a substep
	glm::dvec3 direction {0.0};    		// propagation direction
	Voxel* voxel {nullptr};        		// resident voxel at start of substep (non-owning)
	Voxel* prev_voxel {nullptr};   		// voxel before crossing an interface (non-owning)
	Voxel* last_surface_voxel {nullptr}; // last surface voxel traversed (for emittance recording)

	glm::dvec3 intersect {0.0};    		// voxel boundary intersection
	glm::dvec3 voxel_normal {0.0}; 		// voxel boundary intersection normal

	bool scatters {false};         		// true if photon path scatters at least once
	uint32_t scatter_count {0};			// total number of scattering events for this photon

	// Exit classification (set when photon exits the medium)
	enum class ExitType { NONE, REFLECTED, TRANSMITTED };
	ExitType exit_type {ExitType::NONE};	// how the photon exited the medium

	// Source properties (merged from Source class) - now using Source
	Source source;

	// Path tracking (integrated from PhotonPath)
	std::shared_ptr<PhotonNode> path_head{nullptr};
	std::shared_ptr<PhotonNode> path_last{nullptr};
	uint64_t num_seg_int {1}; // internal segments
	uint64_t num_seg_ext {1}; // emittant segments

	// Constructors
	Photon() = default;
	explicit Photon(uint64_t i) noexcept : id(i) {}
	
	// Constructor with source data merged in
	explicit Photon(uint64_t i, const glm::dvec3& src_origin, const glm::dvec3& src_direction) noexcept;
	
	// Constructor with full source data
	explicit Photon(uint64_t i, const glm::dvec3& src_origin, const glm::dvec3& src_direction,
					const glm::dvec3& spec_direction, const glm::dvec3& src_intersect, 
					const Triangle& src_triangle) noexcept;

	// Path management methods (from PhotonPath)
	void add_internal_vertex(std::shared_ptr<PhotonNode> vert) noexcept;
	void add_external_vertex(std::shared_ptr<PhotonNode> vert) noexcept;
	void initialize_path(const glm::dvec3& start_pos, double path_weight);

	// Source methods
	void set_source_data(const Source& src_data);
	void set_source_data(uint64_t src_id, const glm::dvec3& origin, const glm::dvec3& src_direction);
	
	inline const Source& get_source_data() const noexcept { return source; }

	// Material property accessors with null safety (for backward compatibility)
	[[nodiscard]] inline double g() const noexcept { return (voxel && voxel->material) ? voxel->material->g() : 0.0; }
	[[nodiscard]] inline double eta() const noexcept { return (voxel && voxel->material) ? voxel->material->eta() : 0.0; }
	[[nodiscard]] inline double mu_a() const noexcept { return (voxel && voxel->material) ? voxel->material->mu_a() : 0.0; }
	[[nodiscard]] inline double mu_s() const noexcept { return (voxel && voxel->material) ? voxel->material->mu_s() : 0.0; }

	// Path analysis methods
	[[nodiscard]] glm::dvec3 get_entrance_position() const noexcept;
	[[nodiscard]] glm::dvec3 get_entrance_direction() const noexcept;
	[[nodiscard]] glm::dvec3 get_exit_position() const noexcept;
	[[nodiscard]] glm::dvec3 get_exit_direction() const noexcept;
	[[nodiscard]] glm::dvec3 get_termination_position() const noexcept;
	[[nodiscard]] glm::dvec3 get_termination_direction() const noexcept;
	[[nodiscard]] uint32_t get_path_scatter_count() const noexcept;
	[[nodiscard]] double get_total_path_length() const noexcept;
	[[nodiscard]] double get_total_absorption_deposited() const noexcept;
	[[nodiscard]] bool has_exit() const noexcept;
	[[nodiscard]] bool exited_medium() const noexcept;
	[[nodiscard]] std::string get_termination_reason() const noexcept;
	[[nodiscard]] std::vector<glm::dvec3> get_all_exit_positions() const;

	// Backward compatibility aliases for source access
	[[nodiscard]] inline const glm::dvec3& source_origin() const noexcept { return source.origin; }
	[[nodiscard]] inline const glm::dvec3& source_direction() const noexcept { return source.direction; }
	[[nodiscard]] inline const glm::dvec3& specular_direction() const noexcept { return source.specular_direction; }
	[[nodiscard]] inline const glm::dvec3& source_intersect() const noexcept { return source.intersect; }
	[[nodiscard]] inline const Triangle& source_triangle() const noexcept { return source.triangle; }
};

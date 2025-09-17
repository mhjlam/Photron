# Monte Carlo Simulation Analysis

## Critical Shortcomings

### 1. Energy Conservation Violations

- **Dual Recording**: Energy recorded in both medium records AND voxel data, causing double-counting
- **True Splitting Issues**: In partial reflection mode, same energy recorded 20+ times per photon
- **Photon Termination**: Weight amplification in Russian Roulette not properly tracked in energy budget

### 2. Boundary Detection Complexity

- **Three-Level Detection System**: Voxel boundaries + mesh boundaries + geometry tests create conflicts
- **Interior Exit Bug**: Photons can exit from non-surface voxels (physically impossible)
- **Numerical Instability**: Boundary nudging with inconsistent epsilon values
- **Inconsistent Normals**: Surface normals sometimes incorrectly oriented

### 3. Performance Issues

- **Excessive Ray Marching**: 10× oversampling in `track_voxel_path_and_deposit()`
- **Hot Path Complexity**: Multiple safety checks and fallbacks in critical loops
- **Memory Allocation**: Dynamic vector operations in tight loops
- **Redundant Calculations**: Fresnel equations computed multiple times for same interface

## Essential Improvements

### 1. Fix Energy Conservation

**Recommendation**: Use single voxel-based system, aggregate when needed.

- Use only voxel-based energy recording
- Remove medium record energy tracking
- Add energy conservation assertions

```cpp
// Proposed: Single energy tracking per photon
struct Photon {
    double initial_weight = 1.0;
    double current_weight;
    double absorbed_weight = 0.0;
    double transmitted_weight = 0.0;
    
    bool energy_conserved() const {
        return abs(initial_weight - (current_weight + absorbed_weight + transmitted_weight)) < 1e-12;
    }
};
```

### 2. Simplify Boundary Detection

**Recommendation**: Unified spatial data structure (octree/BVH).

- Unified boundary detection function
- Single ray-casting pass
- Consistent normal orientation

```cpp
// Proposed: Unified boundary system
enum class BoundaryType { None, VoxelEdge, MaterialInterface, GeometryExit };

struct BoundaryIntersection {
    BoundaryType type;
    double distance;
    glm::dvec3 point;
    glm::dvec3 normal;
    Material* new_material;
};
```

### 3. Optimize Critical Paths

- Replace ray marching with DDA voxel traversal
- Cache Fresnel calculations for repeated interfaces
- Remove safety checks from hot loops
- Use fast spatial data structure (octree/BVH)

```cpp
// Current: O(n) ray marching for each photon step
// Proposed: Precomputed voxel traversal using DDA algorithm
void traverse_voxels_DDA(const Ray& ray, double max_distance);
```

### 4. Modernize Random Number Generation

**Recommendation**: The current RNG is adequate but could be improved.

- PCG (Permuted Congruential Generator) for better performance
- SIMD-friendly generators for vectorization
- Better stratified sampling for variance reduction

### 4. Add Validation Tests

- Energy conservation check: `initial_energy == absorbed + transmitted + reflected`
- Limiting cases: zero absorption, zero scattering

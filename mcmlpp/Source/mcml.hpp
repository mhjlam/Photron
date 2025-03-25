/*******************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *  Copyright M.H.J. Lam, 2025.
 *	Monte Carlo simulation of light transport in multi-layered turbid mediums in C++20.
 ****
 *	Version 1.x:    10/1991.
 *	Version 2.0:    02/1996.
 *  Version 3.0:    03/2025.
 *
 *	Lihong Wang, Ph.D.
 *	Bioengineering Program
 *	Texas A&M University
 *	College Station, Texas
 *
 *	Steven L. Jacques, Ph.D.
 *	M.D. Anderson Cancer Center
 *	University of Texas
 *	Houston, Texas
 *
 *	Liqiong Zheng, B.S.
 *	Department of Computer Science
 *	University of Houston
 *	Houston, Texas
 *
 *  M.H.J. Lam, MSc.
 *  Graduate School of Natural Sciences
 *  Utrecht University
 *  Utrecht, Netherlands
 *
 *	This program was based on:
 *	(1) The Pascal code written by Marleen Keijzer and Steven L. Jacques in this
 *  laboratory in 1989, which deals with multi-layered turbid mediums.
 *
 *	(2) Algorithm for semi-infinite turbid name by S.A. Prahl, M. Keijzer,
 *  S.L. Jacques, A.J. Welch, SPIE Institute Series Vol. IS 5 (1989), and by
 *  A.N. Witt, The Astrophysical Journal Supplement Series 35, 1-6 (1977).
 *
 *	Major modifications in version 1.x include:
 *		- Conform to ANSI Standard C.
 *		- Removal of limit on number of array elements, because arrays in this
 *        program are dynamically allocated.
 *        This means that the program can accept any number of layers or grid
 *        lines as long as the memory permits.
 *		- Avoiding global variables whenever possible.
 *		- Grouping variables logically using structures.
 *		- Top-down design, keep each subroutine clear & short.
 *		- Reflectance and transmittance are angularly resolved.
 *	Major modifications in version 2.0 include:
 *		- Allow interactive input of parameters.
 *		- Allow to score various simulation quantities.
 *		- Time-resolved simulation.
 *		- Adjustable source position.
 *		- Support Isotropic source.
 *		- Simulation time control in addition to photon control.
 *		- Compute the standard errors of some physical quantities.
 *		- Allow continuation simulations to reduce standard errors.
 *	Major modifications in version 3.0 include:
 *      - Conform to C++20 and modern C++ coding standards.
 *      - Using standard library containers and algorithms wherever applicable.
 *      - Using standard library random number generation facilities
 *        (Mersenne Twister instead of Delayed Fibonacci generator).
 *      - Reduce usage of pointers and raw arrays.
 *      - Object-oriented design.
 *      - JSON input/output files.
 ****
 *	Dimension of depth: centimeters (cm).
 *  Dimension of angle: steradians (sr).
 *	Dimension of time: picoseconds (ps).
 *
 *  Ballistic quantities relate to photons that travel along a straight path
 *  without scattering.
 *
 ****/


#include "random.hpp"

#include <random>
#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <unordered_set>


template <typename T> using vec1 = std::vector<T>;
template <typename T> using vec2 = std::vector<std::vector<T>>;
template <typename T> using vec3 = std::vector<std::vector<std::vector<T>>>;

using unique_strings = std::unordered_set<std::string>;


// Roulette survival chance
constexpr double ROULETTE_SURVIVAL = 0.1;

// Speed of light in vacuum (C) [cm/ps]
constexpr double SPEED_OF_LIGHT = 0.0299792458;

// Inverse speed of light (1/C) [ps/cm]
constexpr double SPEED_OF_LIGHT_INV = 33.35640952;

// Split photon if true, otherwise statistical reflection.
constexpr bool PARTIAL_REFLECTION = false;

// If 1 - cos(theta) <= COS_ZERO_TOLERANCE, abs(theta) <= 1e-6 rad.
// If 1 + cos(theta) <= COS_ZERO_TOLERANCE, abs(PI - theta) <= 1e-6 rad.
constexpr double COS_0_TOLERANCE = 1.0E-12;

// If cos(theta) <= COS_90_TOLERANCE, theta >= PI/2 - 1e-6 rad.
constexpr double COS_90_TOLERANCE = 1.0E-6;


/********************************* Structures *********************************/

enum class RunType
{
    StartNew,                           // Start a new simulation
    Continue                            // Continue previous simulation
};

enum class BeamType
{
    Pencil,                             // Pencil beam
    Isotropic                           // Isotropic source
};

enum class ControlBit
{
    NumPhotons,                         // Photon number only
    TimeLimit,                          // Time limit only
    Both                                // Both photon number and time limit
};

enum class FileFormat
{
    Ascii,                              // Only ASCII is supported
    Binary
};

enum class IoMode
{
    Read,                               // Read mode
    Write                               // Write mode
};

enum class ScaleMode
{
    Scale,                              // Scale results
    Unscale                             // Unscale results
};

enum class PunchMode
{
    ResetTimer,                         // Reset timer
    TimeElapsed,                        // Return elapsed time since last reset
    TimeElapsedStr,                     // Return time elapsed as string as well
};


/*******************************************************************************
 *  Photon packet.
 ****/
struct Photon
{
    bool    alive{};		            //  Photon alive/terminated
    long    num_scatters{};	            //  Number of scatterings
    short   current_layer{};	        //  Index to layer where photon currently is

    double  x{}, y{}, z{};		        //  Cartesian coordinates [cm]
    double  ux{}, uy{}, uz{};	        //  Directional cosines of a photon
    double  weight{};		            //  Weight
    double  step_size{};		        //  Current step size [cm]
    double  step_size_left{};	        //  Step size left. dimensionless [-]
    double  flight_time{};	            //  Flight time [picosec]
};

/*******************************************************************************
 *  Specify which quantity/quantities should be scored.
 ****/
struct Record
{
    bool Rd_rat{ false };                 // Diffuse reflectance [1/(cm² sr ps)]
    bool Rd_ra{ false };                  // Diffuse reflectance [1/(cm² sr)]
    bool Rd_rt{ false };                  // Diffuse reflectance [1/sr ps]
    bool Rd_at{ false };                  // Diffuse reflectance [1/(cm² ps)]
    bool Rd_r{ false };                   // Diffuse reflectance [1/cm²]
    bool Rd_a{ false };                   // Diffuse reflectance [1/sr]
    bool Rd_t{ false };                   // Diffuse reflectance [1/ps]

    bool Td_rat{ false };                 // Diffuse transmittance [1/(cm² sr ps)]
    bool Td_ra{ false };                  // Diffuse transmittance [1/(cm² sr)]
    bool Td_rt{ false };                  // Diffuse transmittance [1/sr ps]
    bool Td_at{ false };                  // Diffuse transmittance [1/(cm² ps)]
    bool Td_r{ false };                   // Diffuse transmittance [1/cm²]
    bool Td_a{ false };                   // Diffuse transmittance [1/sr]
    bool Td_t{ false };                   // Diffuse transmittance [1/ps]

    bool A_rzt{ false };                  // Absorption [1/(cm² sr ps)]
    bool A_rz{ false };                   // Absorption [1/cm²]
    bool A_zt{ false };                   // Absorption [1/(cm² ps)]
    bool A_z{ false };                    // Absorption [1/cm²]
    bool A_t{ false };                    // Absorption [1/ps]
};

/*******************************************************************************
 *	Structure used to describe the geometry and optical properties of a layer.
 *	Top and bottom are the z coordinates for the upper boundary and lower
 *  boundary respectively.
 *
 *	cos_crit0 and cos_crit1 are the cosines of the critical angle of total
 *  internal reflection for the upper boundary and lower boundary respectively.
 *
 *	They are set to zero if no total internal reflection exists.
 *	They are used for computation speed.
 ****/
struct Layer
{
    std::string name{};                 // Name of the layer

    double      top_z{};                // Top z coordinate [cm] 
    double      bot_z{};                // Bottom z coordinate [cm]

    double      eta{};	                // Refractive index
    double      mua{};	                // Absorption coefficient [1/cm] 
    double      mus{};	                // Scattering coefficient [1/cm] 
    double      ani{};	                // Anisotropy

    double      cos_crit0{};
    double      cos_crit1{};
};

/*******************************************************************************
 *	Input parameters for each independent run.
 *
 *	z and r are for the cylindrical coordinate system. [cm]
 *	alpha denotes the angle between the photon exiting direction and the surface
 *  normal. [radian]
 *
 *	Member layers will point to an array of structures which store parameters of
 *  each layer. This array has (number_layers + 2) elements; one for each layer.
 *	The first and last layers are the top and bottom ambient layers respectively.
 ****/
struct RunParams
{
    std::string output_filename{};      // Output file name
    FileFormat  output_file_format{};   // Output file format
    ControlBit  control_bit{};          // Control of simulation termination.

    long        num_photons{};	        // Number of photons to be traced
    long        add_num_photons{};	    // Additional photon number
    long        time_limit{};	        // Computation time limit [sec]
    long        add_limit{};	        // Additional computation time

    short       seed{};		            // Random number seed (unused)
    short       num_runs{};		        // Number of runs
    double      weight_treshold{};		// Play roulette if photon weight is less than this threshold

    BeamType    source{};	            // Beam type
    double      source_z{};		        // Z coordinate of source
    short       source_layer{};		    // Put source in source_layer
    std::string source_medium_name{};   // Medium name of source layer

    double      grid_z{};		        // Z grid line separation [cm]
    double      grid_r{};		        // R grid line separation [cm]
    double      grid_a{};		        // Alpha grid line separation [rad]
    double      grid_t{};		        // Time grid line separation [ps]

    std::size_t num_z{};		        // Number of Z grid lines
    std::size_t num_r{};		        // Number of R grid lines
    std::size_t num_a{};		        // Number of alpha grid lines
    std::size_t num_t{};		        // Number of time grid lines

    double      max_z{};		        // Maximum z [cm]
    double      max_r{};		        // Maximum r [cm]
    double      max_alpha{};		    // Maximum alpha [rad]
    double      max_time{};		        // Maximum time [ps]

    Record      record;		            // Recorded quantities
    vec1<Layer> layers;                 // Layer parameters
    vec1<Layer> mediums;                // Medium parameters

    std::size_t num_layers{};           // Number of intermediate layers

    unique_strings unique_outputs;      // Unique output file names
};


struct Reflectance
{
    vec3<double> rat;   // Diffuse reflectance per unit area, per unit solid angle, per unit time [1/(cm² sr ps)]

    vec2<double> ra;    // Diffuse reflectance per unit area, per unit solid angle [1/(cm² sr)]
    vec2<double> rt;    // Diffuse reflectance per unit solid angle, per unit time [1/sr ps]
    vec2<double> at;    // Diffuse reflectance per unit area, per unit time [1/cm² ps]

    vec1<double> r;     // Diffuse reflectance per unit area [1/cm²]
    vec1<double> a;     // Diffuse reflectance per unit solid angle [1/sr]
    vec1<double> t;     // Diffuse reflectance per unit time [1/ps]

    double dr{};	    // Total diffuse reflectance
    double de{};	    // Standard error for diffuse reflectance
    double di{};	    // Diffuse reflectance of the i-th photon

    double sp{};	    // Specular reflectance

    double br{};	    // Ballistic reflectance
    double be{};	    // Standard error for ballistic reflectance
    double bi{};	    // Ballistic reflectance of the i-th photon
};

struct Transmittance
{
    vec3<double> rat;   // Diffuse transmittance per unit area, per unit solid angle, per unit time [1/(cm² sr ps)]

    vec2<double> ra;    // Diffuse transmittance per unit area, per unit solid angle [1/(cm² sr)]
    vec2<double> rt;    // Diffuse transmittance per unit solid angle, per unit time [1/sr ps]
    vec2<double> at;    // Diffuse transmittance per unit area, per unit time [1/cm² ps]

    vec1<double> r;	    // Diffuse reflectance per unit area [1/cm²] 
    vec1<double> a;	    // Diffuse reflectance per unit solid angle [1/sr] 
    vec1<double> t;	    // Diffuse reflectance per unit time [1/ps] 

    double dr{};	    // Total diffuse transmittance
    double de{};	    // Standard error for diffuse transmittance
    double di{};	    // Diffuse transmittance of the i-th photon

    double br{};	    // Ballistic transmittance
    double be{};	    // Standard error for ballistic transmission
    double bi{};	    // Ballistic transmittance of the i-th photon
};

struct Absorption
{
    vec3<double> rzt;   // Rate of absorption per unit volume, per unit time [1/(cm³ ps]

    vec2<double> rz;	// Rate of absorption per unit volume [1/cm³]
    vec2<double> zt;	// Rate of absorption per unit time [1/(cm ps)]

    vec1<double> z;	    // Absorption per unit depth [1/cm]
    vec1<double> t;	    // Absorption per unit time [1/ps]

    double ab{};		// Total absorption
    double ae{};		// Standard error for absorption
    double ai{};		// Absorption of the i-th photon

    vec2<double> bzt;	// Ballistic absorption per unit depth, per unit time [1/(cm ps)]
    vec1<double> bz;	// Ballistic absorption per unit depth [1/cm]
};

/*******************************************************************************
 *	Structures for scoring physical quantities.
 ****/
struct Tracer
{
    Absorption A;
    Reflectance R;
    Transmittance T;
};


class Random;
extern Random Rand;

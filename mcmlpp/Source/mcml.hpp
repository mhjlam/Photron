/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
 *	Monte Carlo simulation of light transport in
 *	multi-layered turbid mediums in ANSI Standard C.
 *
 *	This header file is shared by both MCML and CONV.
 ****
 *	Starting Date:	10/1991.
 *	Current Date:	02/1996.
 *
 *	Lihong Wang, Ph.D.
 *	Bioengineering Program
 *	Texas A&M University
 *	College Station, Texas 77843-3120
 *	USA
 *
 *	Steven L. Jacques, Ph.D.
 *	Laser Biology Research Laboratory - 17
 *	M.D. Anderson Cancer Center
 *	University of Texas
 *	1515 Holcombe Blvd.
 *	Houston, Texas 77030
 *	USA
 *
 *	Liqiong Zheng, B.S.
 *	Department of Computer Science
 *	University of Houston
 *	Houston, Texas
 *
 *	This program was based on:
 *	(1) The Pascal code written by Marleen Keijzer and
 *	Steven L. Jacques in this laboratory in 1989, which
 *	deals with multi-layered turbid mediums.
 *
 *	(2) Algorithm for semi-infinite turbid name by
 *	S.A. Prahl, M. Keijzer, S.L. Jacques, A.J. Welch,
 *	SPIE Institute Series Vol. IS 5 (1989), and by
 *	A.N. Witt, The Astrophysical Journal Supplement
 *	Series 35, 1-6 (1977).
 *
 *	Major modifications in version 1.x include:
 *		. Conform to ANSI Standard C.
 *		. Removal of limit on number of array elements,
 *		  because arrays in this program are dynamically
 *		  allocated. This means that the program can accept
 *		  any number of layers or gridlines as long as the
 *		  memory permits.
 *		. Avoiding global variables whenever possible.  This
 *		  program has not used global variables so far.
 *		. Grouping variables logically using structures.
 *		. Top-down design, keep each subroutine clear &
 *		  short.
 *		. Reflectance and transmittance are angularly
 *		  resolved.
 *	Major modifications in version 2.0 include:
 *		. Allow interactive input of parameters.
 *		. Allow to score various simulation quantities.
 *		. Time-resolved simulation.
 *		. Adjustable source position.
 *		. Support Isotropic source.
 *		. Simulation time control in addition to photon control.
 *		. Compute the standard errors of some physical quantities.
 *		. Allow continuation simulations to reduce standard errors.
 ****
 *	General Naming Conventions:
 *	Preprocessor names: all capital letters,
 *		e.aniso. #define PREPROCESSORS
 *	Globals: first letter of each word is capital, no
 *		underscores,
 *		e.aniso. short GlobalVar;
 *	Dummy variables:  first letter of each word is capital,
 *		and words are connected by underscores,
 *		e.aniso. void NiceFunction(char Dummy_Var);
 *	Local variables:  all lower cases, words may be connected
 *		by underscores,
 *		e.aniso. short local_var;
 *	Function names or data types:  same as Globals.
 ****
 *	Dimension of length: cm.
 *	Dimension of time:   ps.
 ****/


#include <random>
#include <string>
#include <vector>
#include <cstdint>


constexpr double    PI = 3.1415926;
constexpr double    CHANCE = 0.1;		        //  Chance of roulette survival. 
constexpr double    C_LIGHT = 0.0299792458;     //  Speed of light in vacuum [cm/ps]. 
constexpr double    ONE_OVER_C = 33.35640952;   //  1/C [ps/cm]. 

template <typename T> constexpr int sign(T x) { return (x >= 0) ? 1 : -1; }


/********************************* Stuctures **********************************/

enum class BeamType
{
    Pencil, Isotropic
};

enum class ControlBit
{
    NumPhotons, TimeLimit, Both
};

enum class FileFormat
{
    ASCII, Binary
};


/****
 *  Photon packet.
 ****/
struct Photon
{
    double  x, y, z;		//  Cartesian coordinates.[cm] 
    double  ux, uy, uz;	    //  Directional cosines of a photon. 
    double  min_weight;		//  Min weight. 
    bool    alive;		    //  Photon alive/terminated. 
    short   current_layer;	//  Index to current_layer where photon is.
    double  step_size;		//  Current step size. [cm]. 
    double  step_size_left;	//  Step size left. dimensionless [-]. 
    long    num_scatters;	//  Number of scatterings. 
    double  flight_time;	//  Flight time [picosec]. 
};

/****
 *   Specify which quantity is to be scored.
 *
 *   Data categories:
 *   Rd_rat                  Td_rat                  A_rzt
 *   Rd_ra   Rd_rt   Rd_at   Td_ra   Td_rt   Rd_at   A_rz    A_zt
 *   Rd_r    Rd_a    Rd_t    Td_r    Td_a    Td_t    A_z     A_t
 ****/
struct Record
{
    //  use bit field to save space.
    int Rd_rat : 1;
    int Rd_ra : 1;
    int Rd_rt : 1;
    int Rd_at : 1;
    int Rd_r : 1;
    int Rd_a : 1;
    int Rd_t : 1;

    int Td_rat : 1;
    int Td_ra : 1;
    int Td_rt : 1;
    int Td_at : 1;
    int Td_r : 1;
    int Td_a : 1;
    int Td_t : 1;

    int A_rzt : 1;
    int A_rz : 1;
    int A_zt : 1;
    int A_z : 1;
    int A_t : 1;
};

/****
 *	Structure used to describe the geometry and optical
 *	properties of a current_layer.
 *	top_z and bot_z are the z coordinates for the upper boundary
 *	and lower boundary respectively.
 *
 *	cos_crit0 and cos_crit1 are the cosines of the
 *	critical angle of total internal reflection for the
 *	upper boundary and lower boundary respectively.
 *	They are set to zero if no total internal reflection
 *	exists.
 *	They are used for computation speed.
 ****/
struct Layer
{
    std::string name    = "";  // Name of the current_layer
    double top_z        = 0.0;  // Top z coordinate  [cm] 
    double bot_z        = 0.0;  // Bottom z coordinate [cm]
    double eta          = 0.0;	// Refractive index of a current_layer
    double mua          = 0.0;	// Absorption coefficient. [1/cm] 
    double mus          = 0.0;	// Scattering coefficient. [1/cm] 
    double aniso        = 0.0;	// Anisotropy
    double cos_crit0    = 0.0;
    double cos_crit1    = 0.0;
};

/****
 *	Input parameters for each independent run.
 *
 *	z and r are for the cylindrical coordinate system. [cm]
 *	a is for the angle alpha between the photon exiting direction and the surface normal. [radian]
 *
 *	The grid line separations in z, r, and alpha directions are grid_z, grid_r, and grid_alpha respectively. 
 *  The numbers of grid lines in z, r, and alpha directions are num_z, num_r, and num_alpha respectively.
 *
 *	The member layers will point to an array of structures which store parameters of each layer.
 *	This array has (number_layers + 2) elements. One element is for a layer.
 *	The layers 0 and (num_layers + 1) are for top ambient name and the bottom ambient name respectively.
 *
 * 	Output-file format:
 *  Use 'A' for ASCII, and 'B' for binary.
 *
 *	Control bit:
 *      1 - photon number only.
 *		2 - time limit only.
 *		3 - both.
 ****/
struct RunParams
{
    std::string         output_filename;    // Output file name. 
    FileFormat          output_file_format; // Output file format. 

    long                num_photons;	    // Number of photons to be traced.
    long                add_num_photons;	// Additional photon number. 
    long                time_limit_seconds;	// Computation time limit. 
    long                add_limit_seconds;	// Additional computation time. 
    ControlBit          control_bit;        // Control of simulation termination. 
    double              min_weight;		    // Play roulette if photon weight < min_weight. 

    BeamType            source;	            // Beam type. 
    double              source_z;		    // Z coordinate of source. 
    short               source_layer;		// Put source in source_layer. 
    std::string         source_medium_name; // Medium name of source_layer. 

    double              grid_z;		        // Z grid separation.[cm] 
    double              grid_r;		        // R grid separation.[cm] 
    double              grid_alpha;		    // Alpha grid separation. [rad] 
    double              grid_time;		    // Time grid separation.[ps] 

    short               num_z;		        // Array range 0..num_z-1. 
    short               num_r;		        // Array range 0..num_r-1. 
    short               num_alpha;		    // Array range 0..num_alpha-1. 
    short               num_time;		    // Array range 0..num_time-1. 

    double              max_z;		        // Maximum z [cm] 
    double              max_r;		        // Maximum r [cm] 
    double              max_alpha;		    // Maximum alpha [rad]
    double              max_time;		    // Maximum time [ps] 

    std::vector<Layer>  layers;             // Layer parameters
    std::vector<Layer>  mediums;            // Medium parameters
    Record              record;		        // Recorded quantities

    short               num_runs;		    // Number of runs
};


struct Reflectance
{
    // Diffuse reflectance. [1/(cm2 sr ps)] 
    std::vector<std::vector<std::vector<double>>> rat;

    std::vector<std::vector<double>> ra;    // [1/(cm2 sr)]
    std::vector<std::vector<double>> rt;    // [1/sr ps]
    std::vector<std::vector<double>> at;    // [1/cm2 ps]

    std::vector<double> r;                  // [1/cm2]
    std::vector<double> a;                  // [1/sr]
    std::vector<double> t;                  // [1/ps]

    double d;	                            // Total diffuse reflectance.
    double de;	                            // Standard error for diffuse reflectance.
    double di;	                            // Diffuse reflectance of the i-th photon.

    double b;	                            // Ballistic reflectance.
    double be;	                            // Standard error for ballistic reflectance.
    double bi;	                            // Ballistic reflectance of the i-th photon.

    double sp;	                            // Specular reflectance.
};

struct Transmittance
{
    // Diffuse transmittance. [1/(cm2 sr ps)] 
    std::vector<std::vector<std::vector<double>>> rat;

    std::vector<std::vector<double>> ra;    // [1/(cm2 sr)]
    std::vector<std::vector<double>> rt;    // [1/sr ps]
    std::vector<std::vector<double>> at;    // [1/cm2 ps]

    std::vector<double> r;	                // [1/cm2] 
    std::vector<double> a;	                // [1/sr] 
    std::vector<double> t;	                // [1/ps] 

    double d;	                            // Total diffuse transmittance.
    double de;	                            // Standard error for Td.
    double di;	                            // Td of the i-th photon.

    double b;	                            // Ballistic transmittance.
    double be;	                            // Standard error for Tb.
    double bi;	                            // Tb of the i-th photon.
};

struct Absorption
{
    // Absorption. [1/(cm3 psabsorption]
    std::vector<std::vector<std::vector<double>>> rzt;

    std::vector<std::vector<double>> rz;	    // [1/cm3] 
    std::vector<std::vector<double>> zt;	    // [1/(cm ps)] 

    std::vector<double> z;	                    // [1/cm] 
    std::vector<double> t;	                    // [1/ps] 

    std::vector<std::vector<double>> bzt;	    // Ballistic absorption. [1/(cm ps)] 
    std::vector<double> bz;	                    // [1/cm] 

    double a;		                            // Total absorption.
    double e;		                            // Standard error for A.
    double i;		                            // A of the i-th photon.
};

/****
 *	Structures for scoring physical quantities.
 *	z and r represent z and r coordinates of the
 *	cylindrical coordinate system. [cm]
 *	a is the angle alpha between the photon exiting
 *	direction and the normal to the surfaces. [radian]
 *	R represents reflectance, T transmission, A absorption.
 *
 *	See comments of RunParams.
 *	See manual for the physical quantities.
 ****/
struct Tracer
{
    Absorption A;
    Reflectance R;
    Transmittance T;
};


extern std::mt19937 RandomEngine;
extern std::uniform_real_distribution<double> Distribution;

static double RandomNumber()
{
    return Distribution(RandomEngine);
}

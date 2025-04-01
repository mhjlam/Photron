#pragma once

#include <string>
#include <memory>
#include <vector>
#include <istream>
#include <variant>
#include <iostream>
#include <optional>

#include "mcml.hpp"


using double_or_string = std::variant<double, std::string>;
using opt_double_or_string = std::optional<std::variant<double, std::string>>;

using double2 = std::tuple<double, double>;
using double3 = std::tuple<double, double, double>;
using double4 = std::tuple<double, double, double, double>;
using double5 = std::tuple<double, double, double, double, double>;

using OutputFile = std::pair<std::string, FileFormat>;


class Random;


class Reader
{
public:
    Reader(std::string filename, std::string_view version = MCI_VERSION);
    ~Reader() = default;

    // Read the input parameters for all runs and count number of runs.
    virtual void ReadParams(std::istream& input, RunParams& params);

    // Read the mediums list.
    virtual vec1<Layer> ReadMediums(std::istream& input);

    // Read the input name and the input format.
    virtual std::string ReadOutput(std::istream& input);

    // Read the parameters of all layers.
    virtual vec1<Layer> ReadLayers(std::istream& input, RunParams& params);

    // Read the beam source type (Pencil or Isotropic) and starting position.
    virtual LightSource ReadSource(std::istream& input, RunParams& params);

    // Read the grid separation parameters (z, r, t) and number of grid lines (z, r, t, and alpha).
    virtual Grid ReadGrid(std::istream& input);

    // Read which quantity is to be scored.
    virtual Record ReadRecord(std::istream& input, RunParams& params);

    // Read the number of photons and computation time limit.
    virtual Target ReadTarget(std::istream& input, RunParams& params, bool add = false);

    // Read the weight threshold.
    virtual double ReadWeight(std::istream& input);


    // Read the seed for random number generator (unused).
    long ReadSeed(std::istream& input);

    // Check whether the input version is the same as version.
    bool ReadVersion(std::istream& input, const std::string_view& version);

    // Read in the input parameters for one run.
    void ReadRunParams(std::istream& input, RunParams& params);

    // Read and restore the status of random number generater from previous output.
    void ReadRandomizer(std::istream& input, std::shared_ptr<Random>& random);

    // Read result back from a output file.
    Radiance ReadRadiance(std::istream& input, RunParams& params, std::shared_ptr<Random>& random);

    void SkipLine(std::istream& input, std::size_t num_lines = 1);


protected:
    // Skip space or comment lines and return a data line.
    std::string readNextLine(std::istream& input);

    std::vector<double_or_string> extract(const std::string& input,
                                          const std::vector<double_or_string>& expected,
                                          std::string parse_err, bool allow_opt = false);

    // Check consistancy of input parameters.
    bool checkInputParams(RunParams& params);

    // Return string as uppercase
    std::string& toUpperCase(std::string& string);

private:
    // Diffuse reflectance per unit area, per unit solid angle, per unit time [1/(cm² sr ps)]
    vec3<double> ReadR_rat(std::istream& input, std::size_t Nr, std::size_t Na, std::size_t Nt);

    // Diffuse reflectance per unit area, per unit solid angle [1/(cm² sr)]
    vec2<double> ReadR_ra(std::istream& input, std::size_t Nr, std::size_t Na);

    // Diffuse reflectance per unit solid angle, per unit time [1/sr ps]
    vec2<double> ReadR_rt(std::istream& input, std::size_t Nr, std::size_t Nt);

    // Diffuse reflectance per unit area, per unit time [1/cm² ps]
    vec2<double> ReadR_at(std::istream& input, std::size_t Na, std::size_t Nt);

    // Diffuse reflectance distribution per unit area [1/cm²]
    vec1<double> ReadR_r(std::istream& input, std::size_t Nr);

    // Diffuse reflectance per unit solid angle [1/sr]
    vec1<double> ReadR_a(std::istream& input, std::size_t Na);

    // Diffuse reflectance per unit time [1/ps]
    vec1<double> ReadR_t(std::istream& input, std::size_t Nt);

    // Diffuse transmittance per unit area, per unit solid angle, per unit time [1/(cm² sr ps)]
    vec3<double> ReadT_rat(std::istream& input, std::size_t Nr, std::size_t Na, std::size_t Nt);

    // Diffuse transmittance per unit area, per unit solid angle [1/(cm² sr)]
    vec2<double> ReadT_ra(std::istream& input, std::size_t Nr, std::size_t Na);

    // Diffuse transmittance per unit solid angle, per unit time [1/sr ps]
    vec2<double> ReadT_rt(std::istream& input, std::size_t Nr, std::size_t Nt);

    // Diffuse transmittance per unit area, per unit time [1/cm² ps]
    vec2<double> ReadT_at(std::istream& input, std::size_t Na, std::size_t Nt);

    // Diffuse reflectance per unit area [1/cm²]
    vec1<double> ReadT_r(std::istream& input, std::size_t Nr);

    // Diffuse reflectance per unit solid angle [1/sr]
    vec1<double> ReadT_a(std::istream& input, std::size_t Na);

    // Diffuse reflectance per unit time [1/ps]
    vec1<double> ReadT_t(std::istream& input, std::size_t Nt);

    // Rate of absorption per unit volume, per unit time [1/(cm³ ps]
    vec3<double> ReadA_rzt(std::istream& input, std::size_t Nr, std::size_t Nz, std::size_t Nt);

    // Rate of absorption per unit volume [1/cm³]
    vec2<double> ReadA_rz(std::istream& input, std::size_t Nr, std::size_t Nz);

    // Rate of absorption per unit time [1/(cm ps)]
    vec2<double> ReadA_zt(std::istream& input, std::size_t Nz, std::size_t Nt);

    // Absorption per unit depth [1/cm]
    vec1<double> ReadA_z(std::istream& input, std::size_t Nz);

    // Absorption per unit time [1/ps]
    vec1<double> ReadA_t(std::istream& input, std::size_t Nt);

    // Ballistic absorption per unit depth, per unit time [1/(cm ps)]
    vec2<double> ReadAb_zt(std::istream& input, std::size_t Nz, std::size_t Nt);

    // Ballistic absorption per unit depth [1/cm]
    vec1<double> ReadAb_z(std::istream& input, std::size_t Nz);


public:
    operator std::istream& () { return *m_input; }

protected:
    std::unique_ptr<std::istream> m_input;
};

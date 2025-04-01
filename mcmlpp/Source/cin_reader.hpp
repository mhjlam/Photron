#pragma once

#include <tuple>
#include <string>

#include "mcml.hpp"
#include "reader.hpp"

class CinReader : public Reader
{
public:
    CinReader();
    ~CinReader() = default;

    // Read the input parameters for all runs and count number of runs.
    void ReadParams(std::istream& input, RunParams& params) override;

    // Read the mediums list.
    vec1<Layer> ReadMediums(std::istream& input) override;

    // Read the input name and the input format.
    std::string ReadOutput(std::istream& input) override;

    // Read the parameters of all layers.
    vec1<Layer> ReadLayers(std::istream& input, RunParams& params) override;

    // Read the beam source type (Pencil or Isotropic) and starting position.
    LightSource ReadSource(std::istream& input, RunParams& params) override;

    // Read the grid separation parameters (z, r, t) and number of grid lines (z, r, t, and alpha).
    Grid ReadGrid(std::istream& input) override;

    // Read which quantity is to be scored.
    Record ReadRecord(std::istream& input, RunParams& params) override;

    // Read the number of photons and computation time limit.
    Target ReadTarget(std::istream& input, RunParams& params, bool add = false) override;

    // Read the weight threshold.
    double ReadWeight(std::istream& input) override;

protected:
    template <typename T> T read(std::istream& input, std::string err_msg);
    template <typename T, typename U> std::tuple<T, U> read(std::istream& input, std::string err_msg);
};

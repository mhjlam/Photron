#include "cin_reader.hpp"

#include <format>
#include <fstream>
#include <numbers>
#include <iostream>
#include <algorithm>


CinReader::CinReader() : Reader{ std::cin }
{
}


void CinReader::ReadParams(std::istream& input, RunParams& params)
{
    params.mediums = ReadMediums(*m_input);
    params.output_filename = ReadOutput(*m_input);
    params.layers = ReadLayers(*m_input, params);
    params.source = ReadSource(*m_input, params);
    params.grid = ReadGrid(*m_input);
    params.record = ReadRecord(*m_input, params);
    params.target = ReadTarget(*m_input, params);
    params.weight_threshold = ReadWeight(*m_input);
}

template <typename T>
T CinReader::read(std::istream& input, std::string err_msg)
{
    T value{};

    do {
        try {
            std::string in = readNextLine(input);
            auto extracted = extract(in, { T{} }, err_msg);
            return std::get<T>(extracted[0]);
        }
        catch (const std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
        }
    } while (true);
}

template <typename T, typename U>
std::tuple<T, U> CinReader::read(std::istream& input, std::string err_msg)
{
    std::tuple<T, U> values{};

    do {
        try {
            std::string in = readNextLine(input);
            auto extracted = extract(in, { T{}, U{} }, err_msg);
            values = std::make_tuple(extracted[0], extracted[1]);
        }
        catch (const std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
        }
    } while (true);
}

vec1<Layer> CinReader::ReadMediums(std::istream& input)
{
    std::cout << "Specify medium list. Total number of mediums: ";

    std::size_t num_media = static_cast<std::size_t>(read<double>(input, "Invalid number of mediums."));

    // Allocate space for the layer parameters.
    vec1<Layer> mediums;
    mediums.resize(num_media);

    for (std::size_t i = 0; i < num_media; i++) {
        std::string name;
        double eta{};
        double mua{};
        double mus{};
        double g{};

        bool duplicate = false;

        std::cout << "Specify medium " << i + 1 << ": " << std::endl << "  Medium name: ";
        do { // TODO: make conditional part of read function
            name = read<std::string>(input, "Invalid medium name or duplicate.");
            duplicate = std::ranges::any_of(mediums, [&](const Layer& l) {
                return l.name == name;
            });
        } while (duplicate);

        std::cout << "  Refractive index n (>= 1.0): ";
        do {
            eta = read<double>(input, "Invalid refractive index. Input again (>= 1.0): ");
        } while (eta < 1.0);

        std::cout << "  Absorption coefficient mua (>= 0.0 /cm): ";
        do {
            mua = read<double>(input, "Invalid absorption coefficient. Input again (>= 0.0): ");
        } while (mua < 0.0);

        std::cout << "  Scattering coefficient mus (>= 0.0 /cm): ";
        do {
            mus = read<double>(input, "Invalid scattering coefficient. Input again (>= 0.0): ");
        } while (mus < 0.0);

        std::cout << "  Anisotropy factor g (0.0 - 1.0): ";
        do {
            g = read<double>(input, "Invalid anisotropy factor. Input again (0.0 - 1.0): ");
        } while (g < 0.0 && g > 1.0);

        mediums[i] = { i, name, eta, mua, mus, g };

        std::cout << std::endl;
    }

    return mediums;
}

std::string CinReader::ReadOutput(std::istream& input)
{
    std::string file_name;
    std::string file_mode;

    do {
        std::cout << "Specify output filename with extension .mco: ";
        std::getline(std::cin, file_name);
        file_mode[0] = 'w';

        // input exists.
        std::ifstream file(file_name, std::ios::in);
        if (!file.is_open()) {
            std::cout << "File " << file_name << " exists, w=overwrite, n=new filename: ";

            // avoid null line.
            do {
                std::getline(std::cin, file_mode);
            } while (file_mode.empty());

            file.close();
        }
    } while (file_mode[0] != 'w');

    std::cout << std::endl;

    return file_name;
}

vec1<Layer> CinReader::ReadLayers(std::istream& input, RunParams& params)
{
    auto make_layer = [&](double z0, double z1) -> Layer {
        std::string name;
        std::size_t index{};
        bool name_exists{ false };

        do {
            name = read<std::string>(input, "Invalid medium name. Input again: ");
            name_exists = std::ranges::any_of(params.mediums, [&](const Layer& medium) {
                if (medium.name == name) {
                    index = medium.index;
                }
                return medium.name == name;
            });
        } while (!name_exists);

        return Layer{ index,
            params.mediums[index].name,
            params.mediums[index].eta,
            params.mediums[index].mu_a,
            params.mediums[index].mu_s,
            params.mediums[index].g,
            z0, z1, 0.0, 0.0
        };
    };

    std::cout << std::endl << "Specify layer list. Available medium types:" << std::endl;
    for (int i = 0, j = 1; i < params.mediums.size(); i++, j++) {
        std::cout << std::format("{:<16}", params.mediums[i].name);
        if (j % 4 == 0) {
            std::cout << std::endl;
        }
    }

    std::size_t num_layers{};
    std::cout << "  Total number of layers: ";
    do {
        num_layers = static_cast<std::size_t>(read<double>(input, "Invalid layer number. Input again: "));
    } while (num_layers < 1);


    vec1<Layer> layers;
    layers.resize(num_layers);

    // z coordinate of the current layer
    double z{ 0.0 };

    // First and last layers are for ambient
    for (std::size_t i = 0; i < num_layers; ++i) {
        Layer layer{};

        if (i == 0) {
            std::cout << std::endl << "  Name of upper ambient medium: ";
            layer = make_layer(0.0, 0.0);
        }
        else if (i == num_layers-1) {
            std::cout << std::endl << "  Name of lower ambient medium: ";
            layer = make_layer(z, z);
        }
        else {
            std::cout << "  Thickness of layer " << i << " (thickness > 0.0 cm): ";
            double thickness = read<double>(input, "Invalid thickness. Input again (thickness > 0.0 cm): ");

            std::cout << std::endl << "  Name of layer " << i << ": ";
            layer = make_layer(z, z + thickness);

            z += thickness;
        }

        layers[i] = layer;
    }

    std::cout << std::endl;

    return layers;
}

LightSource CinReader::ReadSource(std::istream& input, RunParams& params)
{
    // Compute the index to layer according to the z coordinate. 
    // If the z is on an interface between layers, the returned index will point to the upper layer.
    // Layer 0 is the top ambient and layer num_layers+1 is the bottom ambient layer.
    auto layer_index = [&](double z, RunParams& params) -> std::size_t {
        for (std::size_t i = 1; i <= params.num_layers; i++) {
            if (z >= params.layers[i].z0 && z <= params.layers[i].z1) {
                return i;
            }
        }
        return std::numeric_limits<std::size_t>::max();
    };

    LightSource source{};

    std::cout << "Input source type (P = pencil / I = isotropic): ";

    char ch;
    do {
        std::cin.get(ch);
        char ch_upper = std::toupper(ch);

        if (ch_upper == 'P') {
            source.beam = BeamType::Pencil;
            break;
        }
        else if (ch_upper == 'I') {
            source.beam = BeamType::Isotropic;
            break;
        }
        else {
            std::cout << "Invalid type. Input again (P = pencil / I = isotropic): ";
        }
    } while (true);

    std::cout << std::endl;

    std::cout << "Input the z coordinate of source (0.0 - " << params.layers[params.num_layers-1].z1 << " cm) and the medium" << std::endl;
    std::cout << "where the source is if the z is on an interface (e.g. 1.0 [air]): ";

    try {
        std::string in = readNextLine(input);
        auto extracted = extract(in, { double{}, std::string{} }, "Invalid starting position of photon source.", true);

        if (extracted.size() == 1) {
            source.z = std::get<double>(extracted[0]);
            source.layer_index = layer_index(source.z, params);

            if (source.layer_index == std::numeric_limits<std::size_t>::max()) {
                throw std::runtime_error("Invalid starting position of photon source.");
            }
        }
        else if (extracted.size() == 2) {
            source.z = std::get<double>(extracted[0]);
            source.layer_index = layer_index(source.z, params);

            if (source.layer_index == std::numeric_limits<std::size_t>::max()) {
                throw std::runtime_error("Invalid starting position of photon source.");
            }

            source.medium_name = std::get<std::string>(extracted[1]);

            if (params.layers[source.layer_index].name == source.medium_name) {
                if ((std::fabs(source.z - params.layers[source.layer_index].z1) < std::numeric_limits<double>::epsilon()) && 
                    (params.layers[source.layer_index + 1].name == source.medium_name)) {
                    source.layer_index++;
                    if (source.layer_index > params.num_layers) {
                        throw std::runtime_error("Source is outside of the last layer.");
                    }
                }
                else {
                    throw std::runtime_error("Medium name and z coordinate do not match.");
                }
            }
        }

        if (params.source.beam == BeamType::Isotropic && source.z == 0.0) {
            throw std::runtime_error("Can not put isotropic source in upper ambient medium.");
        }
    }
    catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
    }

    std::cout << std::endl;

    return source;
}

Grid CinReader::ReadGrid(std::istream& input)
{
    std::cout << "Specify dz, dr, dt in one line (all > 0.0 cm, e.g., 0.1 0.1 0.1): ";

    std::string buf;

    buf = readNextLine(input);
    auto extracted_d = extract(buf, { double{}, double{}, double{} }, "Error reading dz, dr, dt.");

    double step_z = std::get<double>(extracted_d[0]);
    double step_r = std::get<double>(extracted_d[1]);
    double step_t = std::get<double>(extracted_d[2]);

    if (step_z <= 0) { throw std::runtime_error("Nonpositive dz."); }
    if (step_r <= 0) { throw std::runtime_error("Nonpositive dr."); }
    if (step_t <= 0) { throw std::runtime_error("Nonpositive dt."); }

    std::cout << "Specify nz, nr, nt, na in one line (all > 0, e.g., 100 100 100 100): ";

    buf = readNextLine(input);
    auto extracted_n = extract(buf, { double{}, double{}, double{}, double{} }, "Error reading number of dz, dr, dt, da.");

    std::size_t num_z = static_cast<std::size_t>(std::get<double>(extracted_n[0]));
    std::size_t num_r = static_cast<std::size_t>(std::get<double>(extracted_n[1]));
    std::size_t num_t = static_cast<std::size_t>(std::get<double>(extracted_n[2]));
    std::size_t num_a = static_cast<std::size_t>(std::get<double>(extracted_n[3]));

    if (num_z <= 0) { throw std::runtime_error("Nonpositive number of dz."); }
    if (num_r <= 0) { throw std::runtime_error("Nonpositive number of dr."); }
    if (num_t <= 0) { throw std::runtime_error("Nonpositive number of dt."); }
    if (num_a <= 0) { throw std::runtime_error("Nonpositive number of da."); }

    double step_a = 0.5 * std::numbers::pi / num_a;

    std::cout << std::endl;

    return Grid{
        .step_z = step_z,
        .step_r = step_r,
        .step_a = step_a,
        .step_t = step_t,
        .num_z = num_z,
        .num_r = num_r,
        .num_t = num_t,
        .num_a = num_a,
        .max_z = step_z * num_z,
        .max_r = step_r * num_r,
        .max_a = step_a * num_a,
        .max_t = step_t * num_t
    };
}

Record CinReader::ReadRecord(std::istream& input, RunParams& params)
{
    std::cout << "Select scored quantities from the following data categories:" << std::endl;
    std::cout << "\tR_rat\t\t\tT_rat\t\t\tA_rzt" << std::endl;
    std::cout << "\tR_ra\tR_rt\tR_at\tT_ra\tT_rt\tR_at\tA_rz\tA_zt" << std::endl;
    std::cout << "\tR_r\tR_a\tR_t\tT_r\tT_a\tT_t\tA_z\tA_t" << std::endl;

    Record record = Reader::ReadRecord(input, params);
    
    std::cout << std::endl;
    
    return record;
}

Target CinReader::ReadTarget(std::istream& input, RunParams& params, bool add)
{
    if (add) {
        std::cout << params.target.num_photons << " photons have been traced in the previous simulation.\n";
        std::cout << "Specify additional photons or compution time in hh:mm format,\n";
        std::cout << "or both in one line (e.g. 10000 5:30): ";
    }
    else {
        std::cout << "Specify number of photons or time in hh:mm format,\n";
        std::cout << "or both in one line (e.g. 10000 5:30): ";
    }

    Target target = Reader::ReadTarget(input, params, add);
    std::cout << std::endl;

    return target;
}

double CinReader::ReadWeight(std::istream& input)
{
    std::cout << "Input weight threshold (0 <= W < 1.0, 0.0001 recommended): ";

    bool valid = false;
    double weight{};

    do {
        weight = read<double>(input, "Invalid weight treshold. Input again (0 <= W < 1.0): ");
        valid = weight >= 0.0 && weight < 1.0;
    } while (!valid);

    std::cout << std::endl;

    return weight;
}

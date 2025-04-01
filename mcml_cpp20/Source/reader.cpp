#include "reader.hpp"

#include <format>
#include <limits>
#include <ranges>
#include <fstream>
#include <numbers>
#include <sstream>
#include <algorithm>
#include <string_view>

#include "mcml.hpp"
#include "random.hpp"


using namespace std::literals;


static constexpr bool is_double(const double_or_string& v)
{
    return std::holds_alternative<double>(v);
}

static constexpr bool is_string(const double_or_string& v)
{
    return std::holds_alternative<std::string>(v);
}


Reader::Reader(std::string filename, std::string_view version) : m_input{ std::make_unique<std::istream>(nullptr) }
{
    if (filename.empty()) {
        m_input = std::make_unique<std::istream>(std::cin.rdbuf());
    }
    else {
        auto file_stream = std::make_unique<std::ifstream>(filename);

        if (!file_stream->is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }
        m_input = std::move(file_stream);

        if (!ReadVersion(*m_input, version)) {
            throw std::runtime_error("Invalid file version.");
        }
    }
}


void Reader::ReadParams(std::istream& input, RunParams& params)
{
    auto endOfRuns = [&]() {
        // Found end of runs
        bool end_found = false;

        // Record input position
        std::streampos file_pos = input.tellg();

        std::string buf = readNextLine(input);
        if (buf.empty()) {
            end_found = true;
            std::cout << "Missing end." << std::endl;
        }
        else if (buf.find("end") != std::string::npos) {
            end_found = true;
        }

        // Restore postion
        input.seekg(file_pos, std::ios::beg);
        return end_found;
    };

    // Read list of mediums
    params.mediums = ReadMediums(input);

    // Save current position in input file.
    std::streampos file_pos = input.tellg();

    short run_index = 0;
    do {
        ++run_index;
        std::cout << "Checking input data for run " << run_index << std::endl;

        // Read the input parameters for one run
        ReadRunParams(input, params);

        // Attempt insert and detect duplicate output filenames
        if (!params.unique_output_filenames.insert(params.output_filename).second) {
            throw std::runtime_error("File name " + params.output_filename + " duplicated.");
            std::exit(1);
        }
    } while (!endOfRuns());

    params.num_runs = run_index;

    // Restore file position
    input.seekg(file_pos, std::ios::beg);
}

bool Reader::ReadVersion(std::istream& input, const std::string_view& version)
{
    // Find next data line
    std::string line = readNextLine(input);

    if (line.find(version) == std::string::npos) {
        std::cerr << "Invalid file version.";
        return false;
    }
    return true;
}

vec1<Layer> Reader::ReadMediums(std::istream& input)
{
    // Get current output position
    std::streampos file_pos = input.tellg();

    // Find number of mediums
    std::size_t num_mediums = 0;

    std::string buf;
    while (buf.find("end") == std::string::npos) {
        try {
            buf = readNextLine(input);
        }
        catch (...) {
            throw std::runtime_error("No media found.");
        }

        if (buf.find("end") != std::string::npos) {
            break;
        }
        else if (buf.empty()) {
            throw std::runtime_error("Missing end.");
        }
        else {
            num_mediums++;

            auto extracted = extract(buf, 
                { std::string{}, double{}, double{}, double{}, double{} }, 
                "Error reading medium parameters.");

            std::string name = std::get<std::string>(extracted[0]);
            double n = std::get<double>(extracted[1]);
            double mu_a = std::get<double>(extracted[2]);
            double mu_s = std::get<double>(extracted[3]);
            double g = std::get<double>(extracted[4]);

            // Verify optical parameters
            if (n <= 0.0 || mu_a < 0.0 || mu_s < 0.0 || g < -1.0 || g > 1.0) {
                throw std::runtime_error("Bad optical parameters in " + name);
            }
        }
    }

    if (num_mediums < 1) {
        throw std::runtime_error("No media found.");
    }

    // Seek to previous output position
    input.seekg(file_pos, std::ios::beg);

    vec1<Layer> mediums;
    mediums.resize(num_mediums);

    for (std::size_t i = 0; i < num_mediums; i++) {
        std::string buf = readNextLine(input);
        auto extracted = extract(buf, { std::string{}, double{}, double{}, double{}, double{} }, "Error reading medium parameters.");

        mediums[i].index = i;
        mediums[i].name = std::get<std::string>(extracted[0]);
        mediums[i].eta = std::get<double>(extracted[1]);
        mediums[i].mu_a = std::get<double>(extracted[2]);
        mediums[i].mu_s = std::get<double>(extracted[3]);
        mediums[i].g = std::get<double>(extracted[4]);

        if (mediums[i].eta <= 0.0 || 
            mediums[i].mu_a < 0.0 || 
            mediums[i].mu_s < 0.0 || 
            mediums[i].g < -1.0 || 
            mediums[i].g > 1.0) {
            throw std::runtime_error("Bad optical parameters in " + mediums[i].name);
        }
    }

    // Skip the signal end
    readNextLine(input);

    return mediums;
}

std::string Reader::ReadOutput(std::istream& input)
{
    std::string buf = readNextLine(input);
    auto extracted = extract(buf, { std::string{} }, "Error reading file name.");
    return std::get<std::string>(extracted[0]);
}

Grid Reader::ReadGrid(std::istream& input)
{
    std::string buf = readNextLine(input);
    auto extracted_step = extract(buf, { double{}, double{}, double{} }, "Error reading dz, dr, dt.");

    double step_z = std::get<double>(extracted_step[0]);
    double step_r = std::get<double>(extracted_step[1]);
    double step_t = std::get<double>(extracted_step[2]);

    if (step_z <= 0) { throw std::runtime_error("Nonpositive dz."); }
    if (step_r <= 0) { throw std::runtime_error("Nonpositive dr."); }
    if (step_t <= 0) { throw std::runtime_error("Nonpositive dt."); }

    buf = readNextLine(input);
    auto extracted_d = extract(buf, { double{}, double{}, double{}, double{} }, "Error reading number of dz, dr, dt, da.");

    std::size_t num_z = static_cast<std::size_t>(std::get<double>(extracted_d[0]));
    std::size_t num_r = static_cast<std::size_t>(std::get<double>(extracted_d[1]));
    std::size_t num_t = static_cast<std::size_t>(std::get<double>(extracted_d[2]));
    std::size_t num_a = static_cast<std::size_t>(std::get<double>(extracted_d[3]));

    if (num_z <= 0) { throw std::runtime_error("Nonpositive number of dz."); }
    if (num_r <= 0) { throw std::runtime_error("Nonpositive number of dr."); }
    if (num_t <= 0) { throw std::runtime_error("Nonpositive number of dt."); }
    if (num_a <= 0) { throw std::runtime_error("Nonpositive number of da."); }

    double step_a = 0.5 * std::numbers::pi / num_a;

    return Grid {
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

Record Reader::ReadRecord(std::istream& input, RunParams& params)
{
    Record record{};

    std::string buf = readNextLine(input);
    if (buf.empty()) {
        throw std::runtime_error("Error reading scored quantities.");
    }

    std::string string;
    std::stringstream iss(buf);

    do {
        iss >> string;

        // Stop when comment is found
        if (string.starts_with("#")) {
            break;
        }

        // Trim and uppercase
        string = std::format("{:}", string);
        string = toUpperCase(string);

        if      (string == "R_R"sv)     { record.R_r = true; }
        else if (string == "R_A"sv)     { record.R_a = true; }
        else if (string == "R_RA"sv)    { record.R_ra = true; }
        else if (string == "R_T"sv)     { record.R_t = true; }
        else if (string == "R_RT"sv)    { record.R_rt = true; }
        else if (string == "R_AT"sv)    { record.R_at = true; }
        else if (string == "R_RAT"sv)   { record.R_rat = true; }
        else if (string == "T_R"sv)     { record.T_r = true; }
        else if (string == "T_A"sv)     { record.T_a = true; }
        else if (string == "T_RA"sv)    { record.T_ra = true; }
        else if (string == "T_T"sv)     { record.T_t = true; }
        else if (string == "T_RT"sv)    { record.T_rt = true; }
        else if (string == "T_AT"sv)    { record.T_at = true; }
        else if (string == "T_RAT"sv)   { record.T_rat = true; }
        else if (string == "A_Z"sv)     { record.A_z = true; }
        else if (string == "A_RZ"sv)    { record.A_rz = true; }
        else if (string == "A_T"sv)     { record.A_t = true; }
        else if (string == "A_ZT"sv)    { record.A_zt = true; }
        else if (string == "A_RZT"sv)   { record.A_rzt = true; }

        else {
            throw std::runtime_error("Unknown quantity: "s + string);
        }
    } while (!iss.fail() && !string.empty());

    // Check for conflicting records
    if (record.R_rat) {
        record.R_ra = record.R_rt = record.R_at = record.R_r = record.R_a = record.R_t = false;
    }
    if (record.R_ra) {
        record.R_r = record.R_a = false;
    }
    if (record.R_rt) {
        record.R_r = record.R_t = false;
    }
    if (record.R_at) {
        record.R_a = record.R_t = false;
    }
    if (record.T_rat) {
        record.T_ra = record.T_rt = record.T_at = record.T_r = record.T_a = record.T_t = false;
    }
    if (record.T_ra) {
        record.T_r = record.T_a = false;
    }
    if (record.T_rt) {
        record.T_r = record.T_t = false;
    }
    if (record.T_at) {
        record.T_a = record.T_t = false;
    }
    if (record.A_rzt) {
        record.A_rz = record.A_zt = record.A_z = record.A_t = false;
    }
    if (record.A_rz) {
        record.A_z = false;
    }
    if (record.A_zt) {
        record.A_z = record.A_t = false;
    }
    if (record.A_zt) {
        record.A_z = record.A_t = false;
    }

    return record;
}

double Reader::ReadWeight(std::istream& input)
{
    std::string buf = readNextLine(input);
    auto extracted = extract(buf, { double{} }, "Error reading threshold weight.");

    double weight = std::get<double>(extracted[0]);

    if (weight < 0 || weight >= 1.0) {
        throw std::runtime_error("Threshold weight out of range (0-1).");
    }

    return weight;
}

long Reader::ReadSeed(std::istream& input)
{
    std::string buf = readNextLine(input);
    auto extracted = extract(buf, { double{} }, "Error reading seed value.");

    long seed = static_cast<long>(std::get<double>(extracted[0]));

    if (seed < 0) {
        throw std::runtime_error("Negative seed value.");
    }

    return seed;
}

vec1<Layer> Reader::ReadLayers(std::istream& input, RunParams& params)
{
    std::string name;
    double thickness = 0.0;

    // Z coordinate of the current layer
    double z = 0.0;

    // Save current output position
    std::streampos file_pos = input.tellg();

    std::string buf;
    std::size_t num_layers = 0;

    // While "end" has not been found
    do {
        try {
            buf = readNextLine(input);

            if (buf.find("end") != std::string::npos) {
                break;
            }
        }
        catch (...) {
            throw std::runtime_error("No layers found.");
        }

        if (buf.empty()) {
            throw std::runtime_error("Missing end.");
        }
        else {
            // Read layer name
            extract(buf, { std::string{} }, "Error reading layer name.");
            num_layers++;
        }
    } while (true);

    if (num_layers < 3) {
        throw std::runtime_error("No layers found.");
    }

    // Seek to previous output position 
    input.seekg(file_pos, std::ios::beg);

    // First and last layers are for ambient
    vec1<Layer> layers;
    layers.resize(num_layers);

    for (std::size_t i = 0; i < num_layers; ++i) {
        std::string buf = readNextLine(input);

        // Top and bottom layers (get only name)
        if (i == 0 || i == num_layers-1) {
            // Get name only
            auto extracted = extract(buf, { std::string{} }, "Error reading layer specs.");
            name = std::get<std::string>(extracted[0]);
        }
        else {
            // Get name and thickness
            auto extracted = extract(buf, { std::string{}, double{} }, "Error reading layer specs.");

            name = std::get<std::string>(extracted[0]);
            thickness = std::get<double>(extracted[1]);

            if (thickness <= 0.0) {
                std::cerr << "Nonpositive layer thickness." << std::endl;
                return {};
            }
        }

        // Attempt to find medium by name
        // TODO: Refactor
        bool found = false;
        int medium_i;
        for (short i = 0; i < params.mediums.size(); i++) {
            if (name == params.mediums[i].name) {
                found = true;
                medium_i = i;
                break;
            }
        }

        if (!found) {
            std::cerr << "  Invalid medium name. " << std::endl;
            return {};
        }

        layers[i].name = params.mediums[medium_i].name;
        layers[i].eta = params.mediums[medium_i].eta;
        layers[i].mu_a = params.mediums[medium_i].mu_a;
        layers[i].mu_s = params.mediums[medium_i].mu_s;
        layers[i].g = params.mediums[medium_i].g;

        // Intermediate layers
        if (i != 0 && i != (num_layers + 1)) {
            layers[i].z0 = z;
            z += thickness;
            layers[i].z1 = z;
        }
        // Top and bottom layers
        else {
            layers[i].z0 = z;
            layers[i].z1 = z;
        }
    }

    // Skip the signal "end" of layers.
    readNextLine(input);

    return layers; // NOTE: num_layers = layers.size() - 2
}

Target Reader::ReadTarget(std::istream& input, RunParams& params, bool add)
{
    Target target = params.target;

    std::string buf = readNextLine(input);
    auto extracted = extract(buf, { double{}, std::string{} }, "Error reading number of photons or time limit.", true);

    if (extracted.size() == 1 && is_double(extracted[0])) {
        int num_photons = static_cast<int>(std::get<double>(extracted[0]));

        if (num_photons > 0) {
            target.control_bit = ControlBit::NumPhotons;

            if (add) {
                target.add_num_photons = num_photons;
            }
            else {
                target.num_photons = num_photons;
            }
        }
        else {
            throw std::runtime_error("Nonpositive number of photons.");
        }
    }
    else if (extracted.size() == 1 && is_string(extracted[0])) {
        std::string time = std::get<std::string>(extracted[0]);

        int hours = 0;
        int minutes = 0;
        char separator;

        std::istringstream iss(time);
        iss >> hours >> separator >> minutes;

        if (!iss || separator != ':') {
            throw std::runtime_error("Invalid time limit format.");
        }

        if ((hours * 3600 + minutes * 60) > 0) {
            target.control_bit = ControlBit::TimeLimit;

            if (add) {
                target.add_time_limit = hours * 3600 + minutes * 60;
            }
            else {
                target.time_limit = hours * 3600 + minutes * 60;
            }
        }
        else {
            throw std::runtime_error("Nonpositive time limit.");
        }
    }
    else if (extracted.size() == 2 && is_double(extracted[0]) && is_string(extracted[1])) {
        int num_photons = static_cast<int>(std::get<double>(extracted[0]));
        std::string time = std::get<std::string>(extracted[1]);

        int hours = 0;
        int minutes = 0;
        char separator;

        std::istringstream iss(time);
        iss >> hours >> separator >> minutes;

        if (!iss || separator != ':') {
            throw std::runtime_error("Invalid time limit format.");
        }

        if (num_photons > 0 && (hours * 3600 + minutes * 60) >= 0) {
            target.control_bit = ControlBit::Both;

            if (add) {
                target.add_num_photons = num_photons;
                target.add_time_limit = hours * 3600 + minutes * 60;
            }
            else {
                target.num_photons = num_photons;
                target.time_limit = hours * 3600 + minutes * 60;
            }
        }
        else {
            throw std::runtime_error("Nonpositive number of photons or time limit.");
        }
    }
    else {
        throw std::runtime_error("Invalid number of photons or time limit.");
    }

    return target;
}

LightSource Reader::ReadSource(std::istream& input, RunParams& params)
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

    std::string buf = readNextLine(input);
    auto extracted_source = extract(buf, { std::string{} }, "Error reading photon source type.");

    std::string source_type = std::get<std::string>(extracted_source[0]);

    if (toUpperCase(source_type) == "PENCIL"sv) {
        source.beam = BeamType::Pencil;
    }
    else if (toUpperCase(source_type) == "ISOTROPIC"sv) {
        source.beam = BeamType::Isotropic;
    }
    else {
        throw std::runtime_error("Unknow photon source type.");
    }

    buf = readNextLine(input);
    auto extracted_start = extract(buf, { double{}, std::string{} }, "Invalid starting position of photon source.", true);

    if (extracted_start.size() == 1) {
        source.z = std::get<double>(extracted_start[0]);
        source.layer_index = layer_index(source.z, params);

        if (source.layer_index == std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error("Invalid starting position of photon source.");
        }
    }
    else if (extracted_start.size() == 2) {
        source.z = std::get<double>(extracted_start[0]);
        source.medium_name = std::get<std::string>(extracted_start[1]);

        if (source.medium_name[0] != '#' && source.medium_name[0] != '\n') {
            source.layer_index = layer_index(source.z, params);

            if (source.layer_index == std::numeric_limits<std::size_t>::max()) {
                throw std::runtime_error("Invalid starting position of photon source.");
            }

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
    }

    if (source.beam == BeamType::Isotropic && source.z == 0.0) {
        throw std::runtime_error("Can not put isotropic source in upper ambient medium.");
    }

    return source;
}

void Reader::ReadRunParams(std::istream& input, RunParams& params)
{
    try {
        params.output_filename = ReadOutput(input);
        params.layers = ReadLayers(input, params);
        params.num_layers = params.layers.size() - 2;
        params.source = ReadSource(input, params);
        params.grid = ReadGrid(input);
        params.record = ReadRecord(input, params);
        params.target = ReadTarget(input, params);
        params.weight_threshold = ReadWeight(input);
        params.seed = ReadSeed(input);

        // Compute the critical angles for total internal reflection according to the 
        // relative refractive index of the layer.
        for (short i = 1; i <= params.num_layers; i++) {
            double eta_0 = params.layers[i - 1].eta;
            double eta_1 = params.layers[i].eta;
            double eta_2 = params.layers[i + 1].eta;

            params.layers[i].cos_theta_c0 = eta_1 > eta_0 ? std::sqrt(1.0 - eta_0 * eta_0 / (eta_1 * eta_1)) : 0.0;
            params.layers[i].cos_theta_c1 = eta_1 > eta_2 ? std::sqrt(1.0 - eta_2 * eta_2 / (eta_1 * eta_1)) : 0.0;
        }

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        std::exit(1);
    }
}

void Reader::ReadRandomizer(std::istream& input, std::shared_ptr<Random>& random)
{
    std::string buf;
    std::vector<std::mt19937::result_type> status(624);

    do {
        std::getline(input, buf);
    } while (buf[0] != '#');

    for (int i = 0; i < status.size(); i++) {
        input >> status[i];
    }

    // Restore the status
    random->restore_state(status);
}

Radiance Reader::ReadRadiance(std::istream& input, RunParams& params, std::shared_ptr<Random>& random)
{
    Radiance radiance;
    ReadRandomizer(input, random);

    // skip comment line
    std::string buf = readNextLine(input);

    buf = readNextLine(input);
    std::istringstream iss(buf);
    iss >> radiance.R_spec;

    buf = readNextLine(input);
    iss = std::istringstream(buf);
    iss >> radiance.Rb_total >> radiance.Rb_error;

    buf = readNextLine(input);
    iss = std::istringstream(buf);
    iss >> radiance.R_total >> radiance.R_error;

    buf = readNextLine(input);
    iss = std::istringstream(buf);
    iss >> radiance.A_total >> radiance.A_error;

    buf = readNextLine(input);
    iss = std::istringstream(buf);
    iss >> radiance.Tb_total >> radiance.Tb_error;

    buf = readNextLine(input);
    iss = std::istringstream(buf);
    iss >> radiance.T_total >> radiance.T_error;


    // Reflectance
    if (params.record.R_rat) { radiance.R_rat = ReadR_rat(input, params.grid.num_r, params.grid.num_a, params.grid.num_t); }
    if (params.record.R_ra) { radiance.R_ra = ReadR_ra(input, params.grid.num_r, params.grid.num_a); }
    if (params.record.R_rt) { radiance.R_rt = ReadR_rt(input, params.grid.num_r, params.grid.num_t); }
    if (params.record.R_at) { radiance.R_at = ReadR_at(input, params.grid.num_a, params.grid.num_t); }
    if (params.record.R_r) { radiance.R_r = ReadR_r(input, params.grid.num_r); }
    if (params.record.R_a) { radiance.R_a = ReadR_a(input, params.grid.num_a); }
    if (params.record.R_t) { radiance.R_t = ReadR_t(input, params.grid.num_t); }

    // Transmittance
    if (params.record.T_rat) { radiance.T_rat = ReadT_rat(input, params.grid.num_r, params.grid.num_a, params.grid.num_t); }
    if (params.record.T_ra) { radiance.T_ra = ReadT_ra(input, params.grid.num_r, params.grid.num_a); }
    if (params.record.T_rt) { radiance.T_rt = ReadT_rt(input, params.grid.num_r, params.grid.num_t); }
    if (params.record.T_at) { radiance.T_at = ReadT_at(input, params.grid.num_a, params.grid.num_t); }
    if (params.record.T_r) { radiance.T_r = ReadT_r(input, params.grid.num_r); }
    if (params.record.T_a) { radiance.T_a = ReadT_a(input, params.grid.num_a); }
    if (params.record.T_t) { radiance.T_t = ReadT_t(input, params.grid.num_t); }

    // Absorption
    if (params.record.A_rzt) {
        radiance.Ab_zt = ReadAb_zt(input, params.grid.num_z, params.grid.num_t);
        radiance.A_rzt = ReadA_rzt(input, params.grid.num_r, params.grid.num_z, params.grid.num_t);
    }
    if (params.record.A_rz) {
        radiance.Ab_z = ReadAb_z(input, params.grid.num_z);
        radiance.A_rz = ReadA_rz(input, params.grid.num_r, params.grid.num_z);
    }
    if (params.record.A_zt) { radiance.A_zt = ReadA_zt(input, params.grid.num_z, params.grid.num_t); }
    if (params.record.A_z) { radiance.A_z = ReadA_z(input, params.grid.num_z); }
    if (params.record.A_t) { radiance.A_t = ReadA_t(input, params.grid.num_t); }

    return radiance;
}

void Reader::SkipLine(std::istream& input, std::size_t num_lines)
{
    for (std::size_t i = 0; i < num_lines; i++) {
        readNextLine(input);
    }
}


std::string Reader::readNextLine(std::istream& input)
{
    std::string line;
    while (std::getline(input, line)) {
        // Find first non-whitespace character
        auto it = std::ranges::find_if(line, [](char c) {
            return !std::isspace(c);
        });

        // Skip whitespace-only lines
        if (it == line.end()) {
            continue;
        }

        // Skip comment lines
        if (*it == '#') {
            continue;
        }

        // Return data line
        return line;
    }

    // No datalines found
    throw std::runtime_error("No data line found");
}

bool Reader::checkInputParams(RunParams& params)
{
    for (int i = 0; i <= params.num_layers + 1; i++) {
        // Find index of the medium name in the medium list.
        auto it = std::ranges::find_if(params.mediums, [&](const Layer& m) {
            return std::ranges::any_of(params.layers, [&](const Layer& l) {
                return l.name == m.name;
            });
        });

        if (it == params.mediums.end()) {
            return std::cerr << "Invalid medium name of layer " << i << ".\n", 0;
        }

        std::size_t index = std::distance(params.mediums.begin(), it);
        params.layers[i].eta = params.mediums[index].eta;
        params.layers[i].mu_a = params.mediums[index].mu_a;
        params.layers[i].mu_s = params.mediums[index].mu_s;
        params.layers[i].g = params.mediums[index].g;
    }

    if ((params.source.beam == BeamType::Isotropic) && (params.source.z == 0.0)) {
        return std::cerr << "Can not put isotropic source in upper ambient medium.\n", 0;
    }

    // Find the index of the layer according to the z coordinate.
    if (params.source.z < 0.0) {
        return std::cerr << "Nonpositive z coordinate.\n", false;
    }
    else if (params.source.z > params.layers.back().z1) {
        return std::cerr << "Source is outside of the last layer.\n", false;
    }
    else {
        params.source.layer_index =
            std::ranges::lower_bound(params.layers, params.source.z, {}, &Layer::z1) -
            params.layers.begin();
    }

    // Check the medium name and z coordinate of the source.
    if (params.source.medium_name[0] != '\0') {
        if (params.layers[params.source.layer_index].name == params.source.medium_name) {
            if ((std::abs(params.source.z - params.layers[params.source.layer_index].z1) < std::numeric_limits<double>::epsilon()) &&
                (params.layers[params.source.layer_index + 1].name == params.source.medium_name)) {
                params.source.layer_index++;
            }
            else {
                std::cerr << "Medium name and z coordinate do not match." << std::endl;
                return false;
            }
        }
    }
    return true;
}

std::vector<double_or_string> Reader::extract(const std::string& input, const std::vector<double_or_string>& expected, std::string parse_err, bool allow_opt)
{
    auto parse = [&](const std::string& str, const double_or_string& type) -> opt_double_or_string {
        std::istringstream iss(str);
        if (is_double(type)) {
            double value;
            if (iss >> value) { return value; }
        }
        else if (is_string(type)) {
            return str;
        }
        return std::nullopt;
    };

    std::vector<double_or_string> results;
    std::istringstream iss(input);
    std::string token;

    for (const auto& type : expected) {
        if (!(iss >> token)) { break; }
        opt_double_or_string value = parse(token, type);
        if (value.has_value()) {
            results.push_back(value.value());
        }
        else if (allow_opt) {
            continue;
        }
        else {
            throw std::runtime_error(parse_err);
        }
    }

    if (!allow_opt && results.size() != expected.size()) {
        throw std::runtime_error(parse_err);
    }

    return results;
}

std::string& Reader::toUpperCase(std::string& string)
{
    std::ranges::transform(string, string.begin(), [](unsigned char ch) {
        return std::toupper(ch);
    });
    return string;
}


vec3<double> Reader::ReadR_rat(std::istream& input, std::size_t Nr, std::size_t Na, std::size_t Nt)
{
    readNextLine(input);
    vec3<double> R_rat(Nr, vec2<double>(Na, vec1<double>(Nt)));

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t ia = 0; ia < Na; ia++) {
            for (std::size_t it = 0; it < Nt; it++) {
                input >> R_rat[ir][ia][it];
            }
        }
    }
    return R_rat;
}

vec2<double> Reader::ReadR_ra(std::istream& input, std::size_t Nr, std::size_t Na)
{
    readNextLine(input);
    vec2<double> R_ra(Nr, vec1<double>(Na));

    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t ia = 0; ia < Na; ia++) {
            input >> R_ra[ir][ia];
        }
    }
    return R_ra;
}

vec2<double> Reader::ReadR_rt(std::istream& input, std::size_t Nr, std::size_t Nt)
{
    readNextLine(input);
    vec2<double> R_rt(Nr, vec1<double>(Nt));

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t it = 0; it < Nt; it++) {
            input >> R_rt[ir][it];
        }
    }
    return R_rt;
}

vec2<double> Reader::ReadR_at(std::istream& input, std::size_t Na, std::size_t Nt)
{
    readNextLine(input);
    vec2<double> R_at(Na, vec1<double>(Nt));

    std::size_t i = 0;
    for (std::size_t ia = 0; ia < Na; ia++) {
        for (std::size_t it = 0; it < Nt; it++) {
            input >> R_at[ia][it];
        }
    }
    return R_at;
}

vec1<double> Reader::ReadR_r(std::istream& input, std::size_t Nr)
{
    readNextLine(input);
    vec1<double> R_r(Nr);

    for (std::size_t ir = 0; ir < Nr; ir++) {
        input >> R_r[ir];
    }
    return R_r;
}

vec1<double> Reader::ReadR_a(std::istream& input, std::size_t Na)
{
    readNextLine(input);
    vec1<double> R_a(Na);

    for (std::size_t ia = 0; ia < Na; ia++) {
        input >> R_a[ia];
    }
    return R_a;
}

vec1<double> Reader::ReadR_t(std::istream& input, std::size_t Nt)
{
    readNextLine(input);
    vec1<double> R_t(Nt);

    for (std::size_t it = 0; it < Nt; it++) {
        input >> R_t[it];
    }
    return R_t;
}

vec3<double> Reader::ReadT_rat(std::istream& input, std::size_t Nr, std::size_t Na, std::size_t Nt)
{
    readNextLine(input);
    vec3<double> T_rat(Nr, vec2<double>(Na, vec1<double>(Nt)));

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t ia = 0; ia < Na; ia++) {
            for (std::size_t it = 0; it < Nt; it++) {
                input >> T_rat[ir][ia][it];
            }
        }
    }
    return T_rat;
}

vec2<double> Reader::ReadT_ra(std::istream& input, std::size_t Nr, std::size_t Na)
{
    readNextLine(input);
    vec2<double> T_ra(Nr, vec1<double>(Na));

    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t ia = 0; ia < Na; ia++) {
            input >> T_ra[ir][ia];
        }
    }
    return T_ra;
}

vec2<double> Reader::ReadT_rt(std::istream& input, std::size_t Nr, std::size_t Nt)
{
    readNextLine(input);
    vec2<double> T_rt(Nr, vec1<double>(Nt));

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t it = 0; it < Nt; it++) {
            input >> T_rt[ir][it];
        }
    }
    return T_rt;
}

vec2<double> Reader::ReadT_at(std::istream& input, std::size_t Na, std::size_t Nt)
{
    readNextLine(input);
    vec2<double> T_at(Na, vec1<double>(Nt));

    std::size_t i = 0;
    for (std::size_t ia = 0; ia < Na; ia++) {
        for (std::size_t it = 0; it < Nt; it++) {
            input >> T_at[ia][it];
        }
    }
    return T_at;
}

vec1<double> Reader::ReadT_r(std::istream& input, std::size_t Nr)
{
    readNextLine(input);
    vec1<double> T_r(Nr);

    for (std::size_t ir = 0; ir < Nr; ir++) {
        input >> T_r[ir];
    }
    return T_r;
}

vec1<double> Reader::ReadT_a(std::istream& input, std::size_t Na)
{
    readNextLine(input);
    vec1<double> T_a(Na);

    for (std::size_t ia = 0; ia < Na; ia++) {
        input >> T_a[ia];
    }
    return T_a;
}

vec1<double> Reader::ReadT_t(std::istream& input, std::size_t Nt)
{
    readNextLine(input);
    vec1<double> T_t(Nt);

    for (std::size_t it = 0; it < Nt; it++) {
        input >> T_t[it];
    }
    return T_t;
}

vec3<double> Reader::ReadA_rzt(std::istream& input, std::size_t Nr, std::size_t Nz, std::size_t Nt)
{
    readNextLine(input);
    vec3<double> A_rzt(Nr, vec2<double>(Nz, vec1<double>(Nt)));

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t iz = 0; iz < Nz; iz++) {
            for (std::size_t it = 0; it < Nt; it++) {
                input >> A_rzt[ir][iz][it];
            }
        }
    }
    return A_rzt;
}

vec2<double> Reader::ReadA_rz(std::istream& input, std::size_t Nr, std::size_t Nz)
{
    readNextLine(input);
    vec2<double> A_rz(Nr, vec1<double>(Nz));

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t iz = 0; iz < Nz; iz++) {
            input >> A_rz[ir][iz];
        }
    }
    return A_rz;
}

vec2<double> Reader::ReadA_zt(std::istream& input, std::size_t Nz, std::size_t Nt)
{
    readNextLine(input);
    vec2<double> A_zt(Nz, vec1<double>(Nt));

    std::size_t i = 0;
    for (std::size_t iz = 0; iz < Nz; iz++) {
        for (std::size_t it = 0; it < Nt; it++) {
            input >> A_zt[iz][it];
        }
    }
    return A_zt;
}

vec1<double> Reader::ReadA_z(std::istream& input, std::size_t Nz)
{
    readNextLine(input);
    vec1<double> A_z(Nz);

    for (std::size_t iz = 0; iz < Nz; iz++) {
        input >> A_z[iz];
    }
    return A_z;
}

vec1<double> Reader::ReadA_t(std::istream& input, std::size_t Nt)
{
    readNextLine(input);
    vec1<double> A_t(Nt);

    for (std::size_t it = 0; it < Nt; it++) {
        input >> A_t[it];
    }
    return A_t;
}

vec2<double> Reader::ReadAb_zt(std::istream& input, std::size_t Nz, std::size_t Nt)
{
    readNextLine(input);
    vec2<double> Ab_zt(Nz, vec1<double>(Nt));

    std::size_t i = 0;
    for (std::size_t iz = 0; iz < Nz; iz++) {
        for (std::size_t it = 0; it < Nt; it++) {
            input >> Ab_zt[iz][it];
        }
    }
    return Ab_zt;
}

vec1<double> Reader::ReadAb_z(std::istream& input, std::size_t Nz)
{
    readNextLine(input);
    vec1<double> Ab_z(Nz);

    for (std::size_t iz = 0; iz < Nz; iz++) {
        input >> Ab_z[iz];
    }
    return Ab_z;
}

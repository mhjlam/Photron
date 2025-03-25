/*******************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *  Copyright M.H.J. Lam, 2025.
 *	Input/output of data.
 ****/


#include "mcml.hpp"

#include <format>
#include <ranges>
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string_view>

#include <numbers>
#include <variant>
#include <optional>


using namespace std::literals;

using SpecValue = std::variant<double, std::string>;
using OptSpecValue = std::optional<SpecValue>;


/*******************************************************************************
 *	Parse a line of input to extract a list of double and string values.
 ****/
std::vector<SpecValue> ParseLine(std::string& input, std::string err, std::size_t min_expected = 1)
{
    std::istringstream stream(input);
    std::vector<SpecValue> values;

    std::string token;
    while (stream >> token) {
        // Attempt to parse as double
        try {
            std::size_t pos;
            double val = std::stod(token, &pos);
            if (pos == token.length()) {
                values.emplace_back(val);
                continue;
            }
        }
        catch (...) {}

        // Interpret as a string instead
        if (!token.empty()) {
            values.emplace_back(token);
        }
    }

    if (values.empty() || values.size() < min_expected) {
        std::cerr << err << std::endl;
        return {};
    }

    return values;
}

OptSpecValue ParseInput(std::string err)
{
    std::string input;
    std::getline(std::cin, input);

    // Attempt to parse as double
    try {
        std::size_t pos;
        double val = std::stod(input, &pos);
        if (pos == input.length()) {
            return val;
        }
    }
    catch (...) {}

    // Interpret as a string instead
    if (!input.empty()) {
        return input;
    }

    std::cerr << err << std::endl;
    return {};
}

bool IsString(SpecValue& val)
{
    return std::holds_alternative<std::string>(val);
}

bool IsNumber(SpecValue& val)
{
    return std::holds_alternative<double>(val);
}

int GetIntValue(OptSpecValue& opt_int, int& target_int)
{
    if (std::holds_alternative<double>(opt_int.value())) {
        int val = static_cast<int>(std::get<double>(opt_int.value()));
        target_int = val;
        return val;
    }
    return 0;
}

double GetDoubleValue(OptSpecValue& opt_double, double& target_double)
{
    if (std::holds_alternative<double>(opt_double.value())) {
        double val = std::get<double>(opt_double.value());
        target_double = val;
        return val;
    }
    return 0.0;
}

std::string GetStringValue(OptSpecValue& opt_string, std::string& target_string)
{
    if (std::holds_alternative<std::string>(opt_string.value())) {
        std::string val = std::get<std::string>(opt_string.value());
        target_string = val;
        return val;
    }
    return {};
}

std::string& ToUpper(std::string& string)
{
    std::ranges::transform(string, string.begin(), [](unsigned char ch) {
        return std::toupper(ch);
    });
    return string;
}


/*******************************************************************************
 *	Print messages about MCML.
 ****/
void AboutMCML()
{
    std::cout << "MCML 3.0, Copyright (c) 1992-1996, 2025" << std::endl;
    std::cout << "Monte Carlo Simulation of Light Transport in Multi-Layered Turbid Media" << std::endl;

    std::cout << std::endl;
    std::cout << "Lihong Wang, Ph.D." << std::endl;
    std::cout << "Bioengineering Program, Texas A&M University" << std::endl;
    std::cout << "College Station, Texas, USA" << std::endl;

    std::cout << "Liqiong Zheng, B.S." << std::endl;
    std::cout << "Dept. of Computer Science," << std::endl;
    std::cout << "University of Houston, Texas, USA." << std::endl;

    std::cout << "Steven L. Jacques, Ph.D." << std::endl;
    std::cout << "Oregon Medical Laser Center, Providence/St. Vincent Hospital" << std::endl;
    std::cout << "Portland, Oregon, USA" << std::endl;

    std::cout << "M.H.J. Lam, MSc." << std::endl;
    std::cout << "Utrecht University" << std::endl;
    std::cout << "Utrecht, Netherlands" << std::endl;

    std::cout << std::endl;
    std::cout << "Obtain the original program from omlc.org/software/mc" << std::endl;

    std::cout << std::endl;
    std::cout << "Please cite the following article in your publications:" << std::endl;
    std::cout << "\tL.-H. Wang, S. L. Jacques, and L.-Q. Zheng, MCML - Monte " << std::endl;
    std::cout << "\tCarlo modeling of photon transport in multi-layered" << std::endl;
    std::cout << "\ttissues, Computer Methods and Programs in Biomedicine, 47," << std::endl;
    std::cout << "\t131-146 (1995)" << std::endl;
}

/*******************************************************************************
 *	Eliminate the chars in a string which are not printing chars or spaces.
 *	Spaces include ' ', '\f', '\t' etc.
 *	Return 1 if no non-printing chars found, otherwise return 0.
 ****/
bool CheckChar(std::string& string)
{
    bool found_non_printing = false;

    string.erase(std::remove_if(string.begin(), string.end(), [&](unsigned char ch) {
        // Keep printable and whitespace characters
        if (std::isprint(ch) || std::isspace(ch)) {
            return false;
        }

        // Remove non-printing characters
        found_non_printing = true;
        return true;
    }), string.end());

    return found_non_printing ? 0 : 1;
}

/*******************************************************************************
 *	Return 1 if this line is a comment line in which the first non-space
 *  character is "#" or a space line. Return 0 otherwise.
 ****/
bool CommentLine(std::string& buf)
{
    auto it = std::ranges::find_if(buf, [](char c) {
        return !std::isspace(c);
    });
    return (it == buf.end() || *it == '#');
}

/*******************************************************************************
 *	Skip space or comment lines and return a data line.
 ****/
std::string FindDataLine(std::istream& file)
{
    std::string line;
    while (std::getline(file, line)) {
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
    return {};
}

/*******************************************************************************
 *	Check whether the input version is the same as version.
 ****/
bool CheckFileVersion(std::fstream& file, const std::string& version)
{
    // Skip comment lines.
    std::string line;
    do { std::getline(file, line); } while (line.empty() || CommentLine(line));

    if (line.find(version) == std::string::npos) {
        std::cerr << "Invalid file version.";
        return false;
    }
    return true;
}

/*******************************************************************************
 *  Get a filename and open it for reading, retry until the input can be opened
 *  with a correct version or when a '.' is typed.
 ****/
bool GetFile(std::string& filename, const std::string& version, std::fstream& file)
{
    while (1) {
        std::cout << "Specify filename (or . to quit to main menu):";

        // Read input buffer
        std::getline(std::cin, filename);

        if (!filename.empty()) {
            // Terminate with a period
            if (filename.size() == 1 && filename[0] == '.') {
                // Return if '.' entered
                return false;
            }

            // Open the input & check the version
            file = std::fstream(filename, std::ios::in);
            if (!file.is_open()) {
                // Cannot open the input
                std::cerr << "File does not exist.";
            }
            else {
                if (CheckFileVersion(file, version)) {
                    return true;
                }
                else {
                    file.close();
                }
            }
        }
    }
    return false;
}

/*******************************************************************************
 *  Find number of mediums in the list and check the optical parameters.
 ****/
int FindNumMedia(std::fstream& file)
{
    short num_media = 0;

    while (1) {
        std::string buf = FindDataLine(file);

        if (buf.empty()) {
            std::cerr << "Missing end." << std::endl;
            return 0;
        }
        else if (buf.find("end") != std::string::npos) {
            break;
        }
        else {
            num_media++;

            auto extracted = ParseLine(buf, "Error reading number of mediums.", 5);
            if (extracted.empty()) {
                return 0;
            }

            std::string name = std::get<std::string>(extracted[0]);
            double n = std::get<double>(extracted[1]);
            double mua = std::get<double>(extracted[2]);
            double mus = std::get<double>(extracted[3]);
            double g = std::get<double>(extracted[4]);

            if (n <= 0 || mua < 0 || mus < 0 || g < -1 || g > 1) {
                std::cerr << "Bad optical parameters in " << name << std::endl;
                return 0;
            }
        }
    }

    return num_media;
}

/*******************************************************************************
 *  Read the parameters of one medium, assumming the parameters have been
 *  checked.
 ****/
bool ReadMedium(std::fstream& file, std::vector<Layer>& media)
{
    std::string buf = FindDataLine(file);
    auto extracted = ParseLine(buf, "Error reading medium parameters.", 5);
    if (extracted.empty()) {
        return false;
    }

    std::string name = std::get<std::string>(extracted[0]);
    double eta = std::get<double>(extracted[1]);
    double mua = std::get<double>(extracted[2]);
    double mus = std::get<double>(extracted[3]);
    double ani = std::get<double>(extracted[4]);

    // Add new medium to the list.
    media.push_back(Layer{ .name = name, .eta = eta, .mua = mua, .mus = mus, .ani = ani });
    return true;
}

/************************************************************************************
 *  Read the mediums list.
 ****/
bool ReadMediumListQ(std::fstream& file, RunParams& params)
{
    // Get current output position
    std::streampos file_pos = file.tellg();

    int num_media = FindNumMedia(file);
    if (num_media < 1) {
        file.seekg(file_pos, std::ios::beg);
        return false;
    }

    // Seek to previous output position
    file.seekg(file_pos, std::ios::beg);

    for (short i = 0; i < num_media; i++) {
        ReadMedium(file, params.mediums);
    }

    // skip the signal end.
    FindDataLine(file);

    return true;
}

/*******************************************************************************
 *	Read the input name and the input format.
 ****/
bool ReadFileName(std::fstream& file, RunParams& params)
{
    std::string buf = FindDataLine(file);
    auto extracted = ParseLine(buf, "Error reading file name.", 1);
    if (extracted.empty()) {
        return false;
    }

    params.output_filename = std::get<std::string>(extracted[0]);

    // Only support ASCII format.
    params.output_file_format = FileFormat::Ascii;

    return true;
}

/*******************************************************************************
 *	Read the grid separation parameters for z, r, and t.
 ****/
bool ReadGridParams(std::istream& input, RunParams& params)
{
    std::string buf = FindDataLine(input);
    auto extracted = ParseLine(buf, "Error reading dz, dr, dt.", 3);
    if (extracted.size() != 3) {
        return false;
    }

    double dz = std::get<double>(extracted[0]);
    double dr = std::get<double>(extracted[1]);
    double dt = std::get<double>(extracted[2]);

    if (dz <= 0) {
        std::cerr << "Nonpositive dz." << std::endl;
        return false;
    }
    if (dr <= 0) {
        std::cerr << "Nonpositive dr." << std::endl;
        return false;
    }
    if (dt <= 0) {
        std::cerr << "Nonpositve dt. " << std::endl;
        return false;
    }

    params.grid_z = dz;
    params.grid_r = dr;
    params.grid_t = dt;

    return true;
}

/*******************************************************************************
 *	Read the grid array numbers z, r, t, and alpha.
 ****/
bool ReadGridSize(std::istream& input, RunParams& params)
{
    std::string buf = FindDataLine(input);
    auto extracted = ParseLine(buf, "Reading number of dz, dr, dt, da.", 4);
    if (extracted.size() != 4) {
        return false;
    }

    int nz = static_cast<int>(std::get<double>(extracted[0]));
    int nr = static_cast<int>(std::get<double>(extracted[1]));
    int nt = static_cast<int>(std::get<double>(extracted[2]));
    int na = static_cast<int>(std::get<double>(extracted[3]));

    if (nz <= 0) {
        std::cerr << "Nonpositive number of dz." << std::endl;
        return false;
    }
    if (nr <= 0) {
        std::cerr << "Nonpositive number of dr." << std::endl;
        return false;
    }
    if (nt <= 0) {
        std::cerr << "Nonpositive number of dt." << std::endl;
        return false;
    }
    if (na <= 0) {
        std::cerr << "Nonpositive number of da." << std::endl;
        return false;
    }

    params.num_z = nz;
    params.num_r = nr;
    params.num_t = nt;
    params.num_a = na;
    params.grid_a = 0.5 * std::numbers::pi / params.num_a;

    return true;
}

/*******************************************************************************
 *   Initialize Record.
 ****/
void InitRecord(RunParams& params)
{
    params.record.Rd_r = 0;
    params.record.Rd_a = 0;
    params.record.Rd_ra = 0;
    params.record.Rd_t = 0;
    params.record.Rd_rt = 0;
    params.record.Rd_at = 0;
    params.record.Rd_rat = 0;
    params.record.Td_r = 0;
    params.record.Td_a = 0;
    params.record.Td_ra = 0;
    params.record.Td_t = 0;
    params.record.Td_rt = 0;
    params.record.Td_at = 0;
    params.record.Td_rat = 0;
    params.record.A_z = 0;
    params.record.A_rz = 0;
    params.record.A_t = 0;
    params.record.A_zt = 0;
    params.record.A_rzt = 0;
}

/*******************************************************************************
 *  Read which quantity is to be scored.
 ****/
bool ReadRecord(std::istream& input, RunParams& params)
{
    std::string buf = FindDataLine(input);
    if (buf.empty()) {
        std::cout << "Read scored quantities." << std::endl;
        return false;
    }

    std::string string;
    std::stringstream iss(buf);

    do {
        iss >> string;

        // Trim and uppercase
        string = std::format("{:}", string);
        string = ToUpper(string);

        if (string == "RD_R"sv) {
            params.record.Rd_r = true;
        }
        else if (string == "RD_A"sv) {
            params.record.Rd_a = true;
        }
        else if (string == "RD_RA"sv) {
            params.record.Rd_ra = true;
        }
        else if (string == "RD_T"sv) {
            params.record.Rd_t = true;
        }
        else if (string == "RD_RT"sv) {
            params.record.Rd_rt = true;
        }
        else if (string == "RD_AT"sv) {
            params.record.Rd_at = true;
        }
        else if (string == "RD_RAT"sv) {
            params.record.Rd_rat = true;
        }
        else if (string == "TD_R"sv) {
            params.record.Td_r = true;
        }
        else if (string == "TD_A"sv) {
            params.record.Td_a = true;
        }
        else if (string == "TD_RA"sv) {
            params.record.Td_ra = true;
        }
        else if (string == "TD_T"sv) {
            params.record.Td_t = true;
        }
        else if (string == "TD_RT"sv) {
            params.record.Td_rt = true;
        }
        else if (string == "TD_AT"sv) {
            params.record.Td_at = true;
        }
        else if (string == "TD_RAT"sv) {
            params.record.Td_rat = true;
        }
        else if (string == "A_Z"sv) {
            params.record.A_z = true;
        }
        else if (string == "A_RZ"sv) {
            params.record.A_rz = true;
        }
        else if (string == "A_T"sv) {
            params.record.A_t = true;
        }
        else if (string == "A_ZT"sv) {
            params.record.A_zt = true;
        }
        else if (string == "A_RZT"sv) {
            params.record.A_rzt = true;
        }
        else {
            std::cout << "Unknown quantity: " << string << std::endl;
            return false;
        }
    } while (!iss.fail() && !string.empty());

    return true;
}

/*******************************************************************************
*   Filter Record.
****/
bool FilterRecord(std::istream& input, RunParams& params)
{
    InitRecord(params);

    if (!ReadRecord(input, params)) {
        return false;
    }

    if (params.record.Rd_rat) {
        params.record.Rd_ra = 0;
        params.record.Rd_rt = 0;
        params.record.Rd_at = 0;
        params.record.Rd_r = 0;
        params.record.Rd_a = 0;
        params.record.Rd_t = 0;
    }
    if (params.record.Rd_ra) {
        params.record.Rd_r = 0;
        params.record.Rd_a = 0;
    }
    if (params.record.Rd_rt) {
        params.record.Rd_r = 0;
        params.record.Rd_t = 0;
    }
    if (params.record.Rd_at) {
        params.record.Rd_a = 0;
        params.record.Rd_t = 0;
    }
    if (params.record.Td_rat) {
        params.record.Td_ra = 0;
        params.record.Td_rt = 0;
        params.record.Td_at = 0;
        params.record.Td_r = 0;
        params.record.Td_a = 0;
        params.record.Td_t = 0;
    }
    if (params.record.Td_ra) {
        params.record.Td_r = 0;
        params.record.Td_a = 0;
    }
    if (params.record.Td_rt) {
        params.record.Td_r = 0;
        params.record.Td_t = 0;
    }
    if (params.record.Td_at) {
        params.record.Td_a = 0;
        params.record.Td_t = 0;
    }
    if (params.record.A_rzt) {
        params.record.A_rz = 0;
        params.record.A_zt = 0;
        params.record.A_z = 0;
        params.record.A_t = 0;
    }
    if (params.record.A_rz) {
        params.record.A_z = 0;
    }
    if (params.record.A_zt) {
        params.record.A_z = 0;
        params.record.A_t = 0;
    }
    if (params.record.A_zt) {
        params.record.A_z = 0;
        params.record.A_t = 0;
    }
    return true;
}

/*******************************************************************************
 *  Read the weight threshold.
 ****/
bool ReadWeightTreshold(std::fstream& file, RunParams& params)
{
    std::string buf = FindDataLine(file);
    auto extracted = ParseLine(buf, "Error reading threshold weight.");
    if (extracted.empty()) {
        return false;
    }

    double wth = std::get<double>(extracted[0]);

    if (wth < 0 || wth >= 1.0) {
        std::cerr << "Threshold weight out of range (0-1)." << std::endl;
        return false;
    }

    params.weight_treshold = wth;
    return true;
}

/*******************************************************************************
 *  Read the seed for random number generator.
 ****/
bool ReadSeed(std::fstream& file, RunParams& params)
{
    std::string buf = FindDataLine(file);
    auto extracted = ParseLine(buf, "Error reading seed value.");
    if (extracted.empty()) {
        return false;
    }

    int seed = static_cast<int>(std::get<double>(extracted[0]));

    if (seed < 0) {
        std::cerr << "Negative seed value." << std::endl;
        return false;
    }

    params.seed = seed;
    return true;
}

/*******************************************************************************
 *  Find number of layers.
 ****/
int ReadNumLayers(std::fstream& file)
{
    std::string buf;
    short num_layers = -1;

    // While "end" has not been found
    do {
        buf = FindDataLine(file);

        if (buf.empty()) {
            std::cout << "Missing end." << std::endl;
            return 0;
        }
        else {
            // Read layer name
            auto extracted = ParseLine(buf, "Error reading layer name.");
            if (!extracted.empty() && std::holds_alternative<std::string>(extracted[0])) {
                num_layers++;
            }
        }
    } while (buf.find("end") == std::string::npos);

    return num_layers - 2;
}

/*******************************************************************************
 *  Check whether the medium name is in the mediums list.
 ****/
bool ValidMediumName(std::string& name, int& index, RunParams& params)
{
    for (short i = 0; i < params.mediums.size(); i++) {
        if (name == params.mediums[i].name) {
            index = i;
            return true;
        }
    }
    return false;
}

/*******************************************************************************
 *	Read the parameters of all layers.
 ****/
bool ReadLayers(std::fstream& file, RunParams& params)
{
    std::string name;
    double thickness = 0.0;

    // Z coordinate of the current layer
    double z = 0.0;

    // Save current output position
    std::streampos file_pos = file.tellg();

    int num_layers = ReadNumLayers(file);
    if (num_layers < 1) {
        file.seekg(file_pos, std::ios::beg);
        return false;
    }

    // Seek to previous output position 
    file.seekg(file_pos, std::ios::beg);

    // First and last layers are for ambient
    params.num_layers = num_layers;
    params.layers.resize(num_layers + 2);

    for (short i = 0; i <= num_layers + 1; i++) {
        std::string buf = FindDataLine(file);

        // Top and bottom layers (get only name)
        if (i == 0 || i == num_layers + 1) {
            // Get name only
            auto extracted = ParseLine(buf, "Error reading layer specs.");
            if (extracted.empty()) {
                return false;
            }

            name = std::get<std::string>(extracted[0]);
        }
        else {
            // Get name and thickness
            auto extracted = ParseLine(buf, "Error reading layer specs.", 2);
            if (extracted.empty()) {
                return false;
            }
            name = std::get<std::string>(extracted[0]);
            thickness = std::get<double>(extracted[1]);

            if (thickness <= 0.0) {
                std::cerr << "Nonpositive layer thickness." << std::endl;
                return false;
            }
        }

        int index;
        if (!ValidMediumName(name, index, params)) {
            std::cerr << "  Invalid medium name. " << std::endl;
            return false;
        }

        params.layers[i].name = params.mediums[index].name;
        params.layers[i].eta = params.mediums[index].eta;
        params.layers[i].mua = params.mediums[index].mua;
        params.layers[i].mus = params.mediums[index].mus;
        params.layers[i].ani = params.mediums[index].ani;

        // Intermediate layers
        if (i != 0 && i != (num_layers + 1)) {
            params.layers[i].top_z = z;
            z += thickness;
            params.layers[i].bot_z = z;
        }
        // Top and bottom layers
        else {
            params.layers[i].top_z = z;
            params.layers[i].bot_z = z;
        }
    }

    return true;
}

/*******************************************************************************
 *  Read the number of photons and computation time limit.
 *  addMode: Add a number of new photons or seconds to the time limit.
 ****/
bool ReadEndCriteria(std::istream& input, RunParams& params, bool addMode = false)
{
    std::string buf = FindDataLine(input);
    auto extracted = ParseLine(buf, "Error reading number of photons or time limit.");
    if (extracted.empty()) {
        return false;
    }

    if (extracted.size() == 1 && IsNumber(extracted[0])) {
        int num_photons = static_cast<int>(std::get<double>(extracted[0]));

        if (num_photons > 0) {
            params.num_photons = num_photons;
            params.time_limit = 0;
            params.control_bit = ControlBit::NumPhotons;
        }
        else {
            std::cerr << "Nonpositive number of photons." << std::endl;
            return false;
        }
    }
    else if (extracted.size() == 1 && std::holds_alternative<std::string>(extracted[0])) {
        std::string time = std::get<std::string>(extracted[0]);

        int hours = 0;
        int minutes = 0;
        char separator;

        std::istringstream iss(time);
        iss >> hours >> separator >> minutes;

        if (!iss || separator != ':') {
            std::cerr << "Invalid time limit format." << std::endl;
            return false;
        }

        if ((hours * 3600 + minutes * 60) > 0) {
            params.num_photons = 0;
            params.time_limit = hours * 3600 + minutes * 60;
            params.control_bit = ControlBit::TimeLimit;
        }
        else {
            std::cerr << "Nonpositive time limit." << std::endl;
            return false;
        }

        params.control_bit = ControlBit::TimeLimit;
    }
    else if (extracted.size() == 2 && 
             std::holds_alternative<double>(extracted[0]) && 
             std::holds_alternative<std::string>(extracted[1])) {
        int num_photons = static_cast<int>(std::get<double>(extracted[0]));
        std::string time = std::get<std::string>(extracted[1]);

        int hours = 0;
        int minutes = 0;
        char separator;

        std::istringstream iss(time);
        iss >> hours >> separator >> minutes;

        if (!iss || separator != ':') {
            std::cerr << "Invalid time limit format." << std::endl;
            return false;
        }

        if (num_photons > 0 && (hours * 3600 + minutes * 60) > 0) {
            if (addMode) {
                params.add_num_photons = num_photons;
                params.add_limit = hours * 3600 + minutes * 60;
            }
            else {
                params.num_photons = num_photons;
                params.time_limit = hours * 3600 + minutes * 60;
            }

            params.time_limit = hours * 3600 + minutes * 60;
        }
        else {
            std::cerr << "Nonpositive number of photons or time limit." << std::endl;
            return false;
        }

        params.control_bit = ControlBit::Both;
    }
    else {
        std::cerr << "Invalid number of photons or time limit." << std::endl;
        return false;
    }

    return true;
}

/*******************************************************************************
 *  Read the beam source type (Pencil or Isotropic).
 ****/
bool ReadSource(std::fstream& file, RunParams& params)
{
    std::string buf = FindDataLine(file);
    auto extracted = ParseLine(buf, "Error reading photon source type.");
    if (extracted.empty()) {
        return false;
    }

    std::string source_type = std::get<std::string>(extracted[0]);
    if (ToUpper(source_type) == "PENCIL"sv) {
        params.source = BeamType::Pencil;
    }
    else if (ToUpper(source_type) == "ISOTROPIC"sv) {
        params.source = BeamType::Isotropic;
    }
    else {
        std::cerr << "Unknow photon source type. " << std::endl;
        return false;
    }

    return true;
}

/*******************************************************************************
 *  Compute the index to layer according to the z coordinate. If the z is on an
 *  interface between layers, the returned index will point to the upper layer.
 *	Index 0 is the top ambient name and index num_layers+1 is the bottom one.
 ****/
bool LayerIndex(double z, short& index, RunParams& params)
{
    // Index to layer
    short i = 0;
    std::size_t num_layers = params.num_layers;

    if (z < 0.0) {
        std::cerr << "Nonpositive z coordinate." << std::endl;
        return false;
    }
    else if (z > params.layers[num_layers - 1].bot_z) {
        std::cerr << "Source is outside of the last layer. " << std::endl;
        return false;
    }
    else {
        while (z > params.layers[i].bot_z) { i++; }
        index = i;
        return true;
    }
}

/*******************************************************************************
 *  Read starting position of photon source.
 ****/
bool ReadPhotonSource(std::istream& input, RunParams& params)
{
    std::string buf = FindDataLine(input);
    auto extracted = ParseLine(buf, "Invalid starting position of photon source.", 1);
    if (extracted.empty()) {
        return false;
    }

    double source_z = 0.0;
    std::string medium_name = "";

    if (extracted.size() == 1) {
        source_z = std::get<double>(extracted[0]);

        if (!LayerIndex(source_z, params.source_layer, params)) {
            return false;
        }
    }
    else if (extracted.size() == 2) {
        std::string medium_name = std::get<std::string>(extracted[1]);

        if (medium_name[0] != '#' && medium_name[0] != '\n') {
            short source_layer;
            if (!LayerIndex(source_z, source_layer, params)) {
                return false;
            }

            if (params.layers[source_layer].name == medium_name) {
                if ((std::abs(source_z - params.layers[source_layer].bot_z) < DBL_EPSILON) && (params.layers[source_layer + 1].name == medium_name)) {
                    source_layer++;
                    if (source_layer > params.num_layers) {
                        puts("Source is outside of the last layer.");
                        return false;
                    }

                }
                else {
                    std::cerr << "Medium name and z coordinate do not match." << std::endl;
                    return false;
                }
            }

        }
    }

    if (params.source == BeamType::Isotropic && source_z == 0.0) {
        std::cerr << "Can not put isotropic source in upper ambient medium." << std::endl;
        return false;
    }

    params.source_z = source_z;
    params.source_medium_name = medium_name;
    return true;
}

/*******************************************************************************
 *	Compute the critical angles for total internal reflection according to the
 *  relative refractive index of the layer. All layers are processed.
 ****/
void CriticalAngle(std::vector<Layer>& layers)
{
    for (short i = 1; i <= layers.size() - 2; i++) {
        double n1 = layers[i].eta;
        double n2 = layers[i - 1].eta;

        layers[i].cos_crit0 = n1 > n2 ? std::sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;

        n2 = layers[i + 1].eta;
        layers[i].cos_crit1 = n1 > n2 ? std::sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;
    }
}

/*******************************************************************************
 *	Read in the input parameters for one run.
 ****/
void ReadRunParams(std::fstream& file, RunParams& params)
{
    if (!ReadFileName(file, params)) {
        std::exit(1);
    }

    // Geometry.
    if (!ReadLayers(file, params)) {
        std::exit(1);
    }

    // Skip the signal "end" of layers.
    FindDataLine(file);

    // Source.
    if (!ReadSource(file, params)) {
        std::exit(1);
    }
    if (!ReadPhotonSource(file, params)) {
        std::exit(1);
    }

    // Grids.
    if (!ReadGridParams(file, params)) {
        std::exit(1);
    }
    if (!ReadGridSize(file, params)) {
        std::exit(1);
    }
    params.max_z = params.grid_z * params.num_z;
    params.max_r = params.grid_r * params.num_r;
    params.max_time = params.grid_t * params.num_t;
    params.max_alpha = params.grid_a * params.num_a;

    // Scored data categories.
    if (!FilterRecord(file, params)) {
        std::exit(1);
    }

    // Simulation control.
    if (!ReadEndCriteria(file, params, 0)) {
        std::exit(1);
    }
    if (!ReadWeightTreshold(file, params)) {
        std::exit(1);
    }
    if (!ReadSeed(file, params)) {
        std::exit(1);
    }

    CriticalAngle(params.layers);
}

/*******************************************************************************
 *  Read the mediums list in interactive mode.
 ****/
void InterReadMediumList(RunParams& params)
{
    OptSpecValue input;

    int num_media = 0;
    std::cout << "Specify medium list. Total number of mediums: ";

    do {
        input = ParseInput("Invalid medium number. Input again: ");
    } while (!input.has_value() && GetIntValue(input, num_media) < 1);


    // Allocate space for the layer parameters.
    params.mediums.resize(num_media);

    for (short i = 0; i < params.mediums.size(); i++) {
        std::string name;
        bool exists = true;
        std::cout << "Specify medium " << i + 1 << ": " << std::endl << "  Medium name : ";

        do {
            input = ParseInput("Invalid medium name or duplicate. Input again: ");
            name = GetStringValue(input, name);
            exists = std::ranges::any_of(params.mediums, [&](const Layer& l) {
                return l.name == name;
            });
        } while (!input.has_value() && exists);

        params.mediums[i].name = name;

        std::cout << "  Refractive index n (>= 1.0): ";
        do {
            input = ParseInput("  Invalid refractive index. Input again (>= 1.0): ");
        } while (!input.has_value() && GetDoubleValue(input, params.mediums[i].eta) < 1.0);

        std::cout << "  Absorption coefficient mua (>= 0.0 /cm): ";
        do {
            input = ParseInput("  Invalid absorption coefficient. Input again (>= 0.0): ");
        } while (!input.has_value() && GetDoubleValue(input, params.mediums[i].mua) < 0.0);

        std::cout << "  Scattering coefficient mus (>= 0.0 /cm): ";
        do {
            input = ParseInput("  Invalid scattering coefficient. Input again (>= 0.0): ");
        } while (!input.has_value() && GetDoubleValue(input, params.mediums[i].mus) < 0.0);

        std::cout << "  Anisotropy factor g (0.0 - 1.0): ";
        do {
            input = ParseInput("  Invalid anisotropy factor. Input again (0.0 - 1.0): ");
            params.mediums[i].ani = GetDoubleValue(input, params.mediums[i].ani);
        } while (!input.has_value() && params.mediums[i].ani < 0.0 && params.mediums[i].ani > 1.0);

        std::cout << std::endl;
    }
}

/*******************************************************************************
 *  Read the input name and the input format interactively.
 ****/
void InterReadFilename(RunParams& params)
{
    std::string fname;
    std::string fmode;

    do {
        std::cout << "Specify output filename with extension .mco: ";
        std::getline(std::cin, fname);
        fmode[0] = 'w';

        // input exists.
        std::ifstream file(fname, std::ios::in);
        if (!file.is_open()) {
            std::cout << "File " << fname << " exists, w=overwrite, n=new filename: ";

            // avoid null line.
            do {
                std::getline(std::cin, fmode);
            } while (fmode.empty());

            file.close();
        }
    } while (fmode[0] != 'w');

    params.output_filename = fname;

    // Only support 'A' format.
    params.output_file_format = FileFormat::Ascii;

    std::cout << std::endl;
}

/*******************************************************************************
 *	Read grid_z, grid_r, grid_t interactively.
 ****/
void InterReadGrid(RunParams& params)
{
    do {
        std::cout << "Specify dz, dr, dt in one line" << std::endl;
        std::cout << "(all > 0.0 cm, e.g., 0.1 0.1 0.1): ";
    } while (!ReadGridParams(std::cin, params));

    std::cout << std::endl;
}

/*******************************************************************************
 *      Read the RunParams members num_z, num_r, num_a interactively.
 ****/
void InterReadGridSize(RunParams& params)
{
    do {
        std::cout << "Specify nz, nr, nt, na in one line" << std::endl;
        std::cout << "(all > 0, e.g., 100 100 100 100): ";
    } while (!ReadGridSize(std::cin, params));

    params.grid_a = 0.5 * std::numbers::pi / params.num_a;
    std::cout << std::endl;
}

/*******************************************************************************
 *	Read and filter the quantities to be scored interactively.
 ****/
void InterFilterRecord(RunParams& params)
{
    do {
        std::cout << "Select scored quantities from the following data categories:" << std::endl;
        std::cout << "\tRd_rat\t\t\tTd_rat\t\t\tA_rzt" << std::endl;
        std::cout << "\tRd_ra\tRd_rt\tRd_at\tTd_ra\tTd_rt\tRd_at\tA_rz\tA_zt" << std::endl;
        std::cout << "\tRd_r\tRd_a\tRd_t\tTd_r\tTd_a\tTd_t\tA_z\tA_t" << std::endl;
    } while (!FilterRecord(std::cin, params));

    std::cout << std::endl;
}

/*******************************************************************************
 *	Read the threshold weight interactively.
 ****/
void InterReadWeightTreshold(RunParams& params)
{
    std::cout << "Input weight threshold (0 <= W < 1.0, 0.0001 recommended): ";

    std::string string;
    std::getline(std::cin, string);
    std::istringstream iss(string);

    while (!(iss >> params.weight_treshold) || params.weight_treshold < 0 || params.weight_treshold >= 1) {
        std::cout << "Invalid weight treshold. Input again (0 <= W < 1.0): ";
        std::getline(std::cin, string);
    }

    std::cout << std::endl;
}

/*******************************************************************************
 ****/
void PrintMediumNames(RunParams& params)
{
    std::cout << "Available medium types:" << std::endl;

    int j = 1;
    for (int i = 0; i < params.mediums.size(); i++) {
        std::cout << std::format("{:<16}", params.mediums[i].name);
        if (j % 4 == 0) {
            std::cout << std::endl;
        }
        j++;
    }

    std::cout << std::endl;
}

/*******************************************************************************
 *	Read layer specifications interactively.
 ****/
void InterReadLayers(RunParams& params)
{
    std::string name;
    int num_layers;

    // Z coordinate of the current layer.
    double z = 0.0;

    int index;

    std::cout << std::endl << "Specify layer list. ";
    PrintMediumNames(params);

    std::cout << std::endl << "Total number of layers: ";

    // TODO: use ParseInput
    std::string string;
    std::cin >> num_layers;

    while (num_layers < 1) {
        std::cout << "Invalid layer number. Input again: ";
        std::cin >> num_layers;
    }

    params.num_layers = num_layers;

    // Layer 0 and layer num_layers + 1 are for ambient.
    for (short i = 0; i <= params.num_layers + 1; i++) {
        bool error = 1;
        while (error) {
            error = 0;
            if (i == 0) {
                std::cout << std::endl << "  Name of upper ambient medium: ";
                std::cin >> name;
            }
            else if (i == params.num_layers + 1) {
                std::cout << std::endl << "  Name of lower ambient medium: ";
                std::cin >> name;
            }
            else {
                std::cout << std::endl << "  Medium name of layer " << i << ": ";
                std::cin >> name;
            }

            if (!ValidMediumName(name, index, params)) {
                std::cout << "  Invalid medium name. Input again.";
                error = 1;
            }
        }

        params.layers.push_back(Layer{
            .name = params.mediums[index].name,
            .eta = params.mediums[index].eta,
            .mua = params.mediums[index].mua,
            .mus = params.mediums[index].mus,
            .ani = params.mediums[index].ani });


        if ((i != 0) && (i != num_layers + 1)) {
            std::cout << "  Input the thickness of layer " << i << " (thickness > 0.0 cm) : ";
            std::getline(std::cin, string);

            double thick = 0.0;
            std::cin >> thick;

            while (thick <= 0) {
                std::cout << "  Invalid thickness. Input again (thickness > 0.0 cm): ";
                std::cin >> thick;
            }

            params.layers[i].top_z = z;
            z = z + thick;
            params.layers[i].bot_z = z;
        }
        else if (i == 0) {
            params.layers[i].top_z = 0.0;
            params.layers[i].bot_z = 0.0;
        }
        else if (i == params.num_layers + 1) {
            params.layers[i].top_z = z;
            params.layers[i].bot_z = z;
        }
    }

    std::cout << std::endl;
}

/*******************************************************************************
 *  Read the number of photons, or computation time interactively.
 ****/
void InterReadNumPhotons(RunParams& params)
{
    std::cout << "Specify number of photons or time in hh:mm format," << std::endl;
    std::cout << "or both in one line (e.g. 10000 5:30): ";

    while (!ReadEndCriteria(std::cin, params, 0)) {
        std::cout << "Input again: ";
    }

    std::cout << std::endl;
}

/*******************************************************************************
 *  Read the beam source type (Pencil/Isotropic).
 ****/
void InterReadSourceType(RunParams& params)
{
    std::cout << "Input source type (P = pencil / I = isotropic): ";

    char c; std::cin.get(c);
    while (std::toupper(c) != 'P' || std::toupper(c) != 'I') {
        std::cout << "Invalid type. Input again (P = pencil / I = isotropic): ";
        std::cin.get(c);
    }

    if (std::toupper(c) == 'P') {
        params.source = BeamType::Pencil;
    }
    else {
        params.source = BeamType::Isotropic;
    }

    std::cout << std::endl;
}

/*******************************************************************************
 *  Read starting position of photon source.
 ****/
void InterReadPhotonSource(RunParams& params)
{
    do {
        std::cout << "Input the z coordinate of source (0.0 - " << params.layers[params.num_layers].bot_z << " cm) and the medium" << std::endl;
        std::cout << "where the source is if the z is on an interface (e.g. 1.0 [air]):";
    } while (!ReadPhotonSource(std::cin, params));

    std::cout << std::endl;
}

/*******************************************************************************
 *  If input is stdin, freeze the screen and print a more message on screen every 20 lines. The line is the line index.
 ****/
void More(std::ostream& output, int& line)
{
    if (&output == &std::cout) {
        if (!(line % 20)) {
            std::cout << "--More-- (Press Return to continue)";
            output.flush();

            char c;
            do {
                // Read one character from standard input
                std::cin.get(c);
            } while (c != '\n');
        }
    }
}

/*******************************************************************************
 * Write medium names to the output.
 ****/
void WriteMediums(std::ostream& output, RunParams& params, int& line)
{
    std::string format;

    output << std::format("# Specify media \n");
    line++;
    output << std::format("#\tname\t\tn\tmua\tmus\tg\n");
    line++;

    for (int i = 0; i < params.mediums.size(); i++) {
        More(output, line);
        Layer s = params.mediums[i];
        if (s.name.size() + 1 > 8) {
            output << std::format("\t{}\t{:G}\t{:G}\t{:G}\t{:G}\n", s.name, s.eta, s.mua, s.mus, s.ani);
        }
        else {
            output << std::format("\t{}\t\t{:G}\t{:G}\t{:G}\t{:G}\n", s.name, s.eta, s.mua, s.mus, s.ani);
        }
        line++;
    }
    output << std::format("end #of media\n");
    line++;
}

void WriteFilename(std::ostream& output, RunParams& params, int& line)
{
    More(output, line);
    output << std::format("{} \t{}\t\t\t# output file name, format.\n", params.output_filename, 'A');
    line++;
}

void WriteGridParams(std::ostream& output, RunParams& params, int& line)
{
    More(output, line);
    output << std::format("{:G}\t{:G}\t{:G}\t\t\t# dz, dr, dt.\n", params.grid_z, params.grid_r, params.grid_t);
    line++;
}

void WriteGridSize(std::ostream& output, RunParams& params, int& line)
{
    More(output, line);
    output << std::format("{:d}\t{:d}\t{:d}\t{:d}\t\t# nz, nr, nt, na.\n", params.num_z, params.num_r, params.num_t, params.num_a);
    line++;
}

void WriteScored(std::ostream& output, RunParams& params, int& line)
{
    More(output, line);
    output << std::format("# This simulation will score the following categories:\n");
    line++;

    More(output, line);
    if (params.record.Rd_r) {
        output << std::format("Rd_r \t");
    }
    if (params.record.Rd_a) {
        output << std::format("Rd_a \t");
    }
    if (params.record.Rd_ra) {
        output << std::format("Rd_ra \t");
    }
    if (params.record.Rd_t) {
        output << std::format("Rd_t \t");
    }
    if (params.record.Rd_rt) {
        output << std::format("Rd_rt \t");
    }
    if (params.record.Rd_at) {
        output << std::format("Rd_at \t");
    }
    if (params.record.Rd_rat) {
        output << std::format("Rd_rat \t");
    }

    if (params.record.Td_r) {
        output << std::format("Td_r \t");
    }
    if (params.record.Td_a) {
        output << std::format("Td_a \t");
    }
    if (params.record.Td_ra) {
        output << std::format("Td_ra \t");
    }
    if (params.record.Td_t) {
        output << std::format("Td_t \t");
    }
    if (params.record.Td_rt) {
        output << std::format("Td_rt \t");
    }
    if (params.record.Td_at) {
        output << std::format("Td_at \t");
    }
    if (params.record.Td_rat) {
        output << std::format("Td_rat \t");
    }

    if (params.record.A_z) {
        output << std::format("A_z \t");
    }
    if (params.record.A_rz) {
        output << std::format("A_rz \t");
    }
    if (params.record.A_t) {
        output << std::format("A_t \t");
    }
    if (params.record.A_zt) {
        output << std::format("A_zt \t");
    }
    if (params.record.A_rzt) {
        output << std::format("A_rzt \t");
    }

    output << std::format("\n");
    line++;
}

void WriteWeightTreshold(std::ostream& output, RunParams& params, int& line)
{
    More(output, line);
    output << std::format("{:G}\t\t\t\t\t# threshold weight.\n", params.weight_treshold);
    line++;
}

void WriteRandomSeed(std::ostream& output, RunParams& params, int& line)
{
    More(output, line);
    output << std::format("{}\t\t\t\t\t# random number seed.\n", params.seed);
    line++;
}

void WriteLayers(std::ostream& output, RunParams& params, int& line)
{
    std::string format;

    More(output, line);
    output << std::format("# \tmedium \t\tthickness\n");
    line++;

    for (int i = 0; i <= params.num_layers + 1; i++) {
        Layer s;
        More(output, line);

        s = params.layers[i];
        if (i != 0 && i != params.num_layers + 1) {
            if (s.name.size() + 1 > 8) {
                output << std::format("\t{} \t{:G}\n", s.name, s.bot_z - s.top_z);
            }
            else {
                output << std::format("\t{} \t\t{:G}\n", s.name, s.bot_z - s.top_z);
            }
        }
        else {
            output << std::format("\t{}\n", s.name);
        }
        line++;
    }

    More(output, line);
    output << std::format("end #of layers\n");
    line++;
}

void WriteEndCriteria(std::ostream& output, RunParams& params, int& line)
{
    More(output, line);

    if (params.control_bit == ControlBit::NumPhotons) {
        output << std::format("{}  \t\t\t\t\t# no. of photons | time\n", params.num_photons);
    }
    else if (params.control_bit == ControlBit::TimeLimit) {
        output << std::format("{}:{}\t\t\t\t\t# no. of photons | time\n", params.time_limit / 3600, params.time_limit % 3600 / 60);
    }
    else {
        output << std::format("{}  \t{}:{}\t\t\t\t# no. of photons | time\n", params.num_photons, params.time_limit / 3600, params.time_limit % 3600 / 60);
    }

    line++;
}

void WriteSourceType(std::ostream& output, RunParams& params, int& line)
{
    More(output, line);

    if (params.source == BeamType::Pencil) {
        output << std::format("pencil \t\t\t\t\t# src type: pencil/isotropic.\n");
    }
    else {
        output << std::format("isotropic \t\t\t\t# src type: pencil/isotropic.\n");
    }

    line++;
}

void WritePhotonSource(std::ostream& output, RunParams& params, int& line)
{
    More(output, line);

    if (params.source_medium_name.empty()) {
        output << std::format("{:G}\t\t\t\t\t# starting position of source.\n", params.source_z);
    }
    else if (params.source_medium_name.size() + 1 > 8) {
        output << std::format("{:G}\t{} \t\t\t# starting position of source.\n", params.source_z, params.source_medium_name);
    }
    else {
        output << std::format("{:G}\t{} \t\t\t\t# starting position of source.\n", params.source_z, params.source_medium_name);
    }

    line++;
}

/*******************************************************************************
 *  Write input parameters to the output file.
 ****/
void WriteInputParams(std::ostream& output, RunParams& params)
{
    // line index.
    output << std::format("mcmli2.0 \t\t\t# file version \n\n");
    int line = 2;
    WriteMediums(output, params, line);

    More(output, line);
    output << std::format("\n# Specify data for run 1\n");
    line += 2;

    WriteFilename(output, params, line);

    // geometry.
    More(output, line);
    output << std::format("\n");
    line++;
    WriteLayers(output, params, line);

    // source.
    More(output, line);
    output << std::format("\n");
    line++;
    WriteSourceType(output, params, line);
    WritePhotonSource(output, params, line);

    // grids.
    More(output, line);
    output << std::format("\n");
    line++;
    WriteGridParams(output, params, line);
    WriteGridSize(output, params, line);

    // scored data categories.
    More(output, line);
    output << std::format("\n");
    line++;
    WriteScored(output, params, line);

    // simulation control.
    More(output, line);
    output << std::format("\n");
    line++;

    WriteEndCriteria(output, params, line);
    WriteWeightTreshold(output, params, line);
    WriteRandomSeed(output, params, line);

    More(output, line);
    output << std::format("end #of runs\n\n");
}

/*******************************************************************************
 *  Read in the input parameters for one run in interactive mode.
 ****/
void InterReadParams(RunParams& params)
{
    InterReadMediumList(params);
    InterReadFilename(params);
    InterReadLayers(params);
    InterReadSourceType(params);
    InterReadPhotonSource(params);
    InterReadGrid(params);
    InterReadGridSize(params);
    InterFilterRecord(params);
    InterReadNumPhotons(params);
    InterReadWeightTreshold(params);

    std::cout << "Do you want to save the input to a file? (Y/N)";
    char c; std::cin.get(c);

    if (std::toupper(c) == 'Y') {
        std::cout << "Give the file name to save input: ( .mci): ";

        std::string string;
        std::getline(std::cin, string);

        std::fstream file(string, std::ios::out);
        if (!file.is_open()) {
            std::cout << "Can not open the file to write." << std::endl;
        }
        else {
            WriteInputParams(file, params);
        }
    }
}

/*******************************************************************************
 *  Check consistancy of input parameters for one run. Such as the consistancy
 *  of the medium names, layer list, source starting position and source type.
 ****/
bool CheckInputParams(RunParams& params)
{
    for (int i = 0; i <= params.num_layers + 1; i++) {
        int index;
        if (!ValidMediumName(params.layers[i].name, index, params)) {
            std::cout << "Invalid medium name of layer " << i << "." << std::endl;
            return 0;
        }
        else {
            params.layers[i].eta = params.mediums[index].eta;
            params.layers[i].mua = params.mediums[index].mua;
            params.layers[i].mus = params.mediums[index].mus;
            params.layers[i].ani = params.mediums[index].ani;
        }
    }

    if ((params.source == BeamType::Isotropic) && (params.source_z == 0.0)) {
        std::cout << "Can not put isotropic source in upper ambient medium." << std::endl;
        return 0;
    }
    if (!LayerIndex(params.source_z, params.source_layer, params)) {
        return 0;
    }

    if (params.source_medium_name[0] != '\0') {
        if (params.layers[params.source_layer].name == params.source_medium_name) {
            if ((std::abs(params.source_z - params.layers[params.source_layer].bot_z) < DBL_EPSILON) && (params.layers[params.source_layer + 1].name == params.source_medium_name)) {
                params.source_layer++;
            }
            else {
                std::cerr << "Medium name and z coordinate do not match." << std::endl;
                return false;
            }
        }
    }
    return 1;
}

/*******************************************************************************
 *  Menu for changing input parameters.
 *****/
void ShowChangeMenu()
{
    std::cout << "  o = Print the input on screen." << std::endl;
    std::cout << "  m = Change media list." << std::endl;
    std::cout << "  f = Change output file name and format." << std::endl;
    std::cout << "  d = Change dz, dr, dt." << std::endl;
    std::cout << "  n = Change nz, nr, nt, na." << std::endl;
    std::cout << "  c = Change scored data categories." << std::endl;
    std::cout << "  w = Change threshold weight." << std::endl;
    std::cout << "  r = Change random number seed." << std::endl;
    std::cout << "  l = Change layer specifications." << std::endl;
    std::cout << "  p = Change photon number and computation time limit." << std::endl;
    std::cout << "  s = Change source type." << std::endl;
    std::cout << "  z = Change source starting position." << std::endl;
    std::cout << "  q = Quit from change menu and start simulation." << std::endl;
    std::cout << "  x = Exit to the main menu." << std::endl;
    std::cout << "  * Commands here are not case-sensitive" << std::endl;
}

/*******************************************************************************
 *  Continue to change input parameters or quit.
 *****/
char QuitOrContinue()
{
    std::string string;

    do {
        std::cout << "Do you want to change them? (Y/N): ";
        do {
            std::getline(std::cin, string);
        } while (!string.empty());
    } while (std::toupper(string[0]) != 'Y' && std::toupper(string[0]) != 'N');

    return std::toupper(string[0]);
}

void EditMediums(RunParams& params)
{
    int line = 1;
    std::cout << "Current medium list: " << std::endl;
    WriteMediums(std::cout, params, line);
    std::cout << std::endl;

    if (QuitOrContinue() == 'Y') {
        params.mediums.clear();
        InterReadMediumList(params);
    }
}

void EditFileName(RunParams& params)
{
    int line = 1;
    std::cout << "Current output file name and format: " << std::endl;
    WriteFilename(std::cout, params, line);
    std::cout << std::endl;
    InterReadFilename(params);
}

void EditGridParams(RunParams& params)
{
    int line = 1;
    std::cout << "Current dz, dr, dt: " << std::endl;
    WriteGridParams(std::cout, params, line);
    std::cout << std::endl;
    InterReadGrid(params);
}

void EditGridSize(RunParams& params)
{
    int line = 1;
    std::cout << "Current nz, nr, nt, na: " << std::endl;
    WriteGridSize(std::cout, params, line);
    std::cout << std::endl;
    InterReadGridSize(params);
}

void EditRecord(RunParams& params)
{
    int line = 1;
    WriteScored(std::cout, params, line);
    std::cout << std::endl;

    if (QuitOrContinue() == 'Y') {
        InterFilterRecord(params);
    }
}

void EditWeightTreshold(RunParams& params)
{
    int line = 1;
    std::cout << "Current threshold weight: " << std::endl;
    WriteWeightTreshold(std::cout, params, line);
    std::cout << std::endl;
    InterReadWeightTreshold(params);
}

void EditLayers(RunParams& params)
{
    int line = 1;
    std::cout << "Current layer sepcification: " << std::endl;
    WriteLayers(std::cout, params, line);
    std::cout << std::endl;

    if (QuitOrContinue() == 'Y') {
        InterReadLayers(params);
    }
}

void EditEndCriteria(RunParams& params)
{
    int line = 1;
    std::cout << "Current value: " << std::endl;
    WriteEndCriteria(std::cout, params, line);
    std::cout << std::endl;
    InterReadNumPhotons(params);
}

void EditSourceType(RunParams& params)
{
    int line = 1;
    std::cout << "Current source type: " << std::endl;
    WriteSourceType(std::cout, params, line);
    std::cout << std::endl;
    InterReadSourceType(params);
}

void EditPhotonSource(RunParams& params)
{
    int line = 1;
    std::cout << "Layer Specification: " << std::endl;
    WriteLayers(std::cout, params, line);
    std::cout << "\nCurrent starting position: " << std::endl;
    WritePhotonSource(std::cout, params, line);
    std::cout << std::endl;
    InterReadPhotonSource(params);
}

/*******************************************************************************
 *  If 'Q': Quit change menu; return 1.
 *  If 'X': Quit to main menu; return 2.
 *  Otherwise; return 0.
 ****/
int BranchChangeMenu(std::string& string, RunParams& params)
{
    switch (std::toupper(string[0])) {
        case 'M':
        {
            EditMediums(params);
            break;
        }

        case 'F':
        {
            EditFileName(params);
            break;
        }

        case 'D':
        {
            EditGridParams(params);
            break;
        }

        case 'N':
        {
            EditGridSize(params);
            break;
        }

        case 'C':
        {
            EditRecord(params);
            break;
        }

        case 'W':
        {
            EditWeightTreshold(params);
            break;
        }

        case 'L':
        {
            EditLayers(params);
            break;
        }

        case 'P':
        {
            EditEndCriteria(params);
            break;
        }

        case 'S':
        {
            EditSourceType(params);
            break;
        }

        case 'Z':
        {
            EditPhotonSource(params);
            break;
        }

        case 'O':
        {
            WriteInputParams(std::cout, params);
            break;
        }

        case 'H':
        {
            ShowChangeMenu();
            break;
        }

        case 'Q':
        {
            return 1;
        }

        case 'X':
        {
            return 2;
        }

        default:
        {
            std::cout << "...Unknown command" << std::endl;
        }
    }

    return 0;
}

/********************************************************************************
 *   Return 1 if quit change and start simulation.
 *   Return 0 if exit to main menu.
 ****/
bool RunNewInput(RunParams& params)
{
    std::string string;
    int branch;

    std::cout << "Any changes to the input parameters? (Y/N)";
    do {
        std::getline(std::cin, string);
    } while (string.empty());

    while (std::toupper(string[0]) == 'Y') {
        do {
            do {
                std::cout << "\n> Change menu (h for help) => ";
                std::getline(std::cin, string);
            } while (string.empty());

            // 'X' or 'Q'
            if (branch = BranchChangeMenu(string, params)) {
                break;
            }
        } while (1);

        std::cout << "Do you want to save the input to a file? (Y/N)";
        std::getline(std::cin, string);
        if (std::toupper(string[0]) == 'Y') {
            std::cout << "Give the file name to save input: ( .mci): ";
            std::getline(std::cin, string);

            std::ofstream file(string, std::ios::out);
            if (!file.is_open()) {
                std::cout << "Can not open the file to write." << std::endl;
            }
            else {
                WriteInputParams(file, params);
            }
        }

        // Quit change menu and start simulation
        if (branch == 1) {
            if (!CheckInputParams(params)) {
                do {
                    std::cout << "Change input or exit to main menu (c/x): ";
                    std::getline(std::cin, string);
                } while (!string.empty() || std::toupper(string[0]) != 'X' && std::toupper(string[0]) != 'C');

                if (std::toupper(string[0]) == 'X') {
                    params.mediums.clear();
                    params.layers.clear();
                    return false;
                }
                else {
                    // Continue to change parameters
                    string[0] = 'Y';
                }
            }
            else {
                return true;
            }

        }
        // Exit to menu
        else {
            params.mediums.clear();
            params.layers.clear();
            return false;
        }
    }

    return true;
}

/*******************************************************************************
 *  Check the whether the flag end is found.
 ****/
bool CheckEndOfRuns(std::istream& input)
{
    // Found end of runs
    bool end_found = false;

    // Record input position
    std::streampos file_pos = input.tellg();

    std::string buf = FindDataLine(input);
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
}

/*******************************************************************************
 *  Check the input parameters for all runs.
 *  This function will count number of runs and assign it to params.num_runs.
 ****/
void CheckParamFromFile(std::fstream& input, RunParams& params)
{
    if (!ReadMediumListQ(input, params)) {
        std::exit(1);
    }

    // Save current position in input file.
    std::streampos file_pos = input.tellg();

    short run_index = 0;
    do {
        std::cout << "Checking input data for run " << ++run_index << std::endl;
        ReadRunParams(input, params);

        // Attempt insert and detect duplicates
        if (!params.unique_outputs.insert(params.output_filename).second) {
            std::cout << "File name " << params.output_filename << " duplicated." << std::endl;
            std::exit(1);
        }
    } while (!CheckEndOfRuns(input));

    params.num_runs = run_index;

    // Restore file position
    input.seekg(file_pos, std::ios::beg);
}

/*******************************************************************************
 *	Allocate the Tracer arrays for one run, and array elements are initialized 
 *  to zeros.
 ****/
void InitTracer(RunParams& params, Tracer& tracer)
{
    std::size_t nz = params.num_z;
    std::size_t nr = params.num_r;
    std::size_t na = params.num_a;
    std::size_t nt = params.num_t;

    tracer.R.sp = 0.0;
    tracer.R.br = 0.0;
    tracer.R.dr = 0.0;
    tracer.T.dr = 0.0;
    tracer.T.br = 0.0;
    tracer.A.ab = 0.0;

    tracer.R.be = 0.0;
    tracer.R.de = 0.0;
    tracer.T.de = 0.0;
    tracer.T.be = 0.0;
    tracer.A.ae = 0.0;

    auto alloc3 = [](std::size_t x, std::size_t y, std::size_t z) {
        return vec3<double>(x, vec2<double>(y, vec1<double>(z, 0.0)));
    };

    auto alloc2 = [](std::size_t x, std::size_t y) {
        return vec2<double>(x, vec1<double>(y, 0.0));
    };

    auto alloc1 = [](std::size_t x) {
        return vec1<double>(x, 0.0);
    };

    // Reflectance
    if (params.record.Rd_rat) { tracer.R.rat = alloc3(nr, na, nt); }
    if (params.record.Rd_ra) { tracer.R.ra = alloc2(nr, na); }
    if (params.record.Rd_rt) { tracer.R.rt = alloc2(nr, nt); }
    if (params.record.Rd_at) { tracer.R.at = alloc2(na, nt); }
    if (params.record.Rd_r) { tracer.R.r = alloc1(nr); }
    if (params.record.Rd_a) { tracer.R.a = alloc1(na); }
    if (params.record.Rd_t) { tracer.R.t = alloc1(nt); }

    // Transmittance
    if (params.record.Td_rat) { tracer.T.rat = alloc3(nr, na, nt); }
    if (params.record.Td_ra) { tracer.T.ra = alloc2(nr, na); }
    if (params.record.Td_rt) { tracer.T.rt = alloc2(nr, nt); }
    if (params.record.Td_at) { tracer.T.at = alloc2(na, nt); }
    if (params.record.Td_r) { tracer.T.r = alloc1(nr); }
    if (params.record.Td_a) { tracer.T.a = alloc1(na); }
    if (params.record.Td_t) { tracer.T.t = alloc1(nt); }

    // Absorption
    if (params.record.A_rzt) {
        tracer.A.rzt = alloc3(nr, nz, nt);
        tracer.A.bzt = alloc2(nz, nt);
    }
    if (params.record.A_rz) {
        tracer.A.rz = alloc2(nr, nz);
        tracer.A.bz = alloc1(nz);
    }
    if (params.record.A_zt) { tracer.A.zt = alloc2(nz, nt); }
    if (params.record.A_z) { tracer.A.z = alloc1(nz); }
    if (params.record.A_t) { tracer.A.t = alloc1(nt); }
}

void ClearRun(RunParams& params, Tracer& tracer)
{
    if (params.record.Rd_rat) { tracer.R.rat.clear(); }
    if (params.record.Rd_ra) { tracer.R.ra.clear(); }
    if (params.record.Rd_rt) { tracer.R.rt.clear(); }
    if (params.record.Rd_at) { tracer.R.at.clear(); }
    if (params.record.Rd_r) { tracer.R.r.clear(); }
    if (params.record.Rd_a) { tracer.R.a.clear(); }
    if (params.record.Rd_t) { tracer.R.t.clear(); }

    if (params.record.Td_rat) { tracer.T.rat.clear(); }
    if (params.record.Td_ra) { tracer.T.ra.clear(); }
    if (params.record.Td_rt) { tracer.T.rt.clear(); }
    if (params.record.Td_at) { tracer.T.at.clear(); }
    if (params.record.Td_r) { tracer.T.r.clear(); }
    if (params.record.Td_a) { tracer.T.a.clear(); }
    if (params.record.Td_t) { tracer.T.t.clear(); }

    if (params.record.A_rzt) {
        tracer.A.rzt.clear();
        tracer.A.bzt.clear();
    }
    if (params.record.A_rz) {
        tracer.A.rz.clear();
        tracer.A.bz.clear();
    }
    if (params.record.A_zt) { tracer.A.zt.clear(); }
    if (params.record.A_z) { tracer.A.z.clear(); }
    if (params.record.A_t) { tracer.A.t.clear(); }

    params.layers.clear();
}

/*******************************************************************************
 *	Scale reflectance and transmittance.
 *	"a" stands for angle alpha.
 ****
 *	Scale Rd(r,a) and Td(r,a) by: 
 *      area perpendicular to photon direction * solid angle * num_photons
 *	or
 *		[2*PI*r*grid_r*cos(a)] * [2*PI*sin(a)*grid_a] * [num_photons]
 *	or
 *		[2*PI*PI*grid_r*grid_a*r*sin(2a)] * [num_photons]
 ****
 *	Scale Rd(r) and Td(r) by
 *		area on the surface * num_photons.
 ****
 *	Scale Rd(a) and Td(a) by
 *		solid angle * num_photons.
 ****/
void ScaleReflectanceTransmittance(RunParams& params, Tracer& tracer, ScaleMode mode)
{
    std::size_t nr = params.num_r;
    std::size_t na = params.num_a;
    std::size_t nt = params.num_t;

    double dr = params.grid_r;
    double da = params.grid_a;
    double dt = params.grid_t;

    double scale1 = (double)params.num_photons;
    if (mode == ScaleMode::Scale) {
        tracer.R.de = 1 / scale1 * sqrt(tracer.R.de - tracer.R.dr * tracer.R.dr / scale1);
        tracer.T.de = 1 / scale1 * sqrt(tracer.T.de - tracer.T.dr * tracer.T.dr / scale1);
        tracer.R.be = 1 / scale1 * sqrt(tracer.R.be - tracer.R.br * tracer.R.br / scale1);
        tracer.T.be = 1 / scale1 * sqrt(tracer.T.be - tracer.T.br * tracer.T.br / scale1);

        tracer.R.dr /= scale1;
        tracer.T.dr /= scale1;
        tracer.R.br = tracer.R.br / scale1 + tracer.R.sp;
        tracer.T.br /= scale1;
    }
    else {
        tracer.R.dr *= scale1;
        tracer.T.dr *= scale1;
        tracer.R.br = (tracer.R.br - tracer.R.sp) * scale1;
        tracer.T.br *= scale1;

        tracer.R.de = (scale1 * tracer.R.de) * (scale1 * tracer.R.de) + 1 / scale1 * tracer.R.dr * tracer.R.dr;
        tracer.T.de = (scale1 * tracer.T.de) * (scale1 * tracer.T.de) + 1 / scale1 * tracer.T.dr * tracer.T.dr;
        tracer.R.be = (scale1 * tracer.R.be) * (scale1 * tracer.R.be) + 1 / scale1 * tracer.R.br * tracer.R.br;
        tracer.T.be = (scale1 * tracer.T.be) * (scale1 * tracer.T.be) + 1 / scale1 * tracer.T.br * tracer.T.br;
    }

    scale1 = dt * params.num_photons;
    if (params.record.Rd_t) {
        for (std::size_t it = 0; it < nt; it++) {
            // scale Rd_t
            if (mode == ScaleMode::Scale) {
                tracer.R.t[it] /= scale1;
            }
            // unscale Rd_t
            else {
                tracer.R.t[it] *= scale1;
            }
        }
    }

    if (params.record.Td_t) {
        for (std::size_t it = 0; it < nt; it++) {
            // scale Td_t
            if (mode == ScaleMode::Scale) {
                tracer.T.t[it] /= scale1;
            }
            // unscale Rd_t
            else {
                tracer.T.t[it] *= scale1;
            }
        }
    }

    scale1 = 2.0 * std::numbers::pi * dr * dr * params.num_photons;
    // area is 2*PI*[(ir+0.5)*grid_r]*grid_r.  ir + 0.5 to be added.

    if (params.record.Rd_r) {
        for (std::size_t ir = 0; ir < nr; ir++) {
            double scale2 = 1.0 / ((ir + 0.5) * scale1);
            // scale Rd_r
            if (mode == ScaleMode::Scale) {
                tracer.R.r[ir] *= scale2;
            }
            // unscale Rd_r
            else {
                tracer.R.r[ir] /= scale2;
            }
        }
    }

    if (params.record.Td_r) {
        for (std::size_t ir = 0; ir < nr; ir++) {
            double scale2 = 1.0 / ((ir + 0.5) * scale1);
            // scale Td_r
            if (mode == ScaleMode::Scale) {
                tracer.T.r[ir] *= scale2;
            }
            // unscale Td_r
            else {
                tracer.T.r[ir] /= scale2;
            }
        }
    }

    scale1 *= dt;
    if (params.record.Rd_rt) {
        for (std::size_t ir = 0; ir < nr; ir++) {
            for (std::size_t it = 0; it < nt; it++) {
                double scale2 = 1.0 / ((ir + 0.5) * scale1);
                // scale Rd_rt
                if (mode == ScaleMode::Scale) {
                    tracer.R.rt[ir][it] *= scale2;
                }
                // unscale Rd_rt
                else {
                    tracer.R.rt[ir][it] *= scale2;
                }
            }
        }
    }

    if (params.record.Td_rt) {
        for (std::size_t ir = 0; ir < nr; ir++) {
            for (std::size_t it = 0; it < nt; it++) {
                double scale2 = 1.0 / ((ir + 0.5) * scale1);
                // scale Td_rt
                if (mode == ScaleMode::Scale) {
                    tracer.T.rt[ir][it] *= scale2;
                }
                // unscale Td_rt
                else {
                    tracer.T.rt[ir][it] /= scale2;
                }
            }
        }
    }

    scale1 = std::numbers::pi * da * params.num_photons;

    // Solid angle times cos(a) is PI*sin(2a)*grid_a. sin(2a) to be added.
    if (params.record.Rd_a) {
        for (std::size_t ia = 0; ia < na; ia++) {
            double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
            // scale Rd_a.
            if (mode == ScaleMode::Scale) {
                tracer.R.a[ia] *= scale2;
            }
            // unscale Rd_a.
            else {
                tracer.R.a[ia] /= scale2;
            }
        }
    }

    if (params.record.Td_a) {
        for (std::size_t ia = 0; ia < na; ia++) {
            double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
            // scale Td_a
            if (mode == ScaleMode::Scale) {
                tracer.T.a[ia] *= scale2;
            }
            // unscale Td_a
            else {
                tracer.T.a[ia] /= scale2;
            }
        }
    }

    scale1 *= dt;
    if (params.record.Rd_at) {
        for (std::size_t ia = 0; ia < na; ia++) {
            for (std::size_t it = 0; it < nt; it++) {
                double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
                // scale Rd_at
                if (mode == ScaleMode::Scale) {
                    tracer.R.at[ia][it] *= scale2;
                }
                // unscale Rd_at
                else {
                    tracer.R.at[ia][it] /= scale2;
                }
            }
        }
    }

    if (params.record.Td_at) {
        for (std::size_t ia = 0; ia < na; ia++) {
            for (std::size_t it = 0; it < nt; it++) {
                double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
                // scale Td_at
                if (mode == ScaleMode::Scale) {
                    tracer.T.at[ia][it] *= scale2;
                }
                // unscale Td_at
                else {
                    tracer.T.at[ia][it] /= scale2;
                }
            }
        }
    }

    scale1 = 2.0 * std::numbers::pi * dr * dr * std::numbers::pi * da * params.num_photons;
    if (params.record.Rd_ra) {
        for (std::size_t ir = 0; ir < nr; ir++) {
            for (std::size_t ia = 0; ia < na; ia++) {
                double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                // scale Rd_ra
                if (mode == ScaleMode::Scale) {
                    tracer.R.ra[ir][ia] *= scale2;
                }
                // unscale Rd_ra
                else {
                    tracer.R.ra[ir][ia] /= scale2;
                }
            }
        }
    }

    if (params.record.Td_ra) {
        for (std::size_t ir = 0; ir < nr; ir++) {
            for (std::size_t ia = 0; ia < na; ia++) {
                double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                // scale Td_ra
                if (mode == ScaleMode::Scale) {
                    tracer.T.ra[ir][ia] *= scale2;
                }
                // unscale Td_ra
                else {
                    tracer.T.ra[ir][ia] /= scale2;
                }
            }
        }
    }

    scale1 *= dt;
    if (params.record.Rd_rat) {
        for (std::size_t ir = 0; ir < nr; ir++) {
            for (std::size_t ia = 0; ia < na; ia++) {
                for (std::size_t it = 0; it < nt; it++) {
                    double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                    // scale Rd_rat
                    if (mode == ScaleMode::Scale) {
                        tracer.R.rat[ir][ia][it] *= scale2;
                    }
                    // unscale Rd_rat
                    else {
                        tracer.R.rat[ir][ia][it] /= scale2;
                    }
                }
            }
        }
    }

    if (params.record.Td_rat) {
        for (std::size_t ir = 0; ir < nr; ir++) {
            for (std::size_t ia = 0; ia < na; ia++) {
                for (std::size_t it = 0; it < nt; it++) {
                    double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                    // scale Td_rat
                    if (mode == ScaleMode::Scale) {
                        tracer.T.rat[ir][ia][it] *= scale2;
                    }
                    // unscale Td_rat
                    else {
                        tracer.T.rat[ir][ia][it] /= scale2;
                    }
                }
            }
        }
    }
}

/*******************************************************************************
 *	Scale absorption quantities.
 ****/
void ScaleAbsorption(RunParams& params, Tracer& tracer, ScaleMode mode)
{
    std::size_t nz = params.num_z;
    std::size_t nr = params.num_r;
    std::size_t nt = params.num_t;

    double dz = params.grid_z;
    double dr = params.grid_r;
    double dt = params.grid_t;

    double scale1 = (double)params.num_photons;

    // scale A
    if (mode == ScaleMode::Scale) {
        tracer.A.ae = 1 / scale1 * sqrt(tracer.A.ae - tracer.A.ab * tracer.A.ab / scale1);
        tracer.A.ab /= scale1;
    }
    // unscale A
    else {
        tracer.A.ab *= scale1;
        tracer.A.ae = (scale1 * tracer.A.ae) * (scale1 * tracer.A.ae) + 1 / scale1 * tracer.A.ab * tracer.A.ab;
    }

    double scale2 = scale1 * dt;
    if (params.record.A_t) {
        for (short it = 0; it < nt; it++) {
            // scale A_t
            if (mode == ScaleMode::Scale) {
                tracer.A.t[it] /= scale2;
            }
            // unscale A_t
            else {
                tracer.A.t[it] *= scale2;
            }
        }
    }

    scale1 *= dz;
    if (params.record.A_z) {
        for (short iz = 0; iz < nz; iz++) {
            // scale A_z
            if (mode == ScaleMode::Scale) {
                tracer.A.z[iz] /= scale1;
            }
            // unscale A_z
            else {
                tracer.A.z[iz] *= scale1;
            }
        }
    }

    scale2 = scale1 * dt;
    if (params.record.A_zt) {
        for (short iz = 0; iz < nz; iz++) {
            for (short it = 0; it < nt; it++) {
                // scale A_zt
                if (mode == ScaleMode::Scale) {
                    tracer.A.zt[iz][it] /= scale2;
                }
                // unscale A_zt
                else {
                    tracer.A.zt[iz][it] *= scale2;
                }
            }
        }
    }

    if (params.record.A_rz) {
        for (short iz = 0; iz < nz; iz++) {
            // scale Ab_z
            if (mode == ScaleMode::Scale) {
                tracer.A.bz[iz] /= scale1;
            }
            // unscale Ab_z
            else {
                tracer.A.bz[iz] *= scale1;
            }
        }
    }

    if (params.record.A_rzt) {
        for (short iz = 0; iz < nz; iz++) {
            for (short it = 0; it < nt; it++) {
                // scale Ab_zt
                if (mode == ScaleMode::Scale) {
                    tracer.A.bzt[iz][it] /= scale2;
                }
                // unscale Ab_zt
                else {
                    tracer.A.bzt[iz][it] *= scale2;
                }
            }
        }
    }

    scale1 = 2.0 * std::numbers::pi * dr * dr * dz * params.num_photons;
    if (params.record.A_rz) {
        for (short ir = 0; ir < nr; ir++) {
            for (short iz = 0; iz < nz; iz++) {
                // scale A_rz
                if (mode == ScaleMode::Scale) {
                    tracer.A.rz[ir][iz] /= (ir + 0.5) * scale1;
                }
                // unscale A_rz
                else {
                    tracer.A.rz[ir][iz] *= (ir + 0.5) * scale1;
                }
            }
        }
    }

    scale2 = scale1 * dt;
    if (params.record.A_rzt) {
        for (short ir = 0; ir < nr; ir++) {
            for (short iz = 0; iz < nz; iz++) {
                for (short it = 0; it < nt; it++) {
                    // scale A_rzt
                    if (mode == ScaleMode::Scale) {
                        tracer.A.rzt[ir][iz][it] /= (ir + 0.5) * scale2;
                    }
                    // unscale A_rzt
                    else {
                        tracer.A.rzt[ir][iz][it] *= (ir + 0.5) * scale2;
                    }
                }
            }
        }
    }
}

/*******************************************************************************
 *	Scale results of current run.
 *  Mode 0: scale result.
 *  Mode 1: unscale result.
 ****/
void ScaleResult(RunParams& params, Tracer& tracer, ScaleMode mode)
{
    ScaleReflectanceTransmittance(params, tracer, mode);
    ScaleAbsorption(params, tracer, mode);
}

/*******************************************************************************
 *	Write the version number as the first string in the input.
 *	Use chars only so that they can be read as either ASCII or binary.
 ****/
void WriteVersion(std::fstream& file, const std::string& version)
{
    file << version << " \t# Version number of the file format.\n" << std::endl;
    file << "####\n# Data categories include: " << std::endl;
    file << "# Rd_r\tRd_a\tRd_ra\tRd_t\tRd_rt\tRd_at\tRd_rat" << std::endl;
    file << "# Td_r\tTd_a\tTd_ra\tTd_t\tTd_rt\tTd_at\tTd_rat" << std::endl;
    file << "# A_z\tA_rz\tA_t\t\tA_zt\tA_rzt" << std::endl;
    file << "####\n" << std::endl;
}

/********************************************************************************
 * Save the status of the random number generator to output file.
 ****/
void WriteRandomizer(std::fstream& file)
{
    auto status = Rand.getState();

    file << std::format("# status of the random number generator:") << std::endl;

    for (int i = 0; i < status.size(); i++) {
        if (i % 5) {
            file << std::format("{:14d}", status[i]);
        }
        else {
            file << std::endl << std::format("{:14d}", status[i]);
        }
    }

    file << std::endl << std::endl;
}

/********************************************************************************
 * Read and restore the status of random number generater from previous
 * output input.
 ****/
void ReadRandomizer(std::fstream& file)
{
    std::string buf;
    std::vector<std::mt19937::result_type> status(624);

    do {
        std::getline(file, buf);
    } while (buf[0] != '#');

    for (int i = 0; i < status.size(); i++) {
        file >> status[i];
    }

    // Restore the status
    Rand.setState(status);
}

/*******************************************************************************
 *	Write reflectance, absorption, transmission.
 ****/
void WriteRAT(std::fstream& file, Tracer& tracer)
{
    file << std::format("RAT #Reflectance, Absorption, Transmittance.\n");
    file << std::format("# Average \tStandard Err \tRel Err\n");
    file << std::format("{:<14.6G} \t\t\t\t#Rsp: Specular reflectance.\n", tracer.R.sp);
    file << std::format("{:<14.6G} \t{:<14.6G} {:6.2f}%\t#Rb: Ballistic reflectance.\n", tracer.R.br, tracer.R.be, (tracer.R.br) ? tracer.R.be / tracer.R.br * 100 : 0);
    file << std::format("{:<14.6G} \t{:<14.6G} {:6.2f}%\t#Rd: Diffuse reflectance.\n", tracer.R.dr, tracer.R.de, (tracer.R.dr) ? tracer.R.de / tracer.R.dr * 100 : 0);
    file << std::format("{:<14.6G} \t{:<14.6G} {:6.2f}%\t#A:  Absorbed fraction.\n", tracer.A.ab, tracer.A.ae, (tracer.A.ab) ? tracer.A.ae / tracer.A.ab * 100 : 0);
    file << std::format("{:<14.6G} \t{:<14.6G} {:6.2f}%\t#Tb: Ballistic transmittance.\n", tracer.T.br, tracer.T.be, (tracer.T.br) ? tracer.T.be / tracer.T.br * 100 : 0);
    file << std::format("{:<14.6G} \t{:<14.6G} {:6.2f}%\t#Td: Diffuse transmittance.\n", tracer.T.dr, tracer.T.de, (tracer.T.dr) ? tracer.T.de / tracer.T.dr * 100 : 0);
    file << std::format("\n");
}

/*******************************************************************************
 *	Read reflectance, absorption, transmission.
 ****/
void ReadRAT(std::fstream& file, Tracer& tracer)
{
    // skip RAT line.
    std::string buf = FindDataLine(file);

    buf = FindDataLine(file);
    std::istringstream iss(buf);
    iss >> tracer.R.sp;

    buf = FindDataLine(file);
    iss = std::istringstream(buf);
    iss >> tracer.R.br >> tracer.R.be;

    buf = FindDataLine(file);
    iss = std::istringstream(buf);
    iss >> tracer.R.dr >> tracer.R.de;

    buf = FindDataLine(file);
    iss = std::istringstream(buf);
    iss >> tracer.A.ab >> tracer.A.ae;

    buf = FindDataLine(file);
    iss = std::istringstream(buf);
    iss >> tracer.T.br >> tracer.T.be;

    buf = FindDataLine(file);
    iss = std::istringstream(buf);
    iss >> tracer.T.dr >> tracer.T.de;
}

/*******************************************************************************
 *  Ballistic absorption per unit depth, per unit time [1/(cm ps)]
 ****/
void IOAb_zt(std::fstream& file, std::size_t Nz, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Ab[z][t]. [1/(cm ps)]",
                            "# Ab[0][0], [0][1],..[0][nt-1]",
                            "# Ab[1][0], [1][1],..[1][nt-1]",
                            "# ...",
                            "# Ab[nz-1][0], [nz-1][1],..[nz-1][nt-1]",
                            "Ab_zt");
    }
    else {
        // skip A_z line
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t iz = 0; iz < Nz; iz++) {
        for (std::size_t it = 0; it < Nt; it++) {
            if (mode == IoMode::Write) {
                file << std::format("{:12.4E} ", tracer.A.zt[iz][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.A.zt[iz][it];
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Rate of absorption per unit volume, per unit time [1/(cm³ ps]
 ****/
void IOA_rzt(std::fstream& file, std::size_t Nr, std::size_t Nz, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    IOAb_zt(file, Nz, Nt, tracer, mode);

    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# A[r][z][t]. [1/(cm³ ps)]",
                            "# A[0][0][0], [0][0][1],..[0][0][nt-1]",
                            "# A[0][1][0], [0][1][1],..[0][1][nt-1]",
                            "# ...",
                            "# A[nr-1][nz-1][0], [nr-1][nz-1][1],..[nr-1][nz-1][nt-1]",
                            "A_rzt");
    }
    else {
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t iz = 0; iz < Nz; iz++) {
            for (std::size_t it = 0; it < Nt; it++) {
                if (mode == IoMode::Write) {
                    file << std::format("{:12.4E} ", tracer.A.rzt[ir][iz][it]);
                    if (++i % 5 == 0) {
                        file << std::format("\n");
                    }
                }
                else {
                    file >> tracer.A.rzt[ir][iz][it];
                }
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Ballistic absorption per unit depth [1/cm] 
 ****/
void IOAb_z(std::fstream& file, std::size_t Nz, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("Ab_z #Ab[0], [1],..Ab[nz-1]. [1/cm]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t iz = 0; iz < Nz; iz++) {
        if (mode == IoMode::Write) {
            file << std::format("{:12.4E}\n", tracer.A.bz[iz]);
        }
        else {
            file >> tracer.A.bz[iz];
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n");
    }
}

/*******************************************************************************
 *  Rate of absorption per unit volume [1/cm³] 
 ****/
void IOA_rz(std::fstream& file, std::size_t Nr, std::size_t Nz, Tracer& tracer, IoMode mode)
{
    IOAb_z(file, Nz, tracer, mode);

    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n",
                            "# A[r][z]. [1/cm³]",
                            "# A[0][0], [0][1],..[0][nz-1]",
                            "# ...",
                            "# A[nr-1][0], [nr-1][1],..[nr-1][nz-1]",
                            "A_rz");
    }
    else {
        // skip A_rz line.
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t iz = 0; iz < Nz; iz++) {
            if (mode == IoMode::Write) {
                file << std::format("{:12.4E} ", tracer.A.rz[ir][iz]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.A.rz[ir][iz];
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Rate of absorption per unit time [1/(cm ps)] 
 ****/
void IOA_zt(std::fstream& file, std::size_t Nz, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# A[z][t]. [1/(cm ps)]",
                            "# A[0][0], [0][1],..[0][nt-1]",
                            "# A[1][0], [1][1],..[1][nt-1]",
                            "# ...",
                            "# A[nz-1][0], [nz-1][1],..[nz-1][nt-1]",
                            "A_zt");
    }
    else {
        // skip A_zt line
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t iz = 0; iz < Nz; iz++) {
        for (std::size_t it = 0; it < Nt; it++) {
            if (mode == IoMode::Write) {
                file << std::format("{:12.4E} ", tracer.A.zt[iz][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.A.zt[iz][it];
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Absorption per unit depth [1/cm] 
 ****/
void IOA_z(std::fstream& file, std::size_t Nz, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("A_z #A[0], [1],..A[nz-1]. [1/cm]\n");
    }
    else {
        // skip A_z line
        FindDataLine(file);
    }

    for (std::size_t iz = 0; iz < Nz; iz++) {
        if (mode == IoMode::Write) {
            file << std::format("{:12.4E}\n", tracer.A.z[iz]);
        }
        else {
            file >> tracer.A.z[iz];
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n");
    }
}

/*******************************************************************************
 *  Absorption per unit time [1/ps] 
 ****/
void IOA_t(std::fstream& file, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("A_t #A[0], [1],..A[nt-1]. [1/ps]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t it = 0; it < Nt; it++) {
        if (mode == IoMode::Write) {
            file << std::format("{:12.4E}\n", tracer.A.t[it]);
        }
        else {
            file >> tracer.A.t[it];
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance per unit area, per unit solid angle, per unit time [1/(cm² sr ps)]
 ****/
void IORd_rat(std::fstream& file, std::size_t Nr, std::size_t Na, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Rd[r][a][t]. [1/(cm² sr ps)]",
                            "# Rd[0][0][0], [0][0][1],..[0][0][nt-1]",
                            "# Rd[0][1][0], [0][1][1],..[0][1][nt-1]",
                            "# ...",
                            "# Rd[nr-1][na-1][0], [nr-1][na-1][1],..[nr-1][na-1][nt-1]",
                            "Rd_rat");
    }
    else {
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t ia = 0; ia < Na; ia++) {
            for (std::size_t it = 0; it < Nt; it++) {
                if (mode == IoMode::Write) {
                    file << std::format("{:12.4E} ", tracer.R.rat[ir][ia][it]);
                    if (++i % 5 == 0) {
                        file << std::format("\n");
                    }
                }
                else {
                    file >> tracer.R.rat[ir][ia][it];
                }
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance per unit area, per unit solid angle [1/(cm² sr)]
 ****/
void IORd_ra(std::fstream& file, std::size_t Nr, std::size_t Na, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Rd[r][angle]. [1/(cm² sr)].",
                            "# Rd[0][0], [0][1],..[0][na-1]",
                            "# Rd[1][0], [1][1],..[1][na-1]",
                            "# ...",
                            "# Rd[nr-1][0], [nr-1][1],..[nr-1][na-1]",
                            "Rd_ra");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t ia = 0; ia < Na; ia++) {
            if (mode == IoMode::Write) {
                file << std::format("{:12.4E} ", tracer.R.ra[ir][ia]);
                if ((ir * Na + ia + 1) % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.R.ra[ir][ia];
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance per unit solid angle, per unit time [1/sr ps]
 ****/
void IORd_rt(std::fstream& file, std::size_t Nr, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Rd[r][t]. [1/(cm² ps)]",
                            "# Rd[0][0], [0][1],..[0][nt-1]",
                            "# Rd[0][0], [0][1],..[0][nt-1]",
                            "# ...",
                            "# Rd[nr-1][0], [nr-1][1],..[nr-1][nt-1]",
                            "Rd_rt");
    }
    else {
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t it = 0; it < Nt; it++) {
            if (mode == IoMode::Write) {
                file << std::format("{:12.4E} ", tracer.R.rt[ir][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.R.rt[ir][it];
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance per unit area, per unit time [1/cm² ps]
 ****/
void IORd_at(std::fstream& file, std::size_t Na, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Rd[a][t]. [1/(sr ps)]",
                            "# Rd[0][0], [0][1],..[0][nt-1]",
                            "# Rd[1][0], [1][1],..[1][nt-1]",
                            "# ...",
                            "# Rd[na-1][0], [na-1][1],..[na-1][nt-1]",
                            "Rd_at");
    }
    else {
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t ia = 0; ia < Na; ia++) {
        for (std::size_t it = 0; it < Nt; it++) {
            if (mode == IoMode::Write) {
                file << std::format("{:12.4E} ", tracer.R.at[ia][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.R.at[ia][it];
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance distribution per unit area [1/cm²]
 ****/
void IORd_r(std::fstream& file, std::size_t Nr, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("Rd_r #Rd[0], [1],..Rd[nr-1]. [1/cm²]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ir = 0; ir < Nr; ir++) {
        if (mode == IoMode::Write) {
            file << std::format("{:12.4E}\n", tracer.R.r[ir]);
        }
        else {
            file >> tracer.R.r[ir];
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance per unit solid angle [1/sr]
 ****/
void IORd_a(std::fstream& file, std::size_t Na, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("Rd_a #Rd[0], [1],..Rd[na-1]. [1/sr]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ia = 0; ia < Na; ia++) {
        if (mode == IoMode::Write) {
            file << std::format("{:12.4E}\n", tracer.R.a[ia]);
        }
        else {
            file >> tracer.R.a[ia];
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance per unit time [1/ps]
 ****/
void IORd_t(std::fstream& file, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("Rd_t #Rd[0], [1],..Rd[nt-1]. [1/ps]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t it = 0; it < Nt; it++) {
        if (mode == IoMode::Write) {
            file << std::format("{:12.4E}\n", tracer.R.t[it]);
        }
        else {
            file >> tracer.R.t[it];
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n");
    }
}

/*******************************************************************************
 *  Diffuse transmittance per unit area, per unit solid angle, per unit time [1/(cm² sr ps)]
 ****/
void IOTd_rat(std::fstream& file, std::size_t Nr, std::size_t Na, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Td[r][a][t]. [1/(cm² sr ps)]",
                            "# Td[0][0][0], [0][0][1],..[0][0][nt-1]",
                            "# Td[0][1][0], [0][1][1],..[0][1][nt-1]",
                            "# ...",
                            "# Td[nr-1][na-1][0], [nr-1][na-1][1],..[nr-1][na-1][nt-1]",
                            "Td_rat");
    }
    else {
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t ia = 0; ia < Na; ia++) {
            for (std::size_t it = 0; it < Nt; it++) {
                if (mode == IoMode::Write) {
                    file << std::format("{:12.4E} ", tracer.T.rat[ir][ia][it]);
                    if (++i % 5 == 0) {
                        file << std::format("\n");
                    }
                }
                else {
                    file >> tracer.T.rat[ir][ia][it];
                }
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Diffuse transmittance per unit area, per unit solid angle [1/(cm² sr)]
 ****/
void IOTd_ra(std::fstream& file, std::size_t Nr, std::size_t Na, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Td[r][angle]. [1/(cm² sr)].",
                            "# Td[0][0], [0][1],..[0][na-1]",
                            "# Td[1][0], [1][1],..[1][na-1]",
                            "# ...",
                            "# Td[nr-1][0], [nr-1][1],..[nr-1][na-1]",
                            "Td_ra");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t ia = 0; ia < Na; ia++) {
            if (mode == IoMode::Write) {
                file << std::format("{:12.4E} ", tracer.T.ra[ir][ia]);
                if ((ir * Na + ia + 1) % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.T.ra[ir][ia];
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Diffuse transmittance per unit solid angle, per unit time [1/sr ps]
 ****/
void IOTd_rt(std::fstream& file, std::size_t Nr, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Td[r][t]. [1/(cm² ps)]",
                            "# Td[0][0], [0][1],..[0][nt-1]",
                            "# Td[0][0], [0][1],..[0][nt-1]",
                            "# ...",
                            "# Td[nr-1][0], [nr-1][1],..[nr-1][nt-1]",
                            "Td_rt");
    }
    else {
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t ir = 0; ir < Nr; ir++) {
        for (std::size_t it = 0; it < Nt; it++) {
            if (mode == IoMode::Write) {
                file << std::format("{:12.4E} ", tracer.T.rt[ir][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.T.rt[ir][it];
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Diffuse transmittance per unit area, per unit time [1/cm² ps]
 ****/
void IOTd_at(std::fstream& file, std::size_t Na, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Td[a][t]. [1/(sr ps)]",
                            "# Td[0][0], [0][1],..[0][nt-1]",
                            "# Td[1][0], [1][1],..[1][nt-1]",
                            "# ...",
                            "# Td[na-1][0], [na-1][1],..[na-1][nt-1]",
                            "Td_at");
    }
    else {
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t ia = 0; ia < Na; ia++) {
        for (std::size_t it = 0; it < Nt; it++) {
            if (mode == IoMode::Write) {
                file << std::format("{:12.4E} ", tracer.T.at[ia][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.T.at[ia][it];
            }
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance per unit area [1/cm²] 
 ****/
void IOTd_r(std::fstream& file, std::size_t Nr, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("Td_r #Td[0], [1],..Td[nr-1]. [1/cm²]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ir = 0; ir < Nr; ir++) {
        if (mode == IoMode::Write) {
            file << std::format("{:12.4E}\n", tracer.T.r[ir]);
        }
        else {
            file >> tracer.T.r[ir];
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance per unit solid angle [1/sr] 
 ****/
void IOTd_a(std::fstream& file, std::size_t Na, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("Td_a #Td[0], [1],..Td[na-1]. [1/sr]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ia = 0; ia < Na; ia++) {
        if (mode == IoMode::Write) {
            file << std::format("{:12.4E}\n", tracer.T.a[ia]);
        }
        else {
            file >> tracer.T.a[ia];
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n");
    }
}

/*******************************************************************************
 *  Diffuse reflectance per unit time [1/ps] 
 ****/
void IOTd_t(std::fstream& file, std::size_t Nt, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        file << std::format("Td_t #Rd[0], [1],..Td[nt-1]. [1/ps]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t it = 0; it < Nt; it++) {
        if (mode == IoMode::Write) {
            file << std::format("{:12.4E}\n", tracer.T.t[it]);
        }
        else {
            file >> tracer.T.t[it];
        }
    }

    if (mode == IoMode::Write) {
        file << std::format("\n");
    }
}

/*******************************************************************************
 *  IoMode::Read: read result back from a output input.
 *  IoMode::Write: write result to a output input.
 ****/
void ReadWriteResults(std::fstream& file, RunParams& params, Tracer& tracer, IoMode mode)
{
    if (mode == IoMode::Write) {
        if (params.output_file_format == FileFormat::Ascii) {
            WriteVersion(file, "mcmloA2.0");
        }
        else {
            WriteVersion(file, "mcmloB2.0");
        }

        WriteInputParams(file, params);
        WriteRandomizer(file);
        WriteRAT(file, tracer);
    }
    else {
        ReadRandomizer(file);
        ReadRAT(file, tracer);
    }

    // Absorption
    if (params.record.A_rzt) {
        IOA_rzt(file, params.num_r, params.num_z, params.num_t, tracer, mode);
    }
    if (params.record.A_rz) {
        IOA_rz(file, params.num_r, params.num_z, tracer, mode);
    }
    if (params.record.A_zt) {
        IOA_zt(file, params.num_z, params.num_t, tracer, mode);
    }
    if (params.record.A_z) {
        IOA_z(file, params.num_z, tracer, mode);
    }
    if (params.record.A_t) {
        IOA_t(file, params.num_t, tracer, mode);
    }

    // Reflectance
    if (params.record.Rd_rat) {
        IORd_rat(file, params.num_r, params.num_a, params.num_t, tracer, mode);
    }
    if (params.record.Rd_ra) {
        IORd_ra(file, params.num_r, params.num_a, tracer, mode);
    }
    if (params.record.Rd_rt) {
        IORd_rt(file, params.num_r, params.num_t, tracer, mode);
    }
    if (params.record.Rd_at) {
        IORd_at(file, params.num_a, params.num_t, tracer, mode);
    }
    if (params.record.Rd_r) {
        IORd_r(file, params.num_r, tracer, mode);
    }
    if (params.record.Rd_a) {
        IORd_a(file, params.num_a, tracer, mode);
    }
    if (params.record.Rd_t) {
        IORd_t(file, params.num_t, tracer, mode);
    }

    // Transmittance
    if (params.record.Td_rat) {
        IOTd_rat(file, params.num_r, params.num_a, params.num_t, tracer, mode);
    }
    if (params.record.Td_ra) {
        IOTd_ra(file, params.num_r, params.num_a, tracer, mode);
    }
    if (params.record.Td_rt) {
        IOTd_rt(file, params.num_r, params.num_t, tracer, mode);
    }
    if (params.record.Td_at) {
        IOTd_at(file, params.num_a, params.num_t, tracer, mode);
    }
    if (params.record.Td_r) {
        IOTd_r(file, params.num_r, tracer, mode);
    }
    if (params.record.Td_a) {
        IOTd_a(file, params.num_a, tracer, mode);
    }
    if (params.record.Td_t) {
        IOTd_t(file, params.num_t, tracer, mode);
    }

    file.close();
}

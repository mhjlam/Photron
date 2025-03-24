/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
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

        // Treat as a string instead
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

    // Treat as a string instead
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

std::string& ToUpper(std::string& str)
{
    std::ranges::transform(str, str.begin(), [](unsigned char c) {
        return std::toupper(c);
    });
    return str;
}


double RandomGen(char Type, long Seed, long* Status)
{
#define MBIG    1000000000
#define MSEED   161803398
#define MZ      0
#define FAC     1.0E-9

    // ma[0] is not used
    static long i1;
    static long i2;
    static long ma[56];

    // set seed
    if (Type == 0) {
        long mk = 1;
        long mj = MSEED - (Seed < 0 ? -Seed : Seed) % MBIG;

        ma[55] = mj;

        for (short i = 1; i <= 54; i++) {
            short ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ) {
                mk += MBIG;
            }
            mj = ma[ii];
        }
        for (short ii = 1; ii <= 4; ii++) {
            for (short i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ) {
                    ma[i] += MBIG;
                }
            }
        }
        i1 = 0;
        i2 = 31;
    }
    // get a number
    else if (Type == 1) {
        if (++i1 == 56) {
            i1 = 1;
        }

        if (++i2 == 56) {
            i2 = 1;
        }

        long mj = ma[i1] - ma[i2];
        if (mj < MZ) {
            mj += MBIG;
        }

        ma[i1] = mj;
        return (mj * FAC);
    }
    // get status
    else if (Type == 2) {
        for (short i = 0; i < 55; i++) {
            Status[i] = ma[i + 1];
        }
        Status[55] = i1;
        Status[56] = i2;
    }
    // restore status
    else if (Type == 3) {
        for (short i = 0; i < 55; i++) {
            ma[i + 1] = Status[i];
        }
        i1 = Status[55];
        i2 = Status[56];
    }
    else {
        puts("Wrong parameter to RandomGen().");
    }
    return (0);

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
}


/**************************************************************************
 *	Print messages about MCML.
 ****/
void AboutMCML(void)
{
    puts("MCML 2.0, Copyright (c) 1992-1996");
    puts("Monte Carlo Simulation of Light Transport in Multi-Layered Turbid Media");

    puts(" ");
    puts("Lihong Wang, Ph.D.");
    puts("Bioengineering Program, Texas A&M University");
    puts("College Station, Texas 77843-3120, USA");

    puts("Liqiong Zheng, B.S.");
    puts("Summer student from Dept. of Computer Science,");
    puts("University of Houston, Texas, USA.");

    puts("Steven L. Jacques, Ph.D.");
    puts("Oregon Medical Laser Center, Providence/St. Vincent Hospital");
    puts("9205 SW Barnes Rd., Portland, Oregon 97225, USA");

    puts(" ");
    puts("Obtain the program thru anonymous ftp to laser.mda.uth.tmc.edu");

    puts(" ");
    puts("Please cite the following article in your publications:");
    std::cout << "\tL.-H. Wang, S. L. Jacques, and L.-Q. Zheng, MCML - Monte " << std::endl;
    std::cout << "\tCarlo modeling of photon transport in multi-layered" << std::endl;
    std::cout << "\ttissues, Computer Methods and Programs in Biomedicine, 47," << std::endl;
    std::cout << "\t131-146 (1995)" << std::endl;
}

/**************************************************************************
 *	Kill the ith char (counting from 0), push the following
 *	chars forward by one.
 ****/
void KillChar(size_t i, std::string Str)
{
    size_t sl = Str.size();
    for (; i < sl; i++) {
        Str[i] = Str[i + 1];
    }
}

/**************************************************************************
 *	Eliminate the chars in a string which are not printing
 *	chars or spaces.
 *
 *	Spaces include ' ', '\f', '\t' etc.
 *
 *	Return 1 if no nonprinting chars found, otherwise
 *	return 0.
 ****/
bool CheckCharQ(std::string& Str)
{
    bool found = 0;	// found bad char.
    size_t sl = Str.size();
    size_t i = 0;

    while (i < sl) {
        if (Str[i] < 0 || Str[i] > 255) {
            std::cerr << "Non-ASCII file" << std::endl;
            return false;
        }
        else if (isprint(Str[i]) || isspace(Str[i])) {
            i++;
        }
        else {
            found = 1;
            KillChar(i, Str);
            sl--;
        }
    }

    return (found);
}

/**************************************************************************
 *	Return 1 if this line is a comment line in which the
 *	first non-space character is "#", or a space line.
 *	Return 0 otherwise.
 ****/
bool CommentLineQ(std::string& buf)
{
    auto it = std::ranges::find_if(buf, [](char c) {
        return !std::isspace(c);
    });
    return (it == buf.end() || *it == '#');
}

/**************************************************************************
 *	Skip space or comment lines and return a data line.
 ****/
std::string FindDataLine(std::fstream& file)
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

/**************************************************************************
 *	Skip space or comment lines and return a data line.
 ****/
std::string FindDataLine(std::istream& stream)
{
    std::string line;

    // skip space or comment lines.
    do {
        std::getline(stream, line);

        if (!line.empty()) {
            CheckCharQ(line);
        }
    } while (CommentLineQ(line));

    return line;
}

/**************************************************************************
 *	Check whether the input version is the same as version.
 ****/
bool CheckFileVersionQ(std::fstream& file, const std::string version)
{
    // Skip comment lines.
    std::string line;
    do { std::getline(file, line); } while (line.empty() || CommentLineQ(line));

    if (line.find(version) == std::string::npos) {
        std::cerr << "Invalid file version.";
        return false;
    }
    return true;
}

/**************************************************************************
 *  Get a filename and open it for reading, retry until the input can be
 *  opened with a correct version or a '.' is typed.
 *	Return a NULL pointer if '.' is typed.
 ****/
bool GetFile(std::string& fname, const std::string version, std::fstream& file)
{
    while (1) {
        // prompt.
        std::cout << "Specify filename (or . to quit to main menu):";

        // Read input buffer
        std::getline(std::cin, fname);

        if (!fname.empty()) {
            // terminate with a period.
            if (fname.size() == 1 && fname[0] == '.') {
                // return if '.' entered.
                return false;
            }

            // open the input & check the version.
            file = std::fstream(fname, std::ios::in);
            if (!file.is_open()) {
                // cannot open the input.
                std::cerr << "File does not exist.";
            }
            else {
                if (CheckFileVersionQ(file, version)) {
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
 *  Find number of mediums in the list.
 *  At the same time, check the optical parameters.
 ****/
int FindNumMediaQ(std::fstream& file)
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
 *  Read the parameters of one name, assumming the
 *  parameters have been checked with FindNumMediaQ().
 ****/
bool ReadOneMediumQ(std::fstream& file, std::vector<Layer>& media)
{
    std::string buf = FindDataLine(file);
    auto extracted = ParseLine(buf, "Error reading medium parameters.", 5);
    if (extracted.empty()) {
        return false;
    }

    std::string name = std::get<std::string>(extracted[0]);
    double n = std::get<double>(extracted[1]);
    double mua = std::get<double>(extracted[2]);
    double mus = std::get<double>(extracted[3]);
    double g = std::get<double>(extracted[4]);

    // Add new medium to the list.
    media.push_back(Layer{ .name = name, .eta = n, .mua = mua, .mus = mus, .aniso = g });
    return true;
}

/*******************************************************************************
 *  Read the mediums list.
 ****/
bool ReadMediumListQ(std::fstream& file, RunParams& params)
{
    // Get current output position
    std::streampos file_pos = file.tellg();

    int num_media = FindNumMediaQ(file);
    if (num_media < 1) {
        file.seekg(file_pos, std::ios::beg);
        return false;
    }

    // Seek to previous output position 
    file.seekg(file_pos, std::ios::beg);

    for (short i = 0; i < num_media; i++) {
        ReadOneMediumQ(file, params.mediums);
    }

    // skip the signal end.
    FindDataLine(file);

    return true;
}

/**************************************************************************
 *	Read the input name and the input format.
 *
 *	The input format can be either A for ASCII or B for binary.
 ****/
bool ReadFnameFormatQ(std::fstream& file, RunParams& params)
{
    std::string buf = FindDataLine(file);
    auto extracted = ParseLine(buf, "Error reading file name.", 1);
    if (extracted.empty()) {
        return false;
    }

    params.output_filename = std::get<std::string>(extracted[0]);

    //if (output_file_format == "B") {
    //    params.output_file_format = FileFormat::Binary;
    //}
    //else {
    //    params.output_file_format = FileFormat::ASCII;
    //}

    // Only support ASCII format.
    params.output_file_format = FileFormat::Ascii;

    return true;
}

/**************************************************************************
 *	Read the RunParams members grid_z, grid_r and grid_time.
 ****/
bool ReadDzDrDtQ(std::istream& input, RunParams& params)
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
    params.grid_time = dt;

    return true;
}

/**************************************************************************
 *	Read the RunParams members num_z, num_r, num_time, num_alpha.
 ****/
bool ReadNzNrNtNaQ(std::istream& input, RunParams& params)
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
    params.num_time = nt;
    params.num_alpha = na;
    params.grid_alpha = 0.5 * std::numbers::pi / params.num_alpha;

    return true;
}

/**************************************************************************
 *   Initialize the Record struct.
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

/**************************************************************************
 *  Read which quantity is to be scored.
 ****/
bool ReadRecordQ(std::istream& input, RunParams& params)
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

/**************************************************************************
*   Filter the Record struct.
****/
bool FilterRecordQ(std::istream& input, RunParams& params)
{
    InitRecord(params);

    if (!ReadRecordQ(input, params)) {
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

/**************************************************************************
 *  Read the threshold min_weight.
 ****/
bool ReadWthQ(std::fstream& file, RunParams& params)
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

    params.min_weight = wth;
    return true;
}

/**************************************************************************
 *  Read the seed for random number generator (unused).
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

/**************************************************************************
 *  Find number of layers.
 ****/
int FindNumLayersQ(std::fstream& file)
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

/**************************************************************************
 *  Check whether the name name is in the mediums list.
 ****/
bool ValidMediumNameQ(std::string& name, int& index, RunParams& params)
{
    for (short i = 0; i < params.mediums.size(); i++) {
        if (name == params.mediums[i].name) {
            index = i;
            return true;
        }
    }
    return false;
}

/**************************************************************************
 *	Read the parameters of all layers.
 ****/
bool ReadLayerSpecsQ(std::fstream& file, RunParams& params)
{
    std::string name;
    double thickness = 0.0;

    // z coordinate of the current layer.
    double z = 0.0;

    // Save current output position
    std::streampos file_pos = file.tellg();

    int num_layers = FindNumLayersQ(file);
    if (num_layers < 1) {
        file.seekg(file_pos, std::ios::beg);
        return false;
    }

    // Seek to previous output position 
    file.seekg(file_pos, std::ios::beg);

    // current_layer 0 and current_layer Num_Layers + 1 are for ambient.
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
        if (!ValidMediumNameQ(name, index, params)) {
            std::cerr << "  Invalid medium name. " << std::endl;
            return false;
        }

        params.layers[i].name = params.mediums[index].name;
        params.layers[i].eta = params.mediums[index].eta;
        params.layers[i].mua = params.mediums[index].mua;
        params.layers[i].mus = params.mediums[index].mus;
        params.layers[i].aniso = params.mediums[index].aniso;

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

/**************************************************************************
 *  Read the number of photons.
 *  Read computation time limit.
 *  type = 0, read from a .mci input file;
 *  type = 1, read from a .mco output file.
 ****/
bool ReadNumPhotonsQ(std::istream& input, RunParams& params, char type)
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
    else if (extracted.size() == 1 && IsString(extracted[0])) {
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
    else if (extracted.size() == 2 && IsNumber(extracted[0]) && IsString(extracted[1])) {
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
            if (type == 0) {
                params.num_photons = num_photons;
                params.time_limit = hours * 3600 + minutes * 60;
            }
            else {
                params.add_num_photons = num_photons;
                params.add_limit = hours * 3600 + minutes * 60;
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

/**************************************************************************
 *  Read the beam source type (Pencil/Isotropic).
 ****/
bool ReadSourceTypeQ(std::fstream& file, RunParams& params)
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

/**************************************************************************
 *  Compute the index to layer according to the z coordinate.
 *	If the z is on an interface between layers, the returned index
 *	will point to the upper layer.
 *	Index 0 is the top ambient name and index num_layers+1 is the bottom one.
 ****/
bool ZToLayerQ(double z, short& index, RunParams& params)
{
    // index to current_layer.
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

/**************************************************************************
 *  Read starting position of photon source.
 ****/
bool ReadStartPQ(std::istream& input, RunParams& params)
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

        if (!ZToLayerQ(source_z, params.source_layer, params)) {
            return false;
        }
    }
    else if (extracted.size() == 2) {
        std::string medium_name = std::get<std::string>(extracted[1]);

        if (medium_name[0] != '#' && medium_name[0] != '\n') {
            short source_layer;
            if (!ZToLayerQ(source_z, source_layer, params)) {
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

/*************************************************************************
 *	Compute the critical angles for total internal reflection according to
 *  the relative refractive index of the current_layer.
 *	All layers are processed.
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

/**************************************************************************
 *	Read in the input parameters for one run.
 ****/
void ReadRunParam(std::fstream& file, RunParams& params)
{
    if (!ReadFnameFormatQ(file, params)) {
        exit(1);
    }

    // Geometry.
    if (!ReadLayerSpecsQ(file, params)) {
        exit(1);
    }

    // Skip the signal "end" of layers.
    FindDataLine(file);

    // Source.
    if (!ReadSourceTypeQ(file, params)) {
        exit(1);
    }
    if (!ReadStartPQ(file, params)) {
        exit(1);
    }

    // Grids.
    if (!ReadDzDrDtQ(file, params)) {
        exit(1);
    }
    if (!ReadNzNrNtNaQ(file, params)) {
        exit(1);
    }
    params.max_z = params.grid_z * params.num_z;
    params.max_r = params.grid_r * params.num_r;
    params.max_time = params.grid_time * params.num_time;
    params.max_alpha = params.grid_alpha * params.num_alpha;

    // Scored data categories.
    if (!FilterRecordQ(file, params)) {
        exit(1);
    }

    // Simulation control.
    if (!ReadNumPhotonsQ(file, params, 0)) {
        exit(1);
    }
    if (!ReadWthQ(file, params)) {
        exit(1);
    }
    if (!ReadSeed(file, params)) {
        exit(1);
    }

    CriticalAngle(params.layers);
}

/**************************************************************************
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


    // Allocate space for the current_layer parameters.
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
            params.mediums[i].aniso = GetDoubleValue(input, params.mediums[i].aniso);
        } while (!input.has_value() && params.mediums[i].aniso < 0.0 && params.mediums[i].aniso > 1.0);

        std::cout << std::endl;
    }
}

/**************************************************************************
 *  Read the input name and the input format interactively.
 ****/
void InterReadFnameFormat(RunParams& params)
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

/**************************************************************************
 *	Read grid_z, grid_r, grid_time interactively.
 ****/
void InterReadDzDrDt(RunParams& params)
{
    do {
        std::cout << "Specify dz, dr, dt in one line" << std::endl;
        std::cout << "(all > 0.0 cm, e.g., 0.1 0.1 0.1): ";
    } while (!ReadDzDrDtQ(std::cin, params));

    std::cout << std::endl;
}

/**************************************************************************
 *      Read the RunParams members num_z, num_r, num_alpha interactively.
 ****/
void InterReadNzNrNtNa(RunParams& params)
{
    do {
        std::cout << "Specify nz, nr, nt, na in one line" << std::endl;
        std::cout << "(all > 0, e.g., 100 100 100 100): ";
    } while (!ReadNzNrNtNaQ(std::cin, params));

    params.grid_alpha = 0.5 * std::numbers::pi / params.num_alpha;
    std::cout << std::endl;
}

/**************************************************************************
 *	Read and filter the quantities to be scored interactively.
 ****/
void InterFilterRecord(RunParams& params)
{
    do {
        std::cout << "Select scored quantities from the following data categories:" << std::endl;
        std::cout << "\tRd_rat\t\t\tTd_rat\t\t\tA_rzt" << std::endl;
        std::cout << "\tRd_ra\tRd_rt\tRd_at\tTd_ra\tTd_rt\tRd_at\tA_rz\tA_zt" << std::endl;
        std::cout << "\tRd_r\tRd_a\tRd_t\tTd_r\tTd_a\tTd_t\tA_z\tA_t" << std::endl;
    } while (!FilterRecordQ(std::cin, params));

    std::cout << std::endl;
}

/**************************************************************************
 *	Read the threshold min_weight interactively.
 ****/
void InterReadWth(RunParams& params)
{
    std::cout << "Input threshold weight (0 <= wth < 1.0, 0.0001 recommended): ";

    std::string string;
    std::getline(std::cin, string);
    std::istringstream iss(string);

    // TODO: use ParseInput

    while (!(iss >> params.min_weight) || params.min_weight < 0 || params.min_weight >= 1) {
        std::cout << "Invalid wth. Input again (0 <= wth < 1.0): ";
        std::getline(std::cin, string);
    }

    std::cout << std::endl;
}

/**************************************************************************
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

/**************************************************************************
 *	Read current_layer specifications interactively.
 ****/
void InterReadLayerSpecs(RunParams& params)
{
    std::string name;
    int num_layers;

    // Z coordinate of the current layer.
    double z = 0.0;

    int index;

    std::cout << "\nSpecify layer list. ";
    PrintMediumNames(params);

    std::cout << "\nTotal number of layers: ";

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

            if (!ValidMediumNameQ(name, index, params)) {
                std::cout << "  Invalid medium name. Input again.";
                error = 1;
            }
        }

        params.layers.push_back(Layer{
            .name = params.mediums[index].name,
            .eta = params.mediums[index].eta,
            .mua = params.mediums[index].mua,
            .mus = params.mediums[index].mus,
            .aniso = params.mediums[index].aniso });


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

/**************************************************************************
 *  Read the number of photons, or computation time interactively.
 ****/
void InterReadNumPhotons(RunParams& params)
{
    std::cout << "Specify number of photons or time in hh:mm format," << std::endl;
    std::cout << "or both in one line (e.g. 10000 5:30): ";

    while (!ReadNumPhotonsQ(std::cin, params, 0)) {
        std::cout << "Input again: ";
    }

    std::cout << std::endl;
}

/**************************************************************************
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

/**************************************************************************
 *  Read starting position of photon source.
 ****/
void InterReadStartP(RunParams& params)
{
    do {
        std::cout << "Input the z coordinate of source (0.0 - " << params.layers[params.num_layers].bot_z << " cm) and the medium" << std::endl;
        std::cout << "where the source is if the z is on an interface (e.g. 1.0 [air]):";
    } while (!ReadStartPQ(std::cin, params));

    std::cout << std::endl;
}

/*************************************************************************
 *  If input is stdin, freeze the screen and print a more message on screen
 *  every 20 lines. The Line is the line index.
 ****/
void More(std::ostream& output, int& Line)
{
    if (&output == &std::cout) {
        if (!(Line % 20)) {
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

/*************************************************************************
 * Write name list to the output.
 * If input is stdout, freeze the screen every 20 lines.
 ****/
void PutMediumListToFile(std::ostream& output, RunParams& params, int& Line)
{
    std::string format;

    output << std::format("# Specify media \n");
    Line++;
    output << std::format("#\tname\t\tn\tmua\tmus\tg\n");
    Line++;

    for (int i = 0; i < params.mediums.size(); i++) {
        More(output, Line);
        Layer s = params.mediums[i];
        if (s.name.size() + 1 > 8) {
            output << std::format("\t{}\t{:G}\t{:G}\t{:G}\t{:G}\n", s.name, s.eta, s.mua, s.mus, s.aniso);
        }
        else {
            output << std::format("\t{}\t\t{:G}\t{:G}\t{:G}\t{:G}\n", s.name, s.eta, s.mua, s.mus, s.aniso);
        }
        Line++;
    }
    output << std::format("end #of media\n");
    Line++;
}

void PutFnameFormatToFile(std::ostream& output, RunParams& params, int& Line)
{
    More(output, Line);
    output << std::format("{} \t{}\t\t\t# output file name, format.\n", params.output_filename, 'A');
    Line++;
}

void PutDzDrDtToFile(std::ostream& output, RunParams& params, int& Line)
{
    More(output, Line);
    output << std::format("{:G}\t{:G}\t{:G}\t\t\t# dz, dr, dt.\n", params.grid_z, params.grid_r, params.grid_time);
    Line++;
}

void PutNzNrNtNaToFile(std::ostream& output, RunParams& params, int& Line)
{
    More(output, Line);
    output << std::format("{:d}\t{:d}\t{:d}\t{:d}\t\t# nz, nr, nt, na.\n", params.num_z, params.num_r, params.num_time, params.num_alpha);
    Line++;
}

void PutScoredToFile(std::ostream& output, RunParams& params, int& Line)
{
    More(output, Line);
    output << std::format("# This simulation will score the following categories:\n");
    Line++;

    More(output, Line);
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
    Line++;
}

void PutWthToFile(std::ostream& output, RunParams& params, int& Line)
{
    More(output, Line);
    output << std::format("{:G}\t\t\t\t\t# threshold weight.\n", params.min_weight);
    Line++;
}

void PutSeedToFile(std::ostream& output, RunParams& params, int& Line)
{
    More(output, Line);
    output << std::format("{}\t\t\t\t\t# random number seed.\n", params.seed);
    Line++;
}

void PutLayerSpecsToFile(std::ostream& output, RunParams& params, int& Line)
{
    std::string format;

    More(output, Line);
    output << std::format("# \tmedium \t\tthickness\n");
    Line++;

    for (int i = 0; i <= params.num_layers + 1; i++) {
        Layer s;
        More(output, Line);

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
        Line++;
    }

    More(output, Line);
    output << std::format("end #of layers\n");
    Line++;
}

void PutNumPhotonsToFile(std::ostream& output, RunParams& params, int& Line)
{
    More(output, Line);

    if (params.control_bit == ControlBit::NumPhotons) {
        output << std::format("{}  \t\t\t\t\t# no. of photons | time\n", params.num_photons);
    }
    else if (params.control_bit == ControlBit::TimeLimit) {
        output << std::format("{}:{}\t\t\t\t\t# no. of photons | time\n", params.time_limit / 3600, params.time_limit % 3600 / 60);
    }
    else {
        output << std::format("{}  \t{}:{}\t\t\t\t# no. of photons | time\n", params.num_photons, params.time_limit / 3600, params.time_limit % 3600 / 60);
    }

    Line++;
}

void PutSourceTypeToFile(std::ostream& output, RunParams& params, int& Line)
{
    More(output, Line);

    if (params.source == BeamType::Pencil) {
        output << std::format("pencil \t\t\t\t\t# src type: pencil/isotropic.\n");
    }
    else {
        output << std::format("isotropic \t\t\t\t# src type: pencil/isotropic.\n");
    }

    Line++;
}

void PutStartPToFile(std::ostream& output, RunParams& params, int& Line)
{
    More(output, Line);

    if (params.source_medium_name.empty()) {
        output << std::format("{:G}\t\t\t\t\t# starting position of source.\n", params.source_z);
    }
    else if (params.source_medium_name.size() + 1 > 8) {
        output << std::format("{:G}\t{} \t\t\t# starting position of source.\n", params.source_z, params.source_medium_name);
    }
    else {
        output << std::format("{:G}\t{} \t\t\t\t# starting position of source.\n", params.source_z, params.source_medium_name);
    }

    Line++;
}

/*************************************************************************
 *  Write input parameters to the input output.
 *  If input is stdout, freeze the screen every 20 lines.
 ****/
void PutInputToFile(std::ostream& output, RunParams& params)
{
    // line index.
    output << std::format("mcmli2.0 \t\t\t# file version \n\n");
    int line = 2;
    PutMediumListToFile(output, params, line);

    More(output, line);
    output << std::format("\n# Specify data for run 1\n");
    line += 2;

    PutFnameFormatToFile(output, params, line);

    // geometry.
    More(output, line);
    output << std::format("\n");
    line++;
    PutLayerSpecsToFile(output, params, line);

    // source.
    More(output, line);
    output << std::format("\n");
    line++;
    PutSourceTypeToFile(output, params, line);
    PutStartPToFile(output, params, line);

    // grids.
    More(output, line);
    output << std::format("\n");
    line++;
    PutDzDrDtToFile(output, params, line);
    PutNzNrNtNaToFile(output, params, line);

    // scored data categories.
    More(output, line);
    output << std::format("\n");
    line++;
    PutScoredToFile(output, params, line);

    // simulation control.
    More(output, line);
    output << std::format("\n");
    line++;

    PutNumPhotonsToFile(output, params, line);
    PutWthToFile(output, params, line);
    PutSeedToFile(output, params, line);

    More(output, line);
    output << std::format("end #of runs\n\n");
}

/**************************************************************************
 *  Read in the input parameters for one run in interactive mode.
 ****/
void InterReadParam(RunParams& params)
{
    InterReadMediumList(params);
    InterReadFnameFormat(params);
    InterReadLayerSpecs(params);
    InterReadSourceType(params);
    InterReadStartP(params);
    InterReadDzDrDt(params);
    InterReadNzNrNtNa(params);
    InterFilterRecord(params);
    InterReadNumPhotons(params);
    InterReadWth(params);

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
            PutInputToFile(file, params);
        }
    }
}

/**************************************************************************
 *  Check consistance of input parameters for one run.
 *  Such as: the consistance of name list, current_layer
 *  list, souce starting position and source type.
 ****/
bool CheckInputConsis(RunParams& params)
{
    for (int i = 0; i <= params.num_layers + 1; i++) {
        int index;
        if (!ValidMediumNameQ(params.layers[i].name, index, params)) {
            std::cout << "Invalid medium name of layer " << i << "." << std::endl;
            return 0;
        }
        else {
            params.layers[i].eta = params.mediums[index].eta;
            params.layers[i].mua = params.mediums[index].mua;
            params.layers[i].mus = params.mediums[index].mus;
            params.layers[i].aniso = params.mediums[index].aniso;
        }
    }

    if ((params.source == BeamType::Isotropic) && (params.source_z == 0.0)) {
        std::cout << "Can not put isotropic source in upper ambient medium." << std::endl;
        return 0;
    }
    if (!ZToLayerQ(params.source_z, params.source_layer, params)) {
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

/****************************************************************************
 *    Menu for changing input parameters.
 *****/
void ShowChangeMenu()
{
    puts("  o = Print the input on screen.");
    puts("  m = Change media list.");
    puts("  f = Change output file name and format.");
    puts("  d = Change dz, dr, dt.");
    puts("  n = Change nz, nr, nt, na.");
    puts("  c = Change scored data categories.");
    puts("  w = Change threshold weight.");
    puts("  r = Change random number seed.");
    puts("  l = Change layer specifications.");
    puts("  p = Change photon number and computation time limit.");
    puts("  s = Change source type.");
    puts("  z = Change source starting position.");
    puts("  q = Quit from change menu and start simulation.");
    puts("  x = Exit to the main menu.");
    puts("  * Commands here are not case-sensitive");
}

/***********************************************************************
 *   Continue to change input parameters or quit.
 *****/
char QuitOrContinue()
{
    std::string string;

    do {
        std::cout << "Do you want to change them? (Y/N): ";
        do {
            std::getline(std::cin, string);
        } while (!string.empty());
    } while (toupper(string[0]) != 'Y' && toupper(string[0]) != 'N');

    return (toupper(string[0]));
}

void ChangeMediumList(RunParams& params)
{
    int line = 1;
    std::cout << "Current medium list: " << std::endl;
    PutMediumListToFile(std::cout, params, line);
    std::cout << std::endl;

    if (QuitOrContinue() == 'Y') {
        params.mediums.clear();
        InterReadMediumList(params);
    }
}

void ChangeFnameFormat(RunParams& params)
{
    int line = 1;
    std::cout << "Current output file name and format: " << std::endl;
    PutFnameFormatToFile(std::cout, params, line);
    std::cout << std::endl;
    InterReadFnameFormat(params);
}

void ChangeDzDrDt(RunParams& params)
{
    int line = 1;
    std::cout << "Current dz, dr, dt: " << std::endl;
    PutDzDrDtToFile(std::cout, params, line);
    std::cout << std::endl;
    InterReadDzDrDt(params);
}

void ChangeNzNrNtNa(RunParams& params)
{
    int line = 1;
    std::cout << "Current nz, nr, nt, na: " << std::endl;
    PutNzNrNtNaToFile(std::cout, params, line);
    std::cout << std::endl;
    InterReadNzNrNtNa(params);
}

void ChangeRecord(RunParams& params)
{
    int line = 1;
    PutScoredToFile(std::cout, params, line);
    std::cout << std::endl;

    if (QuitOrContinue() == 'Y') {
        InterFilterRecord(params);
    }
}

void ChangeWth(RunParams& params)
{
    int line = 1;
    std::cout << "Current threshold weight: " << std::endl;
    PutWthToFile(std::cout, params, line);
    std::cout << std::endl;
    InterReadWth(params);
}

void ChangeLayerSpecs(RunParams& params)
{
    int line = 1;
    std::cout << "Current layer sepcification: " << std::endl;
    PutLayerSpecsToFile(std::cout, params, line);
    std::cout << std::endl;

    if (QuitOrContinue() == 'Y') {
        InterReadLayerSpecs(params);
    }
}

void ChangeNumPhotons(RunParams& params)
{
    int line = 1;
    std::cout << "Current value: " << std::endl;
    PutNumPhotonsToFile(std::cout, params, line);
    std::cout << std::endl;
    InterReadNumPhotons(params);
}

void ChangeSourceType(RunParams& params)
{
    int line = 1;
    std::cout << "Current source type: " << std::endl;
    PutSourceTypeToFile(std::cout, params, line);
    std::cout << std::endl;
    InterReadSourceType(params);
}

void ChangeStartP(RunParams& params)
{
    int line = 1;
    std::cout << "Layer Specification: " << std::endl;
    PutLayerSpecsToFile(std::cout, params, line);
    std::cout << "\nCurrent starting position: " << std::endl;
    PutStartPToFile(std::cout, params, line);
    std::cout << std::endl;
    InterReadStartP(params);
}

/************************************************************************
 *  return 1 if string[0] = Q, quit change menu;
 *  return 2 if string[0] = X, quit to the main menu;
 *  return 0 otherwise.
 ****/
int BranchChangeMenu(std::string& string, RunParams& params)
{
    switch (toupper(string[0])) {
        case 'M':
            ChangeMediumList(params);
            break;

        case 'F':
            ChangeFnameFormat(params);
            break;

        case 'D':
            ChangeDzDrDt(params);
            break;

        case 'N':
            ChangeNzNrNtNa(params);
            break;

        case 'C':
            ChangeRecord(params);
            break;

        case 'W':
            ChangeWth(params);
            break;

        case 'L':
            ChangeLayerSpecs(params);
            break;

        case 'P':
            ChangeNumPhotons(params);
            break;

        case 'S':
            ChangeSourceType(params);
            break;

        case 'Z':
            ChangeStartP(params);
            break;

        case 'O':
            PutInputToFile(std::cout, params);
            break;

        case 'H':
            ShowChangeMenu();
            break;

        case 'Q':
            return 1;

        case 'X':
            return 2;

        default:
            puts("...Unknown command");
    }

    return 0;
}

/***************************************************************************
 *   Return 1 if quit change and start simulation;
 *   return 0 if exit to main menu.
 ****/
bool RunChangedInput(RunParams& params)
{
    std::string string;
    int branch;

    std::cout << "Any changes to the input parameters? (Y/N)";
    do {
        std::getline(std::cin, string);
    } while (string.empty());

    while (toupper(string[0]) == 'Y') {
        do {
            do {
                std::cout << "\n> Change menu (h for help) => ";
                std::getline(std::cin, string);
            } while (string.empty());

            // string[0] is 'X' or 'Q'.
            if (branch = BranchChangeMenu(string, params)) {
                break;
            }
        } while (1);

        std::cout << "Do you want to save the input to a file? (Y/N)";
        std::getline(std::cin, string);
        if (toupper(string[0]) == 'Y') {
            std::cout << "Give the file name to save input: ( .mci): ";
            std::getline(std::cin, string);

            std::ofstream file(string, std::ios::out);
            if (!file.is_open()) {
                std::cout << "Can not open the file to write." << std::endl;
            }
            else {
                PutInputToFile(file, params);
            }
        }

        // quit change menu and start simulation.
        if (branch == 1) {
            if (!CheckInputConsis(params)) {
                do {
                    std::cout << "Change input or exit to main menu (c/x): ";
                    std::getline(std::cin, string);
                } while (!string.empty() || toupper(string[0]) != 'X' && toupper(string[0]) != 'C');

                if (toupper(string[0]) == 'X') {
                    params.mediums.clear();
                    params.layers.clear();
                    return false;
                }
                else {
                    string[0] = 'Y';	// continue to change parameters.
                }
            }
            else {
                return true;
            }

        }
        // exit to menu.
        else {
            params.mediums.clear();
            params.layers.clear();
            return false;
        }
    }

    return true;
}

/*******************************************************************************
 *  Check the whether the flag end is met.
 *  The input position is restored to the current position at the end of the inquery.
 ****/
bool EndOfRunsQ(std::istream& input)
{
    // found end of runs.
    bool end_found = false;

    // record input position.
    std::streampos file_pos = input.tellg();

    std::string buf = FindDataLine(input);
    if (buf.empty()) {
        end_found = true;
        std::cout << "Missing end." << std::endl;
    }
    else if (buf.find("end") != std::string::npos) {
        end_found = true;
    }

    // restore postion.
    input.seekg(file_pos, std::ios::beg);
    return end_found;
}

/**************************************************************************
 *  Check the input parameters for all runs.
 *  This function will count number of runs and assign it to params.num_runs.
 ****/
void CheckParamFromFile(std::fstream& input, RunParams& params)
{
    if (!ReadMediumListQ(input, params)) {
        exit(1);
    }

    // Save current position in input file.
    std::streampos file_pos = input.tellg();

    short run_index = 0;
    do {
        std::cout << "Checking input data for run " << ++run_index << std::endl;
        ReadRunParam(input, params);

        // Attempt insert and detect duplicates
        if (!params.unique_outputs.insert(params.output_filename).second) {
            std::cout << "File name " << params.output_filename << " duplicated." << std::endl;
            exit(1);
        }
    } while (!EndOfRunsQ(input));

    params.num_runs = run_index;

    // Restore file position
    input.seekg(file_pos, std::ios::beg);
}

/**************************************************************************
 *	Allocate the arrays in Tracer for one run, and
 *	array elements are automatically initialized to zeros.
 *
 *	Remember that the indices for Rd_r[], Td_r[],
 *	& A_rz[][iz] start from -1 storing the collimated
 *	responses.
 ****/
void InitOutputData(RunParams& params, Tracer& tracer)
{
    std::size_t nz = params.num_z;
    std::size_t nr = params.num_r;
    std::size_t na = params.num_alpha;
    std::size_t nt = params.num_time;

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
        return List<List<List<double>>>(x, List<List<double>>(y, List<double>(z, 0.0)));
    };

    auto alloc2 = [](std::size_t x, std::size_t y) {
        return List<List<double>>(x, List<double>(y, 0.0));
    };

    auto alloc1 = [](std::size_t x) {
        return List<double>(x, 0.0);
    };

    if (params.record.Rd_rat) { tracer.R.rat = alloc3(nr, na, nt); }
    if (params.record.Rd_ra) { tracer.R.ra = alloc2(nr, na); }
    if (params.record.Rd_rt) { tracer.R.rt = alloc2(nr, nt); }
    if (params.record.Rd_at) { tracer.R.at = alloc2(na, nt); }
    if (params.record.Rd_r) { tracer.R.r = alloc1(nr); }
    if (params.record.Rd_a) { tracer.R.a = alloc1(na); }
    if (params.record.Rd_t) { tracer.R.t = alloc1(nt); }

    if (params.record.Td_rat) { tracer.T.rat = alloc3(nr, na, nt); }
    if (params.record.Td_ra) { tracer.T.ra = alloc2(nr, na); }
    if (params.record.Td_rt) { tracer.T.rt = alloc2(nr, nt); }
    if (params.record.Td_at) { tracer.T.at = alloc2(na, nt); }
    if (params.record.Td_r) { tracer.T.r = alloc1(nr); }
    if (params.record.Td_a) { tracer.T.a = alloc1(na); }
    if (params.record.Td_t) { tracer.T.t = alloc1(nt); }

    if (params.record.A_rzt) { tracer.A.rzt = alloc3(nr, nz, nt); }
    if (params.record.A_rzt) { tracer.A.zt = alloc2(nz, nt); }
    if (params.record.A_rz) { tracer.A.rz = alloc2(nr, nz); }
    if (params.record.A_rz) { tracer.A.z = alloc1(nz); }
    if (params.record.A_zt) { tracer.A.zt = alloc2(nz, nt); }
    if (params.record.A_z) { tracer.A.z = alloc1(nz); }
    if (params.record.A_t) { tracer.A.t = alloc1(nt); }
}

void FreeData(RunParams& params, Tracer& tracer)
{
    if (params.record.Rd_rat)   { tracer.R.rat.clear(); }
    if (params.record.Rd_ra)    { tracer.R.ra.clear(); }
    if (params.record.Rd_rt)    { tracer.R.rt.clear(); }
    if (params.record.Rd_at)    { tracer.R.at.clear(); }
    if (params.record.Rd_r)     { tracer.R.r.clear(); }
    if (params.record.Rd_a)     { tracer.R.a.clear(); }
    if (params.record.Rd_t)     { tracer.R.t.clear(); }

    if (params.record.Td_rat)   { tracer.T.rat.clear(); }
    if (params.record.Td_ra)    { tracer.T.ra.clear(); }
    if (params.record.Td_rt)    { tracer.T.rt.clear(); }
    if (params.record.Td_at)    { tracer.T.at.clear(); }
    if (params.record.Td_r)     { tracer.T.r.clear(); }
    if (params.record.Td_a)     { tracer.T.a.clear(); }
    if (params.record.Td_t)     { tracer.T.t.clear(); }

    if (params.record.A_rzt)    { tracer.A.rzt.clear(); }
    if (params.record.A_rzt)    { tracer.A.zt.clear(); }
    if (params.record.A_rz)     { tracer.A.rz.clear(); }
    if (params.record.A_rz)     { tracer.A.z.clear(); }
    if (params.record.A_zt)     { tracer.A.zt.clear(); }
    if (params.record.A_z)      { tracer.A.z.clear(); }
    if (params.record.A_t)      { tracer.A.t.clear(); }

    params.layers.clear();
}

/**************************************************************************
 *	Scale Rd and Td properly.
 *	"a" stands for angle alpha.
 ****
 *	Scale Rd(r,a) and Td(r,a) by
 *      (area perpendicular to photon direction)
 *		x(solid angle)x(No. of photons).
 *	or
 *		[2*PI*r*grid_r*cos(a)]x[2*PI*sin(a)*grid_alpha]x[No. of photons]
 *	or
 *		[2*PI*PI*grid_r*grid_alpha*r*sin(2a)]x[No. of photons]
 ****
 *	Scale Rd(r) and Td(r) by
 *		(area on the surface)x(No. of photons).
 ****
 *	Scale Rd(a) and Td(a) by
 *		(solid angle) cos(a) x(No. of photons).
 ****
 *  Mode = 0, scale Rd and Td; Mode = 1, unscale Rd and Td.
 ****/
void ScaleRdTd(RunParams& params, Tracer& tracer, char Mode)
{
    std::size_t nr = params.num_r;
    std::size_t na = params.num_alpha;
    std::size_t nt = params.num_time;
    double dr = params.grid_r;
    double da = params.grid_alpha;
    double dt = params.grid_time;

    double scale1 = (double)params.num_photons;
    if (Mode == 0) {
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
        for (short it = 0; it < nt; it++) {
            // scale Rd_t.
            if (Mode == 0) {
                tracer.R.t[it] /= scale1;
            }
            // unscale Rd_t.
            else {
                tracer.R.t[it] *= scale1;
            }
        }
    }

    if (params.record.Td_t) {
        for (short it = 0; it < nt; it++) {
            // scale Td_t.
            if (Mode == 0) {
                tracer.T.t[it] /= scale1;
            }
            // unscale Rd_t.
            else {
                tracer.T.t[it] *= scale1;
            }
        }
    }

    scale1 = 2.0 * std::numbers::pi * dr * dr * params.num_photons;
    // area is 2*PI*[(ir+0.5)*grid_r]*grid_r.  ir + 0.5 to be added.

    if (params.record.Rd_r) {
        for (short ir = 0; ir < nr; ir++) {
            double scale2 = 1.0 / ((ir + 0.5) * scale1);
            // scale Rd_r.
            if (Mode == 0) {
                tracer.R.r[ir] *= scale2;
            }
            // unscale Rd_r.
            else {
                tracer.R.r[ir] /= scale2;
            }
        }
    }

    if (params.record.Td_r) {
        for (short ir = 0; ir < nr; ir++) {
            double scale2 = 1.0 / ((ir + 0.5) * scale1);
            // scale Td_r.
            if (Mode == 0) {
                tracer.T.r[ir] *= scale2;
            }
            // unscale Td_r.
            else {
                tracer.T.r[ir] /= scale2;
            }
        }
    }

    scale1 *= dt;
    if (params.record.Rd_rt) {
        for (short ir = 0; ir < nr; ir++) {
            for (short it = 0; it < nt; it++) {
                double scale2 = 1.0 / ((ir + 0.5) * scale1);
                // scale Rd_rt.
                if (Mode == 0) {
                    tracer.R.rt[ir][it] *= scale2;
                }
                // unscale Rd_rt.
                else {
                    tracer.R.rt[ir][it] *= scale2;
                }
            }
        }
    }

    if (params.record.Td_rt) {
        for (short ir = 0; ir < nr; ir++) {
            for (short it = 0; it < nt; it++) {
                double scale2 = 1.0 / ((ir + 0.5) * scale1);
                // scale Td_rt.
                if (Mode == 0) {
                    tracer.T.rt[ir][it] *= scale2;
                }
                // unscale Td_rt.
                else {
                    tracer.T.rt[ir][it] /= scale2;
                }
            }
        }
    }

    scale1 = std::numbers::pi * da * params.num_photons;
    // solid angle times cos(a) is PI*sin(2a)*grid_alpha. sin(2a) to be added.

    if (params.record.Rd_a) {
        for (short ia = 0; ia < na; ia++) {
            double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
            // scale Rd_a.
            if (Mode == 0) {
                tracer.R.a[ia] *= scale2;
            }
            // unscale Rd_a.
            else {
                tracer.R.a[ia] /= scale2;
            }
        }
    }

    if (params.record.Td_a) {
        for (short ia = 0; ia < na; ia++) {
            double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
            // scale Td_a.
            if (Mode == 0) {
                tracer.T.a[ia] *= scale2;
            }
            // unscale Td_a.
            else {
                tracer.T.a[ia] /= scale2;
            }
        }
    }

    scale1 *= dt;
    if (params.record.Rd_at) {
        for (short ia = 0; ia < na; ia++) {
            for (short it = 0; it < nt; it++) {
                double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
                // scale Rd_at.
                if (Mode == 0) {
                    tracer.R.at[ia][it] *= scale2;
                }
                // unscale Rd_at.
                else {
                    tracer.R.at[ia][it] /= scale2;
                }
            }
        }
    }

    if (params.record.Td_at) {
        for (short ia = 0; ia < na; ia++) {
            for (short it = 0; it < nt; it++) {
                double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
                // scale Td_at.
                if (Mode == 0) {
                    tracer.T.at[ia][it] *= scale2;
                }
                // unscale Td_at.
                else {
                    tracer.T.at[ia][it] /= scale2;
                }
            }
        }
    }

    scale1 = 2.0 * std::numbers::pi * dr * dr * std::numbers::pi * da * params.num_photons;
    if (params.record.Rd_ra) {
        for (short ir = 0; ir < nr; ir++) {
            for (short ia = 0; ia < na; ia++) {
                double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                // scale Rd_ra.
                if (Mode == 0) {
                    tracer.R.ra[ir][ia] *= scale2;
                }
                // unscale Rd_ra.
                else {
                    tracer.R.ra[ir][ia] /= scale2;
                }
            }
        }
    }

    if (params.record.Td_ra) {
        for (short ir = 0; ir < nr; ir++) {
            for (short ia = 0; ia < na; ia++) {
                double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                // scale Td_ra.
                if (Mode == 0) {
                    tracer.T.ra[ir][ia] *= scale2;
                }
                // unscale Td_ra.
                else {
                    tracer.T.ra[ir][ia] /= scale2;
                }
            }
        }
    }

    scale1 *= dt;
    if (params.record.Rd_rat) {
        for (short ir = 0; ir < nr; ir++) {
            for (short ia = 0; ia < na; ia++) {
                for (short it = 0; it < nt; it++) {
                    double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                    // scale Rd_rat.
                    if (Mode == 0) {
                        tracer.R.rat[ir][ia][it] *= scale2;
                    }
                    // unscale Rd_rat.
                    else {
                        tracer.R.rat[ir][ia][it] /= scale2;
                    }
                }
            }
        }
    }

    if (params.record.Td_rat) {
        for (short ir = 0; ir < nr; ir++) {
            for (short ia = 0; ia < na; ia++) {
                for (short it = 0; it < nt; it++) {
                    double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                    // scale Td_rat.
                    if (Mode == 0) {
                        tracer.T.rat[ir][ia][it] *= scale2;
                    }
                    // unscale Td_rat.
                    else {
                        tracer.T.rat[ir][ia][it] /= scale2;
                    }
                }
            }
        }
    }
}

/**************************************************************************
 *	Scale absorption arrays properly.
 *  Mode = 0, scale A; Mode = 1, unscale A.
 ****/
void ScaleA(RunParams& params, Tracer& tracer, char Mode)
{
    std::size_t nz = params.num_z;
    std::size_t nr = params.num_r;
    std::size_t nt = params.num_time;
    double dz = params.grid_z;
    double dr = params.grid_r;
    double dt = params.grid_time;
    double scale1 = (double)params.num_photons;

    // scale A.
    if (Mode == 0) {
        tracer.A.ae = 1 / scale1 * sqrt(tracer.A.ae - tracer.A.ab * tracer.A.ab / scale1);
        tracer.A.ab /= scale1;
    }
    // unscale A.
    else {
        tracer.A.ab *= scale1;
        tracer.A.ae = (scale1 * tracer.A.ae) * (scale1 * tracer.A.ae) + 1 / scale1 * tracer.A.ab * tracer.A.ab;
    }

    double scale2 = scale1 * dt;
    if (params.record.A_t) {
        for (short it = 0; it < nt; it++) {
            // scale A_t.
            if (Mode == 0) {
                tracer.A.t[it] /= scale2;
            }
            // unscale A_t.
            else {
                tracer.A.t[it] *= scale2;
            }
        }
    }

    scale1 *= dz;
    if (params.record.A_z) {
        for (short iz = 0; iz < nz; iz++) {
            // scale A_z.
            if (Mode == 0) {
                tracer.A.z[iz] /= scale1;
            }
            // unscale A_z.
            else {
                tracer.A.z[iz] *= scale1;
            }
        }
    }

    scale2 = scale1 * dt;
    if (params.record.A_zt) {
        for (short iz = 0; iz < nz; iz++) {
            for (short it = 0; it < nt; it++) {
                // scale A_zt.
                if (Mode == 0) {
                    tracer.A.zt[iz][it] /= scale2;
                }
                // unscale A_zt.
                else {
                    tracer.A.zt[iz][it] *= scale2;
                }
            }
        }
    }

    if (params.record.A_rz) {
        for (short iz = 0; iz < nz; iz++) {
            // scale Ab_z.
            if (Mode == 0) {
                tracer.A.z[iz] /= scale1;
            }
            // unscale Ab_z.
            else {
                tracer.A.z[iz] *= scale1;
            }
        }
    }

    if (params.record.A_rzt) {
        for (short iz = 0; iz < nz; iz++) {
            for (short it = 0; it < nt; it++) {
                // scale Ab_zt.
                if (Mode == 0) {
                    tracer.A.zt[iz][it] /= scale2;
                }
                // unscale Ab_zt.
                else {
                    tracer.A.zt[iz][it] *= scale2;
                }
            }
        }
    }

    scale1 = 2.0 * std::numbers::pi * dr * dr * dz * params.num_photons;
    if (params.record.A_rz) {
        for (short ir = 0; ir < nr; ir++) {
            for (short iz = 0; iz < nz; iz++) {
                // scale A_rz.
                if (Mode == 0) {
                    tracer.A.rz[ir][iz] /= (ir + 0.5) * scale1;
                }
                // unscale A_rz.
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
                    // scale A_rzt.
                    if (Mode == 0) {
                        tracer.A.rzt[ir][iz][it] /= (ir + 0.5) * scale2;
                    }
                    // unscale A_rzt.
                    else {
                        tracer.A.rzt[ir][iz][it] *= (ir + 0.5) * scale2;
                    }
                }
            }
        }
    }
}

/**************************************************************************
 *	Scale results of current run.
 *  Mode = 0, scale result; Mode = 1, unscale result.
 ****/
void ScaleResult(RunParams& params, Tracer& tracer, char Mode)
{
    ScaleRdTd(params, tracer, Mode);
    ScaleA(params, tracer, Mode);
}

/**************************************************************************
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

/***************************************************************************
 * Save the status of the random number generator to output file.
 ****/
void SaveRandomStatus(std::fstream& file)
{
    //file << std::format("# status of the random number generator:") << std::endl;
    //file << RandomEngine << std::endl;
    //file << std::endl;

    // Get the status
    long status[57];
    RandomGen(2, 0, status);
    file << std::format("# status of the random number generator:") << std::endl;

    for (int i = 0; i < 57; i++) {
        if (i % 5) {
            file << std::format("{:14d}", status[i]);
        }
        else {
            file << std::endl << std::format("{:14d}", status[i]);
        }
    }

    file << std::endl << std::endl;
}

/***************************************************************************
 * Read and restore the status of random number generater from previous
 * output input.
 ****/
void RestoreRandomStatus(std::fstream& file)
{
    //std::string buf;
    //do {
    //    std::getline(file, buf);
    //} while (buf[0] != '#');

    //file >> RandomEngine;

    std::string buf;
    long status[57];

    do {
        std::getline(file, buf);
    } while (buf[0] != '#');

    for (int i = 0; i < 57; i++) {
        file >> status[i];
    }

    // Restore the status
    RandomGen(3, 0, status);
}

/**************************************************************************
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

/**************************************************************************
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

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOAb_zt(std::fstream& file, std::size_t Nz, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Ab[z][t]. [1/(cm ps)]",
                            "# Ab[0][0], [0][1],..[0][nt-1]",
                            "# Ab[1][0], [1][1],..[1][nt-1]",
                            "# ...",
                            "# Ab[nz-1][0], [nz-1][1],..[nz-1][nt-1]",
                            "Ab_zt");
    }
    else {
        // skip A_z line.
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t iz = 0; iz < Nz; iz++) {
        for (std::size_t it = 0; it < Nt; it++) {
            if (Mode == 1) {
                file << std::format("{:12.4E} ", tracer.A.zt[iz][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.A.zt[iz][it];
                //fscanf_s(file, "%lf", &(tracer.A.zt[iz][it]));
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOA_rzt(std::fstream& file, std::size_t Nr, std::size_t Nz, std::size_t Nt, Tracer& tracer, char Mode)
{
    IOAb_zt(file, Nz, Nt, tracer, Mode);

    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# A[r][z][t]. [1/(cm3 ps)]",
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
                if (Mode == 1) {
                    file << std::format("{:12.4E} ", tracer.A.rzt[ir][iz][it]);
                    if (++i % 5 == 0) {
                        file << std::format("\n");
                    }
                }
                else {
                    file >> tracer.A.rzt[ir][iz][it];
                    //fscanf_s(file, "%lf", &(tracer.A.rzt[ir][iz][it]));
                }
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  1 number each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOAb_z(std::fstream& file, std::size_t Nz, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        file << std::format("Ab_z #Ab[0], [1],..Ab[nz-1]. [1/cm]\n");	// flag.
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t iz = 0; iz < Nz; iz++) {
        if (Mode == 1) {
            file << std::format("{:12.4E}\n", tracer.A.z[iz]);
        }
        else {
            file >> tracer.A.z[iz];
            //fscanf_s(file, "%lf", &(tracer.A.z[iz]));
        }
    }

    if (Mode == 1) {
        file << std::format("\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOA_rz(std::fstream& file, std::size_t Nr, std::size_t Nz, Tracer& tracer, char Mode)
{
    IOAb_z(file, Nz, tracer, Mode);

    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n",
                            "# A[r][z]. [1/cm3]",
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
            if (Mode == 1) {
                file << std::format("{:12.4E} ", tracer.A.rz[ir][iz]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.A.rz[ir][iz];
                //fscanf_s(file, "%lf", &(tracer.A.rz[ir][iz]));
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOA_zt(std::fstream& file, std::size_t Nz, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# A[z][t]. [1/(cm ps)]",
                            "# A[0][0], [0][1],..[0][nt-1]",
                            "# A[1][0], [1][1],..[1][nt-1]",
                            "# ...",
                            "# A[nz-1][0], [nz-1][1],..[nz-1][nt-1]",
                            "A_zt");
    }
    else {
        // skip A_zt line.
        FindDataLine(file);
    }

    std::size_t i = 0;
    for (std::size_t iz = 0; iz < Nz; iz++) {
        for (std::size_t it = 0; it < Nt; it++) {
            if (Mode == 1) {
                file << std::format("{:12.4E} ", tracer.A.zt[iz][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.A.zt[iz][it];
                //fscanf_s(file, "%lf", &(tracer.A.zt[iz][it]));
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  1 number each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOA_z(std::fstream& file, std::size_t Nz, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("A_z #A[0], [1],..A[nz-1]. [1/cm]\n");
    }
    else {
        // skip A_z line.
        FindDataLine(file);
    }

    for (std::size_t iz = 0; iz < Nz; iz++) {
        if (Mode == 1) {
            file << std::format("{:12.4E}\n", tracer.A.z[iz]);
        }
        else {
            file >> tracer.A.z[iz];
            //fscanf_s(file, "%lf", &(tracer.A.z[iz]));
        }
    }

    if (Mode == 1) {
        file << std::format("\n");
    }
}

/**************************************************************************
 *  1 number each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOA_t(std::fstream& file, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("A_t #A[0], [1],..A[nt-1]. [1/ps]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t it = 0; it < Nt; it++) {
        if (Mode == 1) {
            file << std::format("{:12.4E}\n", tracer.A.t[it]);
        }
        else {
            file >> tracer.A.t[it];
            //fscanf_s(file, "%lf", &(tracer.A.t[it]));
        }
    }

    if (Mode == 1) {
        file << std::format("\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IORd_rat(std::fstream& file, std::size_t Nr, std::size_t Na, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Rd[r][a][t]. [1/(cm2 sr ps)]",
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
                if (Mode == 1) {
                    file << std::format("{:12.4E} ", tracer.R.rat[ir][ia][it]);
                    if (++i % 5 == 0) {
                        file << std::format("\n");
                    }
                }
                else {
                    file >> tracer.R.rat[ir][ia][it];
                    //fscanf_s(file, "%lf", &(tracer.R.rat[ir][ia][it]));
                }
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IORd_ra(std::fstream& file, std::size_t Nr, std::size_t Na, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Rd[r][angle]. [1/(cm2 sr)].",
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
            if (Mode == 1) {
                file << std::format("{:12.4E} ", tracer.R.ra[ir][ia]);
                if ((ir * Na + ia + 1) % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.R.ra[ir][ia];
                //fscanf_s(file, "%lf", &(tracer.R.ra[ir][ia]));
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IORd_rt(std::fstream& file, std::size_t Nr, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Rd[r][t]. [1/(cm2 ps)]",
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
            if (Mode == 1) {
                file << std::format("{:12.4E} ", tracer.R.rt[ir][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.R.rt[ir][it];
                //fscanf_s(file, "%lf", &(tracer.R.rt[ir][it]));
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IORd_at(std::fstream& file, std::size_t Na, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
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
            if (Mode == 1) {
                file << std::format("{:12.4E} ", tracer.R.at[ia][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.R.at[ia][it];
                //fscanf_s(file, "%lf", &(tracer.R.at[ia][it]));
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  1 number each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IORd_r(std::fstream& file, std::size_t Nr, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("Rd_r #Rd[0], [1],..Rd[nr-1]. [1/cm2]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ir = 0; ir < Nr; ir++) {
        if (Mode == 1) {
            file << std::format("{:12.4E}\n", tracer.R.r[ir]);
        }
        else {
            file >> tracer.R.r[ir];
            //fscanf_s(file, "%lf", &(tracer.R.r[ir]));
        }
    }

    if (Mode == 1) {
        file << std::format("\n");
    }
}

/**************************************************************************
 *  1 number each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IORd_a(std::fstream& file, std::size_t Na, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("Rd_a #Rd[0], [1],..Rd[na-1]. [1/sr]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ia = 0; ia < Na; ia++) {
        if (Mode == 1) {
            file << std::format("{:12.4E}\n", tracer.R.a[ia]);
        }
        else {
            file >> tracer.R.a[ia];
            //fscanf_s(file, "%lf", &(tracer.R.a[ia]));
        }
    }

    if (Mode == 1) {
        file << std::format("\n");
    }
}

/**************************************************************************
 *  1 number each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IORd_t(std::fstream& file, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("Rd_t #Rd[0], [1],..Rd[nt-1]. [1/ps]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t it = 0; it < Nt; it++) {
        if (Mode == 1) {
            file << std::format("{:12.4E}\n", tracer.R.t[it]);
        }
        else {
            file >> tracer.R.t[it];
            //fscanf_s(file, "%lf", &(tracer.R.t[it]));
        }
    }

    if (Mode == 1) {
        file << std::format("\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOTd_rat(std::fstream& file, std::size_t Nr, std::size_t Na, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Td[r][a][t]. [1/(cm2 sr ps)]",
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
                if (Mode == 1) {
                    file << std::format("{:12.4E} ", tracer.T.rat[ir][ia][it]);
                    if (++i % 5 == 0) {
                        file << std::format("\n");
                    }
                }
                else {
                    file >> tracer.T.rat[ir][ia][it];
                    //fscanf_s(file, "%lf", &(tracer.T.rat[ir][ia][it]));
                }
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOTd_ra(std::fstream& file, std::size_t Nr, std::size_t Na, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Td[r][angle]. [1/(cm2 sr)].",
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
            if (Mode == 1) {
                file << std::format("{:12.4E} ", tracer.T.ra[ir][ia]);
                if ((ir * Na + ia + 1) % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.T.ra[ir][ia];
                //fscanf_s(file, "%lf", &(tracer.T.ra[ir][ia]));
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOTd_rt(std::fstream& file, std::size_t Nr, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("{}\n{}\n{}\n{}\n{}\n{}\n",
                            "# Td[r][t]. [1/(cm2 ps)]",
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
            if (Mode == 1) {
                file << std::format("{:12.4E} ", tracer.T.rt[ir][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.T.rt[ir][it];
                //fscanf_s(file, "%lf", &(tracer.T.rt[ir][it]));
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOTd_at(std::fstream& file, std::size_t Na, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
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
            if (Mode == 1) {
                file << std::format("{:12.4E} ", tracer.T.at[ia][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                file >> tracer.T.at[ia][it];
                //fscanf_s(file, "%lf", &(tracer.T.at[ia][it]));
            }
        }
    }

    if (Mode == 1) {
        file << std::format("\n\n");
    }
}

/**************************************************************************
 *  1 number each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOTd_r(std::fstream& file, std::size_t Nr, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("Td_r #Td[0], [1],..Td[nr-1]. [1/cm2]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ir = 0; ir < Nr; ir++) {
        if (Mode == 1) {
            file << std::format("{:12.4E}\n", tracer.T.r[ir]);
        }
        else {
            file >> tracer.T.r[ir];
            //fscanf_s(file, "%lf", &(tracer.T.r[ir]));
        }
    }

    if (Mode == 1) {
        file << std::format("\n");
    }
}

/**************************************************************************
 *  1 number each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOTd_a(std::fstream& file, std::size_t Na, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("Td_a #Td[0], [1],..Td[na-1]. [1/sr]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t ia = 0; ia < Na; ia++) {
        if (Mode == 1) {
            file << std::format("{:12.4E}\n", tracer.T.a[ia]);
        }
        else {
            file >> tracer.T.a[ia];
            //fscanf_s(file, "%lf", &(tracer.T.a[ia]));
        }
    }

    if (Mode == 1) {
        file << std::format("\n");
    }
}

/**************************************************************************
 *  1 number each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOTd_t(std::fstream& file, std::size_t Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        // flag.
        file << std::format("Td_t #Rd[0], [1],..Td[nt-1]. [1/ps]\n");
    }
    else {
        FindDataLine(file);
    }

    for (std::size_t it = 0; it < Nt; it++) {
        if (Mode == 1) {
            file << std::format("{:12.4E}\n", tracer.T.t[it]);
        }
        else {
            file >> tracer.T.t[it];
            //fscanf_s(file, "%lf", &(tracer.T.t[it]));
        }
    }

    if (Mode == 1) {
        file << std::format("\n");
    }
}

/**************************************************************************
 *  Mode = 0, read result back from a output input.
 *  Mode = 1, write result to a output input;
 ****/
void IOResult(std::fstream& file, RunParams& params, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        if (params.output_file_format == FileFormat::Ascii) {
            WriteVersion(file, "mcmloA2.0");
        }
        else {
            WriteVersion(file, "mcmloB2.0");
        }

        PutInputToFile(file, params);
        SaveRandomStatus(file);
        WriteRAT(file, tracer);
    }
    else {
        RestoreRandomStatus(file);
        ReadRAT(file, tracer);
    }

    // reflectance, absorption, transmittance.
    if (params.record.A_rzt) {
        IOA_rzt(file, params.num_r, params.num_z, params.num_time, tracer, Mode);
    }
    if (params.record.A_rz) {
        IOA_rz(file, params.num_r, params.num_z, tracer, Mode);
    }
    if (params.record.A_zt) {
        IOA_zt(file, params.num_z, params.num_time, tracer, Mode);
    }
    if (params.record.A_z) {
        IOA_z(file, params.num_z, tracer, Mode);
    }
    if (params.record.A_t) {
        IOA_t(file, params.num_time, tracer, Mode);
    }

    if (params.record.Rd_rat) {
        IORd_rat(file, params.num_r, params.num_alpha, params.num_time, tracer, Mode);
    }
    if (params.record.Rd_ra) {
        IORd_ra(file, params.num_r, params.num_alpha, tracer, Mode);
    }
    if (params.record.Rd_rt) {
        IORd_rt(file, params.num_r, params.num_time, tracer, Mode);
    }
    if (params.record.Rd_at) {
        IORd_at(file, params.num_alpha, params.num_time, tracer, Mode);
    }
    if (params.record.Rd_r) {
        IORd_r(file, params.num_r, tracer, Mode);
    }
    if (params.record.Rd_a) {
        IORd_a(file, params.num_alpha, tracer, Mode);
    }
    if (params.record.Rd_t) {
        IORd_t(file, params.num_time, tracer, Mode);
    }

    if (params.record.Td_rat) {
        IOTd_rat(file, params.num_r, params.num_alpha, params.num_time, tracer, Mode);
    }
    if (params.record.Td_ra) {
        IOTd_ra(file, params.num_r, params.num_alpha, tracer, Mode);
    }
    if (params.record.Td_rt) {
        IOTd_rt(file, params.num_r, params.num_time, tracer, Mode);
    }
    if (params.record.Td_at) {
        IOTd_at(file, params.num_alpha, params.num_time, tracer, Mode);
    }
    if (params.record.Td_r) {
        IOTd_r(file, params.num_r, tracer, Mode);
    }
    if (params.record.Td_a) {
        IOTd_a(file, params.num_alpha, tracer, Mode);
    }
    if (params.record.Td_t) {
        IOTd_t(file, params.num_time, tracer, Mode);
    }

    file.close();
}

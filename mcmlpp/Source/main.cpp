/*******************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *  Copyright Maurits Lam, 2025.
 *	Monte Carlo simulation of light transport in multi-layered turbid mediums.
 ****/


#include "mcml.hpp"

#include <regex>
#include <chrono>
#include <fstream>
#include <iostream>

#include "random.hpp"


using namespace std::chrono;


void InitTracer(RunParams&, Tracer&);
std::string FindDataLine(std::fstream&);
std::string FindDataLine(std::istream&);
void InterReadParams(RunParams&);
void LaunchPhoton(double, RunParams&, Tracer&, Photon&);
void AboutMCML();
void IOResult(std::fstream&, RunParams&, Tracer&, char);
void CheckParamFromFile(std::fstream&, RunParams&);
void ReadRunParams(std::fstream&, RunParams&);
void ScaleResult(RunParams&, Tracer&, char);
double Rspecular(std::vector<Layer>&);
bool ReadMediumListQ(std::fstream&, RunParams&);
void TracePhoton(RunParams&, Photon&, Tracer&);
bool RunNewInput(RunParams&);
bool GetFile(std::string&, const std::string, std::fstream&);
void ClearRun(RunParams&, Tracer&);
bool ReadEndCriteria(std::istream&, RunParams&, char);
bool CheckFileVersion(std::fstream&, const std::string);

Random Rand;


/*******************************************************************************
 *  If F = 0: reset the clock and return 0.
 *  If F = 1: pass the real time to msg and return the time elapsed since F=0.
 *  If F = 2: same as F = 1 except no printing.
 ****/
long long PunchTime(char F, std::string& msg, RunParams& params)
{
    // Real time reference
    static system_clock::time_point rt0;

    if (F == 0) {
        rt0 = system_clock::now();
        return 0;
    }

    auto now = system_clock::now();
    duration<double> elapsed = now - rt0;

    // Show & pass real time
    auto real_time_secs = duration_cast<seconds>(now - rt0);

    if (F == 1) {
        real_time_secs += seconds(params.add_limit);
        msg = std::format("This took {:%H:%M:%S}.", elapsed);
    }

    return real_time_secs.count();
}

/*******************************************************************************
 *	Print the current time and the estimated finishing time.
 ****/
void PredictDoneTime(long photons_computed, RunParams& params)
{
    auto now = std::chrono::system_clock::now();
    std::cout << std::format("{:%H:%M %a %m/%d/%Y}\t", now);

    if (params.control_bit == ControlBit::TimeLimit || params.control_bit == ControlBit::Both) {
        std::string msg;

        auto done_time = now + (seconds(PunchTime(2, msg, params)) / seconds(photons_computed) * seconds(params.num_photons - photons_computed));

        if (params.control_bit == ControlBit::Both) {
            if (!(done_time < (now + seconds(params.time_limit) - seconds(PunchTime(2, msg, params))))) {
                done_time = (now + seconds(params.time_limit) - seconds(PunchTime(2, msg, params)));
            }
        }

        std::cout << std::format("{:%H:%M %a %Y/%m/%d}", done_time);
    }
    std::cout << std::endl;
}

/*******************************************************************************
 *	Generate a string representing the user-specified done time.
 ****/
std::string FormDateString(RunParams& params)
{
    std::string msg;

    auto now = system_clock::now();
    auto time_limit = seconds(params.time_limit);
    auto punch_time = seconds(PunchTime(2, msg, params));
    
    auto done_time = now + time_limit - punch_time;

    return std::format("{:%H:%M on %Y/%m/%d}", done_time);
}

/*******************************************************************************
 *	Report how and when the simultion will be terminated.
 ****/
void ReportControlInfo(short NumRunsLeft, RunParams& params)
{
    std::string string;

    std::cout << std::endl << "Starting run #" << params.num_runs - NumRunsLeft  << ". ";
    switch (params.control_bit) {
        case ControlBit::NumPhotons: {
            std::cout << "Tracing " << params.num_photons << " photons." << std::endl << std::endl;
            std::cout << "\tPhotons Done\tCurrent Time\t\tEstimated Done Time" << std::endl;
            std::cout << "\t------------\t--------------------\t--------------------" << std::endl;
            break;
        }
        case ControlBit::TimeLimit: {
            std::string date = FormDateString(params);
            std::cout << "The simulation will terminate on " << date << ".\n" << std::endl;
            std::cout << "\tPhotons Done\tCurrent Time" << std::endl;
            std::cout << "\t------------\t--------------------" << std::endl;
            break;
        }
        case ControlBit::Both: {
            std::string date = FormDateString(params);
            std::cout << "Tracing " << params.num_photons << " photons ";
            std::cout << "with a deadline at " << date << ".\n" << std::endl;
            std::cout << "\tPhotons Done\tCurrent Time\t\tEstimated Done Time" << std::endl;
            std::cout << "\t------------\t--------------------\t--------------------" << std::endl;
            break;
        }
    }
}

/*******************************************************************************
 *	Report the estimated time, number of photons and runs left after calculating
 *  10 photons or every 1/10th of total number of photons.
 ****/
void ReportStatus(long num_traced, RunParams& params)
{
    if (params.control_bit == ControlBit::TimeLimit) {
        std::cout << std::format("{:>11} ({:6.2f}%)\t", num_traced, (float)num_traced * 100 / params.num_photons);
    }
    else if (params.control_bit == ControlBit::Both) {
        std::cout << std::format("\t{:>12}\t", num_traced);
    }
    else {
        std::cout << std::format("{:>11} ({:6.2f}%)\t", num_traced, (float)num_traced * 100 / params.num_photons);
    }

    PredictDoneTime(num_traced, params);
}

/*******************************************************************************
 *	Report time, photon number traced and write results.
 ****/
void ReportResult(RunParams& params, Tracer& tracer)
{
    std::string time_report;
    PunchTime(1, time_report, params);
    std::cout << std::endl << "Finished tracing " << params.num_photons << " photons. " << time_report << std::endl;

    ScaleResult(params, tracer, 0);

    std::fstream file(params.output_filename, std::ios::out);

    if (!file.is_open()) {
        std::cout << "Can not open output file to write." << std::endl;
        std::exit(1);
    }
    IOResult(file, params, tracer, 1);
}

/*******************************************************************************
 *	Execute Monte Carlo simulation for one independent run.
 *  Type 0: start a new simulation.
 *  Type 1: continue previous simulation.
 ****/
void DoOneRun(short remaining_runs, RunParams& params, Tracer& tracer, char run_type)
{
    int tens = 10;

    // Start a new simulation
    if (run_type == 0) {
        if (params.source_layer == 0) {
            tracer.R.sp = Rspecular(params.layers);
        }

        // Initialize the number generator
        Rand.seed(1);
    }
    
    std::string msg;
    PunchTime(0, msg, params);
    ReportControlInfo(remaining_runs, params);

    // Photon number traced
    long i_photon = 1;

    // Switch to terminate simulation
    bool exit_switch = 0;

    do {
        Photon photon;
        LaunchPhoton(tracer.R.sp, params, tracer, photon);
        TracePhoton(params, photon, tracer);

        // Report status every ten photons
        if (i_photon == tens) {
            tens *= 10;
            ReportStatus(i_photon, params);
        }
        i_photon++;
        if (params.control_bit == ControlBit::NumPhotons) {
            exit_switch = (i_photon > params.num_photons);
        }
        else if (params.control_bit == ControlBit::TimeLimit) {
            exit_switch = (PunchTime(2, msg, params) >= params.time_limit);
        }
        else {
            exit_switch = (i_photon > params.num_photons) || (PunchTime(2, msg, params) >= params.time_limit);
        }
    }
    while (!exit_switch);

    params.num_photons = params.add_num_photons + i_photon - 1;
    params.time_limit = params.add_limit + (long)PunchTime(2, msg, params);
    params.control_bit = ControlBit::Both;

    ReportResult(params, tracer);
    ClearRun(params, tracer);
}

/*******************************************************************************
 *	In continuation runs, ask for additional number of photons or time.
 ****/
void AddNumPhotons(RunParams& params)
{
    std::cout << std::endl << params.num_photons << " photons have been traced in the previous simulation.";
    std::cout << std::endl << "Specify additional photons or compution time in hh:mm format,";
    std::cout << std::endl << "or both in one line (e.g. 10000 5:30): ";

    while (!ReadEndCriteria(std::cin, params, 1)) {
        std::cout << "Input again: ";
    }

    std::cout << std::endl;
}

/*******************************************************************************
 *	Start a new simulation non-interactively. Read input parameter from file.
 ****/
void Simulate(std::fstream& file, RunParams& params, Tracer& tracer)
{
    CheckParamFromFile(file, params);
    short num_runs_left = params.num_runs;
    while (num_runs_left--) {
        ReadRunParams(file, params);
        InitTracer(params, tracer);
        DoOneRun(num_runs_left, params, tracer, 0);
    }
    file.close();
    std::exit(0);
}

/*******************************************************************************
 *	Read input parameters from a file with interactive change.
 ****/
void FileInterSimu(RunParams& params, Tracer& tracer)
{
    std::string input_filename;
    std::fstream input_file;

    if (GetFile(input_filename, "mcmli2.0", input_file)) {
        if (ReadMediumListQ(input_file, params)) {
            ReadRunParams(input_file, params);
            std::cout << "The parameters of the first run have been read in." << std::endl;

            if (RunNewInput(params)) {
                InitTracer(params, tracer);
                DoOneRun(0, params, tracer, 0);
                input_file.close();
                std::exit(0);
            }
            input_file.close();
        }
    }
}

/*******************************************************************************
 *	Continue a previous simulation.
 ****/
void ContinueSimu(RunParams& params, Tracer& tracer)
{
    std::cout << "Specify the output file name of a previous simulation. " << std::endl;
    
    std::string input_filename;
    std::fstream file;
    
    if (GetFile(input_filename, "mcmloA2.0", file)) {
        return;
    }

    // skip the line of file version.
    FindDataLine(file);
    if (!ReadMediumListQ(file, params)) {
        std::exit(1);
    }
    ReadRunParams(file, params);

    AddNumPhotons(params);
    InitTracer(params, tracer);
    IOResult(file, params, tracer, 0);
    ScaleResult(params, tracer, 1);

    std::swap(params.num_photons, params.add_num_photons);
    std::swap(params.time_limit, params.add_limit);
    DoOneRun(0, params, tracer, 1);
    std::exit(0);
}

/*******************************************************************************
****/
void QuitProgram()
{
    std::cout << "Do you really want to quit MCML? (y/n): ";

    std::string cmd;
    std::getline(std::cin, cmd);

    // Really quit.
    if (std::toupper(cmd[0]) == 'Y') {
        std::exit(0);
    }
}

/*******************************************************************************
 ****/
void ShowMainMenu()
{
    std::cout << "  a = About MCML." << std::endl;
    std::cout << "  r = Run an input file non-interactively." << std::endl;
    std::cout << "  m = Input and modify parameters of a file (the first run only)." << std::endl;
    std::cout << "  i = Input parameters interactively." << std::endl;
    std::cout << "  c = Continue a previous simulation." << std::endl;
    std::cout << "  q = Quit from the program." << std::endl;
    std::cout << "  * Commands here are not case-sensitive." << std::endl;
}

/*******************************************************************************
 ****/
void BranchMainMenu(std::string& string, RunParams& params, Tracer& tracer)
{
    switch (std::toupper(string[0])) {
        case 'A':
        {
            AboutMCML();
            break;
        }

        // Non-interactive
        case 'R':
        {
            std::string input_filename;
            std::fstream input_file;
            
            if (GetFile(input_filename, "mcmli2.0", input_file)) {
                Simulate(input_file, params, tracer);
            }
            break;
        }

        // Read a file with an interactive change
        case 'M':
        {
            FileInterSimu(params, tracer);
            break;
        }

        // Interactive
        case 'I':
        {
            InterReadParams(params);
            if (RunNewInput(params)) {
                InitTracer(params, tracer);
                DoOneRun(0, params, tracer, 0);
                std::exit(0);
            }
            break;
        }

        case 'C':
        {
            ContinueSimu(params, tracer);
            break;
        }

        case 'H':
        {
            ShowMainMenu();
            break;
        }

        case 'Q':
        {
            QuitProgram();
            break;
        }

        default:
        {
            std::cout << "...Unknown command" << std::endl;
        }
    }
}

/*******************************************************************************
 *	The argument to the command line is the input filename, if any.
 ****/
int main(int argc, char* argv[])
{
    std::cout << "MCML Version 3.0, Copyright (c) 1992-1996, 2025" << std::endl << std::endl;

    // Read parameters from input file
    if (argc >= 2) {
        std::string input_filename = (argc >= 2) ? std::string(argv[1]) : std::string();
        std::fstream input_file_ptr(input_filename, std::ios::in);

        if (input_file_ptr.is_open()) {
            if (CheckFileVersion(input_file_ptr, "mcmli2.0")) {
                RunParams params;
                Tracer tracer;
                Simulate(input_file_ptr, params, tracer);
            }
        }
        std::exit(0);
    }
    // Accept commands from console
    else {
        while (1) {
            std::string string;
            do {
                std::cout << std::endl << "> Main menu (h for help) => ";
                std::getline(std::cin, string);
            }
            while (string.empty());

            Tracer tracer;
            RunParams params;
            BranchMainMenu(string, params, tracer);
        }
    }
}

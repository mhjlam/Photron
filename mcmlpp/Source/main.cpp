/*******************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *  Copyright M.H.J. Lam, 2025.
 *	Monte Carlo simulation of light transport in multi-layered turbid mediums.
 ****/


#include "mcml.hpp"

#include <chrono>
#include <format>
#include <fstream>
#include <iostream>

#include "timer.hpp"
#include "random.hpp"


void InitTracer(RunParams&, Tracer&);
std::string NextDataLine(std::istream&);
void InterReadParams(RunParams&);
void LaunchPhoton(RunParams&, Tracer&, Photon&);
void AboutMCML();
void ReadWriteResults(std::fstream&, RunParams&, Tracer&, IoMode);
void CheckParamFromFile(std::fstream&, RunParams&);
void ReadRunParams(std::fstream&, RunParams&);
void ScaleResult(RunParams&, Tracer&, ScaleMode);
bool ReadMedia(std::fstream&, RunParams&);
void TracePhoton(RunParams&, Photon&, Tracer&);
bool RunNewInput(RunParams&);
bool InputFileName(std::string&, const std::string_view&, std::fstream&);
void ClearRun(RunParams&, Tracer&);
bool ReadEndCriteria(std::istream&, RunParams&, bool = false);
bool CheckFileVersion(std::fstream&, const std::string_view&);


// Random number generator
Random g_rand;

// Timer
Timer g_timer;


/*******************************************************************************
 *	Generate a string representing the user-specified done time.
 ****/
std::string FormDateString(RunParams& params)
{
    auto now = std::chrono::system_clock::now();
    auto time_limit = std::chrono::seconds(params.time_limit);
    auto punch_time = std::chrono::seconds(g_timer.punch());
    
    auto done_time = now + time_limit - punch_time;

    return std::format("{:%H:%M on %Y/%m/%d}", done_time);
}

/*******************************************************************************
 *	Report how and when the simultion will be terminated.
 ****/
void ReportControlInfo(RunParams& params, std::size_t runs_remaining)
{
    std::cout << std::endl << "Starting run #" << params.num_runs - runs_remaining  << ". ";
    switch (params.control_bit) {
        case ControlBit::NumPhotons: {
            std::cout << "Tracing " << params.num_photons << " photons." << std::endl << std::endl;
            std::cout << "\tPhotons Done" << std::endl;
            std::cout << "\t------------" << std::endl;
            break;
        }
        case ControlBit::TimeLimit: {
            std::string date = FormDateString(params);
            std::cout << "The simulation will terminate on " << date << "." << std::endl << std::endl;
            std::cout << "\tPhotons Done" << std::endl;
            std::cout << "\t------------" << std::endl;
            break;
        }
        case ControlBit::Both: {
            std::string date = FormDateString(params);
            std::cout << "Tracing " << params.num_photons << " photons ";
            std::cout << "with a deadline at " << date << "." << std::endl << std::endl;
            std::cout << "\tPhotons Done" << std::endl;
            std::cout << "\t------------" << std::endl;
            break;
        }
    }
}

/*******************************************************************************
 *	Report the estimated time, number of photons and runs left after calculating
 *  10 photons or every 1/10th of total number of photons.
 ****/
void ReportStatus(RunParams& params, long photons_done)
{
    if (params.control_bit == ControlBit::TimeLimit) {
        std::cout << std::format("\t{:>12}\t", photons_done) << std::endl;
    }
    else {
        std::cout << std::format("{:>11} ({:6.2f}%)\t", photons_done, (float)photons_done * 100 / params.num_photons) << std::endl;
    }
}

/*******************************************************************************
 *	Report time, photon number traced and write results.
 ****/
void ReportResult(RunParams& params, Tracer& tracer)
{
    std::string time_report = g_timer.hoursMinSec(g_timer.punch());
    std::cout << std::endl << "Finished tracing " << params.num_photons 
              << " photons. This took " << time_report << "." << std::endl;

    ScaleResult(params, tracer, ScaleMode::Scale);

    std::fstream file(params.output_filename, std::ios::out);

    if (!file.is_open()) {
        std::cout << "Can not open output file to write." << std::endl;
        std::exit(1);
    }

    // Write results to file
    ReadWriteResults(file, params, tracer, IoMode::Write);
}

/*******************************************************************************
 *	Execute Monte Carlo simulation for one independent run.
 *  Type 0: start a new simulation.
 *  Type 1: continue previous simulation.
 ****/
void DoOneRun(RunParams& params, Tracer& tracer, RunType run_type, std::size_t runs_remaining = 0)
{
    // Start a new simulation
    if (run_type == RunType::StartNew) {
        if (params.source.layer_index == 0) {
            // Calculate specular reflectance
            double r = (params.layers[0].eta - params.layers[1].eta) / (params.layers[0].eta + params.layers[1].eta);
            tracer.R_spec = (r * r);
        }

        // Initialize the number generator
        g_rand.seed(1);
    }
    
    g_timer.reset();
    ReportControlInfo(params, runs_remaining);

    // Photon number traced
    long photon_index = 1;

    // Switch to terminate simulation
    bool exit_switch = 0;

    int tens = 10;
    do {
        Photon photon;
        LaunchPhoton(params, tracer, photon);
        TracePhoton(params, photon, tracer);

        // Report status every ten photons
        if (photon_index == tens) {
            tens *= 10;
            ReportStatus(params, photon_index);
        }
        photon_index++;
        if (params.control_bit == ControlBit::NumPhotons) {
            exit_switch = (photon_index > params.num_photons);
        }
        else if (params.control_bit == ControlBit::TimeLimit) {
            exit_switch = (g_timer.punch() >= params.time_limit);
        }
        else {
            exit_switch = (photon_index > params.num_photons) || (g_timer.punch() >= params.time_limit);
        }
    }
    while (!exit_switch);

    params.num_photons = params.add_num_photons + photon_index - 1;
    params.time_limit = params.add_time_limit + (long)g_timer.punch();
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
    std::size_t num_runs_left = params.num_runs;

    while (num_runs_left--) {
        ReadRunParams(file, params);
        InitTracer(params, tracer);
        DoOneRun(params, tracer, RunType::StartNew, num_runs_left);
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

    if (InputFileName(input_filename, MCI_VERSION, input_file)) {
        if (ReadMedia(input_file, params)) {
            ReadRunParams(input_file, params);
            std::cout << "The parameters of the first run have been read in." << std::endl;

            if (RunNewInput(params)) {
                InitTracer(params, tracer);
                DoOneRun(params, tracer, RunType::StartNew);
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
void ResumeSimulation(RunParams& params, Tracer& tracer)
{
    std::cout << "Specify the output file name of a previous simulation. " << std::endl;
    
    std::string input_filename;
    std::fstream file;
    
    if (InputFileName(input_filename, MCO_VERSION, file)) {
        return;
    }

    // skip the line of file version.
    NextDataLine(file);
    if (!ReadMedia(file, params)) {
        std::exit(1);
    }
    ReadRunParams(file, params);

    AddNumPhotons(params);
    InitTracer(params, tracer);
    ReadWriteResults(file, params, tracer, IoMode::Read);
    ScaleResult(params, tracer, ScaleMode::Unscale);

    std::swap(params.num_photons, params.add_num_photons);
    std::swap(params.time_limit, params.add_time_limit);
    DoOneRun(params, tracer, RunType::Continue);
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
            
            if (InputFileName(input_filename, MCI_VERSION, input_file)) {
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
                DoOneRun(params, tracer, RunType::StartNew);
                std::exit(0);
            }
            break;
        }

        case 'C':
        {
            ResumeSimulation(params, tracer);
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
            if (CheckFileVersion(input_file_ptr, MCI_VERSION)) {
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

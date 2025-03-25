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

#include "random.hpp"


using str = std::string;

void    InitTracer(RunParams&, Tracer&);
str     FindDataLine(std::istream&);
void    InterReadParams(RunParams&);
void    LaunchPhoton(double, RunParams&, Tracer&, Photon&);
void    AboutMCML();
void    ReadWriteResults(std::fstream&, RunParams&, Tracer&, IoMode);
void    CheckParamFromFile(std::fstream&, RunParams&);
void    ReadRunParams(std::fstream&, RunParams&);
void    ScaleResult(RunParams&, Tracer&, ScaleMode);
double  Rspecular(std::vector<Layer>&);
bool    ReadMediumListQ(std::fstream&, RunParams&);
void    TracePhoton(RunParams&, Photon&, Tracer&);
bool    RunNewInput(RunParams&);
bool    GetFile(std::string&, const std::string&, std::fstream&);
void    ClearRun(RunParams&, Tracer&);
bool    ReadEndCriteria(std::istream&, RunParams&, bool = false);
bool    CheckFileVersion(std::fstream&, const std::string&);

// Random number generator
Random Rand;

// Real time reference
static std::chrono::system_clock::time_point RealTimeRef;


/*******************************************************************************
 *  Mode = ResetTimer:      reset the clock
 *  Mode = TimeElapsed:     return elapsed time since last reset
 *  Mode = TimeElapsedStr:  return elapsed time since last reset and pass 
 *                          the real time to output string
 ****/
long long PunchTime(PunchMode mode, std::string& msg, RunParams& params)
{
    if (mode == PunchMode::ResetTimer) {
        RealTimeRef = std::chrono::system_clock::now();
        return 0;
    }

    auto now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = now - RealTimeRef;

    // Show & pass real time
    auto real_time_secs = std::chrono::duration_cast<std::chrono::seconds>(now - RealTimeRef);

    
    if (mode == PunchMode::TimeElapsedStr) {
        real_time_secs += std::chrono::seconds(params.add_limit);
        msg = std::format("This took {:%H:%M:%S}.", elapsed);
    }

    // Return time elapsed since last reset
    return real_time_secs.count();
}

/*******************************************************************************
 *	Print the current time and the estimated finishing time.
 ****/
void PredictDoneTime(long photonsComputed, RunParams& params)
{
    auto now = std::chrono::system_clock::now();
    std::cout << std::format("{:%H:%M %a %m/%d/%Y}\t", now);

    if (params.control_bit == ControlBit::TimeLimit || params.control_bit == ControlBit::Both) {
        std::string msg;

        auto done_time = now + (std::chrono::seconds(PunchTime(PunchMode::TimeElapsed, msg, params)) / std::chrono::seconds(photonsComputed) * std::chrono::seconds(params.num_photons - photonsComputed));

        if (params.control_bit == ControlBit::Both) {
            if (!(done_time < (now + std::chrono::seconds(params.time_limit) - std::chrono::seconds(PunchTime(PunchMode::TimeElapsed, msg, params))))) {
                done_time = (now + std::chrono::seconds(params.time_limit) - std::chrono::seconds(PunchTime(PunchMode::TimeElapsed, msg, params)));
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

    auto now = std::chrono::system_clock::now();
    auto time_limit = std::chrono::seconds(params.time_limit);
    auto punch_time = std::chrono::seconds(PunchTime(PunchMode::TimeElapsed, msg, params));
    
    auto done_time = now + time_limit - punch_time;

    return std::format("{:%H:%M on %Y/%m/%d}", done_time);
}

/*******************************************************************************
 *	Report how and when the simultion will be terminated.
 ****/
void ReportControlInfo(short NumRunsLeft, RunParams& params)
{
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
void ReportStatus(long numTraced, RunParams& params)
{
    if (params.control_bit == ControlBit::TimeLimit) {
        std::cout << std::format("{:>11} ({:6.2f}%)\t", numTraced, (float)numTraced * 100 / params.num_photons);
    }
    else if (params.control_bit == ControlBit::Both) {
        std::cout << std::format("\t{:>12}\t", numTraced);
    }
    else {
        std::cout << std::format("{:>11} ({:6.2f}%)\t", numTraced, (float)numTraced * 100 / params.num_photons);
    }

    PredictDoneTime(numTraced, params);
}

/*******************************************************************************
 *	Report time, photon number traced and write results.
 ****/
void ReportResult(RunParams& params, Tracer& tracer)
{
    std::string time_report;
    PunchTime(PunchMode::TimeElapsedStr, time_report, params);
    std::cout << std::endl << "Finished tracing " << params.num_photons << " photons. " << time_report << std::endl;

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
void DoOneRun(short runsRemaining, RunParams& params, Tracer& tracer, RunType runType)
{
    // Start a new simulation
    if (runType == RunType::StartNew) {
        if (params.source_layer == 0) {
            tracer.R.sp = Rspecular(params.layers);
        }

        // Initialize the number generator
        Rand.seed(1);
    }
    
    std::string msg;
    PunchTime(PunchMode::ResetTimer, msg, params);
    ReportControlInfo(runsRemaining, params);

    // Photon number traced
    long photon_index = 1;

    // Switch to terminate simulation
    bool exit_switch = 0;

    int tens = 10;
    do {
        Photon photon;
        LaunchPhoton(tracer.R.sp, params, tracer, photon);
        TracePhoton(params, photon, tracer);

        // Report status every ten photons
        if (photon_index == tens) {
            tens *= 10;
            ReportStatus(photon_index, params);
        }
        photon_index++;
        if (params.control_bit == ControlBit::NumPhotons) {
            exit_switch = (photon_index > params.num_photons);
        }
        else if (params.control_bit == ControlBit::TimeLimit) {
            exit_switch = (PunchTime(PunchMode::TimeElapsed, msg, params) >= params.time_limit);
        }
        else {
            exit_switch = (photon_index > params.num_photons) || (PunchTime(PunchMode::TimeElapsed, msg, params) >= params.time_limit);
        }
    }
    while (!exit_switch);

    params.num_photons = params.add_num_photons + photon_index - 1;
    params.time_limit = params.add_limit + (long)PunchTime(PunchMode::TimeElapsed, msg, params);
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
        DoOneRun(num_runs_left, params, tracer, RunType::StartNew);
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
                DoOneRun(0, params, tracer, RunType::StartNew);
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
    ReadWriteResults(file, params, tracer, IoMode::Read);
    ScaleResult(params, tracer, ScaleMode::Unscale);

    std::swap(params.num_photons, params.add_num_photons);
    std::swap(params.time_limit, params.add_limit);
    DoOneRun(0, params, tracer, RunType::Continue);
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
                DoOneRun(0, params, tracer, RunType::StartNew);
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

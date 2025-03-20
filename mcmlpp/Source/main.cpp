/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
 *	The main program for Monte Carlo simulation of light transport
 *	in multi-layered turbid mediums.
 ****/


#include "mcml.hpp"

#include <regex>
#include <chrono>
#include <fstream>
#include <iostream>


using namespace std::chrono;


// Declare before they are used in main().
void InitOutputData(RunParams&, Tracer&);
std::string FindDataLine(std::fstream&);
std::string FindDataLine(std::istream&);
void InterReadParam(RunParams&);
void LaunchPhoton(double, RunParams&, Tracer&, Photon&);
void AboutMCML();
void IOResult(std::fstream&, RunParams&, Tracer&, char);
void CheckParamFromFile(std::fstream&, RunParams&);
void ReadRunParam(std::fstream&, RunParams&);
void ScaleResult(RunParams&, Tracer&, char);
double Rspecular(std::vector<Layer>&);
bool ReadMediumListQ(std::fstream&, RunParams&);
void TracePhoton(RunParams&, Photon&, Tracer&);
bool RunChangedInput(RunParams&);
bool GetFile(std::string&, const std::string, std::fstream&);
bool ReadNumPhotonsQ(std::istream&, RunParams&, char);
bool CheckFileVersionQ(std::fstream&, const std::string);


std::mt19937 RandomEngine;
std::uniform_real_distribution<double> Distribution;


/**************************************************************************
 *  If F = 0, reset the clock and return 0.
 *
 *  If F = 1, pass the real time to Msg, print Msg on the
 *  screen, return the real time elapsed since F=0.
 *
 *  If F = 2, same as F=1 except no printing.
 *
 *  Note that clock() and time() return user time and real
 *  time respectively.  User time is whatever the system allocates
 *  to the running of the program; real time is wall-clock time.
 *  In a time-shared system, they need not be the same.
 *
 *  clock() only hold 16 bit integer, which is about 32768
 *  clock ticks.  Because this fact can lead to wrap-around, it
 *  is dangerous to use clock() in Monte Carlo simulations.
 *  Therefore, we only keep track of the real time.
 *
 *  On UNIX machines, users can use time command to get the
 *  user time.  On personal computers, the user time is equal
 *  to the real time unless the system is multi-tasking during
 *  the Monte Carlo simulation, which usually is not true.
 ****/
long long PunchTime(char F, std::string& msg, RunParams& params)
{
    // real time reference.
    static system_clock::time_point rt0;

    if (F == 0) {
        rt0 = system_clock::now();
        return 0;
    }

    auto now = system_clock::now();
    duration<double> elapsed = now - rt0;

    // show & pass real time.
    auto real_time_secs = duration_cast<seconds>(now - rt0);

    if (F == 1) {
        real_time_secs += seconds(params.add_limit_seconds);
        msg = std::format("This took {:%H:%M:%S}.", elapsed);
    }

    return real_time_secs.count();
}

/**************************************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 ****/
void PredictDoneTime(long P1, RunParams& params)
{
    auto now = std::chrono::system_clock::now();
    std::cout << std::format("{:%H:%M %a %m/%d/%Y}\t", now);

    if (params.control_bit == ControlBit::TimeLimit || params.control_bit == ControlBit::Both) {
        std::string msg;

        auto done_time = now + (seconds(PunchTime(2, msg, params)) / seconds(P1) * seconds(params.num_photons - P1));

        if (params.control_bit == ControlBit::Both) {
            if (!(done_time < (now + seconds(params.time_limit_seconds) - seconds(PunchTime(2, msg, params))))) {
                done_time = (now + seconds(params.time_limit_seconds) - seconds(PunchTime(2, msg, params)));
            }
        }

        std::cout << std::format("{:%H:%M %a %Y/%m/%d}", done_time);
    }
    std::cout << std::endl;
}

/**************************************************************************
 *	Generate a string representing the user-specified done time.
 ****/
std::string FormDateString(RunParams& params)
{
    std::string msg;

    auto now = system_clock::now();
    auto time_limit = seconds(params.time_limit_seconds);
    auto punch_time = seconds(PunchTime(2, msg, params));
    
    auto done_time = now + time_limit - punch_time;

    return std::format("{:%H:%M on %Y/%m/%d}", done_time);
}

/**************************************************************************
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

/**************************************************************************
 *	Report the estimated time, number of photons and runs left
 *	after calculating 10 photons or every 1/10 of total
 *	number of photons.
 *
 *	photons_traced is the number of traced photons so far.
 ****/
void ReportStatus(long photons_traced, RunParams& params)
{
    if (params.control_bit == ControlBit::TimeLimit) {
        std::cout << std::format("%11ld(%6.2f%%)\t", photons_traced, (float)photons_traced * 100 / params.num_photons);
    }
    else if (params.control_bit == ControlBit::Both) {
        std::cout << std::format("\t%12ld\t", photons_traced);
    }
    else {
        std::cout << std::format("%11ld(%6.2f%%)\t", photons_traced, (float)photons_traced * 100 / params.num_photons);
    }

    PredictDoneTime(photons_traced, params);
}

/**************************************************************************
 *	Report time, photon number traced and write results.
 ****/
void ReportResult(RunParams& params, Tracer& tracer)
{
    std::string time_report;
    PunchTime(1, time_report, params);
    std::cout << std::endl << "Finished tracing " << params.num_photons << " photons. " << time_report << std::endl;

    ScaleResult(params, tracer, 0);

    std::fstream file(params.output_filename, std::ios::in | std::ios::out);

    if (!file.is_open()) {
        std::cout << "Can not open output file to write." << std::endl;
        exit(1);
    }
    IOResult(file, params, tracer, 1);
}

/**************************************************************************
 *	Get the file name of the input data file from the
 *	argument to the command line.
 ****/
std::string GetFnameFromArgv(int argc, char* argv[])
{
    return (argc >= 2) ? std::string(argv[1]) : std::string();
}

/**************************************************************************
 *	Execute Monte Carlo simulation for one independent run.
 *  Type = 0: start a new simulation;
 *  Type = 1: continue previous simulation.
 ****/
void DoOneRun(short NumRunsLeft, RunParams& params, Tracer& tracer, char Type)
{
    int tens = 10;

    // start a new simulation.
    if (Type == 0) {
        if (params.source_layer == 0) {
            tracer.R.sp = Rspecular(params.layers);
        }
    }
    
    std::string msg;
    PunchTime(0, msg, params);
    ReportControlInfo(NumRunsLeft, params);

    // photon number traced. 
    long i_photon = 1;

    // switch to terminate simulation.
    bool exit_switch = 0;

    do {
        Photon photon;
        LaunchPhoton(tracer.R.sp, params, tracer, photon);
        TracePhoton(params, photon, tracer);

        // report status every ten photons.
        if (i_photon == tens) {
            tens *= 10;
            ReportStatus(i_photon, params);
        }
        i_photon++;
        if (params.control_bit == ControlBit::TimeLimit) {
            exit_switch = (i_photon > params.num_photons);
        }
        else if (params.control_bit == ControlBit::Both) {
            exit_switch = (PunchTime(2, msg, params) >= params.time_limit_seconds);
        }
        else {
            exit_switch = (i_photon > params.num_photons) || (PunchTime(2, msg, params) >= params.time_limit_seconds);
        }
    }
    while (!exit_switch);

    params.num_photons = params.add_num_photons + i_photon - 1;
    params.time_limit_seconds = params.add_limit_seconds + (long)PunchTime(2, msg, params);
    params.control_bit = ControlBit::Both;

    ReportResult(params, tracer);
}

/**************************************************************************
 *	In continuation runs, ask for additional number of photons or time.
 ****/
void add_num_photons(RunParams& params)
{
    std::cout << std::endl << params.num_photons << (" photons have been traced in the previous simulation.");
    std::cout << std::endl << "Specify additional photons or compution time in hh:mm format,";
    std::cout << std::endl << "or both in one line (e.g. 10000 5:30): ";

    while (!ReadNumPhotonsQ(std::cin, params, 1)) {
        std::cout << "Input again: ";
    }

    std::cout << std::endl;
}

/**************************************************************************
 *	Start a new simulation non-interactively.
 *  Read input parameter from file input_filename.
 ****/
void Simulate(std::fstream& file, RunParams& params, Tracer& tracer)
{
    CheckParamFromFile(file, params);
    short num_runs_left = params.num_runs;
    while (num_runs_left--) {
        ReadRunParam(file, params);
        InitOutputData(params, tracer);
        DoOneRun(num_runs_left, params, tracer, 0);
    }
    file.close();
    exit(0);
}

/**************************************************************************
 *	Read input parameters from a file with interactive change.
 ****/
void FileInterSimu(RunParams& params, Tracer& tracer)
{
    std::string input_filename;
    std::fstream input_file;

    if (GetFile(input_filename, "mcmli2.0", input_file)) {
        if (ReadMediumListQ(input_file, params)) {
            ReadRunParam(input_file, params);
            std::cout << "The parameters of the first run have been read in." << std::endl;

            if (RunChangedInput(params)) {
                InitOutputData(params, tracer);
                DoOneRun(0, params, tracer, 0);
                input_file.close();
                exit(0);
            }
            input_file.close();
        }
    }
}

/**************************************************************************
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
        exit(1);
    }
    ReadRunParam(file, params);

    add_num_photons(params);
    InitOutputData(params, tracer);
    IOResult(file, params, tracer, 0);
    ScaleResult(params, tracer, 1);

    std::swap(params.num_photons, params.add_num_photons);
    std::swap(params.time_limit_seconds, params.add_limit_seconds);
    DoOneRun(0, params, tracer, 1);
    exit(0);
}

/**************************************************************************
****/
void QuitProgram(void)
{
    std::cout << "Do you really want to quit MCML? (y/n): ";

    std::string cmd;
    std::getline(std::cin, cmd);

    // Really quit.
    if (toupper(cmd[0]) == 'Y') {
        exit(0);
    }
}

/**************************************************************************
 ****/
void ShowMainMenu(void)
{
    puts("  a = About MCML.");
    puts("  r = Run an input file non-interactively.");
    puts("  m = Input and modify parameters of a file (the first run only).");
    puts("  i = Input parameters interactively.");
    puts("  c = Continue a previous simulation.");
    puts("  q = Quit from the program.");
    puts("  * Commands here are not case-sensitive.");
}

/**************************************************************************
 ****/
void BranchMainMenu(std::string& string, RunParams& params, Tracer& tracer)
{
    switch (toupper(string[0])) {
        case 'A':
        {
            AboutMCML();
            break;
        }

        // non-interactive.
        case 'R':
        {
            std::string input_filename;
            std::fstream input_file;
            
            if (GetFile(input_filename, "mcmli2.0", input_file)) {
                Simulate(input_file, params, tracer);
            }
            break;
        }

        // read a file with an interactive change.
        case 'M':
        {
            FileInterSimu(params, tracer);
            break;
        }

        // interactive.
        case 'I':
        {
            InterReadParam(params);
            if (RunChangedInput(params)) {
                InitOutputData(params, tracer);
                DoOneRun(0, params, tracer, 0);
                exit(0);
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
            puts("...Unknown command");
        }
    }
}

/**************************************************************************
 *	The argument to the command line is the input filename, if any.
 ****/
int main(int argc, char* argv[])
{
    std::cout << "MCML Version 2.0, Copyright (c) 1992-1996\n" << std::endl;

    RandomEngine.seed(std::random_device{}());
    Distribution = std::uniform_real_distribution<double>(0.0, 1.0);

    // non-interactive.
    if (argc >= 2) {
        std::string input_filename = (argc >= 2) ? std::string(argv[1]) : std::string();
        std::fstream input_file_ptr(input_filename, std::ios::in);

        if (input_file_ptr.is_open()) {
            if (CheckFileVersionQ(input_file_ptr, "mcmli2.0")) {
                RunParams params;
                Tracer tracer;
                Simulate(input_file_ptr, params, tracer);
            }
        }
        exit(0);
    }
    // accept commands from console.
    else {
        while (1) {
            std::string str;
            do {
                std::cout << "\n> Main menu (h for help) => ";
                std::getline(std::cin, str);
            }
            while (str.empty());

            Tracer tracer;
            RunParams params;
            BranchMainMenu(str, params, tracer);
        }
    }
}

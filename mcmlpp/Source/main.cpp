/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
 *	The main program for Monte Carlo simulation of light transport
 *	in multi-layered turbid media.
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
long long PunchTime(char F, std::string& msg, RunParams& run_params)
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
        real_time_secs += seconds(run_params.add_limit_seconds);
        std::string formatted_time = std::format("This took {:%H:%M:%S}.", elapsed);
    }

    return real_time_secs.count();
}

/**************************************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 ****/
void PredictDoneTime(long P1, RunParams& run_params)
{
    auto now = std::chrono::system_clock::now();
    std::cout << std::format("{:%H:%M %a %m/%d/%Y}\t", now);

    if (run_params.control_bit == ControlBit::TimeLimit || run_params.control_bit == ControlBit::Both) {
        std::string msg;

        auto done_time = now + (seconds(PunchTime(2, msg, run_params)) / seconds(P1) * seconds(run_params.num_photons - P1));

        if (run_params.control_bit == ControlBit::Both) {
            if (!(done_time < (now + seconds(run_params.time_limit_seconds) - seconds(PunchTime(2, msg, run_params))))) {
                done_time = (now + seconds(run_params.time_limit_seconds) - seconds(PunchTime(2, msg, run_params)));
            }
        }

        std::cout << std::format("{:%H:%M %a %Y/%m/%d}", done_time);
    }
    std::cout << std::endl;
}

/**************************************************************************
 *	Generate a string representing the user-specified done time.
 ****/
std::string FormDateString(RunParams& run_params)
{
    std::string msg;

    auto now = system_clock::now();
    auto time_limit = seconds(run_params.time_limit_seconds);
    auto punch_time = seconds(PunchTime(2, msg, run_params));
    
    auto done_time = now + time_limit - punch_time;

    return std::format("{:%H:%M on %Y/%m/%d}", done_time);
}

/**************************************************************************
 *	Report how and when the simultion will be terminated.
 ****/
void ReportControlInfo(short NumRunsLeft, RunParams& run_params)
{
    std::string string;

    printf("\nStarting run #%d. ", run_params.num_runs - NumRunsLeft);
    switch (run_params.control_bit) {
        case ControlBit::NumPhotons: {
            printf("Tracing %ld photons.\n\n", run_params.num_photons);
            std::cout << "\tPhotons Done\tCurrent Time\t\tEstimated Done Time\n";
            std::cout << "\t------------\t--------------------\t--------------------\n";
            break;
        }
        case ControlBit::TimeLimit: {
            std::string date = FormDateString(run_params);
            std::cout << "The simulation will terminate on " << date << ".\n\n";
            std::cout << "\tPhotons Done\tCurrent Time\n";
            std::cout << "\t------------\t--------------------\n";
            break;
        }
        case ControlBit::Both: {
            std::string date = FormDateString(run_params);
            std::cout << "Tracing " << run_params.num_photons << " photons ";
            std::cout << "with a deadline at " << date << ".\n\n";
            std::cout << "\tPhotons Done\tCurrent Time\t\tEstimated Done Time\n";
            std::cout << "\t------------\t--------------------\t--------------------\n";
            break;
        }
    }
}

/**************************************************************************
 *	Report the estimated time, number of photons and runs left
 *	after calculating 10 photons or every 1/10 of total
 *	number of photons.
 *
 *	Pi is the number of traced photons so far.
 ****/
void ReportStatus(long Pi, RunParams& run_params)
{
    if (run_params.control_bit == ControlBit::TimeLimit) {
        printf("%11ld(%6.2f%%)\t", Pi, (float)Pi * 100 / run_params.num_photons);
    }
    else if (run_params.control_bit == ControlBit::Both) {
        printf("\t%12ld\t", Pi);
    }
    else {
        printf("%11ld(%6.2f%%)\t", Pi, (float)Pi * 100 / run_params.num_photons);
    }

    PredictDoneTime(Pi, run_params);
}

/**************************************************************************
 *	Report time, photon number traced and write results.
 ****/
void ReportResult(RunParams& run_params, Tracer& tracer)
{
    std::string time_report;
    PunchTime(1, time_report, run_params);
    printf("\nFinished tracing %ld photons. %s\n", run_params.num_photons, time_report.c_str());

    ScaleResult(run_params, tracer, 0);

    std::fstream file(run_params.output_filename, std::ios::in | std::ios::out);

    if (!file.is_open()) {
        std::cout << "Can not open output file to write.\n";
        exit(1);
    }
    IOResult(file, run_params, tracer, 1);
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
void DoOneRun(short NumRunsLeft, RunParams& run_params, Tracer& tracer, char Type)
{
    int tens = 10;

    // start a new simulation.
    if (Type == 0) {
        if (run_params.slayer == 0) {
            tracer.R.sp = Rspecular(run_params.layers);
        }
    }
    
    std::string msg;
    PunchTime(0, msg, run_params);
    ReportControlInfo(NumRunsLeft, run_params);

    // photon number traced. 
    long i_photon = 1;

    // switch to terminate simulation.
    bool exit_switch = 0;

    do {
        Photon photon;
        LaunchPhoton(tracer.R.sp, run_params, tracer, photon);
        TracePhoton(run_params, photon, tracer);

        // report status every ten photons.
        if (i_photon == tens) {
            tens *= 10;
            ReportStatus(i_photon, run_params);
        }
        i_photon++;
        if (run_params.control_bit == ControlBit::TimeLimit) {
            exit_switch = (i_photon > run_params.num_photons);
        }
        else if (run_params.control_bit == ControlBit::Both) {
            exit_switch = (PunchTime(2, msg, run_params) >= run_params.time_limit_seconds);
        }
        else {
            exit_switch = (i_photon > run_params.num_photons) || (PunchTime(2, msg, run_params) >= run_params.time_limit_seconds);
        }
    }
    while (!exit_switch);

    run_params.num_photons = run_params.add_num_photons + i_photon - 1;
    run_params.time_limit_seconds = run_params.add_limit_seconds + (long)PunchTime(2, msg, run_params);
    run_params.control_bit = ControlBit::Both;

    ReportResult(run_params, tracer);
}

/**************************************************************************
 *	In continuation runs, ask for additional number of photons or time.
 ****/
void add_num_photons(RunParams& run_params)
{
    printf("\n%ld photons have been traced in the previous simulation.", run_params.num_photons);
    std::cout << "\nSpecify additional photons or compution time in hh:mm format,";
    std::cout << "\nor both in one line (e.g. 10000 5:30): ";

    while (!ReadNumPhotonsQ(std::cin, run_params, 1)) {
        std::cout << "Input again: ";
    }

    std::cout << "\n";
}

/**************************************************************************
 *	Start a new simulation non-interactively.
 *  Read input parameter from file input_filename.
 ****/
void NonInterSimu(std::fstream& file, RunParams& run_params, Tracer& tracer)
{
    CheckParamFromFile(file, run_params);
    short num_runs_left = run_params.num_runs;
    while (num_runs_left--) {
        ReadRunParam(file, run_params);
        InitOutputData(run_params, tracer);
        DoOneRun(num_runs_left, run_params, tracer, 0);
    }
    file.close();
    exit(0);
}

/**************************************************************************
 *	Read input parameters from a file with interactive change.
 ****/
void FileInterSimu(RunParams& run_params, Tracer& tracer)
{
    std::string input_filename;
    std::fstream input_file;

    if (GetFile(input_filename, "mcmli2.0", input_file)) {
        if (ReadMediumListQ(input_file, run_params)) {
            ReadRunParam(input_file, run_params);
            std::cout << "The parameters of the first run have been read in.\n";

            if (RunChangedInput(run_params)) {
                InitOutputData(run_params, tracer);
                DoOneRun(0, run_params, tracer, 0);
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
void ContinueSimu(RunParams& run_params, Tracer& tracer)
{
    std::cout << "Specify the output file name of a previous simulation. \n";
    
    std::string input_filename;
    std::fstream file;
    
    if (GetFile(input_filename, "mcmloA2.0", file)) {
        return;
    }

    // skip the line of file version.
    FindDataLine(file);
    if (!ReadMediumListQ(file, run_params)) {
        exit(1);
    }
    ReadRunParam(file, run_params);

    add_num_photons(run_params);
    InitOutputData(run_params, tracer);
    IOResult(file, run_params, tracer, 0);
    ScaleResult(run_params, tracer, 1);

    std::swap(run_params.num_photons, run_params.add_num_photons);
    std::swap(run_params.time_limit_seconds, run_params.add_limit_seconds);
    DoOneRun(0, run_params, tracer, 1);
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
void BranchMainMenu(std::string& string, RunParams& run_params, Tracer& tracer)
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
                NonInterSimu(input_file, run_params, tracer);
            }
            break;
        }

        // read a file with an interactive change.
        case 'M':
        {
            FileInterSimu(run_params, tracer);
            break;
        }

        // interactive.
        case 'I':
        {
            InterReadParam(run_params);
            if (RunChangedInput(run_params)) {
                InitOutputData(run_params, tracer);
                DoOneRun(0, run_params, tracer, 0);
                exit(0);
            }
            break;
        }

        case 'C':
        {
            ContinueSimu(run_params, tracer);
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
    std::cout << "MCML Version 2.0, Copyright (c) 1992-1996\n\n";

    // non-interactive.
    if (argc >= 2) {
        std::string input_filename = (argc >= 2) ? std::string(argv[1]) : std::string();
        std::fstream input_file_ptr(input_filename, std::ios::in);

        if (input_file_ptr.is_open()) {
            if (CheckFileVersionQ(input_file_ptr, "mcmli2.0")) {
                RunParams in_param;
                Tracer out_param;
                NonInterSimu(input_file_ptr, in_param, out_param);
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
            RunParams run_params;
            BranchMainMenu(str, run_params, tracer);
        }
    }
}

/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
 *	Input/output of data.
 ****/


#include "mcml.hpp"

#include <format>
#include <fstream>
#include <iostream>


 /**************************************************************************
  *	Structure used to check against duplicated input names.
  ****/
struct NameList
{
    std::string name;
    struct NameList* next;
};

typedef struct NameList NameNode;
typedef NameNode* NameLink;


/**************************************************************************
 *	Allocate an array with index from nl to nh inclusive.
 *
 *	Original matrix and vector from Numerical Recipes in C
 *	Don't initialize the elements to zero. This will be accomplished by the
 *  following functions.
 ****/
double* AllocArray1D(short nl, short nh)
{
    double* v = (double*)malloc((unsigned)(nh - nl + 1) * sizeof(double));
    if (!v) {
        fprintf(stderr, "%s\n", "allocation failure in AllocArray1D()");
        fprintf(stderr, "...now exiting to system...\n");
        exit(1);
    }

    v = v - nl;

    /* init. */
    for (short i = nl; i <= nh; i++) {
        v[i] = 0.0;
    }

    return v;
}

/**************************************************************************
 *	Allocate a matrix with row index from nrl to nrh inclusive, and column
 *  index from ncl to nch inclusive.
 ****/
double** AllocArray2D(short nrl, short nrh, short ncl, short nch)
{
    double** m = (double**)malloc((unsigned)(nrh - nrl + 1) * sizeof(double*));
    if (!m) {
        fprintf(stderr, "%s\n", "allocation failure 1 in AllocArray2D()");
        fprintf(stderr, "...now exiting to system...\n");
        exit(1);
    }
    m -= nrl;

    for (short i = nrl; i <= nrh; i++) {
        m[i] = (double*)malloc((unsigned)(nch - ncl + 1) * sizeof(double));
        if (!m[i]) {
            fprintf(stderr, "%s\n", "allocation failure 2 in AllocArray2D()");
            fprintf(stderr, "...now exiting to system...\n");
            exit(1);
        }
        m[i] -= ncl;
    }

    for (short i = nrl; i <= nrh; i++) {
        for (short j = ncl; j <= nch; j++) {
            m[i][j] = 0.0;
        }
    }
    return m;
}

/**************************************************************************
 *	Allocate a 3D array with row index from nrl to nrh inclusive, column
 *  index from ncl to nch inclusive, and depth index from ndl to ndh inclusive.
 ****/
double*** AllocArray3D(short nrl, short nrh, short ncl, short nch, short ndl, short ndh)
{
    double*** m = (double***)malloc((unsigned)(nrh - nrl + 1) * sizeof(double**));
    if (!m) {
        fprintf(stderr, "%s\n", "allocation failure 1 in AllocArray3D()");
        fprintf(stderr, "...now exiting to system...\n");
        exit(1);
    }
    m -= nrl;

    for (short i = nrl; i <= nrh; i++) {
        m[i] = AllocArray2D(ncl, nch, ndl, ndh);
    }

    return m;
}

/**************************************************************************
 *	Release the memory.
 ****/
void FreeArray1D(double* v, short nl, short nh)
{
    if (v != NULL) {
        free((char*)(v + nl));
    }
}

/**************************************************************************
 *	Release the memory.
 ****/
void FreeArray2D(double** m, short nrl, short nrh, short ncl, short nch)
{
    if (m != NULL) {
        for (short i = nrh; i >= nrl; i--) {
            free((char*)(m[i] + ncl));
        }
        free((char*)(m + nrl));
    }
}

/**************************************************************************
 *  Release the memory.
 ****/
void FreeArray3D(double*** m, short nrl, short nrh, short ncl, short nch, short ndl, short ndh)
{
    if (m != NULL) {
        for (short i = nrh; i >= nrl; i--) {
            FreeArray2D(m[i], ncl, nch, ndl, ndh);
        }
        free((char*)(m + nrl));
    }
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
    std::cout << "\tL.-H. Wang, S. L. Jacques, and L.-Q. Zheng, MCML - Monte \n";
    std::cout << "\tCarlo modeling of photon transport in multi-layered\n";
    std::cout << "\ttissues, Computer Methods and Programs in Biomedicine, 47,\n";
    std::cout << "\t131-146 (1995)\n";
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
    bool found = 0;	/* found bad char. */
    size_t sl = Str.size();
    size_t i = 0;

    while (i < sl) {
        if (Str[i] < 0 || Str[i] > 255) {
            std::cout << "Non-ASCII file\n";
            return (0);
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
bool CommentLineQ(std::string buf)
{
    /* length spanned by space or tab chars. */
    size_t spn = buf.find_first_not_of('\t');

    /* length before the 1st # or return. */
    size_t cspn = buf.find_first_not_of("#\n");

    /* comment line or space line. */
    return spn == cspn;
}

/**************************************************************************
 *	Skip space or comment lines and return a data line.
 ****/
std::string FindDataLine(std::fstream& file)
{
    static std::string buf;

    /* skip space or comment lines. */
    do {
        std::string line;
        std::getline(file, line);

        if (line.empty()) {
            std::cout << "Incomplete data." << std::endl;
            buf[0] = '\0';
            break;
        }
        else {
            CheckCharQ(buf);
        }
    } while (CommentLineQ(buf));

    return (buf);
}

/**************************************************************************
 *	Check whether the input version is the same as version.
 ****/
bool CheckFileVersionQ(std::fstream& file, const std::string Version)
{
    /* line buffer. */
    std::string line;

    /* skip comment line. */
    do {
        std::getline(file, line);
        if (line.empty()) {
            break;
        }
    } while (CommentLineQ(line));

    if (line.empty() || line.find(Version) == std::string::npos) {
        std::cout << "Wrong file version.";
        return false;
    }
    return true;
}

/**************************************************************************
 *  Get a filename and open it for reading, retry until the input can be
 *  opened with a correct version or a '.' is typed.
 *	Return a NULL pointer if '.' is typed.
 ****/
std::fstream GetFile(std::string& Fname, const std::string Version)
{
    std::fstream file;

    while (1) {
        /* prompt. */
        std::cout << "Specify filename (or . to quit to main menu):";

        // Clear the input buffer (consume any leftover characters)
        if (fgets(Fname, STRLEN, stdin) != NULL) {
            /* Replace newline with null terminator */
            size_t len = Fname.size();
            if (len > 0 && Fname[len - 1] == '\n') {
                Fname[len - 1] = '\0';
            }

            /* terminate with a period. */
            if (Fname.size() == 1 && Fname[0] == '.') {
                /* return a NULL pointer if '.' entered. */
                return (NULL);
            }

            /* open the input & check the version. */
            file = std::ifstream(Fname, std::ios::in);
            if (!file.is_open()) {
                /* cannot open the input. */
                std::cout << "File does not exist.";
            }
            else {
                if (CheckFileVersionQ(file, Version)) {
                    return file;
                }
                else {
                    file.close();
                }
            }
        }
    }
}

/*******************************************************************************
 *  Find number of media in the list. At the same time, check the
 *  optical parameters.
 ****/
bool FindNumMediaQ(std::ifstream file, short* NumMediaP)
{
    std::string name;
    short num_media = 0;

    while (1) {
        std::string buf = FindDataLine(file);

        if (buf[0] == '\0') {
            std::cout << "Missing end.\n";
            return (0);
        }
        else if (buf.find("end") != std::string::npos) {
            break;
        }
        else {
            num_media++;

            double n, mua, mus, g;
            sscanf_s(buf, "%s%lf%lf%lf%lf", name, (unsigned)_countof(name), &n, &mua, &mus, &g);
            if (n <= 0 || mua < 0 || mus < 0 || g < -1 || g > 1) {
                printf("Bad optical parameters in %s\n", name);
                return (0);
            }
        }
    }

    *NumMediaP = num_media;
    return (1);
}

/*******************************************************************************
 *  Read the parameters of one medium, assumming the
 *  parameters have been checked with FindNumMediaQ().
 ****/
bool ReadOneMediumQ(std::ifstream file, Layer& MediumP)
{
    std::string buf;

    strcpy_s(buf, sizeof(buf), FindDataLine(file));
    if (buf[0] == '\0') {
        std::cout << "Shouldn't happen here!";
        return (0);
    }
    sscanf_s(buf, "%s%lf%lf%lf%lf", MediumP->medium, (unsigned)_countof(MediumP->medium), &MediumP->n, &MediumP->mua, &MediumP->mus, &MediumP->g);

    return (1);
}

/*******************************************************************************
 *  Read the media list.
 ****/
bool ReadMediumListQ(std::fstream& file, RunParams& run_params)
{
    long file_pos = ftell(Fp);
    if (!FindNumMediaQ(file, &(run_params.num_media))) {
        return (0);
    }
    fseek(file, file_pos, SEEK_SET);

    /* allocate an array for the layer parameters. */
    run_params.medium_list = (Layer*)malloc((run_params.num_media) * sizeof(Layer));

    if (!(run_params.medium_list)) {
        std::cout << "allocation failure in ReadMediumListQ()";
        return (0);
    }
    for (short i = 0; i < run_params.num_media; i++) {
        ReadOneMediumQ(file, &(run_params.medium_list[i]));
    }
    FindDataLine(Fp);		/* skip the signal end. */

    return (1);
}

/**************************************************************************
 *	Read the input name and the input format.
 *
 *	The input format can be either A for ASCII or B for binary.
 ****/
bool ReadFnameFormatQ(std::ifstream file, RunParams& run_params)
{
    std::string buf;
    strcpy_s(buf, sizeof(buf), FindDataLine(Fp));

    if (sscanf_s(buf, "%s %c", run_params.output_filename, (unsigned)_countof(run_params.output_filename), &run_params.output_file_format, 1) != 2) {
        std::cout << "Reading file name and format.\n";
        return (0);
    }
    /* if (toupper(run_params.output_file_format) != 'B') */
    run_params.output_file_format = 'A';	/* now only support 'A' format. */

    return (1);
}

/**************************************************************************
 *	Read the RunParams members dz, dr and dt.
 ****/
bool ReadDzDrDtQ(std::ifstream file, RunParams& run_params)
{
    std::string buf;
    strcpy_s(buf, sizeof(buf), FindDataLine(Fp));

    if (sscanf_s(buf, "%lf%lf%lf", &run_params.dz, &run_params.dr, &run_params.dt) != 3) {
        std::cout << "Reading dz, dr, dt. \n";
        return (0);
    }
    if (run_params.dz <= 0) {
        std::cout << "Nonpositive dz.\n";
        return (0);
    }
    if (run_params.dr <= 0) {
        std::cout << "Nonpositive dr.\n";
        return (0);
    }
    if (run_params.dt <= 0) {
        std::cout << "Nonpositve dt. \n";
        return (0);
    }
    return (1);
}

/**************************************************************************
 *	Read the RunParams members nz, nr, nt, na.
 ****/
bool ReadNzNrNtNaQ(std::ifstream file, RunParams& run_params)
{
    std::string buf;

    /** read in number of dz, dr, da, dt. **/
    strcpy_s(buf, sizeof(buf), FindDataLine(Fp));

    float fz, fr, fa, ft;
    if (sscanf_s(buf, "%f%f%f%f", &fz, &fr, &ft, &fa) != 4) {
        std::cout << "Reading number of dz, dr, dt, da's.\n";
        return (0);
    }
    if (fz <= 0) {
        std::cout << "Nonpositive number of dz's.\n";
        return (0);
    }
    if (fr <= 0) {
        std::cout << "Nonpositive number of dr's.\n";
        return (0);
    }
    if (fa <= 0) {
        std::cout << "Nonpositive number of da's.\n";
        return (0);
    }
    if (ft <= 0) {
        std::cout << "Nonpositive number of dt's.\n";
        return (0);
    }

    run_params.nz = (short)fz;
    run_params.nr = (short)fr;
    run_params.nt = (short)ft;
    run_params.na = (short)fa;
    run_params.da = 0.5 * PI / run_params.na;

    return (1);
}

/**************************************************************************
 *   Initialize the RecordStru.
 ****/
void InitRecord(RunParams& run_params)
{
    run_params.record.Rd_r = 0;
    run_params.record.Rd_a = 0;
    run_params.record.Rd_ra = 0;
    run_params.record.Rd_t = 0;
    run_params.record.Rd_rt = 0;
    run_params.record.Rd_at = 0;
    run_params.record.Rd_rat = 0;
    run_params.record.Td_r = 0;
    run_params.record.Td_a = 0;
    run_params.record.Td_ra = 0;
    run_params.record.Td_t = 0;
    run_params.record.Td_rt = 0;
    run_params.record.Td_at = 0;
    run_params.record.Td_rat = 0;
    run_params.record.A_z = 0;
    run_params.record.A_rz = 0;
    run_params.record.A_t = 0;
    run_params.record.A_zt = 0;
    run_params.record.A_rzt = 0;
}

/**************************************************************************
 *	Change all characters in a string to upper case.
 ****/
char* ToUpperString(char* string)
{
    for (int i = 0; i < (int)strlen(string); i++) {
        string[i] = toupper(string[i]);
    }

    return string;
}

/**************************************************************************
 *  Read which quantity is to be scored.
 ****/
bool ReadRecordQ(std::ifstream file, RunParams& run_params)
{
    std::string buf;
    std::string string;
    bool error = 0;

    strcpy_s(buf, sizeof(buf), FindDataLine(Fp));
    char* index = buf;

    if (index[0] == '\0') {
        std::cout << "Read scored quantities.\n";
        error = 1;
    }
    while (!error && index[0] != '\n' && index[0] != '#') {
        if (index[0] == '\\') {	/* use '\' to continue in next line. */
            strcpy_s(buf, sizeof(buf), FindDataLine(Fp));
            index = buf;
        }
        sscanf_s(index, "%s", string, (unsigned)_countof(string));
        index = index + (char)strlen(string);
        while (index[0] == ' ' || index[0] == '\t') {
            index++;
        }

        if (strcmp(ToUpperString(string), "RD_R") == 0) {
            run_params.record.Rd_r = 1;
        }
        else if (strcmp(ToUpperString(string), "RD_A") == 0) {
            run_params.record.Rd_a = 1;
        }
        else if (strcmp(ToUpperString(string), "RD_RA") == 0) {
            run_params.record.Rd_ra = 1;
        }
        else if (strcmp(ToUpperString(string), "RD_T") == 0) {
            run_params.record.Rd_t = 1;
        }
        else if (strcmp(ToUpperString(string), "RD_RT") == 0) {
            run_params.record.Rd_rt = 1;
        }
        else if (strcmp(ToUpperString(string), "RD_AT") == 0) {
            run_params.record.Rd_at = 1;
        }
        else if (strcmp(ToUpperString(string), "RD_RAT") == 0) {
            run_params.record.Rd_rat = 1;
        }
        else if (strcmp(ToUpperString(string), "TD_R") == 0) {
            run_params.record.Td_r = 1;
        }
        else if (strcmp(ToUpperString(string), "TD_A") == 0) {
            run_params.record.Td_a = 1;
        }
        else if (strcmp(ToUpperString(string), "TD_RA") == 0) {
            run_params.record.Td_ra = 1;
        }
        else if (strcmp(ToUpperString(string), "TD_T") == 0) {
            run_params.record.Td_t = 1;
        }
        else if (strcmp(ToUpperString(string), "TD_RT") == 0) {
            run_params.record.Td_rt = 1;
        }
        else if (strcmp(ToUpperString(string), "TD_AT") == 0) {
            run_params.record.Td_at = 1;
        }
        else if (strcmp(ToUpperString(string), "TD_RAT") == 0) {
            run_params.record.Td_rat = 1;
        }
        else if (strcmp(ToUpperString(string), "A_Z") == 0) {
            run_params.record.A_z = 1;
        }
        else if (strcmp(ToUpperString(string), "A_RZ") == 0) {
            run_params.record.A_rz = 1;
        }
        else if (strcmp(ToUpperString(string), "A_T") == 0) {
            run_params.record.A_t = 1;
        }
        else if (strcmp(ToUpperString(string), "A_ZT") == 0) {
            run_params.record.A_zt = 1;
        }
        else if (strcmp(ToUpperString(string), "A_RZT") == 0) {
            run_params.record.A_rzt = 1;
        }
        else {
            printf("Unknown quantity: %s\n", string);
            error = 1;
        }
    }

    return (!error);
}

/**************************************************************************
*   Filter the RecordStru.
****/
bool FilterRecordQ(std::ifstream file, RunParams& run_params)
{
    InitRecord(run_params);

    if (!ReadRecordQ(file, run_params)) {
        return (0);
    }

    if (run_params.record.Rd_rat) {
        run_params.record.Rd_ra = 0;
        run_params.record.Rd_rt = 0;
        run_params.record.Rd_at = 0;
        run_params.record.Rd_r = 0;
        run_params.record.Rd_a = 0;
        run_params.record.Rd_t = 0;
    }
    if (run_params.record.Rd_ra) {
        run_params.record.Rd_r = 0;
        run_params.record.Rd_a = 0;
    }
    if (run_params.record.Rd_rt) {
        run_params.record.Rd_r = 0;
        run_params.record.Rd_t = 0;
    }
    if (run_params.record.Rd_at) {
        run_params.record.Rd_a = 0;
        run_params.record.Rd_t = 0;
    }
    if (run_params.record.Td_rat) {
        run_params.record.Td_ra = 0;
        run_params.record.Td_rt = 0;
        run_params.record.Td_at = 0;
        run_params.record.Td_r = 0;
        run_params.record.Td_a = 0;
        run_params.record.Td_t = 0;
    }
    if (run_params.record.Td_ra) {
        run_params.record.Td_r = 0;
        run_params.record.Td_a = 0;
    }
    if (run_params.record.Td_rt) {
        run_params.record.Td_r = 0;
        run_params.record.Td_t = 0;
    }
    if (run_params.record.Td_at) {
        run_params.record.Td_a = 0;
        run_params.record.Td_t = 0;
    }
    if (run_params.record.A_rzt) {
        run_params.record.A_rz = 0;
        run_params.record.A_zt = 0;
        run_params.record.A_z = 0;
        run_params.record.A_t = 0;
    }
    if (run_params.record.A_rz) {
        run_params.record.A_z = 0;
    }
    if (run_params.record.A_zt) {
        run_params.record.A_z = 0;
        run_params.record.A_t = 0;
    }
    if (run_params.record.A_zt) {
        run_params.record.A_z = 0;
        run_params.record.A_t = 0;
    }
    return (1);
}

/**************************************************************************
 *  Read the threshold min_weight.
 ****/
bool ReadWthQ(std::ifstream file, RunParams& run_params)
{
    std::string buf;
    strcpy_s(buf, sizeof(buf), FindDataLine(Fp));

    if (sscanf_s(buf, "%lf", &(run_params.min_weight)) != 1) {
        std::cout << "Reading threshold weight.\n";
        return (0);
    }
    if (run_params.min_weight < 0 || run_params.min_weight >= 1.0) {
        std::cout << "Threshold weight out of range (0-1).\n";
        return (0);
    }
    return (1);
}

/**************************************************************************
 *  Read the random number seed.
 ****/
bool ReadRanSeedQ(std::ifstream file, RunParams& run_params)
{
    std::string buf;
    strcpy_s(buf, sizeof(buf), FindDataLine(Fp));

    if (sscanf_s(buf, "%ld", &(run_params.random_seed)) != 1) {
        std::cout << "Reading random number seed.\n";
        return (0);
    }
    if (run_params.random_seed < 0) {
        std::cout << "Nonpositive random number seed.\n";
        return (0);
    }
    return (1);
}

/**************************************************************************
 *  Find number of layers.
 ****/
bool FindNumLayersQ(std::ifstream file, int* NumLayerP)
{
    std::string buf;
    std::string name;
    short num_layers = 0;
    double thick;

    while (1) {
        strcpy_s(buf, sizeof(buf), FindDataLine(Fp));

        if (buf[0] == '\0') {
            std::cout << "Missing end.\n";
            return (0);
        }
        else if (strstr(buf, "end") != NULL) {
            break;
        }
        else if ((sscanf_s(buf, "%s %lf", name, (unsigned)_countof(name), &thick) == 2) || (sscanf_s(buf, "%s", name, (unsigned)_countof(name)) == 1)) {
            num_layers++;
        }
    }

    *NumLayerP = num_layers - 2;
    return (1);
}

/**************************************************************************
 *  Check whether the medium name is in the media list.
 ****/
bool ValidMediumNameQ(char* NameP, int* Index, RunParams& run_params)
{
    for (short i = 0; i < run_params.num_media; i++) {
        if (!strcmp(NameP, run_params.medium_list[i].medium)) {
            *Index = i;
            return (1);
        }
    }
    return (0);
}

/**************************************************************************
 *	Read the parameters of all layers.
 ****/
bool ReadLayerSpecsQ(std::ifstream file, RunParams& run_params)
{
    std::string buf;
    std::string name;
    short i = 0;

    /* z coordinate of the current layer. */
    double z = 0.0;
    double thick = 0.0;

    long file_pos = ftell(Fp);

    int num_layers;
    if (!FindNumLayersQ(file, &num_layers)) {
        return (0);
    }
    fseek(file, file_pos, SEEK_SET);

    run_params.num_layers = num_layers;
    /* Allocate an array for the layer parameters. */
    /* layer 0 and layer Num_Layers + 1 are for ambient. */
    run_params.layer = (Layer*)malloc((unsigned)(num_layers + 2) * sizeof(Layer));
    if (!(run_params.layer)) {
        std::cout << "allocation failure in ReadLayerSpecsQ()";
        return (0);
    }
    for (i = 0; i <= num_layers + 1; i++) {
        strcpy_s(buf, sizeof(buf), FindDataLine(Fp));
        if (i == 0 || i == num_layers + 1) {
            if (sscanf_s(buf, "%s", name, (unsigned)_countof(name)) != 1) {
                std::cout << "  Error in reading ambient layer name.\n";
                return (0);
            }
        }
        else {
            if (sscanf_s(buf, "%s%lf", name, (unsigned)_countof(name), &thick) != 2) {
                std::cout << "  Error in ReadLayerSpecsQ().\n";
                return (0);
            }
            else if (thick <= 0.0) {
                std::cout << "  Nonpositive layer thickness.\n";
                return (0);
            }
        }

        int index;
        if (!ValidMediumNameQ(name, &index, run_params)) {
            std::cout << "  Invalid medium name. \n";
            return (0);
        }
        strcpy_s(run_params.layer[i].medium, sizeof(run_params.layer[i].medium), name);
        run_params.layer[i].n = run_params.medium_list[index].n;
        run_params.layer[i].mua = run_params.medium_list[index].mua;
        run_params.layer[i].mus = run_params.medium_list[index].mus;
        run_params.layer[i].g = run_params.medium_list[index].g;

        if ((i != 0) && (i != num_layers + 1)) {
            run_params.layer[i].z0 = z;
            z = z + thick;
            run_params.layer[i].z1 = z;
        }
        else if (i == 0) {
            run_params.layer[i].z0 = 0.0;
            run_params.layer[i].z1 = 0.0;
        }
        else if (i == (num_layers + 1)) {
            run_params.layer[i].z0 = z;
            run_params.layer[i].z1 = z;
        }
    }

    return (1);
}

/**************************************************************************
 *  Read the number of photons.
 *  Read computation time limit.
 *  Type = 0, read from a .mci input file;
 *  Type = 1, read from a .mco output file.
 ****/
bool ReadNumPhotonsQ(std::istream& input, RunParams& run_params, char Type)
{
    std::string buf = FindDataLine(reinterpret_cast<std::fstream&>(input));

    if (Type == 0) {
        run_params.add_num_photons = 0;
        run_params.add_limit_seconds = 0;
    }

    float temp;
    int hours, minutes;
    if (sscanf_s(buf, "%f %d:%d", &temp, &hours, &minutes) == 3) {
        if (((long)temp > 0) && (hours * 3600 + minutes * 60) >= 0) {
            if (Type == 0) {
                run_params.num_photons = (long)temp;
                run_params.time_limit_seconds = hours * 3600 + minutes * 60;
            }
            else {
                run_params.add_num_photons = (long)temp;
                run_params.add_limit_seconds = hours * 3600 + minutes * 60;
            }

            run_params.control_bit = ControlBit::Both;
        }
        else {
            std::cout << "Nonpositive number of photons or time limit.\n";
            return (0);
        }

    }
    else if (sscanf_s(buf, "%d:%d", &hours, &minutes) == 2) {
        if ((hours * 3600 + minutes * 60) >= 0) {
            if (Type == 0) {
                run_params.time_limit_seconds = hours * 3600 + minutes * 60;
            }
            else {
                run_params.add_limit_seconds = hours * 3600 + minutes * 60;
            }

            run_params.control_bit = ControlBit::TimeLimit;
        }
        else {
            std::cout << "Nonpositive time limit.\n";
            return (0);
        }

    }
    else if (sscanf_s(buf, "%f", &temp) == 1) {
        if ((long)temp > 0) {
            if (Type == 0) {
                run_params.num_photons = (long)temp;
            }
            else {
                run_params.add_num_photons = (long)temp;
            }
            run_params.control_bit = ControlBit::NumPhotons;
        }
        else {
            std::cout << "Nonpositive number of photons.\n";
            return (0);
        }

    }
    else {
        std::cout << "Invalid number of photons or time limit.\n";
        return (0);
    }

    return (1);
}

/**************************************************************************
 *  Read the beam source type (Pencil/Isotropic).
 ****/
bool ReadSourceTypeQ(std::ifstream file, RunParams& run_params)
{
    std::string b_type;
    std::string buf = FindDataLine(Fp);

    if (sscanf_s(buf, "%s", b_type, (unsigned)_countof(b_type)) != 1) {
        std::cout << "Reading photon source type. \n";
        return (0);
    }
    if (strcmp(b_type, "pencil") == 0) {
        run_params.source = BeamType::Pencil;
        return (1);
    }
    else if (strcmp(b_type, "isotropic") == 0) {
        run_params.source = BeamType::Isotropic;
        return (1);
    }
    else {
        std::cout << "Unknow photon source type. \n";
        return (0);
    }
}

/**************************************************************************
 *      Compute the index to layer according to the z coordinate.
 *	If the z is on an interface between layers, the returned index
 *	will point to the upper layer.
 *	Index 0 is the top ambient medium and index num_layers+1 is the
 *	bottom one.
 ****/
bool ZToLayerQ(double z, short* index, RunParams& run_params)
{
    /* index to layer. */
    short i = 0;
    short num_layers = run_params.num_layers;

    if (z < 0.0) {
        std::cout << "Nonpositive z coordinate.\n";
        return (0);
    }
    else if (z > run_params.layer[num_layers].z1) {
        std::cout << "Source is outside of the last layer. \n";
        return (0);
    }
    else {
        while (z > run_params.layer[i].z1) { i++; }
        *index = i;
        return (1);
    }
}

/**************************************************************************
 *  Read starting position of photon source.
 ****/
bool ReadStartPQ(std::ifstream file, RunParams& run_params)
{
    std::string buf;
    strcpy_s(buf, sizeof(buf), FindDataLine(Fp));

    double source_z;
    short slayer;
    std::string medium_name;

    /* z and medium. */
    if ((sscanf_s(buf, "%lf %s", &source_z, medium_name, (unsigned)_countof(medium_name)) == 2) && (medium_name[0] != '#' && medium_name[0] != '\n')) {
        if (!ZToLayerQ(source_z, &slayer, run_params)) {
            return (0);
        }

        if (strcmp(run_params.layer[slayer].medium, medium_name) != 0) {
            if ((fabs(source_z - run_params.layer[slayer].z1) < DBL_EPSILON) && (strcmp(run_params.layer[slayer + 1].medium, medium_name) == 0)) {
                slayer++;
                if (slayer > run_params.num_layers) {
                    puts("Source is outside of the last layer.");
                    return (0);
                }

            }
            else {
                std::cout << "Medium name and z coordinate do not match.\n";
                return (0);
            }
        }

    }
    /* z only. */
    else if (sscanf_s(buf, "%lf", &(source_z)) == 1) {
        if (!ZToLayerQ(source_z, &run_params.slayer, run_params)) {
            return (0);
        }
        strcpy_s(medium_name, sizeof(medium_name), "");
    }
    else {
        std::cout << "Invalid starting position of photon source.\n";
        return (0);
    }

    if ((run_params.source == BeamType::Isotropic) && (source_z == 0.0)) {
        std::cout << "Can not put isotropic source in upper ambient medium.\n";
        return 0;
    }

    run_params.source_z = source_z;
    strcpy_s(run_params.medium_name, sizeof(run_params.medium_name), medium_name);
    return (1);
}

/*************************************************************************
 *	Compute the critical angles for total internal
 *	reflection according to the relative refractive index
 *	of the layer.
 *	All layers are processed.
 ****/
void CriticalAngle(short Num_Layers, Layer** Layerspecs_PP)
{
    for (short i = 1; i <= Num_Layers; i++) {
        double n1 = (*Layerspecs_PP)[i].n;
        double n2 = (*Layerspecs_PP)[i - 1].n;
        (*Layerspecs_PP)[i].cos_crit0 = n1 > n2 ? sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;

        n2 = (*Layerspecs_PP)[i + 1].n;
        (*Layerspecs_PP)[i].cos_crit1 = n1 > n2 ? sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;
    }
}

/**************************************************************************
 *	Read in the input parameters for one run.
 ****/
void ReadRunParam(std::fstream& file, RunParams& run_params)
{
    if (!ReadFnameFormatQ(file, run_params)) {
        exit(1);
    }

    /* geometry. */
    if (!ReadLayerSpecsQ(file, run_params)) {
        exit(1);
    }
    FindDataLine(Fp);		/* skip the signal "end" of layers. */

    /* source. */
    if (!ReadSourceTypeQ(file, run_params)) {
        exit(1);
    }
    if (!ReadStartPQ(file, run_params)) {
        exit(1);
    }

    /* grids. */
    if (!ReadDzDrDtQ(file, run_params)) {
        exit(1);
    }
    if (!ReadNzNrNtNaQ(file, run_params)) {
        exit(1);
    }
    run_params.zm = run_params.dz * run_params.nz;
    run_params.rm = run_params.dr * run_params.nr;
    run_params.tm = run_params.dt * run_params.nt;
    run_params.am = run_params.da * run_params.na;

    /* scored data categories. */
    if (!FilterRecordQ(file, run_params)) {
        exit(1);
    }

    /* simulation control. */
    if (!ReadNumPhotonsQ(file, run_params, 0)) {
        exit(1);
    }
    if (!ReadWthQ(file, run_params)) {
        exit(1);
    }
    if (!ReadRanSeedQ(file, run_params)) {
        exit(1);
    }

    CriticalAngle(run_params.num_layers, &run_params.layer);
}

/**************************************************************************
 *  Read the media list in interactive mode.
 ****/
void InterReadMediumList(RunParams& run_params)
{
    std::string string;
    std::string medium;

    std::cout << "Specify medium list. Total number of mediums: ";

    fgets(string, STRLEN, stdin);
    while (sscanf_s(string, "%hd", &run_params.num_media) != 1 || run_params.num_media <= 0) {
        std::cout << "Invalid medium number. Input again: ";
        fgets(string, STRLEN, stdin);
    }

    /* allocate an array for the layer parameters. */
    run_params.medium_list = (Layer*)malloc((run_params.num_media) * sizeof(Layer));

    if (!(run_params.medium_list)) {
        std::cout << "allocation failure in ReadMediumListQ()";
        exit(1);
    }

    for (short i = 0; i < run_params.num_media; i++) {
        printf("Specify medium %d: \n  Medium name: ", i + 1);
        fgets(string, STRLEN, stdin);

        bool name_taken = 0;
        do {
            sscanf_s(string, "%s", medium, (unsigned)_countof(medium));
            for (short j = 0; j < i; j++) {
                if (strcmp(run_params.medium_list[j].medium, medium) == 0) {
                    name_taken = 1;
                    std::cout << "  Duplicate medium. Input again: ";
                    fgets(string, STRLEN, stdin);
                    break;
                }
            }
        } while (name_taken);

        strcpy_s(run_params.medium_list[i].medium, sizeof(run_params.medium_list[i].medium), medium);

        std::cout << "  Refractive index n (>= 1.0): ";
        fgets(string, STRLEN, stdin);
        while (sscanf_s(string, "%lf", &run_params.medium_list[i].n) != 1 || run_params.medium_list[i].n < 1.0) {
            std::cout << "  Invalid refractive index. Input again (>= 1.0): ";
            fgets(string, STRLEN, stdin);
        }

        std::cout << "  Absorption coefficient mua (>= 0.0 /cm): ";
        fgets(string, STRLEN, stdin);
        while (sscanf_s(string, "%lf", &run_params.medium_list[i].mua) != 1 || run_params.medium_list[i].mua < 0.0) {
            std::cout << "  Invalid absorption coefficient. Input again (>= 0.0): ";
            fgets(string, STRLEN, stdin);
        }

        std::cout << "  Scattering coefficient mus (>= 0.0 /cm): ";
        fgets(string, STRLEN, stdin);
        while (sscanf_s(string, "%lf", &run_params.medium_list[i].mus) != 1 || run_params.medium_list[i].mus < 0.0) {
            std::cout << "  Invalid scattering coefficient. Input again (>= 0.0): ";
            fgets(string, STRLEN, stdin);
        }

        std::cout << "  Anisotropy factor g (0.0 - 1.0): ";
        fgets(string, STRLEN, stdin);
        while (sscanf_s(string, "%lf", &run_params.medium_list[i].g) != 1 || run_params.medium_list[i].g < 0.0 || run_params.medium_list[i].g > 1.0) {
            std::cout << "  Invalid anisotropy factor. Input again (0.0 - 1.0): ";
            fgets(string, STRLEN, stdin);
        }
        std::cout << "\n";
    }
}

/**************************************************************************
 *  Read the input name and the input format interactively.
 ****/
void InterReadFnameFormat(RunParams& run_params)
{
    FILE* file;
    std::string fname;
    std::string fmode;

    do {
        std::cout << "Specify output filename with extension .mco: ";
        fgets(fname, STRLEN, stdin);
        fmode[0] = 'w';

        /* input exists. */
        if (!fopen_s(&file, fname, "r")) {
            std::cout << "File %s exists, %s", fname, "w=overwrite, n=new filename: ";

            /* avoid null line. */
            do { fgets(fmode, STRLEN, stdin); } while (!strlen(fmode));
            fclose(file);
        }
    } while (fmode[0] != 'w');

    strcpy_s(run_params.output_filename, sizeof(run_params.output_filename), fname);

    //std::cout << "Output input format (A/B): ";
    //fgets(fname, STRLEN, stdin);
    //while (sscanf_s(fname, "%c", &run_params.output_file_format, 1) != 1) {
    //    std::cout << "Error occured. Output input format (A/B): ";
    //    fgets(fname, STRLEN, stdin);
    //}

    //if (toupper(run_params.output_file_format) != 'B') {
    //    run_params.output_file_format = 'A';
    //}

    /* now only support 'A' format. */
    run_params.output_file_format = 'A';

    std::cout << "\n";
}

/**************************************************************************
 *	Read dz, dr, dt interactively.
 ****/
void InterReadDzDrDt(RunParams& run_params)
{
    do {
        std::cout << "Specify dz, dr, dt in one line\n";
        std::cout << "(all > 0.0 cm, e.g., 0.1 0.1 0.1): ";
    } while (!ReadDzDrDtQ(stdin, run_params));

    std::cout << "\n";
}

/**************************************************************************
 *      Read the RunParams members nz, nr, na interactively.
 ****/
void InterReadNzNrNtNa(RunParams& run_params)
{
    do {
        std::cout << "Specify nz, nr, nt, na in one line\n";
        std::cout << "(all > 0, e.g., 100 100 100 100): ";
    } while (!ReadNzNrNtNaQ(stdin, run_params));

    run_params.da = 0.5 * PI / run_params.na;
    std::cout << "\n";
}

/**************************************************************************
 *	Read and filter the quantities to be scored interactively.
 ****/
void InterFilterRecord(RunParams& run_params)
{
    do {
        std::cout << "Select scored quantities from the following data categories:\n";
        std::cout << "\tRd_rat\t\t\tTd_rat\t\t\tA_rzt\n";
        std::cout << "\tRd_ra\tRd_rt\tRd_at\tTd_ra\tTd_rt\tRd_at\tA_rz\tA_zt\n";
        std::cout << "\tRd_r\tRd_a\tRd_t\tTd_r\tTd_a\tTd_t\tA_z\tA_t\n";
    } while (!FilterRecordQ(stdin, run_params));

    std::cout << "\n";
}

/**************************************************************************
 *	Read the threshold min_weight interactively.
 ****/
void InterReadWth(RunParams& run_params)
{
    std::string string;

    std::cout << "Input threshold weight (0 <= wth < 1.0, 0.0001 recommended): ";
    fgets(string, STRLEN, stdin);
    while (sscanf_s(string, "%lf", &run_params.min_weight) != 1 || run_params.min_weight < 0 || run_params.min_weight >= 1) {
        std::cout << "Invalid wth. Input again (0 <= wth < 1.0): ";
        fgets(string, STRLEN, stdin);
    }

    std::cout << "\n";
}

/**************************************************************************
 *	Read the random seed interactively.
 ****/
void InterReadRanSeed(RunParams& run_params)
{
    std::string string;

    std::cout << "Input random number seed (1 <= ran_seed <= 32000): ";
    fgets(string, STRLEN, stdin);
    while (sscanf_s(string, "%ld", &run_params.random_seed) != 1 || run_params.random_seed < 1 || run_params.random_seed > 32000) {
        std::cout << "Invalid ran_seed. Input again (1 <= ran_seed <= 32000): ";
        fgets(string, STRLEN, stdin);
    }

    std::cout << "\n";
}

/**************************************************************************
 ****/
void PrintMediumNames(RunParams& run_params)
{
    std::cout << "Available medium types:\n";

    int j = 1;
    for (int i = 0; i < run_params.num_media; i++) {
        printf("%-16s", run_params.medium_list[i].medium);
        if (j % 4 == 0) {
            std::cout << "\n";
        }
        j++;
    }

    std::cout << "\n";
}

/**************************************************************************
 *	Read layer specifications interactively.
 ****/
void InterReadLayerSpecs(RunParams& run_params)
{
    std::string string;
    std::string name;

    /* z coordinate of the current layer. */
    double z = 0.0;

    int index;

    std::cout << "\nSpecify layer list. ";
    PrintMediumNames(run_params);
    std::cout << "\nTotal number of layers: ";

    fgets(string, STRLEN, stdin);
    while (sscanf_s(string, "%hd", &run_params.num_layers) != 1 || run_params.num_layers <= 0) {
        std::cout << "Invalid layer number. Input again: ";
        fgets(string, STRLEN, stdin);
    }

    /* Allocate an array for the layer parameters. */
    /* layer 0 and layer Num_Layers + 1 are for ambient. */
    run_params.layer = (Layer*)malloc((unsigned)(run_params.num_layers + 2) * sizeof(Layer));

    if (!(run_params.layer)) {
        std::cout << "allocation failure in ReadLayerSpecsQ()";
        exit(1);
    }

    for (short i = 0; i <= run_params.num_layers + 1; i++) {
        bool error = 1;
        while (error) {
            error = 0;
            if (i == 0) {
                std::cout << "\n  Name of upper ambient medium: ";
                fgets(string, STRLEN, stdin);
                sscanf_s(string, "%s", name, (unsigned)_countof(name));
            }
            else if (i == run_params.num_layers + 1) {
                std::cout << "\n  Name of lower ambient medium: ";
                fgets(string, STRLEN, stdin);
                sscanf_s(string, "%s", name, (unsigned)_countof(name));
            }
            else {
                printf("\n  Medium name of layer %d: ", i);
                fgets(string, STRLEN, stdin);
                sscanf_s(string, "%s", name, (unsigned)_countof(name));
            }

            if (!ValidMediumNameQ(name, &index, run_params)) {
                std::cout << "  Invalid medium name. Input again.";
                error = 1;
            }
        }

        strcpy_s(run_params.layer[i].medium, sizeof(run_params.layer[i].medium), name);
        run_params.layer[i].n = run_params.medium_list[index].n;
        run_params.layer[i].mua = run_params.medium_list[index].mua;
        run_params.layer[i].mus = run_params.medium_list[index].mus;
        run_params.layer[i].g = run_params.medium_list[index].g;

        if ((i != 0) && (i != run_params.num_layers + 1)) {
            printf("  Input the thickness of layer %d (thickness > 0.0 cm): ", i);
            fgets(string, STRLEN, stdin);

            double thick = 0.0;
            while (sscanf_s(string, "%lf", &thick) != 1 || thick <= 0) {
                std::cout << "  Invalid thickness. Input again (thickness > 0.0 cm): ";
                fgets(string, STRLEN, stdin);
            }
            run_params.layer[i].z0 = z;
            z = z + thick;
            run_params.layer[i].z1 = z;
        }
        else if (i == 0) {
            run_params.layer[i].z0 = 0.0;
            run_params.layer[i].z1 = 0.0;
        }
        else if (i == run_params.num_layers + 1) {
            run_params.layer[i].z0 = z;
            run_params.layer[i].z1 = z;
        }
    }

    std::cout << "\n";
}

/**************************************************************************
 *  Read the number of photons, or computation time interactively.
 ****/
void InterReadNumPhotons(RunParams& run_params)
{
    std::cout << "Specify number of photons or time in hh:mm format,\n";
    std::cout << "or both in one line (e.g. 10000 5:30): ";

    while (!ReadNumPhotonsQ(stdin, run_params, 0)) {
        std::cout << "Input again: ";
    }

    std::cout << "\n";
}

/**************************************************************************
 *  Read the beam source type (Pencil/Isotropic).
 ****/
void InterReadSourceType(RunParams& run_params)
{
    std::string string;

    std::cout << "Input source type (p = pencil / i = isotropic): ";
    fgets(string, STRLEN, stdin);

    char c;
    while (sscanf_s(string, "%c", &c, 1) != 1 || !(toupper(c) == 'P' || toupper(c) == 'I')) {
        std::cout << "Invalid type. Input again (p = pencil / i = isotropic): ";
        fgets(string, STRLEN, stdin);
    }

    if (toupper(c) == 'P') {
        run_params.source = BeamType::Pencil;
    }
    else {
        run_params.source = BeamType::Isotropic;
    }

    std::cout << "\n";
}

/**************************************************************************
 *  Read starting position of photon source.
 ****/
void InterReadStartP(RunParams& run_params)
{
    do {
        printf("Input the z coordinate of source (0.0 - %f cm) and the medium\n", run_params.layers[run_params.layers.size()-1].z1);
        std::cout << "where the source is if the z is on an interface (e.g. 1.0 [air]):";
    } while (!ReadStartPQ(stdin, run_params));

    std::cout << "\n";
}

/*************************************************************************
 *  If input is stdout, freeze the screen and print a more message on screen
 *  every 20 lines. The Line is the line index.
 ****/
void More(std::fstream& file, int* Line)
{
    if (Fp == stdout) {
        if (!((*Line) % 20)) {
            std::cout << "--More-- (Press Return to continue)";
            fflush(Fp);

            char c;
            do {
                fread(&c, 1, 1, stdin);
            } while (c != '\n');
        }
    }
}

/*************************************************************************
 *     Write medium list to the input input.
 *     if input is stdout, freeze the screen every 20 lines.
 ****/
void PutMediumListToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    std::string format;

    file << std::format("# Specify media \n");
    (*Line)++;
    file << std::format("#\tname\t\tn\tmua\tmus\tg\n");
    (*Line)++;

    for (int i = 0; i < run_params.media.size(); i++) {
        Layer s;

        More(file, Line);
        s = run_params.media[i];
        if (strlen(s.medium) + 1 > 8) {
            strcpy_s(format, sizeof(format), "\t%s \t%G\t%G\t%G\t%G\n");
        }
        else {
            strcpy_s(format, sizeof(format), "\t%s \t\t%G\t%G\t%G\t%G\n");
        }

        file << std::format(format, s.medium, s.n, s.mua, s.mus, s.g);
        (*Line)++;
    }
    file << std::format("end #of media\n");
    (*Line)++;
}

void PutFnameFormatToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    More(file, Line);
    file << std::format("%s \t%c\t\t\t# output file name, format.\n", run_params.output_filename, run_params.output_file_format);
    (*Line)++;
}

void PutDzDrDtToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    More(file, Line);
    file << std::format("%G\t%G\t%G\t\t\t# dz, dr, dt.\n", run_params.dz, run_params.dr, run_params.dt);
    (*Line)++;
}

void PutNzNrNtNaToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    More(file, Line);
    file << std::format("%d\t%d\t%d\t%d\t\t# nz, nr, nt, na.\n", run_params.nz, run_params.nr, run_params.nt, run_params.na);
    (*Line)++;
}

void PutScoredToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    More(file, Line);
    file << std::format("# This simulation will score the following categories:\n");
    (*Line)++;

    More(file, Line);
    if (run_params.record.Rd_r) {
        file << std::format("Rd_r \t");
    }
    if (run_params.record.Rd_a) {
        file << std::format("Rd_a \t");
    }
    if (run_params.record.Rd_ra) {
        file << std::format("Rd_ra \t");
    }
    if (run_params.record.Rd_t) {
        file << std::format("Rd_t \t");
    }
    if (run_params.record.Rd_rt) {
        file << std::format("Rd_rt \t");
    }
    if (run_params.record.Rd_at) {
        file << std::format("Rd_at \t");
    }
    if (run_params.record.Rd_rat) {
        file << std::format("Rd_rat \t");
    }

    if (run_params.record.Td_r) {
        file << std::format("Td_r \t");
    }
    if (run_params.record.Td_a) {
        file << std::format("Td_a \t");
    }
    if (run_params.record.Td_ra) {
        file << std::format("Td_ra \t");
    }
    if (run_params.record.Td_t) {
        file << std::format("Td_t \t");
    }
    if (run_params.record.Td_rt) {
        file << std::format("Td_rt \t");
    }
    if (run_params.record.Td_at) {
        file << std::format("Td_at \t");
    }
    if (run_params.record.Td_rat) {
        file << std::format("Td_rat \t");
    }

    if (run_params.record.A_z) {
        file << std::format("A_z \t");
    }
    if (run_params.record.A_rz) {
        file << std::format("A_rz \t");
    }
    if (run_params.record.A_t) {
        file << std::format("A_t \t");
    }
    if (run_params.record.A_zt) {
        file << std::format("A_zt \t");
    }
    if (run_params.record.A_rzt) {
        file << std::format("A_rzt \t");
    }

    file << std::format("\n");
    (*Line)++;
}

void PutWthToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    More(file, Line);
    file << std::format("%G\t\t\t\t\t# threshold weight.\n", run_params.min_weight);
    (*Line)++;
}

void PutRanSeedToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    More(file, Line);
    file << std::format("%ld\t\t\t\t\t# random number seed.\n", run_params.random_seed);
    (*Line)++;
}

void PutLayerSpecsToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    std::string format;

    More(file, Line);
    file << std::format("# \tmedium \t\tthickness\n");
    (*Line)++;

    for (int i = 0; i <= run_params.num_layers + 1; i++) {
        Layer s;
        More(file, Line);

        s = run_params.layer[i];
        if (i != 0 && i != run_params.num_layers + 1) {
            if (strlen(s.medium) + 1 > 8) {
                strcpy_s(format, sizeof(format), "\t%s \t%G\n");
            }
            else {
                strcpy_s(format, sizeof(format), "\t%s \t\t%G\n");
            }
            file << std::format(format, s.medium, s.z1 - s.z0);
        }
        else {
            file << std::format("\t%s\n", s.medium);
        }
        (*Line)++;
    }

    More(file, Line);
    file << std::format("end #of layers\n");
    (*Line)++;
}

void PutNumPhotonsToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    More(file, Line);

    if (run_params.control_bit == 1) {
        file << std::format("%ld  \t\t\t\t\t# no. of photons | time\n", run_params.num_photons);
    }
    else if (run_params.control_bit == 2) {
        file << std::format("%ld:%ld\t\t\t\t\t# no. of photons | time\n", run_params.time_limit_seconds / 3600, run_params.time_limit_seconds % 3600 / 60);
    }
    else {
        file << std::format("%ld  \t%ld:%ld\t\t\t\t# no. of photons | time\n", run_params.num_photons, run_params.time_limit_seconds / 3600, run_params.time_limit_seconds % 3600 / 60);
    }

    (*Line)++;
}

void PutSourceTypeToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    More(file, Line);

    if (run_params.source == BeamType::Pencil) {
        file << std::format("pencil \t\t\t\t\t# src type: pencil/isotropic.\n");
    }
    else {
        file << std::format("isotropic \t\t\t\t# src type: pencil/isotropic.\n");
    }

    (*Line)++;
}

void PutStartPToFile(std::fstream& file, RunParams& run_params, int* Line)
{
    More(file, Line);

    if (strlen(run_params.medium_name) == 0) {
        file << std::format("%G\t\t\t\t\t# starting position of source.\n", run_params.source_z);
    }
    else if (strlen(run_params.medium_name) + 1 > 8) {
        file << std::format("%G\t%s \t\t\t# starting position of source.\n", run_params.source_z, run_params.medium_name);
    }
    else {
        file << std::format("%G\t%s \t\t\t\t# starting position of source.\n", run_params.source_z, run_params.medium_name);
    }

    (*Line)++;
}

/*************************************************************************
 *  Write input parameters to the input input.
 *  If input is stdout, freeze the screen every 20 lines.
 ****/
void PutInputToFile(std::fstream& file, RunParams& run_params)
{
    /* line index. */
    file << std::format("mcmli2.0 \t\t\t# file version \n\n");
    int line = 2;
    PutMediumListToFile(file, run_params, &line);

    More(file, &line);
    file << std::format("\n# Specify data for run 1\n");
    line += 2;

    PutFnameFormatToFile(file, run_params, &line);

    /* geometry. */
    More(file, &line);
    file << std::format("\n");
    line++;
    PutLayerSpecsToFile(file, run_params, &line);

    /* source. */
    More(file, &line);
    file << std::format("\n");
    line++;
    PutSourceTypeToFile(file, run_params, &line);
    PutStartPToFile(file, run_params, &line);

    /* grids. */
    More(file, &line);
    file << std::format("\n");
    line++;
    PutDzDrDtToFile(file, run_params, &line);
    PutNzNrNtNaToFile(file, run_params, &line);

    /* scored data categories. */
    More(file, &line);
    file << std::format("\n");
    line++;
    PutScoredToFile(file, run_params, &line);

    /* simulation control. */
    More(file, &line);
    file << std::format("\n");
    line++;

    PutNumPhotonsToFile(file, run_params, &line);
    PutWthToFile(file, run_params, &line);
    PutRanSeedToFile(file, run_params, &line);

    More(file, &line);
    file << std::format("end #of runs\n\n");
}

/**************************************************************************
 *  Read in the input parameters for one run in interactive mode.
 ****/
void InterReadParam(RunParams& run_params)
{
    std::string string;

    InterReadMediumList(run_params);
    InterReadFnameFormat(run_params);
    InterReadLayerSpecs(run_params);
    InterReadSourceType(run_params);
    InterReadStartP(run_params);
    InterReadDzDrDt(run_params);
    InterReadNzNrNtNa(run_params);
    InterFilterRecord(run_params);
    InterReadNumPhotons(run_params);
    InterReadWth(run_params);
    InterReadRanSeed(run_params);

    std::cout << "Do you want to save the input to a file? (Y/N)";
    fgets(string, STRLEN, stdin);
    if (toupper(string[0]) == 'Y') {
        std::cout << "Give the file name to save input: ( .mci): ";
        fgets(string, STRLEN, stdin);

        FILE* fp;
        if (fopen_s(&fp, string, "w")) {
            puts("Can not open the file to write.");
        }
        else {
            PutInputToFile(fp, run_params);
        }
    }
}

/**************************************************************************
 *  Check consistance of input parameters for one run.
 *  Such as: the consistance of medium list, layer
 *  list, souce starting position and source type.
 ****/
bool CheckInputConsis(RunParams& run_params)
{
    for (int i = 0; i <= run_params.num_layers + 1; i++) {
        int index;
        if (!ValidMediumNameQ(run_params.layer[i].medium, &index, run_params)) {
            printf("Invalid medium name of layer %d.\n", i);
            return 0;
        }
        else {
            run_params.layer[i].n = run_params.medium_list[index].n;
            run_params.layer[i].mua = run_params.medium_list[index].mua;
            run_params.layer[i].mus = run_params.medium_list[index].mus;
            run_params.layer[i].g = run_params.medium_list[index].g;
        }
    }

    if ((run_params.source == BeamType::Isotropic) && (run_params.source_z == 0.0)) {
        std::cout << "Can not put isotropic source in upper ambient medium.\n";
        return 0;
    }
    if (!ZToLayerQ(run_params.source_z, &run_params.slayer, run_params)) {
        return 0;
    }

    if (run_params.medium_name[0] != '\0') {
        if (strcmp(run_params.layer[run_params.slayer].medium, run_params.medium_name) != 0) {
            if ((fabs(run_params.source_z - run_params.layer[run_params.slayer].z1) < DBL_EPSILON) &&
                (strcmp(run_params.layer[run_params.slayer + 1].medium, run_params.medium_name) == 0)) {
                run_params.slayer++;
            }
            else {
                std::cout << "Medium name and z coordinate do not match.\n";
                return (0);
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
        do { fgets(string, STRLEN, stdin); } while (!strlen(string));
    } while (toupper(string[0]) != 'Y' && toupper(string[0]) != 'N');

    return (toupper(string[0]));
}

void ChangeMediumList(RunParams& run_params)
{
    int line = 1;
    std::cout << "Current medium list: \n";
    PutMediumListToFile(stdout, run_params, &line);
    std::cout << "\n";

    if (QuitOrContinue() == 'Y') {
        free(run_params.medium_list);
        InterReadMediumList(run_params);
    }
}

void ChangeFnameFormat(RunParams& run_params)
{
    int line = 1;
    std::cout << "Current output file name and format: \n";
    PutFnameFormatToFile(stdout, run_params, &line);
    std::cout << "\n";
    InterReadFnameFormat(run_params);
}

void ChangeDzDrDt(RunParams& run_params)
{
    int line = 1;
    std::cout << "Current dz, dr, dt: \n";
    PutDzDrDtToFile(stdout, run_params, &line);
    std::cout << "\n";
    InterReadDzDrDt(run_params);
}

void ChangeNzNrNtNa(RunParams& run_params)
{
    int line = 1;
    std::cout << "Current nz, nr, nt, na: \n";
    PutNzNrNtNaToFile(stdout, run_params, &line);
    std::cout << "\n";
    InterReadNzNrNtNa(run_params);
}

void ChangeRecord(RunParams& run_params)
{
    int line = 1;
    PutScoredToFile(stdout, run_params, &line);
    std::cout << "\n";

    if (QuitOrContinue() == 'Y') {
        InterFilterRecord(run_params);
    }
}

void ChangeWth(RunParams& run_params)
{
    int line = 1;
    std::cout << "Current threshold weight: \n";
    PutWthToFile(stdout, run_params, &line);
    std::cout << "\n";
    InterReadWth(run_params);
}

void ChangeRanSeed(RunParams& run_params)
{
    int line = 1;
    std::cout << "Current random number seed: \n";
    PutRanSeedToFile(stdout, run_params, &line);
    std::cout << "\n";
    InterReadRanSeed(run_params);
}

void ChangeLayerSpecs(RunParams& run_params)
{
    int line = 1;
    std::cout << "Current layer sepcification: \n";
    PutLayerSpecsToFile(stdout, run_params, &line);
    std::cout << "\n";

    if (QuitOrContinue() == 'Y') {
        InterReadLayerSpecs(run_params);
    }
}

void ChangeNumPhotons(RunParams& run_params)
{
    int line = 1;
    std::cout << "Current value: \n";
    PutNumPhotonsToFile(stdout, run_params, &line);
    std::cout << "\n";
    InterReadNumPhotons(run_params);
}

void ChangeSourceType(RunParams& run_params)
{
    int line = 1;
    std::cout << "Current source type: \n";
    PutSourceTypeToFile(stdout, run_params, &line);
    std::cout << "\n";
    InterReadSourceType(run_params);
}

void ChangeStartP(RunParams& run_params)
{
    int line = 1;
    std::cout << "Layer Specification: \n";
    PutLayerSpecsToFile(stdout, run_params, &line);
    std::cout << "\nCurrent starting position: \n";
    PutStartPToFile(stdout, run_params, &line);
    std::cout << "\n";
    InterReadStartP(run_params);
}

/************************************************************************
 *  return 1 if string[0] = Q, quit change menu;
 *  return 2 if string[0] = X, quit to the main menu;
 *  return 0 otherwise.
 ****/
int BranchChangeMenu(std::string& string, RunParams& run_params)
{
    switch (toupper(string[0])) {
        case 'M':
            ChangeMediumList(run_params);
            break;

        case 'F':
            ChangeFnameFormat(run_params);
            break;

        case 'D':
            ChangeDzDrDt(run_params);
            break;

        case 'N':
            ChangeNzNrNtNa(run_params);
            break;

        case 'C':
            ChangeRecord(run_params);
            break;

        case 'W':
            ChangeWth(run_params);
            break;

        case 'R':
            ChangeRanSeed(run_params);
            break;

        case 'L':
            ChangeLayerSpecs(run_params);
            break;

        case 'P':
            ChangeNumPhotons(run_params);
            break;

        case 'S':
            ChangeSourceType(run_params);
            break;

        case 'Z':
            ChangeStartP(run_params);
            break;

        case 'O':
            PutInputToFile(stdout, run_params);
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
bool RunChangedInput(RunParams& run_params)
{
    std::string string;
    int branch;

    std::cout << "Any changes to the input parameters? (Y/N)";
    do { 
        fgets(string, STRLEN, stdin);
    } while (string.empty());

    while (toupper(string[0]) == 'Y') {
        do {
            do {
                std::cout << "\n> Change menu (h for help) => ";
                fgets(string, STRLEN, stdin);
            } while (string.empty());

            /* string[0] is 'X' or 'Q'. */
            if (branch = BranchChangeMenu(string, run_params)) {
                break;
            }
        } while (1);

        std::cout << "Do you want to save the input to a file? (Y/N)";
        fgets(string, STRLEN, stdin);
        if (toupper(string[0]) == 'Y') {
            std::cout << "Give the file name to save input: ( .mci): ";
            fgets(string, STRLEN, stdin);

            std::ofstream file;
            if (fopen_s(&fp, string, "w")) {
                puts("Can not open the file to write.");
            }
            else {
                PutInputToFile(file, run_params);
            }
        }

        /* quit change menu and start simulation. */
        if (branch == 1) {
            if (!CheckInputConsis(run_params)) {
                do {
                    std::cout << "Change input or exit to main menu (c/x): ";
                    fgets(string, STRLEN, stdin);
                } while (!string.empty() || toupper(string[0]) != 'X' && toupper(string[0]) != 'C');

                if (toupper(string[0]) == 'X') {
                    run_params.media.clear();
                    run_params.layers.clear()
                    return false;
                }
                else {
                    string[0] = 'Y';	/* continue to change parameters. */
                }
            }
            else {
                return true;
            }

        }
        /* exit to menu. */
        else {
            run_params.media.clear();
            run_params.layers.clear();
            return false;
        }
    }

    return true;
}

/**************************************************************************
 *	Return 1, if the name in the name list.
 *	Return 0, otherwise.
 ****/
bool NameInList(char* Name, NameLink List)
{
    while (List != NULL) {
        if (strcmp(Name, List->name) == 0) {
            return (1);
        }
        List = List->next;
    };
    return (0);
}

/**************************************************************************
 *	Add the name to the name list.
 ****/
void AddNameToList(char* Name, NameLink* List_Ptr)
{
    //NameLink list = *List_Ptr;

    /* first node. */
    if (*List_Ptr == NULL) {
        /* Allocate if list is null */
        NameLink list = (NameLink)calloc(1, sizeof(NameNode));

        *List_Ptr = list;
        strcpy_s(list->name, STRLEN, Name);
        list->next = NULL;
    }

    /* subsequent nodes. */
    else {
        NameLink list = *List_Ptr;

        /* Move to the last node. */
        while (list->next != NULL) {
            list = list->next;
        }

        /* Append a node to the list. */
        list->next = (NameLink)malloc(sizeof(NameNode));
        list = list->next;
        if (list != NULL) {
            strcpy_s(list->name, STRLEN, Name);
            list->next = NULL;
        }
    }
}

/**************************************************************************
 *	Check against duplicated input names.
 *
 *	A linked list is set up to store the input names used
 *	in this input data input.
 ****/
bool FnameTaken(char* fname, NameLink* List_Ptr)
{
    if (NameInList(fname, *List_Ptr)) {
        return (1);
    }
    else {
        AddNameToList(fname, List_Ptr);
        return (0);
    }
}

/**************************************************************************
 *	Free each node in the input name list.
 ****/
void FreeFnameList(NameLink List)
{
    NameLink next;
    while (List != NULL) {
        next = List->next;
        free(List);
        List = next;
    }
}

/*******************************************************************************
 *  Check the whether the flag end is met.
 *  The input position is restored to the current position
 *  at the end of the inquery.
 ****/
bool EndOfRunsQ(FILE** FilePP)
{
    std::string buf;

    /* found end of runs. */
    bool found = 1;

    /* record input position. */
    long file_pos = ftell(*FilePP);
    strcpy_s(buf, sizeof(buf), FindDataLine(*FilePP));
    if (buf[0] == '\0') {
        found = 0;
        std::cout << "Missing end.\n";
    }
    else if (strstr(buf, "end") == NULL) {
        found = 0;
    }

    fseek(*FilePP, file_pos, SEEK_SET);	/* restore postion. */
    return (found);
}

/**************************************************************************
 *  Check the input parameters for all runs.
 *  This function will count number of runs and assign it to
 *  run_params.num_runs.
 ****/
void CheckParamFromFile(std::fstream& file, RunParams& run_params)
{
    short i_run = 0;
    NameLink head = NULL;

    if (!ReadMediumListQ(file, run_params)) {
        exit(1);
    }

    long file_pos = ftell(Fp);
    do {
        printf("Checking input data for run %d\n", ++i_run);
        ReadRunParam(file, run_params);

        /* output files share the same input name. */
        bool name_taken = FnameTaken(run_params.output_filename, &head);
        if (name_taken) {
            printf("file name %s duplicated.\n", run_params.output_filename);
            exit(1);
        }
        free(run_params.layer);
    } while (!EndOfRunsQ(&Fp));

    run_params.num_runs = i_run;
    FreeFnameList(head);
    fseek(file, file_pos, SEEK_SET);
}

/**************************************************************************
 *	Allocate the arrays in Tracer for one run, and
 *	array elements are automatically initialized to zeros.
 *
 *	Remember that the indices for Rd_r[], Td_r[],
 *	& A_rz[][iz] start from -1 storing the collimated
 *	responses.
 ****/
void InitOutputData(RunParams& run_params, Tracer& tracer)
{
    short nz = run_params.nz;
    short nr = run_params.nr;
    short na = run_params.na;
    short nt = run_params.nt;

    /* remember to use nl+2 because of 2 for ambient. */
    short nl = run_params.layers.size();

    if (nz <= 0 || nr <= 0 || na <= 0 || nl <= 0) {
        std::cout << "Invalid grid parameters.\n";
        exit(1);
    }

    /* Init pure numbers. */
    tracer.Rsp = 0.0;
    tracer.Rb = 0.0;
    tracer.Rd = 0.0;
    tracer.Td = 0.0;
    tracer.Tb = 0.0;
    tracer.A = 0.0;

    tracer.Rbe = 0.0;
    tracer.Rde = 0.0;
    tracer.Tde = 0.0;
    tracer.Tbe = 0.0;
    tracer.Ae = 0.0;

    /* Allocate the 1D, 2D and 3D arrays. */
    tracer.Rd_rat = (run_params.record.Rd_rat) ? AllocArray3D(0, nr - 1, 0, na - 1, 0, nt - 1) : NULL;
    tracer.Rd_ra = (run_params.record.Rd_ra) ? AllocArray2D(0, nr - 1, 0, na - 1) : NULL;
    tracer.Rd_rt = (run_params.record.Rd_rt) ? AllocArray2D(0, nr - 1, 0, nt - 1) : NULL;
    tracer.Rd_at = (run_params.record.Rd_at) ? AllocArray2D(0, na - 1, 0, nt - 1) : NULL;
    tracer.Rd_r = (run_params.record.Rd_r) ? AllocArray1D(0, nr - 1) : NULL;
    tracer.Rd_a = (run_params.record.Rd_a) ? AllocArray1D(0, na - 1) : NULL;
    tracer.Rd_t = (run_params.record.Rd_t) ? AllocArray1D(0, nt - 1) : NULL;

    tracer.Td_rat = (run_params.record.Td_rat) ? AllocArray3D(0, nr - 1, 0, na - 1, 0, nt - 1) : NULL;
    tracer.Td_ra = (run_params.record.Td_ra) ? AllocArray2D(0, nr - 1, 0, na - 1) : NULL;
    tracer.Td_rt = (run_params.record.Td_rt) ? AllocArray2D(0, nr - 1, 0, nt - 1) : NULL;
    tracer.Td_at = (run_params.record.Td_at) ? AllocArray2D(0, na - 1, 0, nt - 1) : NULL;
    tracer.Td_r = (run_params.record.Td_r) ? AllocArray1D(0, nr - 1) : NULL;
    tracer.Td_a = (run_params.record.Td_a) ? AllocArray1D(0, na - 1) : NULL;
    tracer.Td_t = (run_params.record.Td_t) ? AllocArray1D(0, nt - 1) : NULL;

    tracer.A_rzt = (run_params.record.A_rzt) ? AllocArray3D(0, nr - 1, 0, nz - 1, 0, nt - 1) : NULL;
    tracer.Ab_zt = (run_params.record.A_rzt) ? AllocArray2D(0, nz - 1, 0, nt - 1) : NULL;
    tracer.A_rz = (run_params.record.A_rz) ? AllocArray2D(0, nr - 1, 0, nz - 1) : NULL;
    tracer.Ab_z = (run_params.record.A_rz) ? AllocArray1D(0, nz - 1) : NULL;
    tracer.A_zt = (run_params.record.A_zt) ? AllocArray2D(0, nz - 1, 0, nt - 1) : NULL;
    tracer.A_z = (run_params.record.A_z) ? AllocArray1D(0, nz - 1) : NULL;
    tracer.A_t = (run_params.record.A_t) ? AllocArray1D(0, nt - 1) : NULL;
}

/**************************************************************************
 *	Undo what InitOutputData did.
 *  i.e. free the data allocations.
 ****/
void FreeData(RunParams& run_params, Tracer& tracer)
{
    short nz = run_params.nz;
    short nr = run_params.nr;
    short na = run_params.na;
    short nt = run_params.nt;

    free(run_params.layer);

    FreeArray3D(tracer.Rd_rat, 0, nr - 1, 0, na - 1, 0, nt - 1);
    FreeArray2D(tracer.Rd_ra, 0, nr - 1, 0, na - 1);
    FreeArray2D(tracer.Rd_rt, 0, nr - 1, 0, nt - 1);
    FreeArray2D(tracer.Rd_at, 0, na - 1, 0, nt - 1);
    FreeArray1D(tracer.Rd_r, 0, nr - 1);
    FreeArray1D(tracer.Rd_a, 0, na - 1);
    FreeArray1D(tracer.Rd_t, 0, nt - 1);

    FreeArray3D(tracer.Td_rat, 0, nr - 1, 0, na - 1, 0, nt - 1);
    FreeArray2D(tracer.Td_ra, 0, nr - 1, 0, na - 1);
    FreeArray2D(tracer.Td_rt, 0, nr - 1, 0, nt - 1);
    FreeArray2D(tracer.Td_at, 0, na - 1, 0, nt - 1);
    FreeArray1D(tracer.Td_r, 0, nr - 1);
    FreeArray1D(tracer.Td_a, 0, na - 1);
    FreeArray1D(tracer.Td_t, 0, nt - 1);

    FreeArray3D(tracer.A_rzt, 0, nr - 1, 0, nz - 1, 0, nt - 1);
    FreeArray2D(tracer.Ab_zt, 0, nz - 1, 0, nt - 1);
    FreeArray2D(tracer.A_rz, 0, nr - 1, 0, nz - 1);
    FreeArray1D(tracer.Ab_z, 0, nz - 1);
    FreeArray2D(tracer.A_zt, 0, nz - 1, 0, nt - 1);
    FreeArray1D(tracer.A_z, 0, nz - 1);
    FreeArray1D(tracer.A_t, 0, nt - 1);
}

/**************************************************************************
 *	Scale Rd and Td properly.
 *	"a" stands for angle alpha.
 ****
 *	Scale Rd(r,a) and Td(r,a) by
 *      (area perpendicular to photon direction)
 *		x(solid angle)x(No. of photons).
 *	or
 *		[2*PI*r*dr*cos(a)]x[2*PI*sin(a)*da]x[No. of photons]
 *	or
 *		[2*PI*PI*dr*da*r*sin(2a)]x[No. of photons]
 ****
 *	Scale Rd(r) and Td(r) by
 *		(area on the surface)x(No. of photons).
 ****
 *	Scale Rd(a) and Td(a) by
 *		(solid angle) cos(a) x(No. of photons).
 ****
 *  Mode = 0, scale Rd and Td; Mode = 1, unscale Rd and Td.
 ****/
void ScaleRdTd(RunParams& run_params, Tracer& tracer, char Mode)
{
    short nr = run_params.nr;
    short na = run_params.na;
    short nt = run_params.nt;
    double dr = run_params.dr;
    double da = run_params.da;
    double dt = run_params.dt;

    double scale1 = (double)run_params.num_photons;
    if (Mode == 0) {
        tracer.Rde = 1 / scale1 * sqrt(tracer.Rde - tracer.Rd * tracer.Rd / scale1);
        tracer.Tde = 1 / scale1 * sqrt(tracer.Tde - tracer.Td * tracer.Td / scale1);
        tracer.Rbe = 1 / scale1 * sqrt(tracer.Rbe - tracer.Rb * tracer.Rb / scale1);
        tracer.Tbe = 1 / scale1 * sqrt(tracer.Tbe - tracer.Tb * tracer.Tb / scale1);

        tracer.Rd /= scale1;
        tracer.Td /= scale1;
        tracer.Rb = tracer.Rb / scale1 + tracer.Rsp;
        tracer.Tb /= scale1;
    }
    else {
        tracer.Rd *= scale1;
        tracer.Td *= scale1;
        tracer.Rb = (tracer.Rb - tracer.Rsp) * scale1;
        tracer.Tb *= scale1;

        tracer.Rde = (scale1 * tracer.Rde) * (scale1 * tracer.Rde) + 1 / scale1 * tracer.Rd * tracer.Rd;
        tracer.Tde = (scale1 * tracer.Tde) * (scale1 * tracer.Tde) + 1 / scale1 * tracer.Td * tracer.Td;
        tracer.Rbe = (scale1 * tracer.Rbe) * (scale1 * tracer.Rbe) + 1 / scale1 * tracer.Rb * tracer.Rb;
        tracer.Tbe = (scale1 * tracer.Tbe) * (scale1 * tracer.Tbe) + 1 / scale1 * tracer.Tb * tracer.Tb;
    }

    scale1 = dt * run_params.num_photons;
    if (run_params.record.Rd_t) {
        for (short it = 0; it < nt; it++) {
            /* scale Rd_t. */
            if (Mode == 0) {
                tracer.Rd_t[it] /= scale1;
            }
            /* unscale Rd_t. */
            else {
                tracer.Rd_t[it] *= scale1;
            }
        }
    }

    if (run_params.record.Td_t) {
        for (short it = 0; it < nt; it++) {
            /* scale Td_t. */
            if (Mode == 0) {
                tracer.Td_t[it] /= scale1;
            }
            /* unscale Rd_t. */
            else {
                tracer.Td_t[it] *= scale1;
            }
        }
    }

    scale1 = 2.0 * PI * dr * dr * run_params.num_photons;
    /* area is 2*PI*[(ir+0.5)*dr]*dr.  ir + 0.5 to be added. */

    if (run_params.record.Rd_r) {
        for (short ir = 0; ir < nr; ir++) {
            double scale2 = 1.0 / ((ir + 0.5) * scale1);
            /* scale Rd_r. */
            if (Mode == 0) {
                tracer.Rd_r[ir] *= scale2;
            }
            /* unscale Rd_r. */
            else {
                tracer.Rd_r[ir] /= scale2;
            }
        }
    }

    if (run_params.record.Td_r) {
        for (short ir = 0; ir < nr; ir++) {
            double scale2 = 1.0 / ((ir + 0.5) * scale1);
            /* scale Td_r. */
            if (Mode == 0) {
                tracer.Td_r[ir] *= scale2;
            }
            /* unscale Td_r. */
            else {
                tracer.Td_r[ir] /= scale2;
            }
        }
    }

    scale1 *= dt;
    if (run_params.record.Rd_rt) {
        for (short ir = 0; ir < nr; ir++) {
            for (short it = 0; it < nt; it++) {
                double scale2 = 1.0 / ((ir + 0.5) * scale1);
                /* scale Rd_rt. */
                if (Mode == 0) {
                    tracer.Rd_rt[ir][it] *= scale2;
                }
                /* unscale Rd_rt. */
                else {
                    tracer.Rd_rt[ir][it] *= scale2;
                }
            }
        }
    }

    if (run_params.record.Td_rt) {
        for (short ir = 0; ir < nr; ir++) {
            for (short it = 0; it < nt; it++) {
                double scale2 = 1.0 / ((ir + 0.5) * scale1);
                /* scale Td_rt. */
                if (Mode == 0) {
                    tracer.Td_rt[ir][it] *= scale2;
                }
                /* unscale Td_rt. */
                else {
                    tracer.Td_rt[ir][it] /= scale2;
                }
            }
        }
    }

    scale1 = PI * da * run_params.num_photons;
    /* solid angle times cos(a) is PI*sin(2a)*da. sin(2a) to be added. */

    if (run_params.record.Rd_a) {
        for (short ia = 0; ia < na; ia++) {
            double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
            /* scale Rd_a. */
            if (Mode == 0) {
                tracer.Rd_a[ia] *= scale2;
            }
            /* unscale Rd_a. */
            else {
                tracer.Rd_a[ia] /= scale2;
            }
        }
    }

    if (run_params.record.Td_a) {
        for (short ia = 0; ia < na; ia++) {
            double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
            /* scale Td_a. */
            if (Mode == 0) {
                tracer.Td_a[ia] *= scale2;
            }
            /* unscale Td_a. */
            else {
                tracer.Td_a[ia] /= scale2;
            }
        }
    }

    scale1 *= dt;
    if (run_params.record.Rd_at) {
        for (short ia = 0; ia < na; ia++) {
            for (short it = 0; it < nt; it++) {
                double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
                /* scale Rd_at. */
                if (Mode == 0) {
                    tracer.Rd_at[ia][it] *= scale2;
                }
                /* unscale Rd_at. */
                else {
                    tracer.Rd_at[ia][it] /= scale2;
                }
            }
        }
    }

    if (run_params.record.Td_at) {
        for (short ia = 0; ia < na; ia++) {
            for (short it = 0; it < nt; it++) {
                double scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
                /* scale Td_at. */
                if (Mode == 0) {
                    tracer.Td_at[ia][it] *= scale2;
                }
                /* unscale Td_at. */
                else {
                    tracer.Td_at[ia][it] /= scale2;
                }
            }
        }
    }

    scale1 = 2.0 * PI * dr * dr * PI * da * run_params.num_photons;
    if (run_params.record.Rd_ra) {
        for (short ir = 0; ir < nr; ir++) {
            for (short ia = 0; ia < na; ia++) {
                double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                /* scale Rd_ra. */
                if (Mode == 0) {
                    tracer.Rd_ra[ir][ia] *= scale2;
                }
                /* unscale Rd_ra. */
                else {
                    tracer.Rd_ra[ir][ia] /= scale2;
                }
            }
        }
    }

    if (run_params.record.Td_ra) {
        for (short ir = 0; ir < nr; ir++) {
            for (short ia = 0; ia < na; ia++) {
                double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                /* scale Td_ra. */
                if (Mode == 0) {
                    tracer.Td_ra[ir][ia] *= scale2;
                }
                /* unscale Td_ra. */
                else {
                    tracer.Td_ra[ir][ia] /= scale2;
                }
            }
        }
    }

    scale1 *= dt;
    if (run_params.record.Rd_rat) {
        for (short ir = 0; ir < nr; ir++) {
            for (short ia = 0; ia < na; ia++) {
                for (short it = 0; it < nt; it++) {
                    double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                    /* scale Rd_rat. */
                    if (Mode == 0) {
                        tracer.Rd_rat[ir][ia][it] *= scale2;
                    }
                    /* unscale Rd_rat. */
                    else {
                        tracer.Rd_rat[ir][ia][it] /= scale2;
                    }
                }
            }
        }
    }

    if (run_params.record.Td_rat) {
        for (short ir = 0; ir < nr; ir++) {
            for (short ia = 0; ia < na; ia++) {
                for (short it = 0; it < nt; it++) {
                    double scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
                    /* scale Td_rat. */
                    if (Mode == 0) {
                        tracer.Td_rat[ir][ia][it] *= scale2;
                    }
                    /* unscale Td_rat. */
                    else {
                        tracer.Td_rat[ir][ia][it] /= scale2;
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
void ScaleA(RunParams& run_params, Tracer& tracer, char Mode)
{
    short nz = run_params.nz;
    short nr = run_params.nr;
    short nt = run_params.nt;
    double dz = run_params.dz;
    double dr = run_params.dr;
    double dt = run_params.dt;
    double scale1 = (double)run_params.num_photons;

    /* scale A. */
    if (Mode == 0) {
        tracer.Ae = 1 / scale1 * sqrt(tracer.Ae - tracer.A * tracer.A / scale1);
        tracer.A /= scale1;
    }
    /* unscale A. */
    else {
        tracer.A *= scale1;
        tracer.Ae = (scale1 * tracer.Ae) * (scale1 * tracer.Ae) + 1 / scale1 * tracer.A * tracer.A;
    }

    double scale2 = scale1 * dt;
    if (run_params.record.A_t) {
        for (short it = 0; it < nt; it++) {
            /* scale A_t. */
            if (Mode == 0) {
                tracer.A_t[it] /= scale2;
            }
            /* unscale A_t. */
            else {
                tracer.A_t[it] *= scale2;
            }
        }
    }

    scale1 *= dz;
    if (run_params.record.A_z) {
        for (short iz = 0; iz < nz; iz++) {
            /* scale A_z. */
            if (Mode == 0) {
                tracer.A_z[iz] /= scale1;
            }
            /* unscale A_z. */
            else {
                tracer.A_z[iz] *= scale1;
            }
        }
    }

    scale2 = scale1 * dt;
    if (run_params.record.A_zt) {
        for (short iz = 0; iz < nz; iz++) {
            for (short it = 0; it < nt; it++) {
                /* scale A_zt. */
                if (Mode == 0) {
                    tracer.A_zt[iz][it] /= scale2;
                }
                /* unscale A_zt. */
                else {
                    tracer.A_zt[iz][it] *= scale2;
                }
            }
        }
    }

    if (run_params.record.A_rz) {
        for (short iz = 0; iz < nz; iz++) {
            /* scale Ab_z. */
            if (Mode == 0) {
                tracer.Ab_z[iz] /= scale1;
            }
            /* unscale Ab_z. */
            else {
                tracer.Ab_z[iz] *= scale1;
            }
        }
    }

    if (run_params.record.A_rzt) {
        for (short iz = 0; iz < nz; iz++) {
            for (short it = 0; it < nt; it++) {
                /* scale Ab_zt. */
                if (Mode == 0) {
                    tracer.Ab_zt[iz][it] /= scale2;
                }
                /* unscale Ab_zt. */
                else {
                    tracer.Ab_zt[iz][it] *= scale2;
                }
            }
        }
    }

    scale1 = 2.0 * PI * dr * dr * dz * run_params.num_photons;
    if (run_params.record.A_rz) {
        for (short ir = 0; ir < nr; ir++) {
            for (short iz = 0; iz < nz; iz++) {
                /* scale A_rz. */
                if (Mode == 0) {
                    tracer.A_rz[ir][iz] /= (ir + 0.5) * scale1;
                }
                /* unscale A_rz. */
                else {
                    tracer.A_rz[ir][iz] *= (ir + 0.5) * scale1;
                }
            }
        }
    }

    scale2 = scale1 * dt;
    if (run_params.record.A_rzt) {
        for (short ir = 0; ir < nr; ir++) {
            for (short iz = 0; iz < nz; iz++) {
                for (short it = 0; it < nt; it++) {
                    /* scale A_rzt. */
                    if (Mode == 0) {
                        tracer.A_rzt[ir][iz][it] /= (ir + 0.5) * scale2;
                    }
                    /* unscale A_rzt. */
                    else {
                        tracer.A_rzt[ir][iz][it] *= (ir + 0.5) * scale2;
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
void ScaleResult(RunParams& run_params, Tracer& tracer, char Mode)
{
    ScaleRdTd(run_params, tracer, Mode);
    ScaleA(run_params, tracer, Mode);
}

/**************************************************************************
 *	Write the version number as the first string in the
 *	input.
 *	Use chars only so that they can be read as either
 *	ASCII or binary.
 ****/
void WriteVersion(std::fstream& file, const std::string& version)
{
    file << version << " \t# Version number of the file format.\n\n";
    file << "####\n# Data categories include: \n";
    file << "# InParam, RAT, \n";
    file << "# Rd_r\tRd_a\tRd_ra\tRd_t\tRd_rt\tRd_at\tRd_rat\n";
    file << "# Td_r\tTd_a\tTd_ra\tTd_t\tTd_rt\tTd_at\tTd_rat\n";
    file << "# A_z\tA_rz\tA_t\tA_zt\tA_rzt\n";
    file << "####\n\n";
}

/***************************************************************************
 * Save the status of the random number generater to output input.
 ****/
void SaveRandomStatus(std::ofstream file)
{
    /* get the status. */
    long status[57];
    RandomGen();
    file << std::format("# status of the random number generator:");

    for (int i = 0; i < 57; i++) {
        if (i % 5) {
            file << std::format("%14ld", status[i]);
        }
        else {
            file << std::format("\n%14ld ", status[i]);
        }
    }

    file << std::format("\n\n");
}

/***************************************************************************
 * Read and restore the status of random number generater from previous
 * output input.
 ****/
void RestoreRandomStatus(std::ifstream file)
{
    std::string buf;
    long status[57];

    do {
        fgets(buf, sizeof(buf), Fp);
    } while (buf[0] != '#');

    for (int i = 0; i < 57; i++) {
        fscanf_s(file, "%ld", &status[i]);
    }

    /* restore the status. */
    RandomGen(3, 0, status);
}

/**************************************************************************
 *	Write reflectance, absorption, transmission.
 ****/
void WriteRAT(std::ofstream file, Tracer& tracer)
{
    file << std::format("RAT #Reflectance, Absorption, Transmittance.\n");
    file << std::format("# Average \tStandard Err \tRel Err\n");
    file << std::format("%-14.6G \t\t\t\t#Rsp: Specular reflectance.\n", tracer.R.s);
    file << std::format("%-14.6G \t%-14.6G %6.2f%%\t#Rb: Ballistic reflectance.\n",   tracer.R.b, tracer.R.se, (tracer.R.b) ? tracer.R.se / tracer.R.b * 100 : 0);
    file << std::format("%-14.6G \t%-14.6G %6.2f%%\t#Rd: Diffuse reflectance.\n",     tracer.R.d, tracer.R.de, (tracer.R.d) ? tracer.R.de / tracer.R.d * 100 : 0);
    file << std::format("%-14.6G \t%-14.6G %6.2f%%\t#A:  Absorbed fraction.\n",       tracer.A.a, tracer.A.e,  (tracer.A.a) ? tracer.A.e  / tracer.A.a * 100 : 0);
    file << std::format("%-14.6G \t%-14.6G %6.2f%%\t#Tb: Ballistic transmittance.\n", tracer.T.b, tracer.T.be, (tracer.T.b) ? tracer.T.be / tracer.T.b * 100 : 0);
    file << std::format("%-14.6G \t%-14.6G %6.2f%%\t#Td: Diffuse transmittance.\n",   tracer.T.d, tracer.T.de, (tracer.T.d) ? tracer.T.de / tracer.T.d * 100 : 0);
    file << std::format("\n");
}

/**************************************************************************
 *	Read reflectance, absorption, transmission.
 ****/
void ReadRAT(std::ifstream& file, Tracer& tracer)
{
    /* skip RAT line. */
    std::string buf = FindDataLine(file);

    strcpy_s(buf, sizeof(buf), FindDataLine(file));
    sscanf_s(buf, "%lf", &(tracer.R.s));

    strcpy_s(buf, sizeof(buf), FindDataLine(file));
    sscanf_s(buf, "%lf %lf", &(tracer.R.b), &(tracer.R.be));

    strcpy_s(buf, sizeof(buf), FindDataLine(file));
    sscanf_s(buf, "%lf %lf", &(tracer.R.d), &(tracer.R.de));

    strcpy_s(buf, sizeof(buf), FindDataLine(file));
    sscanf_s(buf, "%lf %lf", &(tracer.A.a), &(tracer.A.e));

    strcpy_s(buf, sizeof(buf), FindDataLine(file));
    sscanf_s(buf, "%lf %lf", &(tracer.T.b), &(tracer.T.be));

    strcpy_s(buf, sizeof(buf), FindDataLine(file));
    sscanf_s(buf, "%lf %lf", &(tracer.T.d), &(tracer.T.de));
}

/**************************************************************************
 *  5 numbers each line.
 *  Mode = 0, read; Mode = 1, write.
 ****/
void IOAb_zt(std::ifstream file, short Nz, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# Ab[z][t]. [1/(cm ps)]",
                "# Ab[0][0], [0][1],..[0][nt-1]",
                "# Ab[1][0], [1][1],..[1][nt-1]",
                "# ...",
                "# Ab[nz-1][0], [nz-1][1],..[nz-1][nt-1]",
                "Ab_zt");
    }
    else {
        /* skip A_z line. */
        FindDataLine(Fp);
    }

    short i = 0;
    for (short iz = 0; iz < Nz; iz++) {
        for (short it = 0; it < Nt; it++) {
            if (Mode == 1) {
                file << std::format("%12.4E ", tracer.Ab_zt[iz][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                fscanf_s(file, "%lf", &(tracer.Ab_zt[iz][it]));
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
void IOA_rzt(std::ifstream file, short Nr, short Nz, short Nt, Tracer& tracer, char Mode)
{
    IOAb_zt(file, Nz, Nt, tracer, Mode);

    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# A[r][z][t]. [1/(cm3 ps)]",
                "# A[0][0][0], [0][0][1],..[0][0][nt-1]",
                "# A[0][1][0], [0][1][1],..[0][1][nt-1]",
                "# ...",
                "# A[nr-1][nz-1][0], [nr-1][nz-1][1],..[nr-1][nz-1][nt-1]",
                "A_rzt");
    }
    else {
        FindDataLine(Fp);
    }

    short i = 0;
    for (short ir = 0; ir < Nr; ir++) {
        for (short iz = 0; iz < Nz; iz++) {
            for (short it = 0; it < Nt; it++) {
                if (Mode == 1) {
                    file << std::format("%12.4E ", tracer.A_rzt[ir][iz][it]);
                    if (++i % 5 == 0) {
                        file << std::format("\n");
                    }
                }
                else {
                    fscanf_s(file, "%lf", &(tracer.A_rzt[ir][iz][it]));
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
void IOAb_z(std::ifstream file, short Nz, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        file << std::format("Ab_z #Ab[0], [1],..Ab[nz-1]. [1/cm]\n");	/* flag. */
    }
    else {
        FindDataLine(Fp);
    }

    for (short iz = 0; iz < Nz; iz++) {
        if (Mode == 1) {
            file << std::format("%12.4E\n", tracer.Ab_z[iz]);
        }
        else {
            fscanf_s(file, "%lf", &(tracer.Ab_z[iz]));
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
void IOA_rz(std::ifstream file, short Nr, short Nz, Tracer& tracer, char Mode)
{
    IOAb_z(file, Nz, tracer, Mode);

    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n",
                "# A[r][z]. [1/cm3]",
                "# A[0][0], [0][1],..[0][nz-1]",
                "# ...",
                "# A[nr-1][0], [nr-1][1],..[nr-1][nz-1]",
                "A_rz");
    }
    else {
        /* skip A_rz line. */
        FindDataLine(Fp);
    }

    short i = 0;
    for (short ir = 0; ir < Nr; ir++) {
        for (short iz = 0; iz < Nz; iz++) {
            if (Mode == 1) {
                file << std::format("%12.4E ", tracer.A_rz[ir][iz]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                fscanf_s(file, "%lf", &(tracer.A_rz[ir][iz]));
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
void IOA_zt(std::ifstream file, short Nz, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# A[z][t]. [1/(cm ps)]",
                "# A[0][0], [0][1],..[0][nt-1]",
                "# A[1][0], [1][1],..[1][nt-1]",
                "# ...",
                "# A[nz-1][0], [nz-1][1],..[nz-1][nt-1]",
                "A_zt");
    }
    else {
        /* skip A_zt line. */
        FindDataLine(Fp);
    }

    short i = 0;
    for (short iz = 0; iz < Nz; iz++) {
        for (short it = 0; it < Nt; it++) {
            if (Mode == 1) {
                file << std::format("%12.4E ", tracer.A_zt[iz][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                fscanf_s(file, "%lf", &(tracer.A_zt[iz][it]));
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
void IOA_z(std::ifstream file, short Nz, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        file << std::format("A_z #A[0], [1],..A[nz-1]. [1/cm]\n");
    }
    else {
        /* skip A_z line. */
        FindDataLine(Fp);
    }

    for (short iz = 0; iz < Nz; iz++) {
        if (Mode == 1) {
            file << std::format("%12.4E\n", tracer.A_z[iz]);
        }
        else {
            fscanf_s(file, "%lf", &(tracer.A_z[iz]));
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
void IOA_t(std::ifstream file, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        file << std::format("A_t #A[0], [1],..A[nt-1]. [1/ps]\n");
    }
    else {
        FindDataLine(Fp);
    }

    for (short it = 0; it < Nt; it++) {
        if (Mode == 1) {
            file << std::format("%12.4E\n", tracer.A_t[it]);
        }
        else {
            fscanf_s(file, "%lf", &(tracer.A_t[it]));
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
void IORd_rat(std::ifstream file, short Nr, short Na, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# Rd[r][a][t]. [1/(cm2 sr ps)]",
                "# Rd[0][0][0], [0][0][1],..[0][0][nt-1]",
                "# Rd[0][1][0], [0][1][1],..[0][1][nt-1]",
                "# ...",
                "# Rd[nr-1][na-1][0], [nr-1][na-1][1],..[nr-1][na-1][nt-1]",
                "Rd_rat");
    }
    else {
        FindDataLine(Fp);
    }

    short i = 0;
    for (short ir = 0; ir < Nr; ir++) {
        for (short ia = 0; ia < Na; ia++) {
            for (short it = 0; it < Nt; it++) {
                if (Mode == 1) {
                    file << std::format("%12.4E ", tracer.Rd_rat[ir][ia][it]);
                    if (++i % 5 == 0) {
                        file << std::format("\n");
                    }
                }
                else {
                    fscanf_s(file, "%lf", &(tracer.Rd_rat[ir][ia][it]));
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
void IORd_ra(std::ifstream file, short Nr, short Na, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# Rd[r][angle]. [1/(cm2 sr)].",
                "# Rd[0][0], [0][1],..[0][na-1]",
                "# Rd[1][0], [1][1],..[1][na-1]",
                "# ...",
                "# Rd[nr-1][0], [nr-1][1],..[nr-1][na-1]",
                "Rd_ra");
    }
    else {
        FindDataLine(Fp);
    }

    for (short ir = 0; ir < Nr; ir++) {
        for (short ia = 0; ia < Na; ia++) {
            if (Mode == 1) {
                file << std::format("%12.4E ", tracer.Rd_ra[ir][ia]);
                if ((ir * Na + ia + 1) % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                fscanf_s(file, "%lf", &(tracer.Rd_ra[ir][ia]));
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
void IORd_rt(std::ifstream file, short Nr, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# Rd[r][t]. [1/(cm2 ps)]",
                "# Rd[0][0], [0][1],..[0][nt-1]",
                "# Rd[0][0], [0][1],..[0][nt-1]",
                "# ...",
                "# Rd[nr-1][0], [nr-1][1],..[nr-1][nt-1]",
                "Rd_rt");
    }
    else {
        FindDataLine(Fp);
    }

    short i = 0; 
    for (short ir = 0; ir < Nr; ir++) {
        for (short it = 0; it < Nt; it++) {
            if (Mode == 1) {
                file << std::format("%12.4E ", tracer.Rd_rt[ir][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                fscanf_s(file, "%lf", &(tracer.Rd_rt[ir][it]));
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
void IORd_at(std::ifstream file, short Na, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# Rd[a][t]. [1/(sr ps)]",
                "# Rd[0][0], [0][1],..[0][nt-1]",
                "# Rd[1][0], [1][1],..[1][nt-1]",
                "# ...",
                "# Rd[na-1][0], [na-1][1],..[na-1][nt-1]",
                "Rd_at");
    }
    else {
        FindDataLine(Fp);
    }

    short i = 0;
    for (short ia = 0; ia < Na; ia++) {
        for (short it = 0; it < Nt; it++) {
            if (Mode == 1) {
                file << std::format("%12.4E ", tracer.Rd_at[ia][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                fscanf_s(file, "%lf", &(tracer.Rd_at[ia][it]));
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
void IORd_r(std::ifstream file, short Nr, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        file << std::format("Rd_r #Rd[0], [1],..Rd[nr-1]. [1/cm2]\n");
    }
    else {
        FindDataLine(Fp);
    }

    for (short ir = 0; ir < Nr; ir++) {
        if (Mode == 1) {
            file << std::format("%12.4E\n", tracer.Rd_r[ir]);
        }
        else {
            fscanf_s(file, "%lf", &(tracer.Rd_r[ir]));
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
void IORd_a(std::ifstream file, short Na, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        file << std::format("Rd_a #Rd[0], [1],..Rd[na-1]. [1/sr]\n");
    }
    else {
        FindDataLine(Fp);
    }

    for (short ia = 0; ia < Na; ia++) {
        if (Mode == 1) {
            file << std::format("%12.4E\n", tracer.Rd_a[ia]);
        }
        else {
            fscanf_s(file, "%lf", &(tracer.Rd_a[ia]));
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
void IORd_t(std::ifstream file, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        file << std::format("Rd_t #Rd[0], [1],..Rd[nt-1]. [1/ps]\n");
    }
    else {
        FindDataLine(Fp);
    }

    for (short it = 0; it < Nt; it++) {
        if (Mode == 1) {
            file << std::format("%12.4E\n", tracer.Rd_t[it]);
        }
        else {
            fscanf_s(file, "%lf", &(tracer.Rd_t[it]));
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
void IOTd_rat(std::ifstream file, short Nr, short Na, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# Td[r][a][t]. [1/(cm2 sr ps)]",
                "# Td[0][0][0], [0][0][1],..[0][0][nt-1]",
                "# Td[0][1][0], [0][1][1],..[0][1][nt-1]",
                "# ...",
                "# Td[nr-1][na-1][0], [nr-1][na-1][1],..[nr-1][na-1][nt-1]",
                "Td_rat");
    }
    else {
        FindDataLine(Fp);
    }

    short i = 0;
    for (short ir = 0; ir < Nr; ir++) {
        for (short ia = 0; ia < Na; ia++) {
            for (short it = 0; it < Nt; it++) {
                if (Mode == 1) {
                    file << std::format("%12.4E ", tracer.Td_rat[ir][ia][it]);
                    if (++i % 5 == 0) {
                        file << std::format("\n");
                    }
                }
                else {
                    fscanf_s(file, "%lf", &(tracer.Td_rat[ir][ia][it]));
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
void IOTd_ra(std::ifstream file, short Nr, short Na, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# Td[r][angle]. [1/(cm2 sr)].",
                "# Td[0][0], [0][1],..[0][na-1]",
                "# Td[1][0], [1][1],..[1][na-1]",
                "# ...",
                "# Td[nr-1][0], [nr-1][1],..[nr-1][na-1]",
                "Td_ra");
    }
    else {
        FindDataLine(Fp);
    }

    for (short ir = 0; ir < Nr; ir++) {
        for (short ia = 0; ia < Na; ia++) {
            if (Mode == 1) {
                file << std::format("%12.4E ", tracer.Td_ra[ir][ia]);
                if ((ir * Na + ia + 1) % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                fscanf_s(file, "%lf", &(tracer.Td_ra[ir][ia]));
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
void IOTd_rt(std::ifstream file, short Nr, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# Td[r][t]. [1/(cm2 ps)]",
                "# Td[0][0], [0][1],..[0][nt-1]",
                "# Td[0][0], [0][1],..[0][nt-1]",
                "# ...",
                "# Td[nr-1][0], [nr-1][1],..[nr-1][nt-1]",
                "Td_rt");
    }
    else {
        FindDataLine(Fp);
    }

    short i = 0;
    for (short ir = 0; ir < Nr; ir++) {
        for (short it = 0; it < Nt; it++) {
            if (Mode == 1) {
                file << std::format("%12.4E ", tracer.Td_rt[ir][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                fscanf_s(file, "%lf", &(tracer.Td_rt[ir][it]));
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
void IOTd_at(std::ifstream file, short Na, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        fprintf(file, "%s\n%s\n%s\n%s\n%s\n%s\n",
                "# Td[a][t]. [1/(sr ps)]",
                "# Td[0][0], [0][1],..[0][nt-1]",
                "# Td[1][0], [1][1],..[1][nt-1]",
                "# ...",
                "# Td[na-1][0], [na-1][1],..[na-1][nt-1]",
                "Td_at");
    }
    else {
        FindDataLine(Fp);
    }

    short i = 0;
    for (short ia = 0; ia < Na; ia++) {
        for (short it = 0; it < Nt; it++) {
            if (Mode == 1) {
                file << std::format("%12.4E ", tracer.Td_at[ia][it]);
                if (++i % 5 == 0) {
                    file << std::format("\n");
                }
            }
            else {
                fscanf_s(file, "%lf", &(tracer.Td_at[ia][it]));
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
void IOTd_r(std::ifstream file, short Nr, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        file << std::format("Td_r #Td[0], [1],..Td[nr-1]. [1/cm2]\n");
    }
    else {
        FindDataLine(Fp);
    }

    for (short ir = 0; ir < Nr; ir++) {
        if (Mode == 1) {
            file << std::format("%12.4E\n", tracer.Td_r[ir]);
        }
        else {
            fscanf_s(file, "%lf", &(tracer.Td_r[ir]));
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
void IOTd_a(std::ifstream file, short Na, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        file << std::format("Td_a #Td[0], [1],..Td[na-1]. [1/sr]\n");
    }
    else {
        FindDataLine(Fp);
    }

    for (short ia = 0; ia < Na; ia++) {
        if (Mode == 1) {
            file << std::format("%12.4E\n", tracer.Td_a[ia]);
        }
        else {
            fscanf_s(file, "%lf", &(tracer.Td_a[ia]));
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
void IOTd_t(std::ifstream file, short Nt, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        /* flag. */
        file << std::format("Td_t #Rd[0], [1],..Td[nt-1]. [1/ps]\n");
    }
    else {
        FindDataLine(Fp);
    }

    for (short it = 0; it < Nt; it++) {
        if (Mode == 1) {
            file << std::format("%12.4E\n", tracer.Td_t[it]);
        }
        else {
            fscanf_s(file, "%lf", &(tracer.Td_t[it]));
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
void IOResult(std::fstream& file, RunParams& run_params, Tracer& tracer, char Mode)
{
    if (Mode == 1) {
        if (run_params.output_file_format == FileFormat::ASCII) {
            WriteVersion(file, "mcmloA2.0");
        }
        else {
            WriteVersion(file, "mcmloB2.0");
        }

        PutInputToFile(file, run_params);
        SaveRandomStatus(file);
        WriteRAT(file, tracer);
    }
    else {
        RestoreRandomStatus(file);
        ReadRAT(file, tracer);
    }

    /* reflectance, absorption, transmittance. */
    if (run_params.record.A_rzt) {
        IOA_rzt(file, run_params.nr, run_params.nz, run_params.nt, tracer, Mode);
    }
    if (run_params.record.A_rz) {
        IOA_rz(file, run_params.nr, run_params.nz, tracer, Mode);
    }
    if (run_params.record.A_zt) {
        IOA_zt(file, run_params.nz, run_params.nt, tracer, Mode);
    }
    if (run_params.record.A_z) {
        IOA_z(file, run_params.nz, tracer, Mode);
    }
    if (run_params.record.A_t) {
        IOA_t(file, run_params.nt, tracer, Mode);
    }

    if (run_params.record.Rd_rat) {
        IORd_rat(file, run_params.nr, run_params.na, run_params.nt, tracer, Mode);
    }
    if (run_params.record.Rd_ra) {
        IORd_ra(file, run_params.nr, run_params.na, tracer, Mode);
    }
    if (run_params.record.Rd_rt) {
        IORd_rt(file, run_params.nr, run_params.nt, tracer, Mode);
    }
    if (run_params.record.Rd_at) {
        IORd_at(file, run_params.na, run_params.nt, tracer, Mode);
    }
    if (run_params.record.Rd_r) {
        IORd_r(file, run_params.nr, tracer, Mode);
    }
    if (run_params.record.Rd_a) {
        IORd_a(file, run_params.na, tracer, Mode);
    }
    if (run_params.record.Rd_t) {
        IORd_t(file, run_params.nt, tracer, Mode);
    }

    if (run_params.record.Td_rat) {
        IOTd_rat(file, run_params.nr, run_params.na, run_params.nt, tracer, Mode);
    }
    if (run_params.record.Td_ra) {
        IOTd_ra(file, run_params.nr, run_params.na, tracer, Mode);
    }
    if (run_params.record.Td_rt) {
        IOTd_rt(file, run_params.nr, run_params.nt, tracer, Mode);
    }
    if (run_params.record.Td_at) {
        IOTd_at(file, run_params.na, run_params.nt, tracer, Mode);
    }
    if (run_params.record.Td_r) {
        IOTd_r(file, run_params.nr, tracer, Mode);
    }
    if (run_params.record.Td_a) {
        IOTd_a(file, run_params.na, tracer, Mode);
    }
    if (run_params.record.Td_t) {
        IOTd_t(file, run_params.nt, tracer, Mode);
    }

    fclose(file);
}

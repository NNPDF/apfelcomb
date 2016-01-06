/*
 *  fk_utils.h
 *  Utilities for applgrid combination
 * *  n.p.hartland@ed.ac.uk 03/12
 */

#pragma once

#include <iostream>
#include <iterator>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>


#ifndef RESULTS_PATH
#define RESULTS_PATH run/
#endif

#ifndef DATA_PATH
#define DATA_PATH data/
#endif

#ifndef APPL_PATH
#define APPL_PATH ./
#endif

#ifndef DB_PATH
#define DB_PATH ./
#endif

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

using namespace std;

namespace appl{
  class grid;
}

namespace Colour {
    enum Code {
        FG_RED      = 31,
        FG_GREEN    = 32,
        FG_YELLOW   = 33,
        FG_BLUE     = 34,
        FG_DEFAULT  = 39,
        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_YELLOW   = 43,
        BG_BLUE     = 44,
        BG_DEFAULT  = 49
    };
    
    std::ostream& operator<<(std::ostream& os, Code code);
}

const double PI=3.1415926535897932384626433832795;

// Display helpers
void DisplayHR();
void Splash();

std::string applPath();
std::string dataPath();
std::string resultsPath();
std::string databasePath();

vector<double> dsplit(string input);
vector<string> ssplit(string input);

// Directory helper functions
void setupDir(int const& theoryID, std::string const& setname);
std::string getOutputFilename(int const& theoryID, std::string const& setname, std::string const& gridname);

// Tostring function
template < typename T > std::string to_string( const T& n )
{
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
}

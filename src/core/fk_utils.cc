/*
 *  appl_utils.cc
 *  Utilities for applgrid combination
 * *  n.p.hartland@ed.ac.uk 03/12
 */


#include "apfelcomb/fk_utils.h"
#include "NNPDF/common.h"

#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"

#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>


using namespace std;

// http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c-using-posix
std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != NULL)
            result += buffer.data();
    }
    return result;
}

std::string applPath()
{
  std::string applDir(STR(APPL_PATH));
  return applDir;
};

std::string dataPath()
{
  std::string dataDir(STR(DATA_PATH));
  return dataDir;
};

std::string resultsPath()
{
  std::string resultsDir(STR(RESULTS_PATH));
  return resultsDir;
};

std::string databasePath()
{
  std::string dbDir(STR(DB_PATH));
  return dbDir;
};

std::string applCommit()
{
  std::stringstream cmdstream;
  cmdstream << "git --git-dir " << applPath() <<".git rev-parse HEAD";
  std::string retval = exec(cmdstream.str().c_str());
  retval.erase(std::remove(retval.begin(), retval.end(), '\n'), retval.end());
  return retval;
}

// Colour stream
namespace Colour{
      std::ostream& operator<<(std::ostream& os, Code code) {
        return os << "\033[" << static_cast<int>(code) << "m";
    }
}

// ************* SPLITTERS *************************

vector<string> ssplit(string in)
{
	stringstream sstr(in);
	istream_iterator<string> itr(sstr);
	istream_iterator<string> end;
	vector<string> results(itr, end);
	return results;
}

vector<double> dsplit(string in)
{
	stringstream sstr(in);
	istream_iterator<double> itr(sstr);
	istream_iterator<double> end;
	vector<double> results(itr, end);
	return results;
}

// ************ Theory Dir **************************

void setupDir(int const& theoryID, std::string const& setname, vector<std::string> const& reqgrids)
{
  // Setup required directories
  stringstream theoryDir;
  theoryDir << resultsPath()<<"theory_" << theoryID<<"/";

  mkdir(theoryDir.str().c_str(),0777);
  mkdir((theoryDir.str() + "apfelcomb/").c_str(),0777);
  mkdir((theoryDir.str() + "apfelcomb/" + setname).c_str(),0777);

  char mode[] = "0777";
  const int imode = strtol(mode, 0, 8);

  // Setup GenFK Scripts
  const std::string genFK  = setname + "_genFK.sh";
  const std::string source = "genFK_scripts/" + genFK;

  struct stat buf;
  if (stat(source.c_str(), &buf) == 0)
  {
    const std::string dest   = theoryDir.str()  + "apfelcomb/" + setname + "/genFK.sh";
    std::ifstream  src(source.c_str(), std::ios::binary);
    std::ofstream  dst(dest.c_str(),   std::ios::binary);
    dst << src.rdbuf(); chmod( dest.c_str() , imode);
  }

  // Write set inventory to file
  const std::string invPath = theoryDir.str()  + "apfelcomb/" + setname + "/inventory";
  std::ofstream invStream(invPath.c_str());
  for (auto grid : reqgrids)
    invStream << grid << std::endl;
  invStream.close();

  // Setup GenAll script
  std::ifstream  src("genFK_scripts/genAll.sh", std::ios::binary);
  const std::string dest = theoryDir.str() + "/apfelcomb/genAll.sh";
  std::ofstream  dst(dest.c_str(),   std::ios::binary);
  dst << src.rdbuf(); chmod( dest.c_str() , imode);
}

std::string getOutputFilename(int const& theoryID, std::string const& setname, std::string const& gridname)
{
  // Final output filename
  stringstream fname;
  fname << resultsPath()<<"theory_" << theoryID<<"/apfelcomb/" << setname<<"/";
  const std::string outname = fname.str() + "FK_"+gridname+".dat";
  return outname;
}

void DisplayHR()
{
  cout << Colour::FG_YELLOW;
  cout << "************************************************************" <<endl;
  cout << Colour::FG_DEFAULT;
}

void Splash()
{
  cout << Colour::FG_BLUE << endl;
  cout << "  ███╗   ██╗███╗   ██╗██████╗ ██████╗ ███████╗ " << endl;
  cout << "  ████╗  ██║████╗  ██║██╔══██╗██╔══██╗██╔════╝ " << endl;
  cout << "  ██╔██╗ ██║██╔██╗ ██║██████╔╝██║  ██║█████╗   " << endl;
  cout << "  ██║╚██╗██║██║╚██╗██║██╔═══╝ ██║  ██║██╔══╝ " << endl;
  cout << "  ██║ ╚████║██║ ╚████║██║     ██████╔╝██║ " << endl;
  cout << "  ╚═╝  ╚═══╝╚═╝  ╚═══╝╚═╝     ╚═════╝ ╚═╝ 2012-2015" << Colour::FG_DEFAULT <<endl;
  cout << "  ____ver____ : "<<NNPDF::getVersion()<<", __coredevs__ : N.H., S.C.\n" << endl;
  DisplayHR();
  cout << "                          APFELcomb                        " <<endl;
  DisplayHR();
}



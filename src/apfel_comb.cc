// applconv.cc
// Code to generate FastKernel-Like grids
// from applgrid files and evolution grids.
//
// Input: appl_comb parameter file
// Output: Combined FastKernel grid
//
// n.p.hartland@ed.ac.uk  03/12

#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <cstdio>

#include "apfelcomb/fk_appl.h"
#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_pdf.h"
#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_grids.h"

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"

using namespace std;

extern "C" {
    extern struct {
      double ad, au, as, ac;
    } eft_;
}

void SetupAPFELeftCoeff(int theoryid)
{
  double Lambda2=pow(1000.0,2.0);
    // ad = as; au = ac
    if (theoryid == 181) // BP1
      {
	eft_.ad = -0.26/Lambda2;
	eft_.au = +0.22/Lambda2;
	eft_.as = eft_.ad;
	eft_.ac = eft_.au;
      } 
    else if (theoryid == 182) // BP2
      {
	eft_.ad = +0.06/Lambda2;
	eft_.au = -0.08/Lambda2;
	eft_.as = eft_.ad;
	eft_.ac = eft_.au;
      }
    else if (theoryid == 183) // BP3
      {
	eft_.ad = 0.0;
	eft_.au = +0.16/Lambda2;
	eft_.as = eft_.ad;
	eft_.ac = eft_.au;
      }
    else if (theoryid == 184) // BP4
      {
	eft_.ad = -0.19/Lambda2;
	eft_.au = -0.04/Lambda2;
	eft_.as = eft_.ad;
	eft_.ac = eft_.au;
      }
    else if (theoryid == 191) // BP1b
      {
	eft_.ad = +0.10/Lambda2;
	eft_.au = -0.28/Lambda2;
	eft_.as = eft_.ad;
	eft_.ac = eft_.au;
      }
    else if (theoryid == 192) // BP2b
      {
	eft_.ad = -0.01/Lambda2;
	eft_.au = +0.10/Lambda2;
	eft_.as = eft_.ad;
	eft_.ac = eft_.au;
      }
    else if (theoryid == 193) // BP3b
      {
	eft_.ad = +0.18/Lambda2;
	eft_.au = -0.04/Lambda2;
	eft_.as = eft_.ad;
	eft_.ac = eft_.au;
      }
  // ad = -as; au = -ac
    else if (theoryid == 185) // BP5
      {
	eft_.ad = -0.48/Lambda2;
	eft_.au = +0.57/Lambda2;
	eft_.as = - eft_.ad;
	eft_.ac = - eft_.au;
      }
    else if (theoryid == 186) // BP6
      {
	eft_.ad = +0.43/Lambda2;
	eft_.au = -0.61/Lambda2;
	eft_.as = - eft_.ad;
	eft_.ac = - eft_.au;
      }
    else if (theoryid == 187) // BP7
      {
	eft_.ad = +0.09/Lambda2;
	eft_.au = +0.15/Lambda2;
	eft_.as = - eft_.ad;
	eft_.ac = - eft_.au;
      }
    else if (theoryid == 188) // BP8
      {
	eft_.ad = +0.11/Lambda2;
	eft_.au = 0.0;
	eft_.as = - eft_.ad;
	eft_.ac = - eft_.au;
      }
    else if (theoryid == 194) // BP11
      {
	eft_.ad = -0.75/Lambda2;
	eft_.au = +0.95/Lambda2;
	eft_.as = - eft_.ad;
	eft_.ac = - eft_.au;
      }
    else if (theoryid == 195) // BP12
      {
	eft_.ad = +0.7/Lambda2;
	eft_.au = -1.0/Lambda2;
	eft_.as = - eft_.ad;
	eft_.ac = - eft_.au;
      }
  // au = as = -ad = -ac
    else if (theoryid == 196) // BP9b
      {
	eft_.ad = +0.51/Lambda2;
	eft_.au = - eft_.ad;
	eft_.as = - eft_.ad;
	eft_.ac = eft_.ad;
      }
    else if (theoryid == 189) // BP9
      {
	eft_.ad = +0.33/Lambda2;
	eft_.au = - eft_.ad;
	eft_.as = - eft_.ad;
	eft_.ac = eft_.ad;
      }
    else if (theoryid == 197) // BP10b
      {
	eft_.ad = -0.75/Lambda2;
	eft_.au = - eft_.ad;
	eft_.as = - eft_.ad;
	eft_.ac = eft_.ad;
      }
    else if (theoryid == 190) // BP10
      {
	eft_.ad = -0.495/Lambda2;
	eft_.au = - eft_.ad;
	eft_.as = - eft_.ad;
	eft_.ac = eft_.ad;
      }
  // ad = 0, au = 0
    else if (theoryid == 198) // BP13
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 1.45/Lambda2;
	eft_.ac = 2.8/Lambda2;
      }
    else if (theoryid == 199) // BP14
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = -3.2/Lambda2;
	eft_.ac = 0.5/Lambda2;
      }
    else if (theoryid == 201) // BP16
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 1.2/Lambda2;
	eft_.ac = 1.0/Lambda2;
      }
  // au = -0.0672 - 1.32 ad
    else if (theoryid == 202) // BP17
      {
	eft_.ad = -1.7/Lambda2;
	eft_.au = - 0.0672 - 1.32*eft_.ad;
	eft_.as = +1.7/Lambda2;
	eft_.ac = 0.0;
      }
    else if (theoryid == 203) // BP18
      {
	eft_.ad = -0.3/Lambda2;
	eft_.au = -0.0672 - 1.32*eft_.ad;
	eft_.as = -3.2/Lambda2;
	eft_.ac = 0.0;
      }
    else if (theoryid == 204) // BP19
      {
	eft_.ad = +0.5/Lambda2;
	eft_.au = - 0.0672 - 1.32*eft_.ad;
	eft_.as = 0.0;
	eft_.ac = +3.0/Lambda2;
      }
    else if (theoryid == 205) // BP20
      {
	eft_.ad = -1.8/Lambda2;
	eft_.au = - 0.0672 - 1.32*eft_.ad;
	eft_.as = 0.0;
	eft_.ac = -1.2/Lambda2;
      }
  // ad = as = ac = 0
    else if (theoryid == 206) // BP21
      {
	eft_.ad = 0.0;
	eft_.au = -0.18/Lambda2;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 207) // BP22
      {
	eft_.ad = 0.0;
	eft_.au = +0.17/Lambda2;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 214)  // BP29
      {
	eft_.ad = 0.0;
	eft_.au = -1.0/Lambda2;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 215) // BP30
      {
	eft_.ad = 0.0;
	eft_.au = +1.0/Lambda2;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 216) // BP31
      {
	eft_.ad = 0.0;
	eft_.au = -0.5/Lambda2;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 217) // BP32
      {
	eft_.ad = 0.0;
	eft_.au = +0.5/Lambda2;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 218) // BP33
      {
	eft_.ad = 0.0;
	eft_.au = -0.3/Lambda2;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 219) // BP34
      {
	eft_.ad = 0.0;
	eft_.au = +0.3/Lambda2;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
  //  au = as = ac = 0
    else if (theoryid == 208) // BP23
      {
	eft_.ad = -0.22/Lambda2;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 209) // BP24
      {
	eft_.ad = +0.10/Lambda2;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 220) // BP35
      {
	eft_.ad = -1.0/Lambda2;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 221) // BP36
      {
	eft_.ad = +1.0/Lambda2;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 222) // BP37
      {
	eft_.ad = -0.5/Lambda2;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 223) // BP38
      {
	eft_.ad = +0.5/Lambda2;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 224) // BP39
      {
	eft_.ad = -0.3/Lambda2;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
    else if (theoryid == 225) // BP40
      {
	eft_.ad = +0.3/Lambda2;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = 0.0;
      }
  //  au = ad = ac = 0
    else if (theoryid == 210) // BP25
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = -2.6/Lambda2;
	eft_.ac = 0.0;
      }
    else if (theoryid == 211) // BP26
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = +0.4/Lambda2;
	eft_.ac = 0.0;
      }
    else if (theoryid == 226) // BP41
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = -0.4/Lambda2;
	eft_.ac = 0.0;
      }
    else if (theoryid == 227) // BP42
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = +2.6/Lambda2;
	eft_.ac = 0.0;
      }
    else if (theoryid == 228) // BP43
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = -1.0/Lambda2;
	eft_.ac = 0.0;
      }
    else if (theoryid == 229) // BP44
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = +1.0/Lambda2;
	eft_.ac = 0.0;
      }
    else if (theoryid == 230) // BP45
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = -0.7/Lambda2;
	eft_.ac = 0.0;
      }
    else if (theoryid == 231) // BP46
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.7/Lambda2;
	eft_.ac = 0.0;
      }
  // au = ad = as = 0
    else if (theoryid == 212) // BP27
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = -0.2/Lambda2;
      }
    else if (theoryid == 213) // BP28
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = +0.8/Lambda2;
      }
    else if (theoryid == 232) // BP47
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = -1.0/Lambda2;
      }
    else if (theoryid == 233) // BP48
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = +1.0/Lambda2;
      }
    else if (theoryid == 234) // BP49
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = -2.0/Lambda2;
      }
    else if (theoryid == 235) // BP50
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = +2.0/Lambda2;
      }
    else if (theoryid == 200) // BP15
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = -0.5/Lambda2;
      }
    else if (theoryid == 236) // BP51
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = +0.5/Lambda2;
      }
    else if (theoryid == 237) // BP52
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = -0.8/Lambda2;
      }
    else if (theoryid == 238) // BP53
      {
	eft_.ad = 0.0;
	eft_.au = 0.0;
	eft_.as = 0.0;
	eft_.ac = +0.2/Lambda2;
      }
    else 
      {
	cerr << "Wrong theory ID!" << endl;
      }
}

int main(int argc, char* argv[]) {

  if (argc!=4)
  {
    cout << "Usage: "<<argv[0]<<" <source=app/dis/dyp> <database id> <theory id>"<<endl;
    exit(1);
  }

  // APPLgrid and theory indices
  const string source = argv[1];
  const int iDB = atoi(argv[2]);
  const int iTh = atoi(argv[3]);
  setupDir(iTh);

  // Assign eft coeff
  SetupAPFELeftCoeff(iTh);

  // Parse parameters
  Splash(); QCD::qcd_param par; QCD::parse_input(iTh, par);

  NNPDF::SetVerbosity(0);
  NNPDF::IndexDB grid_db(databasePath()+"apfelcomb.db", "grids");
  NNPDF::IndexDB subgrid_db(databasePath()+"apfelcomb.db", source+"_subgrids");

  // Read grid information
  const std::string fktarget = NNPDF::dbquery<string>(subgrid_db,iDB,"fktarget");
  const std::vector<int> grid_matches = NNPDF::dbmatch(grid_db, "name", fktarget);

  if ( fktarget == "" ) // A bit clumsy, doing better would require modifying dbquery to throw an error when no match is found
    throw std::runtime_error("Cannot find subgrid "+to_string(iDB)+" in apfelcomb database");
  if ( grid_matches.size() == 0 )
    throw std::runtime_error("Cannot find FK target "+fktarget+" in apfelcomb database");
  const int target = grid_matches[0];
  FKTarget table(grid_db, target, par.global_nx); table.ReadSubGrids(subgrid_db);

  // // Initialise QCD
  if (table.GetSource() == FKTarget::DYP) QCD::setFTDYmode(true);
  QCD::initQCD(par, table.GetPositivity(), table.GetQ2max());
  QCD::initEvolgrid(table.GetNX(), table.GetXmin());

  // Initialise empty mFK table
  NNPDF::FKHeader FKhead; table.SetFKHeader(FKhead); QCD::set_params(par, FKhead);
  std::stringstream IO; FKhead.Print(IO);
  NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );

  // Compute FK table
  DisplayHR();  cout << "                        Combination                  "<<endl;
  table.GetSubgrid(iDB)->Combine(par, FK); FK->Finalise();

  DisplayHR();  cout << "                        Verification                  "<<endl;
  APFELPDFSet apfelPDF;
  const vector<double> xsec = table.Compute(par);
  const NNPDF::ThPredictions theory(&apfelPDF, FK);

  cout  << setw(10) << left << "<idat>"
        << setw(15) << left << "<FK>"
        << setw(15) << left << "<Source>"
        << setw(15) << left << "<Rel.Err>"
        << endl;

  double max_relerr = 0;
  for (auto targets : table.GetSubgrid(iDB)->GetDataMap())
    for (int i : targets)
    {
      const double applpred = xsec[i];
      const double FKpred  = theory.GetObsCV(i);
      const double rel_err = abs((FKpred-applpred)/applpred);
      max_relerr = max(max_relerr, rel_err);
      cout  << setw(10) << left << i
            << setw(15) << left << FKpred
            << setw(15) << left << applpred
            << setw(15) << left << rel_err
            << endl;
    }

  if ( max_relerr > 1E-2 && table.GetSource() != FKTarget::DYP && table.GetPositivity() == false )
  {
    cerr << "Error: Relative error exceeds expectations - something has gone wrong in the combination"<<endl;
    exit(1);
  }

  // // Print to file
  const std::string outFile = resultsPath()+"theory_" + to_string(iTh) + "/subgrids/FK_"+table.GetTargetName()+"_"+to_string(iDB) + ".subgrid.dat";
  FK->Print(outFile, true);

  DisplayHR();
  NNPDF::FKTable *impFK = new NNPDF::FKTable(outFile);
  NNPDF::ThPredictions rat = (NNPDF::ThPredictions(&apfelPDF, impFK) - NNPDF::ThPredictions(&apfelPDF, FK))
                           /  NNPDF::ThPredictions(&apfelPDF, FK);
  for (int i=0; i<theory.GetNData(); i++)
  {
    const double rel_err = abs(rat.GetObsCV(i));
    if (rel_err > 1E-5)
    {
      cerr << "Error: FK Table Export Verification failed"<<endl;
      rat.Print(cout);
      remove(outFile.c_str());
      exit(1);
    }
  }

  //  // --  Cleanup ************************************************************

  delete FK;

  cout << "                      APPLComb Complete "<<endl;
  DisplayHR();
  exit(0);
}


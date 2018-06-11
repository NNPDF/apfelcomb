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

#include "apfelcomb/fk_evol.h"
#include "apfelcomb/fk_appl.h"
#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_pdf.h"
#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_grids.h"

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"

using namespace std;

// PDF evolution tables generate their own CommonData rather than reading them from file
NNPDF::CommonData ReadCommonData(QCD::qcd_param const& par, string const& setname, string const& source)
{
    if (source == "PDF")
    {
        const int    nfl   = 14; // All for now
        const int    nQ2   = 50;
        const double Q2Min = pow(par.Q0,2);
        const double Q2Max = 1E5*1E5;

        const int nx_out  = 100;
        const double xmin = 1E-9;
        const double xmed = 0.1;
        const double xmax = 1.0;

        // Generate target x-grid
        const vector<double> xg_out = generate_xgrid(nx_out, xmin, xmed, xmax);
        const vector< vector<double> > subgrids = generate_q2subgrids(qcdPar, nQ2, Q2Min, Q2Max)[0];

        vector<double> q2_out;
        switch(setname):
        {
            case "EVOLSUBGRID1":
                q2_out = subgrids[0];
            case "EVOLSUBGRID2":
                q2_out = subgrids[1];
            case "EVOLSUBGRID3":
                q2_out = subgrids[2];
        }

         const int npoints = q2_out.size()*nx_out*nfl;
         const dataInfoRaw subgrid_info = {npoints, 0, "EVOL", "PDF"};
         const EVL::EvolutionData ed(subgrid_info, xg_out, q2_out);
         std::cout << std::endl;
         for (auto j : q2_out)
             std::cout <<  "  "<< sqrt(j) <<"  " << std::endl;
        return ed;
    }
    return NNPDF::CommonData::ReadFile(dataPath() + "commondata/DATA_" + setname + ".dat",
            dataPath() + "commondata/systypes/SYSTYPE_" + setname + "_DEFAULT.dat");
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
    const string setname = NNPDF::dbquery<string>(grid_db,target,"name");
    const NNPDF::CommonData cd = ReadCommonData(setname, source);
    FKTarget table(cd, grid_db, target, par.global_nx);
    table.ReadSubGrids(subgrid_db);

    // // Initialise QCD
    if (table.GetSource() == FKTarget::DYP) QCD::setFTDYmode(true);
    if (table.GetSource() == FKTarget::DIS) QCD::setDISmode(true);
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


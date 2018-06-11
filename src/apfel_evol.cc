// apfel_evol.cc
// Code to generate evolution-grid FK tables

#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <cstdio>
#include <map>

#include "apfelcomb/fk_evol.h"

#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_pdf.h"
#include "apfelcomb/fk_qcd.h"

using namespace std;


int main(int argc, char* argv[]) { if ( argc != 2 )
    {
        std::cerr << "Usage: " <<argv[0] <<" <target ThID>" <<std::endl;
        exit(-1);
    }

    Splash();

    // target theoryIDs
    const int th = atoi(argv[1]);
    QCD::qcd_param qcdPar;
    QCD::parse_input(th, qcdPar);
    setupDir(th);

    const double Q2Min = pow(qcdPar.Q0,2);
    const double Q2Max = 1E5*1E5;

    const int nx_in   = 195;
    const int nx_out  = 100;

    const double xmin = 1E-9;
    const double xmed = 0.1;
    const double xmax = 1.0;

    // NNPDF3.1 standard grid uses a = 30 (fished out from an old email)
    QCD::initQCD(qcdPar, false, Q2Max);
    QCD::initEvolgrid(nx_in, xmin, 30);

    // Generate target x-grid
    const vector<double> xg_out = EVL::generate_xgrid(nx_out, xmin, xmed, xmax);
    const vector< vector<double> > q2_subgrids = EVL::generate_q2subgrids(qcdPar, 50, Q2Min, Q2Max);

    // Generate and loop over subgrids
    for (size_t i=0; i< q2_subgrids.size(); i++)
    {
        const vector<double> q2_out = q2_subgrids[i];
        const int nfl     = 14; // All for now
        const int npoints = q2_out.size()*xg_out.size()*nfl;
        const dataInfoRaw subgrid_info = {npoints, 0, "EVOL", "PDF"};
        const EVL::EvolutionData ed(subgrid_info, xg_out, q2_out);

        // Initialise empty mFK table
        NNPDF::FKHeader FKhead;
        //
        //table.SetFKHeader(FKhead);
        QCD::set_params(qcdPar, FKhead);

        const string setname     = "PDFEVOL"+std::to_string(i);
        const string description = "PDF Evolution test grid";
        FKhead.AddTag(NNPDF::FKHeader::BLOB,     "GridDesc", description);
        FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "SETNAME", setname);
        FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "NDATA",   ed.GetNData());
        FKhead.AddTag(NNPDF::FKHeader::VERSIONS, "APPLrepo", applCommit() );
        FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "HADRONIC", false);
        FKhead.ResetFlavourMap();

        std::stringstream IO; FKhead.Print(IO);
        NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );

        for (int i=0; i<ed.GetNData(); i++)
        {
            const double xo = ed.GetKinematics(i,0);
            const double Q  = sqrt(ed.GetKinematics(i,1));
            const double fl = ed.GetKinematics(i,2);
            std::array<double, 14> pdf;

            for(int xi=0; xi<FK->GetNx(); xi++)
                for(int ifl=0; ifl<14; ifl++)
                {
                    // Including photon
                    QCD::EvolutionOperator(false, true, xi, xo, ifl, Q, pdf.data());
                    FK->Fill(i, xi, ifl, pdf[fl]);
                }
            std::cout << i <<"/"<<ed.GetNData() <<"\r";
            std::cout.flush();
        }

        // Write to file
        const std::string outFile = resultsPath()+"theory_" + to_string(th) + "/evolution/FK_"+setname+ ".dat";
        FK->Print(outFile, true);
        delete FK;
    }

    exit(0);
}


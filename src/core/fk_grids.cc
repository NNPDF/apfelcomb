// Classes for A matrices (combined evolution-rotation)
// and Sigma matrices (combined applgrid - A matrices)
// n.p.hartland@ed.ac.uk  - 03/12

#include "apfelcomb/fk_grids.h"
#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_utils.h"

#include "NNPDF/common.h"
#include "NNPDF/nnpdfdb.h"
#include "NNPDF/commondata.h"

#include "apfelcomb/fk_appl.h"
#include "apfelcomb/fk_dis.h"
#include "apfelcomb/fk_ftdy.h"

#include <chrono>
#include <ctime>

using namespace std;
using namespace NNPDF;

FKTarget::FKTarget(NNPDF::IndexDB const& db, int const& targetID):
id(targetID),
name(dbquery<string>(db,id,"name")),
setname(dbquery<string>(db,id,"setname")),
description(dbquery<string>(db,id,"description")),
subgrid_source(parse_source(dbquery<string>(db,id,"source"))),
positivity(dbquery<bool>(db,id,"positivity")),
nx(dbquery<int>(db,id,"nx")),
data(CommonData::ReadFile(dataPath() + "commondata/DATA_" + setname + ".dat", dataPath() + "commondata/systypes/SYSTYPE_" + setname + "_DEFAULT.dat"))
{
};

void FKTarget::Splash(ostream& o) const
{
  o << description              << endl
    << "FKTarget: " << name     << endl
    << "Setname:  " << setname  << endl
    << "NX: "       << nx       << endl
    << "Xmin: "		<< GetXmin()<< endl
    << "Qmax: "<<sqrt(GetQ2max())<< endl
    << endl;
  for (auto i : components)
  	i.second->Splash(o);
}

vector<double> FKTarget::Compute(qcd_param const& par) const
{
	vector<double> results(data.GetNData(), 0);
	for (auto subgrid : components) subgrid.second->Compute(par, results);
	return results;
}

void FKTarget::Combine(QCD::qcd_param const& par, NNPDF::FKGenerator* FK) const
{
	for (auto subgrid : components) 
		subgrid.second->Combine(par, FK);
} 


FKTarget::source FKTarget::parse_source(string const& ss)
{
  if (ss.compare("APP") == 0) return FKTarget::APP;
  if (ss.compare("DIS") == 0) return FKTarget::DIS;
  if (ss.compare("DYP") == 0) return FKTarget::DYP;
  return FKTarget::NSR;
}

void FKTarget::ReadSubGrids(NNPDF::IndexDB const& db)
{
	const std::vector<int> subgridIDs = NNPDF::dbmatch(db, "fktarget", name);
	for (int i : subgridIDs) 
		switch(subgrid_source)
		{
			case APP:
				components.insert(pair<int, FKSubGrid*>(i, new APP::SubGrid(*this, db, i)));
				break;
			case DIS:
				components.insert(pair<int, FKSubGrid*>(i, new DIS::SubGrid(*this, db, i)));
				break;
			case DYP:
				components.insert(pair<int, FKSubGrid*>(i, new FTDY::SubGrid(*this, db, i)));
				break;
			case NSR:
				cerr << "ERROR: Cannot handle source" <<endl;
				exit(-1);
		}

	// Compute subgrid maps
	int lastNdat = 0; 
    for (auto imap: components)
    {
    	FKSubGrid* subgrid = imap.second;
		lastNdat += subgrid->incdat;
		subgrid->datamap = vector< vector<int> >(subgrid->GetNdat(), vector<int>(subgrid->muldat, -1));
        for (size_t j=0; j<subgrid->GetNdat(); j++)
	    for (size_t k=0; k<subgrid->muldat; k++)
			subgrid->datamap[j][k] = lastNdat +  j*subgrid->muldat + k;
		lastNdat += subgrid->muldat*subgrid->GetNdat();
    }

	// Verify basic subgrid properties:
	if (lastNdat != data.GetNData())
	{
		std::cerr << "Parsing Error: Total of subgrid datapoints (" << lastNdat <<" ) does not equal commondata points (" <<data.GetNData() << ")"<<std::endl;
		exit(-1);
	}
}

void FKTarget::SetFKHeader(NNPDF::FKHeader& FK) const
{
	FK.AddTag(FKHeader::BLOB, "GridDesc", description);
	// FK.AddTag(FKHeader::BLOB, "Readme", par.readme);
	FK.AddTag(FKHeader::GRIDINFO, "SETNAME", setname);
	FK.AddTag(FKHeader::GRIDINFO, "NDATA", data.GetNData());
	FK.AddTag(FKHeader::VERSIONS, "APPLrepo", applCommit() );
	FK.AddTag(FKHeader::GRIDINFO, "HADRONIC", subgrid_source != DIS);

	// Full flavourmap
	stringstream fMapHeader;

	if (subgrid_source == DIS)
	{
		for (int i=0; i<14; i++)
			fMapHeader <<"1 ";
		fMapHeader<<endl;
	}
	else
	{
		for (int i=0; i<14; i++)
		{
			for (int i=0; i<14; i++)
			  fMapHeader << "1 ";
			fMapHeader<<endl;
		}
	}

	FK.AddTag(FKHeader::BLOB, "FlavourMap", fMapHeader.str());
}

double FKTarget::GetQ2max() const
{
	double q2max = 0;
	for (auto subgrid: components)
	{
		std::cout << subgrid.first <<"  " << subgrid.second->GetQ2max() <<endl;
		q2max = max(q2max, subgrid.second->GetQ2max());
	}
	return q2max;
}

double FKTarget::GetXmin() const
{
	double xmin = 1;
	for (auto subgrid: components)
		xmin = min(xmin, subgrid.second->GetXmin());
	return xmin;
}

double FKTarget::GetComputeXmin() const
{
	double xmin = 1;
	for (auto subgrid: components)
		xmin = min(xmin, subgrid.second->GetComputeXmin());
	return xmin;
}

// **************************************************************************************************

void FKSubGrid::Splash(ostream& o) const
{
	o 	<< "- Subgrid ID: " << id << endl
	    << "- NData: "  	  << GetNdat() << endl
		<< "- Operators: +:"<<incdat <<" *:"<<muldat<<" N:"<<nrmdat <<endl
	    << "- xmin: "     << GetXmin()<< endl
		<< "- Qmax: "<<sqrt(GetQ2max())<< endl;

}

void FKSubGrid::StatusUpdate( time_point const& t1, double const& complete_frac, ostream& o) const
{
      // Elapsed time update
	const time_point t2 = std::chrono::system_clock::now();
	const double percomp= 100.0*complete_frac;
	const time_point t3 = t2 + std::chrono::duration_cast<time_span>( (t2-t1) * ( 100.0 / percomp - 1.0 ));
	const std::time_t end_time = std::chrono::system_clock::to_time_t(t3);

	if (complete_frac < 1.0)
	{
		char eta[80]; strftime (eta,80,"ETA: %R %x.", localtime(&end_time));
		o 	<< "-- "<< setw(6) << setprecision(4)  << percomp << "\% complete. "
			<< eta <<"\r";
		o.flush();
	} else
	{
		o << "-- Completed in " << std::chrono::duration_cast<std::chrono::minutes>(t2 - t1).count() << " minutes"<< endl;
	}
}

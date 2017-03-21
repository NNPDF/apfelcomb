#pragma once
/*
 *  fk_grids.h
 *  FK grid common information
 * *  nph 03/17
 */

#include <stdlib.h>
#include <vector>

// NNPDF
#include "NNPDF/commondata.h"

#include "fk_utils.h"
#include "fk_qcd.h"

using std::string;

class FKSubGrid;
namespace NNPDF{
	class IndexDB;
}

// Common information for a target FK table
class FKTarget
{
public:
	FKTarget(NNPDF::IndexDB const& db, int const& targetID);
	~FKTarget(){};
	enum Source { APP, DIS, DYP, NSR };

	virtual void Splash( ostream& ); 				//!< Write table information to stream
	void ReadSubGrids(NNPDF::IndexDB const& db);	//!< Read information on subgrids from database

	double GetQ2max();								//!< Return maximum Q2-value used in this FK grid
	double GetXmin();								//!< Return minimum x-value used in this FK grid
	double GetComputeXmin();						//!< Return minimal x-value used in computation of this observable (if different to above)

	vector<double> Compute(QCD::qcd_param const&);		//!< Compute the full FK table predictions

	string GetSetName() 		const {return setname;};
	string GetDescription() 	const {return description;};
	NNPDF::CommonData const& GetCommonData() const {return data;};

private:
	const int    id;			//!< Target ID
	const string name;			//!< Name of the FK target
	const string setname;		//!< Parent dataset name
	const string description;	//!< FK table description
	const Source subgrid_source;//!< Source of subgrids (APP/DIS/DYP)

	const int 		nx;				//!< Number of x-points in the interpolation grid
	const NNPDF::CommonData data;	// Reference data file

	map<int,FKSubGrid*> components;	//!< Subgrid components

	Source parse_source(std::string const&); //!< Returns the appropriate enum for the strings "APP/DIS/DYP"
};


class FKSubGrid
{
protected:
	friend class FKTarget;
	FKSubGrid(FKTarget const& _parent, int const& _id, std::string const& operators):
	parent(_parent),
	id(_id),
	incdat(parse_operator<size_t>(operators, "+")),
	muldat(parse_operator<size_t>(operators, "*") == 0 ? 1 : parse_operator<size_t>(operators, "*")),
	nrmdat(parse_operator<double>(operators, "N") == 0 ? 1 : parse_operator<double>(operators, "N"))
	{}

	virtual void Splash( ostream& ); 								//!< Write subgrid information to stream
	virtual	void Compute(QCD::qcd_param const&, vector<double>&) = 0;	//!< Compute the full FK table predictions

	virtual size_t GetNdat()  = 0;									//!< Return number of datapoints in a subgrid
	virtual double GetQ2max() = 0;									//!< Return maximum scale used in a subgrid
	virtual double GetXmin()  = 0;									//!< Return minimum x-value used in this subgrid
	virtual double GetComputeXmin() {return GetXmin();};			//!< Return minimal x-value used in computation of this subgrid (if different to above)

protected:
	FKTarget const& parent;											//!< Parent FKTarget
	const int id;													//!< Subgrid ID
	vector< vector<int> > datamap;									//!< Map of subgrid point to data point

	// Operators
	const size_t incdat;    										//!< Datapoint increment operator
	const size_t muldat;											//!< Datapoint multiplicative operator
	const double nrmdat;											//!< Datapoint normalisation operator

	template< typename T>
	T parse_operator(std::string const& operators, string const& key)
	{
		T retval = 0;
		stringstream str;
		const vector<string> tokens = gsplit(operators, ",");
		if (operators != "")
		for (auto s:tokens) 
		{
			const vector<string> subtoken = gsplit(s, ":");
			if (subtoken.size() != 2)
			{
			  std::cerr << "Cannot parse operator: " <<s <<std::endl;
			  exit(-1);
			}

			if (subtoken[0] == key)
			{
				str << subtoken[1];
				str >> retval;
				break;
			}
		}
		return retval;
	}
};

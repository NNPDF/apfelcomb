// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"

#include "mcgrid/mcgrid.hh"

namespace Rivet {


  class D0_2013_I1253555 : public Analysis {
  public:

    /// Constructor
    D0_2013_I1253555()
      : Analysis("D0_2013_I1253555")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {

      const double minMass = 50.0*GeV;
      const double maxMass = 1.0*TeV;
      const double missingET = 25.0*GeV;

      FinalState fs;
      WFinder wfinder(fs, Cuts::abseta < 2.0 && Cuts::pT > 25*GeV, PID::MUON, minMass, maxMass, missingET, 0, WFinder::NOCLUSTER, WFinder::NOTRACK, WFinder::TRANSMASS);
      addProjection(wfinder, "WFinder");
      
      // Book Histograms
      _h_asym = bookHisto1D(1, 1, 1);

      // A bit of a hack - but the binning is the same
      _tmp_h_plus = bookHisto1D(1, 1, 2);
      _tmp_h_minus = bookHisto1D(2, 1, 2);


      // MCgrid
      const std::string configname = "D0_2013_I1253555.config";
      MCgrid::bookPDF(configname, histoDir(), MCgrid::BEAM_PROTON, MCgrid::BEAM_ANTIPROTON);

      MCgrid::gridArch arch(50,1,5,0);
      _a_plus = MCgrid::bookGrid(_tmp_h_plus, histoDir(), configname, 0, 1E-5, 1, 8315.18, 8315.18, arch);
      _a_minus = MCgrid::bookGrid(_tmp_h_minus, histoDir(), configname, 0, 1E-5, 1, 8315.18, 8315.18, arch);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      MCgrid::PDFHandler::HandleEvent( event, histoDir());
      const WFinder& wfinder = applyProjection<WFinder>(event, "WFinder");

      if (wfinder.bosons().size() != 1) 
        vetoEvent;
      
      // Get lepton constituent
      Particle l=wfinder.constituentLeptons()[0];

      const bool isMinus = l.pid() > 0;
      const bool isPosEta = l.eta() > 0;

      // Includes -A(-eta) = A(eta);
      if (isMinus == isPosEta)
      {
        _tmp_h_minus->fill(fabs(l.eta()), event.weight());
        _a_minus->fill(fabs(l.eta()), event);
      }
      else
      {
        _tmp_h_plus->fill(fabs(l.eta()), event.weight());
        _a_plus->fill(fabs(l.eta()), event);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = 100.0 * 0.2; // Data scaled by 100 * binWidth
      scale ( _tmp_h_plus , norm ) ;
      scale ( _tmp_h_minus , norm ) ;
      _a_plus->scale( norm ) ;
      _a_minus->scale( norm ) ;

      for (size_t i = 0; i < _tmp_h_plus->numBins(); ++i) {
        const double num   = _tmp_h_plus->bin(i).sumW() - _tmp_h_minus->bin(i).sumW();
        const double denom = _tmp_h_plus->bin(i).sumW() + _tmp_h_minus->bin(i).sumW();
        const double asym = (num != 0 && denom != 0) ? num / denom : 0;
        _h_asym->fill(_tmp_h_plus->bin(i).midpoint(), norm*asym);
      }

      // Normalise and export APPLgrids

      _a_plus->exportgrid();
      _a_minus->exportgrid();

      // Clear MCgrid event counter
      MCgrid::PDFHandler::CheckOutAnalysis(histoDir());
    }



  private:

    Histo1DPtr _h_asym;
    Histo1DPtr  _tmp_h_plus, _tmp_h_minus;

    MCgrid::gridPtr _a_plus;
    MCgrid::gridPtr _a_minus;


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2013_I1253555);


}

#include "fjClustering.h"
#include <bits/hhc-process.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"


#include "scales.h"

fjClustering::fjClustering()
{
  double R = 1.0;
  *this = fjClustering(fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);
}

fjClustering::fjClustering(fastjet::JetAlgorithm ja, 
			   double rparam,
			   fastjet::RecombinationScheme rs,
			   fastjet::Strategy s)
{
  rParameter = rparam;

  fjAlgorithm = ja;
  fjStrategy = s;
  fjRecombScheme = rs;
   
  fjJetDefinition = fastjet::JetDefinition(fjAlgorithm, rParameter, fjRecombScheme, fjStrategy);
  std::cout<<fjJetDefinition.description()<<std::endl;
}

fjClustering::~fjClustering()
{
  eventParticles.clear();
  fjJets.clear();
  outputJets.clear();
 }


void fjClustering::doClustering()
{

  fastjet::ClusterSequence cluster_seq(eventParticles, fjJetDefinition);
  fjJets = cluster_seq.inclusive_jets( pT_min * .8 );
  
}

void fjClustering::convert2nlojet()
{
  outputJets.clear();
  outputJets.resize(1, 0);
  for (int  i = 0; i < fjJets.size(); i++)
    {
      outputJets.push_back(
			   _Lv(
			       fjJets[i].px(),
			       fjJets[i].py(),
			       fjJets[i].pz(),
			       fjJets[i].e()));
    };
}

void fjClustering::convert2fj(const nlo::event_hhc& event)
{
  eventParticles.clear();
  for (int i = 1; i <= event.upper(); i++)
    {
      eventParticles.push_back(
			      fastjet::PseudoJet(event[i].X(), 
						 event[i].Y(), 
						 event[i].Z(), 
						 event[i].T()));
    }

}

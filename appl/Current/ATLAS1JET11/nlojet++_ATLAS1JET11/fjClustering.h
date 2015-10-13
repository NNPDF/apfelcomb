#ifndef FJCLUSTERING__HH
#define FJCLUSTERING__HH 1

#include <string>
#include <bits/hhc-process.h>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>
#include <cstdio>

typedef nlo::lorentzvector<double> _Lv;

class fjClustering
{

 public:
  fjClustering();
  fjClustering(fastjet::JetAlgorithm, double, fastjet::RecombinationScheme, fastjet::Strategy);
  ~fjClustering();

  nlo::bounded_vector<_Lv>& operator()(const nlo::event_hhc& p)
    {
      convert2fj(p);
      doClustering();
      convert2nlojet();
      return outputJets;
    };

 private:

  std::vector<fastjet::PseudoJet> eventParticles;
  std::vector<fastjet::PseudoJet> fjJets;
  nlo::bounded_vector<_Lv> outputJets;
  void convert2fj(const nlo::event_hhc&);
  void convert2nlojet();
  void doClustering();

 private:
  fastjet::JetAlgorithm fjAlgorithm;

  fastjet::JetDefinition fjJetDefinition;
  double rParameter;
  fastjet::Strategy fjStrategy;
  fastjet::RecombinationScheme fjRecombScheme;
};

#endif

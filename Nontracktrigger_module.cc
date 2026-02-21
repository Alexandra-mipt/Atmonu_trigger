////////////////////////////////////////////////////////////////////////
// Class:       Nontracktrigger
// Plugin Type: filter (art v2_12_01)
// File:        Nontracktrigger_module.cc
//
// Generated at Tue Feb 10 14:24:26 2026 by alexandra_mipt using cetskelgen
// from cetlib version v3_06_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "DDTBaseDataProducts/HitList.h"
#include "DDTBaseDataProducts/DAQHit.h"

#include "Nontracktrigger.h"

#include <memory>

namespace novaddt {
  class Nontracktrigger;
}

class novaddt::Nontracktrigger : public art::EDFilter {
public:
  explicit Nontracktrigger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Nontracktrigger(Nontracktrigger const &) = delete;
  Nontracktrigger(Nontracktrigger &&) = delete;
  Nontracktrigger & operator = (Nontracktrigger const &) = delete;
  Nontracktrigger & operator = (Nontracktrigger &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:
  std::string _SliceModuleTag;
  TriggerNonTrack _trigger;
  unsigned int _trigger_counts = 0;
};


novaddt::Nontracktrigger::Nontracktrigger(fhicl::ParameterSet const & p)
  : _SliceModuleTag(p.get<std::string>("SliceModuleTag")),
    _trigger(p)
{
  // Call appropriate produces<>() functions here.
}

bool novaddt::Nontracktrigger::filter(art::Event & e)
{
  std::vector<Hit> hits;
  
  auto hitHandle = e.getValidHandle<std::vector<std::vector<novaddt::DAQHit>>>(_SliceModuleTag);

  for (size_t hitSet_id = 0; hitSet_id < hitHandle->size(); ++hitSet_id) {
    for (const auto& h : hitHandle->at(hitSet_id)) {
      hits.emplace_back(h, hitSet_id);
    }
  }

  bool passed = _trigger.run_algorithm(hits, e);
  std::cout<<"LOOOOK"<<passed<<"\n";
  _trigger_counts++;
  
  return passed;
}

DEFINE_ART_MODULE(novaddt::Nontracktrigger)

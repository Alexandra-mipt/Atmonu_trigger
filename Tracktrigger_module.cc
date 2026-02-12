////////////////////////////////////////////////////////////////////////
// Class:       Tracktrigger
// Plugin Type: filter (art v2_12_01)
// File:        Tracktrigger_module.cc
//
// Generated at Tue Feb 10 14:24:05 2026 by alexandra_mipt using cetskelgen
// from cetlib version v3_06_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Tracktrigger.h"
#include "DDTBaseDataProducts/Track3D.h"
#include "DDTBaseDataProducts/HitList.h"
#include "DDTBaseDataProducts/DAQHit.h"
#include "DDTBaseDataProducts/TriggerDecision.h"

#include <memory>
#include <vector>
#include <utility>

namespace novaddt {
  class Tracktrigger;
}


Hit::Hit(const novaddt::DAQHit& h, int slice_id):
    hitSet_id(slice_id),
    plane(h.Plane().val),
    cell(h.Cell().val),
    adc(h.ADC().val),
    tdc(h.TDC().val),
    view(h.View().val),
    used(false)
    {}


Track::Track(const novaddt::Track3D& t, int slice_id):
    fView(t.View()),
    StartX(t.Start().X()),
    StartY(t.Start().Y()),
    StartZ(t.Start().Z()),
    EndX(t.End().X()),
    EndY(t.End().Y()),
    EndZ(t.End().Z()),
    sliceID(slice_id)
    { }


class novaddt::Tracktrigger : public art::EDFilter {
public:
  explicit Tracktrigger(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  Tracktrigger(Tracktrigger const &) = delete;
  Tracktrigger(Tracktrigger &&) = delete;
  Tracktrigger & operator = (Tracktrigger const &) = delete;
  Tracktrigger & operator = (Tracktrigger &&) = delete;

  bool filter(art::Event & e) override;

private:

      std::string _TrackModuleTag;
      std::string _AssnsTag;
      
      Trigger _trigger;

      unsigned int _trigger_counts = 0;

};


novaddt::Tracktrigger::Tracktrigger(fhicl::ParameterSet const & p)
 : _TrackModuleTag (p.get<std::string>("TrackModuleTag")),
   _AssnsTag       (p.get<std::string>("AssnsTag")),
   _trigger(p)

// Initialize member data here.
{
  //  produces< std::vector<novaddt::TriggerDecision>>();
  //  produces< art::Assns<novaddt::TriggerDecision, novaddt::HitList> >();
}

bool novaddt::Tracktrigger::filter(art::Event & e)
{
  // Implementation of required member function here.
    // Tracks and Hits
    std::vector<Hit> hits;
    std::vector<Track> tracks;

    // Get from the event
    auto trackHandle = e.getValidHandle<std::vector<novaddt::Track3D>>(_TrackModuleTag);

    // Connect track with HitList
    //art::FindOneP<novaddt::HitList> fohl(trackHandle, e, _TrackModuleLabel);
    art::FindOneP<novaddt::HitList> fohl(trackHandle, e, _AssnsTag);

    unsigned long hitSet_number = 0;

    for (size_t i=0; i<trackHandle->size(); ++i) {
        art::Ptr<novaddt::HitList> hitList = fohl.at(i);
 
        tracks.emplace_back(trackHandle->at(i), hitSet_number);

       for (auto const& h : *hitList) 
            hits.emplace_back(h, hitSet_number);
        
        hitSet_number++;
    }


    // Run trigger algorithm
    bool passed = _trigger.run_algorithm(hits, tracks);
    _trigger_counts++;

    // Trigger decision


    return passed;
}

DEFINE_ART_MODULE(novaddt::Tracktrigger)

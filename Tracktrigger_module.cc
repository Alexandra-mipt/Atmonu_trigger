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
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "DDTBaseDataProducts/Track3D.h"
#include "DDTBaseDataProducts/HitList.h"
#include "DDTBaseDataProducts/DAQHit.h"
#include "DDTBaseDataProducts/TriggerDecision.h"

#include "Tracktrigger.h"

#include <memory>
#include <vector>
#include <utility>
#include <map>

namespace novaddt {
  class Tracktrigger;
}

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
  std::string _SlicesTag;     
  Trigger _trigger;

  unsigned int _prescale;
  unsigned int _trigger_counts;

};


novaddt::Tracktrigger::Tracktrigger(fhicl::ParameterSet const & p)
 : _TrackModuleTag(p.get<std::string>("TrackModuleTag")),
   _AssnsTag(p.get<std::string>("AssnsTag")),
   _SlicesTag(p.get<std::string>("SlicesTag")),
   _trigger(p),
   _prescale(p.get<unsigned>("prescale")),
   _trigger_counts(0)
{
   produces< std::vector<novaddt::TriggerDecision>>();
  //  produces< art::Assns<novaddt::TriggerDecision, novaddt::HitList> >();
}

bool novaddt::Tracktrigger::filter(art::Event & e)
{
    // Tracks and Hits
    std::vector<Hit> hits;
    std::vector<Track> tracks;

    // Get from the event
    auto slices = e.getValidHandle<std::vector<novaddt::HitList>>(_SlicesTag);
    art::FindMany<novaddt::Track3D> fm_tracks_in_slice(slices, e, _AssnsTag);

    for (size_t slice_id=0; slice_id<fm_tracks_in_slice.size(); slice_id++){ 
        for(auto const& track : fm_tracks_in_slice.at(slice_id)) 
        {
            tracks.emplace_back(*track, slice_id);    
        }
        for(auto const& hit : slices -> at(slice_id))
            hits.emplace_back(hit, slice_id);
     
    }

    for(const auto& track : tracks)
        std::cout<<"Track"<<track.sliceID<<" "<< track.StartX <<"\n";

    for(const auto& hit : hits)  
        std::cout<<"Hits"<<hit.hitSet_id<<" "<< hit.adc<<"\n";

   
/* auto trackHandle = e.getValidHandle<std::vector<novaddt::Track3D>>(_TrackModuleTag);

    // Connect track with HitList
    //art::FindOneP<novaddt::HitList> fohl(trackHandle, e, _TrackModuleLabel);
    art::FindOneP<novaddt::HitList> fohl(trackHandle, e, _AssnsTag);

    // Для отслеживания уникальных слайсов и их ID
    std::map<art::Ptr<novaddt::HitList>, size_t> slice_map;
    size_t next_slice_id = 0;

    // Проходим по всем трекам
    for (size_t i = 0; i < trackHandle->size(); ++i) {
        // Получаем указатель на HitList (слайс) для этого трека
        art::Ptr<novaddt::HitList> hitListPtr = fohl.at(i);
        
        // Определяем slice_id для этого слайса
        size_t slice_id;
        auto it = slice_map.find(hitListPtr);
        if (it == slice_map.end()) {
            // Новый слайс - присваиваем новый ID и запоминаем
            slice_id = next_slice_id++;
            slice_map[hitListPtr] = slice_id;

            // Добавляем ВСЕ хиты из этого слайса (только один раз!)
            for (const auto& h : *hitListPtr) {
                hits.emplace_back(h, slice_id);
            }

        } else {
            // Уже известный слайс
            slice_id = it->second;
        }

        // Добавляем трек с правильным slice_id
        tracks.emplace_back(trackHandle->at(i), slice_id); 
    }    
*/
//  std::unique_ptr<std::vector<novaddt::TriggerDecision> > 
  //  trigger_decisions(new std::vector<novaddt::TriggerDecision>());
  
    std::unique_ptr<std::vector<novaddt::TriggerDecision>> trigger_decisions(new std::vector<novaddt::TriggerDecision>);

  // Run trigger algorithm
    bool passed = _trigger.run_algorithm(hits, tracks, e);


  if (passed) {
    for (auto td : _trigger.TriggerDecisions()) {
      _trigger_counts++;

      if (_trigger_counts % _prescale) {
        td.setPrescale(_prescale);
        trigger_decisions->push_back(td);
      }
    }
  }
    e.put(std::move(trigger_decisions));
    // Trigger decision
    return passed;
}

DEFINE_ART_MODULE(novaddt::Tracktrigger)

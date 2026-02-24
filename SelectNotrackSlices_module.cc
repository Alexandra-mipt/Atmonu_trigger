////////////////////////////////////////////////////////////////////////
// Class:       SelectNotrackSlices
// Plugin Type: filter (art v2_12_01)
// File:        SelectNotrackSlices_module.cc
//
// Generated at Fri Feb  6 17:57:44 2026 by Andrey Sheshukov using cetskelgen
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
#include "canvas/Persistency/Common/FindMany.h"

#include <memory>

#include "DDTBaseDataProducts/Track3D.h"
#include "DDTBaseDataProducts/HitList.h"
#include "DDTBaseDataProducts/GroupedHitList.h"
#include "DDTUtilities/DetectorUtils.h"

namespace novaddt {
  class SelectNotrackSlices;
}


class novaddt::SelectNotrackSlices : public art::EDFilter {
public:
  explicit SelectNotrackSlices(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SelectNotrackSlices(SelectNotrackSlices const &) = delete;
  SelectNotrackSlices(SelectNotrackSlices &&) = delete;
  SelectNotrackSlices & operator = (SelectNotrackSlices const &) = delete;
  SelectNotrackSlices & operator = (SelectNotrackSlices &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;
  size_t _size = 0;
private:
  art::InputTag fSlicesTag;
  art::InputTag fAssnsTag;
};


novaddt::SelectNotrackSlices::SelectNotrackSlices(fhicl::ParameterSet const & p):
  fSlicesTag(p.get<std::string>("slices_tag")),
  fAssnsTag(p.get<std::string>("assns_tag"))
{
  produces<std::vector<novaddt::HitList>>();
}

bool novaddt::SelectNotrackSlices::filter(art::Event & e)
{
  auto slices = e.getValidHandle<std::vector<novaddt::HitList>>(fSlicesTag);
  art::FindMany<novaddt::Track3D> fm_tracks_in_slice(slices, e, fAssnsTag); //// assosiations
  //prepare resulting vector
  auto result=std::make_unique<std::vector<novaddt::HitList>>();
  for (size_t slice_id=0; slice_id<fm_tracks_in_slice.size(); slice_id++){
    const auto & tracks = fm_tracks_in_slice.at(slice_id);
    novaddt::HitList product = slices->at(slice_id);
    _size+=slices->at(slice_id).size();
    if(tracks.size()==0) result->push_back(product);
  }
  std::cout << "slices-size: " << _size << "\n";
  mf::LogDebug("SelectNotrackSlices")<<" selected "<<result->size()<<" slices";
  e.put(std::move(result));
  return !slices->empty();
}

DEFINE_ART_MODULE(novaddt::SelectNotrackSlices)

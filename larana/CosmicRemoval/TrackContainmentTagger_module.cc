////////////////////////////////////////////////////////////////////////
// Class:       TrackContainmentTagger
// Module Type: producer
// File:        TrackContainmentTagger_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "TTree.h"

#include "TrackContainment/TrackContainmentAlg.hh"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"

namespace trk {
  class TrackContainmentTagger;
}

class trk::TrackContainmentTagger : public art::EDProducer {
public:
  explicit TrackContainmentTagger(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackContainmentTagger(TrackContainmentTagger const &) = delete;
  TrackContainmentTagger(TrackContainmentTagger &&) = delete;
  TrackContainmentTagger & operator = (TrackContainmentTagger const &) = delete;
  TrackContainmentTagger & operator = (TrackContainmentTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  // Declare member data here.
  trk::TrackContainmentAlg fAlg;

  std::vector<std::string> fTrackModuleLabels;
  std::vector<bool>        fApplyTags;
};


trk::TrackContainmentTagger::TrackContainmentTagger(fhicl::ParameterSet const & p)
  : EDProducer{p}
{
  art::ServiceHandle<art::TFileService const> tfs;
  fAlg.SetupOutputTree(tfs->make<TTree>("myanatree","MyAnalysis Tree"));

  fAlg.Configure(p.get<fhicl::ParameterSet>("TrackContainmentAlg"));
  fTrackModuleLabels = p.get< std::vector<std::string> >("TrackModuleLabels");
  fApplyTags = p.get< std::vector<bool> >("ApplyTags",std::vector<bool>(fTrackModuleLabels.size(),true));

  if(fApplyTags.size()!=fTrackModuleLabels.size())
    throw cet::exception("TrackContainmentTagger::reconfigure")
      << "ApplyTags not same size as TrackModuleLabels. ABORT!!!";

  fAlg.setMakeCosmicTags();

  produces< std::vector<anab::CosmicTag> >();
  produces< art::Assns<recob::Track, anab::CosmicTag> >();
}

void trk::TrackContainmentTagger::produce(art::Event & e)
{

  std::unique_ptr< std::vector< anab::CosmicTag > >              cosmicTagTrackVector( new std::vector<anab::CosmicTag> );
  std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag > >  assnOutCosmicTagTrack( new art::Assns<recob::Track, anab::CosmicTag>);

  fAlg.SetRunEvent(e.run(),e.event());

  std::vector< std::vector<recob::Track> > trackVectors;
  std::vector< art::Handle< std::vector<recob::Track> > > trackHandles;
  for(size_t i_l=0; i_l<fTrackModuleLabels.size(); ++i_l){
    art::Handle< std::vector<recob::Track> >  trackHandle;
    e.getByLabel(fTrackModuleLabels[i_l],trackHandle);
    trackVectors.push_back(*trackHandle);
    trackHandles.push_back(trackHandle);
  }

  art::ServiceHandle<geo::Geometry const> geoHandle;
  fAlg.ProcessTracks(trackVectors,*geoHandle);

  auto const& cosmicTags = fAlg.GetTrackCosmicTags();

  for(size_t i_tc=0; i_tc<cosmicTags.size(); ++i_tc){
    if(!fApplyTags[i_tc]) continue;
    for(size_t i_t=0; i_t<fAlg.GetTrackCosmicTags()[i_tc].size(); ++i_t){
      cosmicTagTrackVector->emplace_back(fAlg.GetTrackCosmicTags()[i_tc][i_t]);
      util::CreateAssn(*this, e, *cosmicTagTrackVector, art::Ptr<recob::Track>(trackHandles[i_tc],i_t), *assnOutCosmicTagTrack );
    }
  }

  e.put(std::move(cosmicTagTrackVector));
  e.put(std::move(assnOutCosmicTagTrack));

}

DEFINE_ART_MODULE(trk::TrackContainmentTagger)

/////////////////////////////////////////////////////////////////
//  \fileMVAAlg.h
//  m.haigh@warwick.ac.uk
////////////////////////////////////////////////////////////////////
#ifndef MVAAlg_H
#define MVAAlg_H

#include <map>
#include <vector>

#include "art/Framework/Core/Frameworkfwd.h"
#include "art/Framework/Principal/fwd.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/Vector3Dfwd.h"
#include "TGraph2D.h"
#include "TLorentzVector.h"
#include "TMVA/Reader.h"
#include "TVector3.h"

namespace mvapid {

  //---------------------------------------------------------------
  class MVAAlg {
  public:
    struct SortedObj {
      TVector3 start, end, dir;
      double length;
      std::map<double, const art::Ptr<recob::Hit>> hitMap;
    };

    struct SumDistance2 {
      // the TGraph is a data member of the object
      TGraph2D* fGraph;

      SumDistance2(TGraph2D* g) : fGraph(g) {}

      // implementation of the function to be minimized
      double
      operator()(const double* p)
      {

        ROOT::Math::XYZVector x0(p[0], p[2], p[4]);
        ROOT::Math::XYZVector u(p[1], p[3], p[5]);

        u = u.Unit();
        double* x = fGraph->GetX();
        double* y = fGraph->GetY();
        double* z = fGraph->GetZ();
        int npoints = fGraph->GetN();
        double sum = 0;
        for (int i = 0; i < npoints; ++i) {
          ROOT::Math::XYZVector xp(x[i], y[i], z[i]);
          sum += ((xp - x0).Cross(u)).Mag2();
        }
        return sum;
      }
    };

    MVAAlg(fhicl::ParameterSet const& pset);

    void GetDetectorEdges();

    void GetWireNormals();

    void RunPID(art::Event& evt,
                std::vector<anab::MVAPIDResult>& result,
                art::Assns<recob::Track, anab::MVAPIDResult, void>& trackAssns,
                art::Assns<recob::Shower, anab::MVAPIDResult, void>& showerAssns);

  private:
    int IsInActiveVol(const TVector3& pos);

    void PrepareEvent(const art::Event& event, const detinfo::DetectorClocksData& clockData);

    void FitAndSortTrack(art::Ptr<recob::Track> track, int& isStoppingReco, SortedObj& sortedObj);

    //void SortShower(art::Ptr<recob::Shower> shower,TVector3 dir,int& isStoppingReco,
    //		    mvapid::MVAAlg::SortedObj& sortedShower);
    void SortShower(art::Ptr<recob::Shower> shower,
                    int& isStoppingReco,
                    mvapid::MVAAlg::SortedObj& sortedShower);

    void RunPCA(std::vector<art::Ptr<recob::Hit>>& hits,
                std::vector<double>& eVals,
                std::vector<double>& eVecs);

    void _Var_Shape(const SortedObj& track,
                    double& coreHaloRatio,
                    double& concentration,
                    double& conicalness);

    double CalcSegmentdEdxFrac(const detinfo::DetectorClocksData& clock_data,
                               const detinfo::DetectorPropertiesData& det_prop,
                               const SortedObj& track,
                               double start,
                               double end);

    double CalcSegmentdEdxDist(const detinfo::DetectorClocksData& clock_data,
                               const detinfo::DetectorPropertiesData& det_prop,
                               const SortedObj& track,
                               double start,
                               double end);

    double CalcSegmentdEdxDistAtEnd(const detinfo::DetectorClocksData& clock_data,
                                    const detinfo::DetectorPropertiesData& det_prop,
                                    const mvapid::MVAAlg::SortedObj& track,
                                    double distAtEnd);

    int LinFit(const art::Ptr<recob::Track> track, TVector3& trackPoint, TVector3& trackDir);

    int LinFitShower(const art::Ptr<recob::Shower> shower,
                     TVector3& showerPoint,
                     TVector3& showerDir);

    const calo::CalorimetryAlg fCaloAlg;

    double fEventT0;

    double fDetMinX, fDetMaxX, fDetMinY, fDetMaxY, fDetMinZ, fDetMaxZ;

    std::map<int, double> fNormToWiresY;
    std::map<int, double> fNormToWiresZ;

    std::string fTrackLabel;
    std::string fShowerLabel;
    std::string fHitLabel;
    std::string fSpacePointLabel;
    std::string fTrackingLabel;

    std::vector<art::Ptr<recob::Track>> fTracks;
    std::vector<art::Ptr<recob::Shower>> fShowers;
    std::vector<art::Ptr<recob::SpacePoint>> fSpacePoints;
    std::vector<art::Ptr<recob::Hit>> fHits;

    std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<recob::Hit>>> fTracksToHits;
    std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<recob::SpacePoint>>> fTracksToSpacePoints;
    std::map<art::Ptr<recob::Shower>, std::vector<art::Ptr<recob::Hit>>> fShowersToHits;
    std::map<art::Ptr<recob::Shower>, std::vector<art::Ptr<recob::SpacePoint>>>
      fShowersToSpacePoints;
    std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> fHitsToSpacePoints;
    std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> fSpacePointsToHits;

    anab::MVAPIDResult fResHolder;

    TMVA::Reader fReader;

    std::vector<std::string> fMVAMethods;
    std::vector<std::string> fWeightFiles;

    bool fCheatVertex;

    TLorentzVector fVertex4Vect;

  }; // class MVAAlg

} // namespace mvapid

#endif // ifndef MVAAlg_H

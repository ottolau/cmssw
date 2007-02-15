
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "RecoBTag/SoftLepton/interface/SoftElectronProducer.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/Vector/interface/GlobalVector.h"

#include <vector>
#include <iostream>
#include <boost/regex.hpp>

#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

#include "ElectronIdMLP.h"


//------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//------------------------------------------------------------------------------

SoftElectronProducer::SoftElectronProducer(const edm::ParameterSet &iConf) :
  theConf(iConf), theTrackAssociator(0), theElecNN(0)
{
  theTrackTag = theConf.getParameter<InputTag>("TrackTag");

  theBasicClusterTag = theConf.getParameter<InputTag>("BasicClusterTag");
  theBasicClusterShapeTag = theConf.getParameter<InputTag>("BasicClusterShapeTag");

  theHOverEConeSize = theConf.getParameter<double>("HOverEConeSize");

  theTrackAssociator = new TrackDetectorAssociator();
  theTrackAssociator->useDefaultPropagator();

  // fill data labels
  theTrackAssociator->theEBRecHitCollectionLabel = theConf.getParameter<InputTag>("EBRecHitCollectionLabel");
  theTrackAssociator->theEERecHitCollectionLabel = theConf.getParameter<InputTag>("EERecHitCollectionLabel");

  theDiscriminatorCut = theConf.getParameter<double>("DiscriminatorCut");

  theElecNN = new ElectronIdMLP;

  string theCollectionName = theConf.getParameter<string>("ElectronCollection");

  // register the product
  produces<reco::ElectronCollection>(theCollectionName);
}

//------------------------------------------------------------------------------

SoftElectronProducer::~SoftElectronProducer()
{
  if(theElecNN) delete theElecNN;
  if(theTrackAssociator) delete theTrackAssociator;
}

//------------------------------------------------------------------------------

void SoftElectronProducer::beginJob (edm::EventSetup const & theEventSetup)
{

}

//------------------------------------------------------------------------------

void SoftElectronProducer::produce(edm::Event &iEvent,
                                      const edm::EventSetup &iSetup)
{

  auto_ptr<reco::ElectronCollection> candidates(new reco::ElectronCollection());

  Handle<reco::TrackCollection> handleTrack;
  reco::TrackCollection::const_iterator itTrack;

  Handle<reco::BasicClusterCollection> handleCluster;
  reco::BasicClusterCollection::const_iterator itCluster;

  Handle<reco::ClusterShapeCollection> handleShape;
  reco::ClusterShapeCollection::const_iterator itShape;

  Handle<HBHERecHitCollection> handleRecHit;
  CaloRecHitMetaCollectionV::const_iterator itRecHit;

  ESHandle<CaloGeometry> handleCaloGeom;
  
  double hcalEnergy, clusEnergy, trkP;
  double x[2], y[2], z[2], eta[2], phi[2];
  double covEtaEta, covEtaPhi, covPhiPhi, emFraction, deltaE;
  double eMax, e2x2, e3x3, e5x5, v1, v2, v3, v4;
  double value, dist, distMin;

  const reco::BasicCluster *matchedCluster;
  const reco::ClusterShape *matchedShape;
  const reco::Track *track;

  // get basic clusters
  iEvent.getByLabel(theBasicClusterTag, handleCluster);

  // get basic cluster shapes
  iEvent.getByLabel(theBasicClusterShapeTag, handleShape);

  // get rec. hits
  iEvent.getByType(handleRecHit);
  HBHERecHitMetaCollection metaRecHit(*handleRecHit);

  // get calorimeter geometry
  iSetup.get<IdealGeometryRecord>().get(handleCaloGeom);

  CaloConeSelector selectorRecHit(theHOverEConeSize, handleCaloGeom.product(), DetId::Hcal);

  // get tracks
  iEvent.getByLabel(theTrackTag, handleTrack);

  FreeTrajectoryState tmpFTS;
  TrackDetMatchInfo info;
  TrackDetectorAssociator::AssociatorParameters parameters;
  parameters.useEcal = true ;
  parameters.useHcal = false ;
  parameters.useHO = false ;
  parameters.useCalo = false ;
  parameters.useMuon = false ;
  parameters.dREcal = 0.03;

  unsigned int counterTrack;

  // loop over tracks
  for(itTrack = handleTrack->begin(), counterTrack = 0;
      itTrack != handleTrack->end();
      ++itTrack, ++counterTrack)
  {
    track = &(*itTrack);

    tmpFTS = theTrackAssociator->getFreeTrajectoryState(iSetup, *track);
    info = theTrackAssociator->associate(iEvent, iSetup, tmpFTS, parameters);

    x[0] = info.trkGlobPosAtEcal.x();
    y[0] = info.trkGlobPosAtEcal.y();
    z[0] = info.trkGlobPosAtEcal.z();
    eta[0] = info.trkGlobPosAtEcal.eta();
    phi[0] = info.trkGlobPosAtEcal.phi();

    // analyse only tracks passing quality cuts
    if(track->numberOfValidHits() >= 8 && track->pt() > 2.0 &&
       abs(track->eta()) < 1.2 && info.isGoodEcal)
    {
      distMin = 1.0e6;
      matchedCluster = 0;
      matchedShape = 0;

      // loop over basic clusters
      for(itCluster = handleCluster->begin(), itShape = handleShape->begin();
          itCluster != handleCluster->end() && itShape != handleShape->end();
          ++itCluster, ++itShape)
      {
        x[1] = itCluster->x();
        y[1] = itCluster->y();
        z[1] = itCluster->z();

        eta[1] = itCluster->eta();
        phi[1] = itCluster->phi();


        dist = hypot(x[0] - x[1], y[0] - y[1]);
        dist = hypot(dist, z[0] - z[1]);

        if(dist < distMin)
        {
          distMin = dist;
          matchedCluster = &(*itCluster);
          matchedShape = &(*itShape);
        }
      }

      // identify electrons based on cluster properties
      if(matchedCluster && matchedShape && distMin < 5.0)
      {

        GlobalPoint position(matchedCluster->x(), matchedCluster->y(), matchedCluster->z());
        auto_ptr<CaloRecHitMetaCollectionV> chosen = selectorRecHit.select(position, metaRecHit);
        hcalEnergy = 0.0;
        for(itRecHit = chosen->begin(); itRecHit != chosen->end(); ++itRecHit)
        {
          hcalEnergy += itRecHit->energy();
        }

        clusEnergy = matchedCluster->energy();
        trkP = track->p();

        deltaE = (clusEnergy - trkP)/(clusEnergy + trkP);
        emFraction =  clusEnergy/(clusEnergy + hcalEnergy);

        eMax = matchedShape->eMax();
        e2x2 = matchedShape->e2x2();
        e3x3 = matchedShape->e3x3();
        e5x5 = matchedShape->e5x5();
        v1 = eMax/e3x3;
        v2 = eMax/e2x2;
        v3 = e2x2/e5x5;
        v4 = (e3x3 - eMax)/(e5x5 - eMax);

        covEtaEta = matchedShape->covEtaEta();
        covEtaPhi = matchedShape->covEtaPhi();
        covPhiPhi = matchedShape->covPhiPhi();


        value = theElecNN->value(0, covEtaEta, covEtaPhi, covPhiPhi,
                                 v1, v2, v3, v4, emFraction, deltaE);

        if (value > theDiscriminatorCut)
        {
          const reco::Particle::LorentzVector  p4(0.0, 0.0, 0.0, clusEnergy);
          const reco::Particle::Point vtx(0.0, 0.0, 0.0);

          reco::Electron newCandidate(0, p4, vtx);
          reco::TrackRef refTrack(handleTrack, counterTrack);
  
          newCandidate.setTrack(refTrack);

          candidates->push_back(newCandidate);
        }
      }
    }
  }
  
  // put the product in the event
  iEvent.put(candidates, theCollectionName);

}

//------------------------------------------------------------------------------


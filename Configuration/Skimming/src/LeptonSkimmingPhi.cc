
// -*- C++ -*-
//
// Package:    SkimmingForB/LeptonSkimmingPhi
// Class:      LeptonSkimmingPhi
// 
/**\class LeptonSkimmingPhi LeptonSkimmingPhi.cc SkimmingForB/LeptonSkimming/plugins/LeptonSkimmingPhi.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ka Tung Lau ka.tung.lau@cern.ch
//         Created:  Tue, 29 Jan 2019 15:23:09 GMT
//
//


#include "Configuration/Skimming/interface/LeptonSkimmingPhi.h"


using namespace edm;
using namespace reco;
using namespace std;

LeptonSkimmingPhi::LeptonSkimmingPhi(const edm::ParameterSet& iConfig):
  electronsToken_(consumes<std::vector<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>  ("electrons"))),
  muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  Tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  conversionsToken_(consumes< reco::ConversionCollection > (iConfig.getParameter<edm::InputTag> ("conversions"))),

  trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag> ("triggerresults"))),
  trigobjectsToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag> ("triggerobjects"))),
  HLTFilter_(iConfig.getParameter<vector<string> >("HLTFilter")),
  HLTPath_(iConfig.getParameter<vector<string> >("HLTPath"))
{
  edm::ParameterSet runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");     
  PtTrack_Cut=runParameters.getParameter<double>("PtTrack_Cut");
  EtaTrack_Cut=runParameters.getParameter<double>("EtaTrack_Cut");
  //cout<<PtTrack_Cut<<endl;
  MinChi2Track_Cut=runParameters.getParameter<double>("MinChi2Track_Cut");
  MaxChi2Track_Cut=runParameters.getParameter<double>("MaxChi2Track_Cut");
  MuTrkMinDR_Cut=runParameters.getParameter<double>("MuTrkMinDR_Cut");
  MaxMee_Cut=runParameters.getParameter<double>("MaxMee_Cut");
  MinMee_Cut=runParameters.getParameter<double>("MinMee_Cut");
  Probee_Cut=runParameters.getParameter<double>("Probee_Cut");
  Cosee_Cut=runParameters.getParameter<double>("Cosee_Cut");
 
  PtKleadTrack_Cut=runParameters.getParameter<double>("PtKleadTrack_Cut");
  PtKsubleadTrack_Cut=runParameters.getParameter<double>("PtKsubleadTrack_Cut");
  MaxMB_Cut=runParameters.getParameter<double>("MaxMB_Cut");
  MinMB_Cut=runParameters.getParameter<double>("MinMB_Cut");
  TrkTrkMinDR_Cut=runParameters.getParameter<double>("TrkTrkMinDR_Cut");

  MaxMphi_Cut=runParameters.getParameter<double>("MaxMphi_Cut");
  MinMphi_Cut=runParameters.getParameter<double>("MinMphi_Cut");

  TrackSdxy_Cut=runParameters.getParameter<double>("TrackSdxy_Cut");
 
  MuTrgMatchCone=runParameters.getParameter<double>("MuTrgMatchCone");
  SkipIfNoMuMatch=runParameters.getParameter<bool>("SkipIfNoMuMatch");
  EpairZvtx_Cut=runParameters.getParameter<double>("EpairZvtx_Cut");
  Ksdxy_Cut=runParameters.getParameter<double>("Ksdxy_Cut");
  ProbeeKK_Cut=runParameters.getParameter<double>("ProbeeKK_Cut");
  CoseeKK_Cut=runParameters.getParameter<double>("CoseeKK_Cut");
  TrackMuDz_Cut=runParameters.getParameter<double>("TrackMuDz_Cut");
  TrgExclusionCone=runParameters.getParameter<double>("TrgExclusionCone");
  SLxy_Cut=runParameters.getParameter<double>("SLxy_Cut");
  PtB_Cut=runParameters.getParameter<double>("PtB_Cut");
  PtMu_Cut=runParameters.getParameter<double>("PtMu_Cut");
  PtEl_Cut=runParameters.getParameter<double>("PtEl_Cut");
  QualMu_Cut=runParameters.getParameter<double>("QualMu_Cut");
  MuTrgExclusionCone=runParameters.getParameter<double>("MuTrgExclusionCone");
  ElTrgExclusionCone=runParameters.getParameter<double>("ElTrgExclusionCone");
  TrkObjExclusionCone=runParameters.getParameter<double>("TrkObjExclusionCone");
  MuTrgMuDz_Cut=runParameters.getParameter<double>("MuTrgMuDz_Cut");
  ElTrgMuDz_Cut=runParameters.getParameter<double>("ElTrgMuDz_Cut");
  ObjPtLargerThanTrack=runParameters.getParameter<bool>("ObjPtLargerThanTrack");

}


LeptonSkimmingPhi::~LeptonSkimmingPhi()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


bool LeptonSkimmingPhi::hltFired(const edm::Event& iEvent, const edm::EventSetup& iSetup,std::vector<string> HLTPath ){
  using namespace std;  using namespace edm;  using namespace reco;
  using namespace trigger;
 
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  
  bool fire = false;
  if( trigResults.failedToGet() ) return false;
  for (unsigned int ip=0; ip<HLTPath.size(); ip++){
    int N_Triggers = trigResults->size();
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);    
    //cout << "new" << endl;
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (!trigResults->accept(i_Trig))   continue;
      const std::string & TrigPath = trigName.triggerName(i_Trig); 
      if (TrigPath.find(HLTPath[ip]) != std::string::npos) fire = true;
    } 
  }
  return fire;
}

std::array<float,5> LeptonSkimmingPhi::hltObject(const edm::Event& iEvent, const edm::EventSetup& iSetup,std::vector<string> Seed ){
  using namespace std;  using namespace edm;  using namespace reco;
  using namespace trigger;
 
  edm::Handle<trigger::TriggerEvent> triggerObjectsSummary;
  iEvent.getByToken(trigobjectsToken_ ,triggerObjectsSummary);
  trigger::TriggerObjectCollection selectedObjects;

  std::vector<std::array<float,5>> max_per_trigger; 

  for (unsigned int ipath=0; ipath<Seed.size(); ipath++){ 
    std::vector<std::array<float, 5>> tot_tr_obj_pt_eta_phi;
    if (!triggerObjectsSummary.isValid()) continue;
    size_t filterIndex = (*triggerObjectsSummary).filterIndex(InputTag(Seed[ipath],"","HLT"));
    trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();     
    if (filterIndex < (*triggerObjectsSummary).sizeFilters()) { 
      const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
      for (size_t j = 0; j < keys.size(); j++) {
        const trigger::TriggerObject & foundObject = (allTriggerObjects)[keys[j]];
        std::array<float, 5> tr_obj_pt_eta_phi;
        if (fabs(foundObject.id())!=13) continue;
        tr_obj_pt_eta_phi[0]=foundObject.pt();
        tr_obj_pt_eta_phi[1]=foundObject.eta();
        tr_obj_pt_eta_phi[2]=foundObject.phi();
        tr_obj_pt_eta_phi[3]=foundObject.id()/fabs(foundObject.id());
        tot_tr_obj_pt_eta_phi.push_back(tr_obj_pt_eta_phi);
      }  
    }
    //take the max per hlt
    if (!tot_tr_obj_pt_eta_phi.empty()){
      std::sort(tot_tr_obj_pt_eta_phi.begin(),tot_tr_obj_pt_eta_phi.end(),
               [](const std::array<float, 5>& a, const std::array<float, 5>& b) {
               return a[0] > b[0];
                });
       max_per_trigger.push_back(tot_tr_obj_pt_eta_phi.at(0));   }
    } 
//we know that at least a trigger fired
//find the total max
  std::sort( max_per_trigger.begin(), max_per_trigger.end(),
         [](const std::array<float,5>& a, const std::array<float,5>& b) {
               return a[0] > b[0];
                });
  return  max_per_trigger.at(0);

}




// ------------ method called on each new Event  ------------
bool
LeptonSkimmingPhi::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
  // using namespace PhysicsTools;

  test_ev++;
  //Get a few collections to apply basic electron ID
  //Get electrons
  edm::Handle<std::vector<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronsToken_, electrons);

 
  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_,muons);

  //Get conversions
  edm::Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(conversionsToken_, conversions);    
  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  
  //Get vertices 
  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  //continue if there are no vertices
  if (vertices->empty()) return false;
  edm::Handle<vector<reco::Track>> tracks;
  iEvent.getByToken(Tracks_, tracks);
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  KalmanVertexFitter theKalmanFitter(false);
  TransientVertex LLvertex;
  
 
  // trigger1=0; trigger2=0; trigger3=0; trigger4=0; trigger5=0; trigger6=0;
  nmuons=0; nel=0; ntracks=0;
    
  SelectedMu_index=-1;
  SelectedMu_DR=std::numeric_limits<float>::max(); muon_pt.clear(); muon_eta.clear(); muon_phi.clear();
  Result=false;  el_pt.clear(); el_eta.clear(); el_phi.clear();
  Trk_container.clear();  MuTracks.clear();  ElTracks.clear(); 
  object_container.clear(); object_id.clear();  cleanedTracks.clear(); 
  Epair_ObjectId.clear(); muon_soft.clear(); muon_medium.clear(); muon_tight.clear();
  Epair_ObjectIndex.clear(); cleanedObjTracks.clear(); cleanedPairTracks.clear();
  Epair_ObjectIndex.clear(); Epair_ObjectId.clear(); Epair_TrkIndex.clear();
  //internal stuff
  ZvertexTrg=-1*std::numeric_limits<float>::max();
  std::array<float,5> SelectedTrgObj_PtEtaPhiCharge{ {-999,-999,-999,-999,-999}} ;
  vertex_point.SetCoordinates(-1*std::numeric_limits<float>::max(),-1*std::numeric_limits<float>::max(),-1*std::numeric_limits<float>::max());
  for (const reco::Vertex& vtx : *vertices) {
    bool isFake = vtx.isFake();
    if (isFake) continue;
    vertex_point.SetCoordinates(vtx.x(),vtx.y(),vtx.z());
    if (vertex_point.x()!=-1*std::numeric_limits<float>::max()) break;
  }
  if (vertex_point.x()==-1*std::numeric_limits<float>::max())
    return false;

  beam_x= theBeamSpot->x0(); beam_y= theBeamSpot->y0(); beam_z= theBeamSpot->z0();
  if(!hltFired(iEvent,iSetup,HLTPath_)) 
    return false;   
      
  SelectedTrgObj_PtEtaPhiCharge=hltObject(iEvent,iSetup,HLTFilter_);
//  std::cout<<"pt "<<SelectedTrgObj_PtEtaPhiCharge[0]<<" eta "<<SelectedTrgObj_PtEtaPhiCharge[1]<<" phi "<<SelectedTrgObj_PtEtaPhiCharge[2]<<endl; 
  SelectedMu_DR=std::numeric_limits<float>::max(); 
  MuTracks.clear();  object_container.clear(); object_id.clear(); nmuons=0;
  for (const reco::Muon  &mu: *muons){
    if (fabs(mu.eta())>EtaTrack_Cut) continue;
    bool tight=false,soft=false;
    if(vertices.isValid()){
      tight=isTightMuonCustom(*(&mu),(*vertices)[0]);
      soft=muon::isSoftMuon(*(&mu),(*vertices)[0]);
    }    
    const Track * mutrack= mu.bestTrack();
    muon_medium.push_back(isMediumMuonCustom(*(&mu)));    
    muon_tight.push_back(tight); muon_soft.push_back(soft);
    muon_pt.push_back(mu.pt()); muon_eta.push_back(mu.eta()); muon_phi.push_back(mu.phi()); 
    auto muTrack=std::make_shared<reco::Track>(*mutrack);
    MuTracks.push_back(muTrack);
    object_container.push_back(nmuons);
    object_id.push_back(13); 
   
    if ( deltaR(mu.eta(),mu.phi(), SelectedTrgObj_PtEtaPhiCharge[1], SelectedTrgObj_PtEtaPhiCharge[2])<MuTrgMatchCone &&  SelectedMu_DR>deltaR(mu.eta(),mu.phi(), SelectedTrgObj_PtEtaPhiCharge[1], SelectedTrgObj_PtEtaPhiCharge[2]) ){
      SelectedMu_DR=deltaR(mu.eta(),mu.phi(), SelectedTrgObj_PtEtaPhiCharge[1], SelectedTrgObj_PtEtaPhiCharge[2]);
      ZvertexTrg=mu.vz();}
    nmuons++;
    //delete mutrack;
  }
  
  if (SelectedMu_DR==std::numeric_limits<float>::max() && SkipIfNoMuMatch){
      return false;
   }

  ElTracks.clear();
  for(const reco::GsfElectron &el : *electrons){  
    bool passConvVeto = !ConversionTools::hasMatchedConversion(*(&el), conversions, theBeamSpot->position());
    if (!passConvVeto) continue;
    if (fabs(el.eta())>EtaTrack_Cut) continue;
    if (el.pt()<PtEl_Cut) continue;
    const Track * eltrack= el.bestTrack();
    auto ElTrack=std::make_shared<reco::Track>(*eltrack);
    ElTracks.push_back(ElTrack); object_container.push_back(nel);
    el_pt.push_back(el.pt()); el_eta.push_back(el.eta()); el_phi.push_back(el.phi());
    //     cout<<"el "<<nel<<" pt  "<<el.pt()<<endl;
    nel++; object_id.push_back(11);
  }
//  cout<<nmuons<<" el "<<nel<<endl;
  cleanedTracks.clear(); 
  trk_index=0;
  
  for (const reco::Track& trk : *tracks){
    if (!trk.quality(Track::highPurity)) continue;
    if (trk.pt()<PtTrack_Cut) continue;
    if (fabs(trk.eta())>EtaTrack_Cut) continue;
    if(trk.charge()==0) continue;
    if(trk.normalizedChi2()>MaxChi2Track_Cut || trk.normalizedChi2()<MinChi2Track_Cut) continue;
    if (fabs(trk.dxy())/trk.dxyError()<TrackSdxy_Cut) continue;
    double minDR=1000;
    for (const reco::Muon& mu: *muons){
      double tempDR=deltaR(mu.eta(),mu.phi(),trk.eta(),trk.phi());
      if (minDR<tempDR) continue;
      minDR=tempDR;
    }
    if (minDR<MuTrkMinDR_Cut) continue;
    if (SelectedMu_DR<std::numeric_limits<float>::max() ){
      if (fabs(ZvertexTrg-trk.vz())>TrackMuDz_Cut ) continue;
      if ( deltaR(trk.eta(),trk.phi(),SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2])<TrgExclusionCone) continue;
    }
    //assignments   
    auto cleanTrack=std::make_shared<reco::Track>(trk);
    cleanedTracks.push_back(cleanTrack);
    Trk_container.push_back(trk_index);
    trk_index++;
  }
 
  //create mother ee combination
 
  // fit track pairs
  cleanedObjTracks.clear(); cleanedPairTracks.clear();     
  TLorentzVector vel1,vel2;
  std::vector<std::shared_ptr<reco::Track>> cleanedObjects; 
  //objects    
  for(auto & vec: MuTracks) cleanedObjects.push_back(vec);
  for(auto & vec: ElTracks) cleanedObjects.push_back(vec);  

  if (cleanedObjects.empty())
    return false;

  for(auto& obj: cleanedObjects){
    //auto obj=cleanedObjects.at(iobj);
    auto tranobj=std::make_shared<reco::TransientTrack>(reco::TransientTrack(*obj,&(*bFieldHandle)));
    unsigned int iobj=&obj-&cleanedObjects[0];
    unsigned int index=object_container.at(iobj);
    if ( object_id.at(iobj)==13 && QualMu_Cut==1 && !muon_soft.at(index)) continue;
    if ( object_id.at(iobj)==13 && QualMu_Cut==2 && !muon_medium.at(index)) continue;
    if ( object_id.at(iobj)==13 && QualMu_Cut==3 && !muon_tight.at(index)) continue;
    if (object_id.at(iobj)==13) vel1.SetPtEtaPhiM(muon_pt.at(index),muon_eta.at(index),muon_phi.at(index),0.0005);
       
    else vel1.SetPtEtaPhiM(el_pt.at(index),el_eta.at(index),el_phi.at(index),0.0005);
    if (object_id.at(iobj)==13 && vel1.Pt()<PtMu_Cut) continue;
    for(auto& trk2: cleanedTracks){
      unsigned int itrk2=&trk2-&cleanedTracks[0];
      //auto trk2=cleanedTracks.at(itrk2);
      if (obj->charge()*trk2->charge()==1) continue;
      if (ObjPtLargerThanTrack && vel1.Pt()<trk2->pt()) continue;
      vel2.SetPtEtaPhiM(trk2->pt(),trk2->eta(),trk2->phi(),0.0005);
      if (object_id.at(iobj)==13 && deltaR(obj->eta(),obj->phi(),SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2])<MuTrgExclusionCone) continue;
      if (object_id.at(iobj)==11 && deltaR(obj->eta(),obj->phi(),SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2])<ElTrgExclusionCone) continue;
      if (SelectedMu_DR<std::numeric_limits<float>::max() ){
        if (object_id.at(iobj)==13 && fabs(ZvertexTrg- obj->vz())>MuTrgMuDz_Cut ) continue;
        if (object_id.at(iobj)==11 && fabs(ZvertexTrg-obj->vz())>ElTrgMuDz_Cut ) continue;
      }
      if ((vel1+vel2).M()>MaxMee_Cut || (vel1+vel2).M()<MinMee_Cut ) continue;        
      auto trantrk2=std::make_shared<reco::TransientTrack>(reco::TransientTrack(*trk2,&(*bFieldHandle)));
      std::vector<reco::TransientTrack> tempTracks;
      tempTracks.reserve(2);
      tempTracks.push_back(*tranobj); tempTracks.push_back(*trantrk2);
      LLvertex = theKalmanFitter.vertex(tempTracks);
      if (!LLvertex.isValid()) continue;
      if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom())<Probee_Cut)  continue;
      if (ZvertexTrg>-1*std::numeric_limits<float>::max() && fabs(ZvertexTrg-LLvertex.position().z())>EpairZvtx_Cut ) continue;
      GlobalError err =LLvertex.positionError();
      GlobalPoint Dispbeamspot(-1*((theBeamSpot->x0()-LLvertex.position().x())+(LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),-1*((theBeamSpot->y0()-LLvertex.position().y())+ (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);
      math::XYZVector pperp((vel1+vel2).Px(),(vel1+vel2).Py(),0);
      math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);
      float tempCos=vperp.Dot(pperp)/(vperp.R()*pperp.R());
      if (tempCos<Cosee_Cut) continue;
      cleanedObjTracks.push_back(obj);
      cleanedPairTracks.push_back(trk2);
      Epair_ObjectIndex.push_back(object_container.at(iobj));
      Epair_ObjectId.push_back(object_id.at(iobj));
      Epair_TrkIndex.push_back(Trk_container.at(itrk2));
    }
  }
//    cout<<trk_index<<" comb "<<Epair_ObjectId.size()<<endl;   
    
  // Bs recontrsuvtion
  TLorentzVector vKlead,vKsublead; 
  for(unsigned int iobj=0; iobj<cleanedObjTracks.size(); iobj++){
    auto objtrk=cleanedObjTracks.at(iobj);
    auto pairtrk=cleanedPairTracks.at(iobj);
    auto tranobj=std::make_shared<reco::TransientTrack>(reco::TransientTrack(*objtrk,&(*bFieldHandle)));
    auto tranpair=std::make_shared<reco::TransientTrack>(reco::TransientTrack(*pairtrk,&(*bFieldHandle)));
    unsigned int index=Epair_ObjectIndex.at(iobj);
    if ( Epair_ObjectId.at(iobj)==13) vel1.SetPtEtaPhiM(muon_pt.at(index),muon_eta.at(index),muon_phi.at(index),0.0005);
      
    else vel1.SetPtEtaPhiM(el_pt.at(index),el_eta.at(index),el_phi.at(index),0.0005);
       
    for (auto itrk1 = cleanedTracks.begin(); itrk1 != cleanedTracks.end(); ++itrk1) {
      auto & trk1 = *itrk1;
      if (deltaR(objtrk->eta(),objtrk->phi(),trk1->eta(),trk1->phi())<TrkObjExclusionCone) continue;
      if (fabs(trk1->dxy(vertex_point))/trk1->dxyError()<Ksdxy_Cut) continue;
      if (deltaR(pairtrk->eta(),pairtrk->phi(),trk1->eta(),trk1->phi())<TrkTrkMinDR_Cut) continue;
        for (auto itrk2 = itrk1 + 1; itrk2 != cleanedTracks.end(); ++itrk2) {
	auto & trk2 = *itrk2;
	if (deltaR(objtrk->eta(),objtrk->phi(),trk2->eta(),trk2->phi())<TrkObjExclusionCone) continue;
        if (fabs(trk2->dxy(vertex_point))/trk2->dxyError()<Ksdxy_Cut) continue;
        if (deltaR(pairtrk->eta(),pairtrk->phi(),trk2->eta(),trk2->phi())<TrkTrkMinDR_Cut) continue;
	if (trk1->charge() == trk2->charge()) continue;

	auto leadtrk = trk1->pt() > trk2->pt() ? trk1 : trk2;
	auto subleadtrk = trk1->pt() > trk2->pt() ? trk2 : trk1;

	vel2.SetPtEtaPhiM(pairtrk->pt(),pairtrk->eta(),pairtrk->phi(),0.0005);
	vKlead.SetPtEtaPhiM(leadtrk->pt(),leadtrk->eta(),leadtrk->phi(),0.493);
	vKsublead.SetPtEtaPhiM(subleadtrk->pt(),subleadtrk->eta(),subleadtrk->phi(),0.493);

	if (vKlead.Pt() < PtKleadTrack_Cut) continue;
        if (vKsublead.Pt() < PtKsubleadTrack_Cut) continue;

	if ((vKlead+vKsublead).M() > MaxMphi_Cut || (vKlead+vKsublead).M() < MinMphi_Cut) continue;

	if ((vel1+vel2+vKlead+vKsublead).M()> MaxMB_Cut || (vel1+vel2+vKlead+vKsublead).M()< MinMB_Cut) continue;
	if ((vel1+vel2+vKlead+vKsublead).Pt()<PtB_Cut) continue;
	auto trantrklead=std::make_shared<reco::TransientTrack>(reco::TransientTrack(*leadtrk,&(*bFieldHandle)));
	auto trantrksublead=std::make_shared<reco::TransientTrack>(reco::TransientTrack(*subleadtrk,&(*bFieldHandle)));

	std::vector<reco::TransientTrack> tempTracks;
	tempTracks.reserve(4);
	tempTracks.push_back(*tranobj); tempTracks.push_back(*tranpair);
	tempTracks.push_back(*trantrklead); tempTracks.push_back(*trantrksublead);
	 
	LLvertex = theKalmanFitter.vertex(tempTracks);
	if (!LLvertex.isValid()) continue;
	if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom())<ProbeeKK_Cut) continue;
	GlobalError err = LLvertex.positionError();
	GlobalPoint Dispbeamspot(-1*((theBeamSpot->x0()-LLvertex.position().x())+(LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),-1*((theBeamSpot->y0()-LLvertex.position().y())+ (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);
      
	math::XYZVector pperp((vel1+vel2+vKlead+vKsublead).Px(),(vel1+vel2+vKlead+vKsublead).Py(),0);
	math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);
	float tempCos=vperp.Dot(pperp)/(vperp.R()*pperp.R());
	if (tempCos<CoseeKK_Cut) continue;
	if (SLxy_Cut>Dispbeamspot.perp()/sqrt(err.rerr(Dispbeamspot))) continue;
	Result=true;
//        std::cout<<"fired "<<test_ev<<std::endl;
	break;
     }  
   }
   if (Result) break;
 }
  
 return Result;
   
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
LeptonSkimmingPhi::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
LeptonSkimmingPhi::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
LeptonSkimmingPhi::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
LeptonSkimmingPhi::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
LeptonSkimmingPhi::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
LeptonSkimmingPhi::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonSkimmingPhi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(LeptonSkimmingPhi);

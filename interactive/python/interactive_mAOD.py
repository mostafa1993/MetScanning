''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
from math import cos, sin, atan2, sqrt, pi

#RootTools
from RootTools.core.standard import *

def toDict(p):
    return {'pt':p.pt(), 'eta':p.eta(), 'phi':p.phi(), 'pdgId':p.pdgId()}

def bold(s):
    return '\033[1m'+s+'\033[0m'

#s0 = FWLiteSample.fromFiles("evangelos", files = ["./data/pickEvents_evangelos_265523163.root"] )

#s0 = FWLiteSample.fromFiles("bobak", files = ["./data/bobak_pickevents.root"] )
#products = {
#    'pfJets':{'type':'vector<reco::PFJet>', 'label':( "ak4PFJets" ) },
#    'gen':{'type':'vector<reco::GenParticle>', 'label':( "genParticles" ) },
#    'pfMet':{'type':'vector<reco::PFMET>', 'label':( "pfMet" )},
#    'pf':{'type':'vector<reco::PFCandidate>', 'label':( "particleFlow" ) },
#    'muons':{'type':'vector<reco::Muon>', 'label': ("muons") }
#}

#s0 = FWLiteSample.fromFiles("bobak", files = ["./data/bobak_pickevents_miniAOD.root"] )
#s0 = FWLiteSample.fromFiles("bobak", files = ["./data/pickevents_maria_badMuon.root"] )
s0 = FWLiteSample.fromFiles("bobak", files = ["./data/maria_badevents_mAOD.root "] )
#s0 = FWLiteSample.fromFiles("bobak", files = ["./data/raman_pickevents_BadMuFilterFailed.root "] )
#s0 = FWLiteSample.fromFiles("bobak", files = ["root://eoscms.cern.ch//store/data/Run2016B/JetHT/MINIAOD/PromptReco-v2/000/273/150/00000/66051AAF-D819-E611-BD3D-02163E011D55.root"] )
#s0 = FWLiteSample.fromFiles("T1tttt", files = ["root://eoscms.cern.ch//store/mc/RunIISpring16MiniAODv1/SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/60000/26C76166-0FFE-E511-BA96-0025905D1D60.root"] )
#s0 = FWLiteSample.fromFiles("ece", files = ["/afs/cern.ch/user/e/easilar/public/forMETScanners/pickevents_basChargedFail_10TeVPfMET.root"] )

products = {
    'muons':{'type':'vector<pat::Muon>', 'label': "slimmedMuons"},
    'pf':{'type':'vector<pat::PackedCandidate>', 'label':( "packedPFCandidates" ) },
#    'pf':{'type':'vector<reco::PFCandidate>', 'label':( "particleFlow" ) },
#    'muons':{'type':'vector<reco::Muon>', 'label': ("muons") }
    'pfMet':{'type':'vector<pat::MET>', 'label':( "slimmedMETs" )},
}

r = s0.fwliteReader( products = products )
r.start()

def printHisto( p, prefix = ""):
    print "%s gen %10s status %i nMothers %i pt/eta/phi %5.1f/%5.1f/%5.1f" % ( prefix, \
        pdgToName(p.pdgId()), p.status(), p.numberOfMothers(), p.pt(), p.eta(), p.phi(),
    )
    for i in range( p.numberOfMothers() ):
        printHisto( p.mother(i), "  " + prefix )

def deltaPhi(phi1, phi2):
    dphi = phi2-phi1
    if  dphi > pi:
        dphi -= 2.0*pi
    if dphi <= -pi:
        dphi += 2.0*pi
    return abs(dphi)

def deltaR(eta1, phi1, eta2, phi2):
    return sqrt( (eta1-eta2)**2 + deltaPhi(phi1, phi2)**2 )
 

def badChargedHadronFilter( muons, pf, discardHPTracks = True):
    maxDR = 0.001
    minMuonTrackRelErr = 0.5
    minPtDiffRel = -0.5
    minMuPt = 100
    flagged = False
    for muon in muons:
        if muon.pt()<minMuPt : continue
        if muon.innerTrack().isNonnull():
            it = muon.innerTrack()
            if it.quality(it.highPurity): 
                if discardHPTracks:continue
            # All events had a drastically high pt error on the inner muon track (fac. ~10). Require at least 0.5
            if it.ptError()/it.pt() < minMuonTrackRelErr: continue
            for c in pf:
                if abs(c.pdgId()) == 211:
                    # Require very loose similarity in pt (one-sided). 
                    dPtRel =  ( c.pt() - it.pt() )/(0.5*(c.pt() + it.pt()))
                    # Flag the event bad if dR is tiny
                    if deltaR( it.eta(), it.phi(), c.eta(), c.phi() ) < maxDR and dPtRel > minPtDiffRel:
                        flagged = True
                        break
    return not flagged


def badPFMuonFilter( pf ):
    minMuPt = 100
    minMuPtError = -1
    minDZ = 0.1
    pfMuons = filter(lambda p: abs(p.pdgId())==13 and p.pt()>minMuPt, pf)
    if len(pfMuons)==0: return True
    for pfMu in pfMuons:
        tref = pfMu.muonRef().innerTrack()
        if tref.isNull() or tref.quality(tref.highPurity):
            continue 
        if abs(tref.dz()) < minDZ: continue
        if tref.ptError()<minMuPtError: continue
        #print "PFMu NON HP IT pdgId %i pt %3.2f eta %3.2f phi %3.2f track pt %3.2f +/- %3.2f dxy %3.2f dz %3.2f" %( pfMu.pdgId(), pfMu.pt(), pfMu.eta(), pfMu.phi(), tref.pt(), tref.ptError(), tref.dxy(), tref.dz() )
        return False
    return True

def badPFMuonFilter2( muons, pf, verbose = False, discardHPTracks = False):
    minMuPt = 100
    minMuPtRelError = 0.3
    maxDR = 0.001
   
    #minDZ = 0.1
    if len(muons)==0: return True
    for i_muon, muon in enumerate(muons):
        #if verbose: print "Testing ", i_muon
        it = muon.innerTrack()
        if it.isNull(): continue
        if discardHPTracks and ( it.quality(it.highPurity) ):
            if verbose: print "muon pt %3.2f Inner track is HP -> continue"%muon.pt()
            continue 

        reco_muon_is_bad = ( (it.pt()> minMuPt) and it.ptError()/it.pt() > minMuPtRelError ) 
        gt = muon.globalTrack()
        if (not reco_muon_is_bad) and ( (not gt.isNull()) and (gt.pt()> minMuPt) and gt.ptError()/gt.pt()>minMuPtRelError ):
            reco_muon_is_bad = True
        
        if verbose: 
            print "reco::muon bad?  %r IT-HP? %r pdgId %i pt %3.2f eta %3.2f phi %3.2f track pt %3.2f +/- %3.2f " %( reco_muon_is_bad, it.quality(it.highPurity), muon.pdgId(), muon.pt(), muon.eta(), muon.phi(), it.pt(), it.ptError() )
            if not gt.isNull(): print "global pt %3.2f +/- %3.2f"% (gt.pt(), gt.ptError())

        pfMuons = filter(lambda p: abs(p.pdgId())==13 and p.pt()>minMuPt, pf )

        if len(pfMuons)>0:
            minDR = min([deltaR( muon.eta(), muon.phi(), c.eta(), c.phi() ) for c in pfMuons] )
            if verbose: print "minDR", minDR
            flagged = minDR<maxDR
        else:
            flagged = False
            if verbose: print "no PF muon found"
        if flagged: return False

    return True

verbose = False

c_total = 0
c_inf = 0
c_badChargedHadronFilter = 0
c_badChargedHadronFilter_HP = 0
c_badPFMuonFilter = 0
c_badPFMuonFilter_HP = 0
while r.run():
    if 23672899 in r.evt: break

    c_total += 1
    pf = list(r.products['pf'])
    pfMuons = filter(lambda p: abs(p.pdgId())==13, pf)
    #pf.sort(key = lambda p: -p.pt())
    muons = list(r.products['muons'])
    #if len(muons)>0:
    #    muon = muons[0]
    #    print "reco::muon pdgId %i pt %3.2f eta %3.2f phi %3.2f inner track pt %3.2f +/- %3.2f. Flagged? %r" %( muon.pdgId(), muon.pt(), muon.eta(), muon.phi(), muon.innerTrack().pt(), muon.innerTrack().ptError(), badChargedHadronFilter(muons, pf) )
    # print "badPFMuonFilter? %r badPFMuonFilter2 %r"%( badPFMuonFilter( pf ), badPFMuonFilter2( muons) )
    print "%30s"%("%i:%i:%i"% r.evt), 
    if not badChargedHadronFilter(muons, pf, discardHPTracks = True): 
        #print "Failed badChargedHadronFilter"
        c_badChargedHadronFilter += 1
        print "badCh: 0",
    else:
        print "badCh: 1",
      
    if not badChargedHadronFilter(muons, pf, discardHPTracks = False):
        # print "Failed badChargedHadronFilter if HP tracks are also considered!"
        c_badChargedHadronFilter_HP += 1
        print "badCh(HP): 0",
    else:
        print "badCh(HP): 1",

    if not badPFMuonFilter2( muons, pf, verbose = verbose, discardHPTracks = True):
        c_badPFMuonFilter += 1
        # print "Flagged %i:%i:%i"% r.evt,"MET %3.2f"%r.products['pfMet'][0].pt()
        print "bad PF Mu: 0",
    else:
        print "bad PF Mu: 1",
    if not badPFMuonFilter2( muons, pf, verbose = verbose, discardHPTracks = False): 
        c_badPFMuonFilter_HP += 1
        print "bad PF Mu (HP): 0",
    else:
        print "bad PF Mu (HP): 1",
        # print "Failed badPFMuonFilter if HP tracks are also considered!"
    #else:
    #    if verbose: print "Not Flagged %i:%i:%i"% r.evt,"MET %3.2f"%r.products['pfMet'][0].pt()
    pfInf = filter(lambda p: not p.pt()<float('Inf'), pf)
    if len(pfInf)>0:
        #for p in pfInf:
            # print "Found INF opt particle pdgId %i eta %3.2f phi %3.2f" % (p.pdgId(), p.eta(), p.phi())
        c_inf += 1
        print "Inf Cand: 0"
    else:
        print "Inf Cand: 1"

print "total %i inf-cand %i badChargedHadronFilter %i badChargedHadronFilter_HP %i badPFMuonFilter %i badPFMuonFilter_HP %i"\
        %( c_total, c_inf, c_badChargedHadronFilter, c_badChargedHadronFilter_HP, c_badPFMuonFilter, c_badPFMuonFilter_HP) 
         
    

#Type                                  Module                      Label             Process   
#----------------------------------------------------------------------------------------------
#LHEEventProduct                       "externalLHEProducer"       ""                "LHE"     
#GenEventInfoProduct                   "generator"                 ""                "SIM"     
#edm::TriggerResults                   "TriggerResults"            ""                "SIM"     
#edm::TriggerResults                   "TriggerResults"            ""                "HLT"     
#int                                   "addPileupInfo"             "bunchSpacing"    "HLT"     
#vector<PileupSummaryInfo>             "addPileupInfo"             ""                "HLT"     
#vector<int>                           "genParticles"              ""                "HLT"     
#vector<reco::GenJet>                  "ak4GenJets"                ""                "HLT"     
#vector<reco::GenJet>                  "ak4GenJetsNoNu"            ""                "HLT"     
#vector<reco::GenJet>                  "ak8GenJets"                ""                "HLT"     
#vector<reco::GenJet>                  "ak8GenJetsNoNu"            ""                "HLT"     
#vector<reco::GenMET>                  "genMetCalo"                ""                "HLT"     
#vector<reco::GenMET>                  "genMetTrue"                ""                "HLT"     
#vector<reco::GenParticle>             "genParticles"              ""                "HLT"     
#trigger::TriggerEvent                 "hltTriggerSummaryAOD"      ""                "HLT"     
#ClusterSummary                        "clusterSummaryProducer"    ""                "RECO"    
#EBDigiCollection                      "selectDigi"                "selectedEcalEBDigiCollection"   "RECO"    
#EEDigiCollection                      "selectDigi"                "selectedEcalEEDigiCollection"   "RECO"    
#HcalNoiseSummary                      "hcalnoise"                 ""                "RECO"    
#HcalUnpackerReport                    "castorDigis"               ""                "RECO"    
#HcalUnpackerReport                    "hcalDigis"                 ""                "RECO"    
#L1GlobalTriggerObjectMaps             "l1L1GtObjectMap"           ""                "RECO"    
#L1GlobalTriggerReadoutRecord          "gtDigis"                   ""                "RECO"    
#double                                "fixedGridRhoAll"           ""                "RECO"    
#double                                "fixedGridRhoFastjetAll"    ""                "RECO"    
#double                                "fixedGridRhoFastjetAllCalo"   ""                "RECO"    
#double                                "fixedGridRhoFastjetAllTmp"   ""                "RECO"    
#double                                "fixedGridRhoFastjetCentral"   ""                "RECO"    
#double                                "fixedGridRhoFastjetCentralCalo"   ""                "RECO"    
#double                                "fixedGridRhoFastjetCentralChargedPileUp"   ""                "RECO"    
#double                                "fixedGridRhoFastjetCentralNeutral"   ""                "RECO"    
#double                                "ak4CaloJets"               "rho"             "RECO"    
#double                                "ak4PFJets"                 "rho"             "RECO"    
#double                                "ak4PFJetsCHS"              "rho"             "RECO"    
#double                                "ak4TrackJets"              "rho"             "RECO"    
#double                                "ak5CastorJets"             "rho"             "RECO"    
#double                                "ak7CastorJets"             "rho"             "RECO"    
#double                                "ak8PFJetsCHS"              "rho"             "RECO"    
#double                                "ak8PFJetsCHSSoftDrop"      "rho"             "RECO"    
#double                                "cmsTopTagPFJetsCHS"        "rho"             "RECO"    
#double                                "ak4CaloJets"               "sigma"           "RECO"    
#double                                "ak4PFJets"                 "sigma"           "RECO"    
#double                                "ak4PFJetsCHS"              "sigma"           "RECO"    
#double                                "ak4TrackJets"              "sigma"           "RECO"    
#double                                "ak5CastorJets"             "sigma"           "RECO"    
#double                                "ak7CastorJets"             "sigma"           "RECO"    
#double                                "ak8PFJetsCHS"              "sigma"           "RECO"    
#double                                "ak8PFJetsCHSSoftDrop"      "sigma"           "RECO"    
#double                                "cmsTopTagPFJetsCHS"        "sigma"           "RECO"    
#edm::Association<vector<reco::DeDxHitInfo> >    "dedxHitInfo"               ""                "RECO"    
#edm::Association<vector<reco::Vertex> >    "offlinePrimaryVertices"    ""                "RECO"    
#edm::Association<vector<reco::Vertex> >    "offlinePrimaryVerticesWithBS"   ""                "RECO"    
#edm::AssociationMap<edm::OneToOne<vector<reco::SuperCluster>,vector<reco::HFEMClusterShape>,unsigned int> >    "hfEMClusters"              ""                "RECO"    
#edm::AssociationMap<edm::OneToOne<vector<reco::Track>,vector<reco::Track>,unsigned int> >    "tevMuons"                  "default"         "RECO"    
#edm::AssociationMap<edm::OneToOne<vector<reco::Track>,vector<reco::Track>,unsigned int> >    "tevMuons"                  "dyt"             "RECO"    
#edm::AssociationMap<edm::OneToOne<vector<reco::Track>,vector<reco::Track>,unsigned int> >    "tevMuons"                  "firstHit"        "RECO"    
#edm::AssociationMap<edm::OneToOne<vector<reco::Track>,vector<reco::Track>,unsigned int> >    "tevMuons"                  "picky"           "RECO"    
#edm::AssociationVector<edm::RefProd<vector<reco::PFTau> >,vector<edm::Ref<vector<reco::PFTauTransverseImpactParameter>,reco::PFTauTransverseImpactParameter,edm::refhelper::FindUsingAdvance<vector<reco::PFTauTransverseImpactParameter>,reco::PFTauTransverseImpactParameter> > >,edm::Ref<vector<reco::PFTau>,reco::PFTau,edm::refhelper::FindUsingAdvance<vector<reco::PFTau>,reco::PFTau> >,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "hpsPFTauTransverseImpactParameters"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<edm::RefVector<vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<vector<reco::Track>,reco::Track> > >,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "ak4JetTracksAssociatorAtVertex"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<edm::RefVector<vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<vector<reco::Track>,reco::Track> > >,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "ak4JetTracksAssociatorAtVertexPF"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<edm::RefVector<vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<vector<reco::Track>,reco::Track> > >,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "ak4JetTracksAssociatorExplicit"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedCvsBJetTags"     ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedCvsLJetTags"     ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedInclusiveSecondaryVertexV2BJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedInclusiveSecondaryVertexV2BJetTagsEI"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedMVABJetTags"     ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedMVAV2BJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedSecondaryVertexBJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedSecondaryVertexSoftLeptonBJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedSecondaryVertexV2BJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfJetBProbabilityBJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfJetProbabilityBJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfSimpleSecondaryVertexHighEffBJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfSimpleSecondaryVertexHighPurBJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfTrackCountingHighEffBJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfTrackCountingHighPurBJetTags"   ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "softPFElectronBJetTags"    ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "softPFMuonBJetTags"        ""                "RECO"    
#edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<reco::JetExtendedAssociation::JetExtendedData>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "ak4JetExtender"            ""                "RECO"    
#edm::ConditionsInEventBlock           "conditionsInEdm"           ""                "RECO"    
#edm::OwnVector<TrackingRecHit,edm::ClonePolicy<TrackingRecHit> >    "displacedGlobalMuons"      ""                "RECO"    
#edm::OwnVector<TrackingRecHit,edm::ClonePolicy<TrackingRecHit> >    "displacedStandAloneMuons"   ""                "RECO"    
#edm::OwnVector<TrackingRecHit,edm::ClonePolicy<TrackingRecHit> >    "refittedStandAloneMuons"   ""                "RECO"    
#edm::OwnVector<TrackingRecHit,edm::ClonePolicy<TrackingRecHit> >    "standAloneMuons"           ""                "RECO"    
#edm::RangeMap<CSCDetId,edm::OwnVector<CSCSegment,edm::ClonePolicy<CSCSegment> >,edm::ClonePolicy<CSCSegment> >    "cscSegments"               ""                "RECO"    
#edm::RangeMap<DTChamberId,edm::OwnVector<DTRecSegment4D,edm::ClonePolicy<DTRecSegment4D> >,edm::ClonePolicy<DTRecSegment4D> >    "dt4DCosmicSegments"        ""                "RECO"    
#edm::RangeMap<DTChamberId,edm::OwnVector<DTRecSegment4D,edm::ClonePolicy<DTRecSegment4D> >,edm::ClonePolicy<DTRecSegment4D> >    "dt4DSegments"              ""                "RECO"    
#edm::RangeMap<RPCDetId,edm::OwnVector<RPCRecHit,edm::ClonePolicy<RPCRecHit> >,edm::ClonePolicy<RPCRecHit> >    "rpcRecHits"                ""                "RECO"    
#edm::SortedCollection<CastorRecHit,edm::StrictWeakOrdering<CastorRecHit> >    "castorreco"                ""                "RECO"    
#edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "reducedEcalRecHitsEB"      ""                "RECO"    
#edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "reducedEcalRecHitsEE"      ""                "RECO"    
#edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "reducedEcalRecHitsES"      ""                "RECO"    
#edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >    "reducedHcalRecHits"        "hbhereco"        "RECO"    
#edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> >    "reducedHcalRecHits"        "hfreco"          "RECO"    
#edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> >    "reducedHcalRecHits"        "horeco"          "RECO"    
#edm::TriggerResults                   "TriggerResults"            ""                "RECO"    
#edm::ValueMap<bool>                   "PhotonIDProd"              "PhotonCutBasedIDLoose"   "RECO"    
#edm::ValueMap<bool>                   "PhotonIDProdGED"           "PhotonCutBasedIDLoose"   "RECO"    
#edm::ValueMap<bool>                   "PhotonIDProd"              "PhotonCutBasedIDLooseEM"   "RECO"    
#edm::ValueMap<bool>                   "PhotonIDProdGED"           "PhotonCutBasedIDLooseEM"   "RECO"    
#edm::ValueMap<bool>                   "PhotonIDProd"              "PhotonCutBasedIDTight"   "RECO"    
#edm::ValueMap<bool>                   "PhotonIDProdGED"           "PhotonCutBasedIDTight"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidAllArbitrated"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidGMStaChiCompatibility"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidGMTkChiCompatibility"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidGMTkKinkTight"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidGlobalMuonPromptTight"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidRPCMuLoose"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTM2DCompatibilityLoose"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTM2DCompatibilityTight"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMLastStationAngLoose"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMLastStationAngTight"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMLastStationLoose"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMLastStationOptimizedLowPtLoose"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMLastStationOptimizedLowPtTight"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMLastStationTight"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMOneStationAngLoose"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMOneStationAngTight"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMOneStationLoose"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTMOneStationTight"   "RECO"    
#edm::ValueMap<bool>                   "muons"                     "muidTrackerMuonArbitrated"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueCharged03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueCharged04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueChargedAll03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueChargedAll04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueGamma03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueGamma04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueGammaHighThreshold03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueGammaHighThreshold04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueNeutral03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueNeutral04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueNeutralHighThreshold03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValueNeutralHighThreshold04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValuePU03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFMeanDRIsoValuePU04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueCharged03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueCharged04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueChargedAll03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueChargedAll04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueGamma03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueGamma04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueGammaHighThreshold03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueGammaHighThreshold04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueNeutral03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueNeutral04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueNeutralHighThreshold03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValueNeutralHighThreshold04"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValuePU03"   "RECO"    
#edm::ValueMap<double>                 "muons"                     "muPFSumDRIsoValuePU04"   "RECO"    
#edm::ValueMap<edm::Ptr<reco::PFCandidate> >    "particleFlow"              "electrons"       "RECO"    
#edm::ValueMap<edm::Ptr<reco::PFCandidate> >    "particleFlow"              "muons"           "RECO"    
#edm::ValueMap<edm::Ptr<reco::PFCandidate> >    "particleFlow"              "photons"         "RECO"    
#edm::ValueMap<float>                  "ak8PFJetsCHSSoftDropMass"   ""                "RECO"    
#edm::ValueMap<float>                  "eidLoose"                  ""                "RECO"    
#edm::ValueMap<float>                  "eidRobustHighEnergy"       ""                "RECO"    
#edm::ValueMap<float>                  "eidRobustLoose"            ""                "RECO"    
#edm::ValueMap<float>                  "eidRobustTight"            ""                "RECO"    
#edm::ValueMap<float>                  "eidTight"                  ""                "RECO"    
#edm::ValueMap<float>                  "electronEcalPFClusterIsolationProducer"   ""                "RECO"    
#edm::ValueMap<float>                  "electronHcalPFClusterIsolationProducer"   ""                "RECO"    
#edm::ValueMap<float>                  "offlinePrimaryVertices"    ""                "RECO"    
#edm::ValueMap<float>                  "offlinePrimaryVerticesWithBS"   ""                "RECO"    
#edm::ValueMap<float>                  "photonEcalPFClusterIsolationProducer"   ""                "RECO"    
#edm::ValueMap<float>                  "photonHcalPFClusterIsolationProducer"   ""                "RECO"    
#edm::ValueMap<int>                    "offlinePrimaryVertices"    ""                "RECO"    
#edm::ValueMap<int>                    "offlinePrimaryVerticesWithBS"   ""                "RECO"    
#edm::ValueMap<reco::CastorJetID>      "ak5CastorJetID"            ""                "RECO"    
#edm::ValueMap<reco::CastorJetID>      "ak7CastorJetID"            ""                "RECO"    
#edm::ValueMap<reco::DeDxData>         "dedxHarmonic2"             ""                "RECO"    
#edm::ValueMap<reco::JetID>            "ak4JetID"                  ""                "RECO"    
#edm::ValueMap<reco::MuonCosmicCompatibility>    "muons"                     "cosmicsVeto"     "RECO"    
#edm::ValueMap<reco::MuonMETCorrectionData>    "muonMETValueMapProducer"   "muCorrData"      "RECO"    
#edm::ValueMap<reco::MuonShower>       "muons"                     "muonShowerInformation"   "RECO"    
#edm::ValueMap<reco::MuonTimeExtra>    "muons"                     "combined"        "RECO"    
#edm::ValueMap<reco::MuonTimeExtra>    "muons"                     "csc"             "RECO"    
#edm::ValueMap<reco::MuonTimeExtra>    "muons"                     "dt"              "RECO"    
#edm::ValueMap<vector<edm::Ref<vector<reco::PFCandidate>,reco::PFCandidate,edm::refhelper::FindUsingAdvance<vector<reco::PFCandidate>,reco::PFCandidate> > > >    "particleBasedIsolation"    "gedGsfElectrons"   "RECO"    
#edm::ValueMap<vector<edm::Ref<vector<reco::PFCandidate>,reco::PFCandidate,edm::refhelper::FindUsingAdvance<vector<reco::PFCandidate>,reco::PFCandidate> > > >    "particleBasedIsolation"    "gedPhotons"      "RECO"    
#edm::ValueMap<unsigned int>           "muons"                     "cosmicsVeto"     "RECO"    
#int                                   "tcdsDigis"                 "nibble"          "RECO"    
#reco::BeamHaloSummary                 "BeamHaloSummary"           ""                "RECO"    
#reco::BeamSpot                        "offlineBeamSpot"           ""                "RECO"    
#reco::GlobalHaloData                  "GlobalHaloData"            ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauChargedIsoPtSum"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByDeadECALElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByDecayModeFinding"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByDecayModeFindingNewDMs"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByDecayModeFindingOldDMs"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByIsolationMVA3newDMwLTraw"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByIsolationMVA3newDMwoLTraw"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByIsolationMVA3oldDMwLTraw"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByIsolationMVA3oldDMwoLTraw"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseChargedIsolation"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseIsolation"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseIsolationMVA3newDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseIsolationMVA3newDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseMuonRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseMuonRejection2"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLooseMuonRejection3"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByLoosePileupWeightedIsolation3Hits"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVA5LooseElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVA5MediumElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVA5TightElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVA5VLooseElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVA5VTightElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVA5rawElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVALooseMuonRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVAMediumMuonRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVATightMuonRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVArawMuonRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumChargedIsolation"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumIsolation"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumIsolationMVA3newDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumIsolationMVA3newDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumMuonRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumMuonRejection2"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMediumPileupWeightedIsolation3Hits"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByRawChargedIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByRawGammaIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByRawPileupWeightedIsolation3Hits"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightChargedIsolation"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightElectronRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightIsolation"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightIsolationMVA3newDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightIsolationMVA3newDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightIsolationMVA3oldDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightIsolationMVA3oldDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightMuonRejection"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightMuonRejection2"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightMuonRejection3"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByTightPileupWeightedIsolation3Hits"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVLooseChargedIsolation"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVLooseIsolation"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVTightIsolationMVA3newDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVTightIsolationMVA3newDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwoLT"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauNeutralIsoPtSum"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauPUcorrPtSum"       ""                "RECO"    
#reco::PFTauDiscriminator              "pfTausDiscriminationByDecayModeFinding"   ""                "RECO"    
#reco::PFTauDiscriminator              "pfTausDiscriminationByIsolation"   ""                "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByIsolationMVA3newDMwLTraw"   "category"        "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByIsolationMVA3newDMwoLTraw"   "category"        "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByIsolationMVA3oldDMwLTraw"   "category"        "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByIsolationMVA3oldDMwoLTraw"   "category"        "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVA5rawElectronRejection"   "category"        "RECO"    
#reco::PFTauDiscriminator              "hpsPFTauDiscriminationByMVArawMuonRejection"   "category"        "RECO"    
#vector<BeamSpotOnline>                "scalersRawToDigi"          ""                "RECO"    
#vector<DcsStatus>                     "scalersRawToDigi"          ""                "RECO"    
#vector<L1AcceptBunchCrossing>         "scalersRawToDigi"          ""                "RECO"    
#vector<L1TriggerScalers>              "scalersRawToDigi"          ""                "RECO"    
#vector<Level1TriggerScalers>          "scalersRawToDigi"          ""                "RECO"    
#vector<LumiScalers>                   "scalersRawToDigi"          ""                "RECO"    
#vector<double>                        "ak4PFJetsCHS"              "rhos"            "RECO"    
#vector<double>                        "ak8PFJetsCHS"              "rhos"            "RECO"    
#vector<double>                        "ak8PFJetsCHSSoftDrop"      "rhos"            "RECO"    
#vector<double>                        "cmsTopTagPFJetsCHS"        "rhos"            "RECO"    
#vector<double>                        "ak4PFJetsCHS"              "sigmas"          "RECO"    
#vector<double>                        "ak8PFJetsCHS"              "sigmas"          "RECO"    
#vector<double>                        "ak8PFJetsCHSSoftDrop"      "sigmas"          "RECO"    
#vector<double>                        "cmsTopTagPFJetsCHS"        "sigmas"          "RECO"    
#vector<edm::ErrorSummaryEntry>        "logErrorHarvester"         ""                "RECO"    
#vector<edm::FwdPtr<reco::PFCandidate> >    "particleFlowPtrs"          ""                "RECO"    
#vector<edm::FwdPtr<reco::PFCandidate> >    "particleFlowTmpPtrs"       ""                "RECO"    
#vector<edm::FwdPtr<reco::PFCandidate> >    "pfIsolatedElectronsEI"     ""                "RECO"    
#vector<edm::FwdPtr<reco::PFCandidate> >    "pfIsolatedMuonsEI"         ""                "RECO"    
#vector<float>                         "generalTracks"             "MVAValues"       "RECO"    
#vector<l1extra::L1EmParticle>         "l1extraParticles"          "Isolated"        "RECO"    
#vector<l1extra::L1EmParticle>         "l1extraParticles"          "NonIsolated"     "RECO"    
#vector<l1extra::L1EtMissParticle>     "l1extraParticles"          "MET"             "RECO"    
#vector<l1extra::L1EtMissParticle>     "l1extraParticles"          "MHT"             "RECO"    
#vector<l1extra::L1HFRings>            "l1extraParticles"          ""                "RECO"    
#vector<l1extra::L1JetParticle>        "l1extraParticles"          "Central"         "RECO"    
#vector<l1extra::L1JetParticle>        "l1extraParticles"          "Forward"         "RECO"    
#vector<l1extra::L1JetParticle>        "l1extraParticles"          "IsoTau"          "RECO"    
#vector<l1extra::L1JetParticle>        "l1extraParticles"          "Tau"             "RECO"    
#vector<l1extra::L1MuonParticle>       "l1extraParticles"          ""                "RECO"    
#vector<reco::BasicJet>                "ak5CastorJets"             ""                "RECO"    
#vector<reco::BasicJet>                "ak7CastorJets"             ""                "RECO"    
#vector<reco::BasicJet>                "ak8PFJetsCHSSoftDrop"      ""                "RECO"    
#vector<reco::BasicJet>                "cmsTopTagPFJetsCHS"        ""                "RECO"    
#vector<reco::CaloCluster>             "hfEMClusters"              ""                "RECO"    
#vector<reco::CaloCluster>             "particleFlowEGamma"        "EBEEClusters"    "RECO"    
#vector<reco::CaloCluster>             "particleFlowEGamma"        "ESClusters"      "RECO"    
#vector<reco::CaloCluster>             "hybridSuperClusters"       "hybridBarrelBasicClusters"   "RECO"    
#vector<reco::CaloCluster>             "multi5x5SuperClusters"     "multi5x5EndcapBasicClusters"   "RECO"    
#vector<reco::CaloCluster>             "particleFlowSuperClusterECAL"   "particleFlowBasicClusterECALBarrel"   "RECO"    
#vector<reco::CaloCluster>             "particleFlowSuperClusterECAL"   "particleFlowBasicClusterECALEndcap"   "RECO"    
#vector<reco::CaloCluster>             "particleFlowSuperClusterECAL"   "particleFlowBasicClusterECALPreshower"   "RECO"    
#vector<reco::CaloCluster>             "hybridSuperClusters"       "uncleanOnlyHybridBarrelBasicClusters"   "RECO"    
#vector<reco::CaloJet>                 "ak4CaloJets"               ""                "RECO"    
#vector<reco::CaloMET>                 "caloMet"                   ""                "RECO"    
#vector<reco::CaloMET>                 "caloMetBE"                 ""                "RECO"    
#vector<reco::CaloMET>                 "caloMetBEFO"               ""                "RECO"    
#vector<reco::CaloMET>                 "caloMetM"                  ""                "RECO"    
#vector<reco::CastorTower>             "CastorTowerReco"           ""                "RECO"    
#vector<reco::Conversion>              "allConversions"            ""                "RECO"    
#vector<reco::Conversion>              "conversions"               ""                "RECO"    
#vector<reco::Conversion>              "particleFlowEGamma"        ""                "RECO"    
#vector<reco::Conversion>              "uncleanedOnlyAllConversions"   ""                "RECO"    
#vector<reco::DeDxHitInfo>             "dedxHitInfo"               ""                "RECO"    
#vector<reco::GsfElectron>             "gedGsfElectrons"           ""                "RECO"    
#vector<reco::GsfElectron>             "uncleanedOnlyGsfElectrons"   ""                "RECO"    
#vector<reco::GsfElectronCore>         "gedGsfElectronCores"       ""                "RECO"    
#vector<reco::GsfElectronCore>         "uncleanedOnlyGsfElectronCores"   ""                "RECO"    
#vector<reco::GsfTrack>                "electronGsfTracks"         ""                "RECO"    
#vector<reco::HFEMClusterShape>        "hfEMClusters"              ""                "RECO"    
#vector<reco::JPTJet>                  "JetPlusTrackZSPCorJetAntiKt4"   ""                "RECO"    
#vector<reco::Muon>                    "muons"                     ""                "RECO"    
#vector<reco::Muon>                    "muonsFromCosmics"          ""                "RECO"    
#vector<reco::Muon>                    "muonsFromCosmics1Leg"      ""                "RECO"    
#vector<reco::PFCandidate>             "particleFlow"              ""                "RECO"    
#vector<reco::PFCandidate>             "pfIsolatedElectronsEI"     ""                "RECO"    
#vector<reco::PFCandidate>             "pfIsolatedMuonsEI"         ""                "RECO"    
#vector<reco::PFCandidate>             "particleFlowTmp"           "AddedMuonsAndHadrons"   "RECO"    
#vector<reco::PFCandidate>             "particleFlowTmp"           "CleanedCosmicsMuons"   "RECO"    
#vector<reco::PFCandidate>             "particleFlowTmp"           "CleanedFakeMuons"   "RECO"    
#vector<reco::PFCandidate>             "particleFlowTmp"           "CleanedHF"       "RECO"    
#vector<reco::PFCandidate>             "particleFlowTmp"           "CleanedPunchThroughMuons"   "RECO"    
#vector<reco::PFCandidate>             "particleFlowTmp"           "CleanedPunchThroughNeutralHadrons"   "RECO"    
#vector<reco::PFCandidate>             "particleFlowTmp"           "CleanedTrackerAndGlobalMuons"   "RECO"    
#vector<reco::PFJet>                   "ak4PFJets"                 ""                "RECO"    
#vector<reco::PFJet>                   "ak4PFJetsCHS"              ""                "RECO"    
#vector<reco::PFJet>                   "ak8PFJetsCHS"              ""                "RECO"    
#vector<reco::PFJet>                   "pfJetsEI"                  ""                "RECO"    
#vector<reco::PFJet>                   "ak8PFJetsCHSSoftDrop"      "SubJets"         "RECO"    
#vector<reco::PFJet>                   "cmsTopTagPFJetsCHS"        "caTopSubJets"    "RECO"    
#vector<reco::PFMET>                   "pfChMet"                   ""                "RECO"    
#vector<reco::PFMET>                   "pfMet"                     ""                "RECO"    
#vector<reco::PFMET>                   "pfMetEI"                   ""                "RECO"    
#vector<reco::PFRecHit>                "particleFlowRecHitECAL"    "Cleaned"         "RECO"    
#vector<reco::PFRecHit>                "particleFlowRecHitHBHE"    "Cleaned"         "RECO"    
#vector<reco::PFRecHit>                "particleFlowRecHitHF"      "Cleaned"         "RECO"    
#vector<reco::PFRecHit>                "particleFlowRecHitHO"      "Cleaned"         "RECO"    
#vector<reco::PFRecHit>                "particleFlowRecHitPS"      "Cleaned"         "RECO"    
#vector<reco::PFTau>                   "hpsPFTauProducer"          ""                "RECO"    
#vector<reco::PFTau>                   "pfTausEI"                  ""                "RECO"    
#vector<reco::PFTauTransverseImpactParameter>    "hpsPFTauTransverseImpactParameters"   "PFTauTIP"        "RECO"    
#vector<reco::Photon>                  "gedPhotons"                ""                "RECO"    
#vector<reco::Photon>                  "photons"                   ""                "RECO"    
#vector<reco::PhotonCore>              "gedPhotonCore"             ""                "RECO"    
#vector<reco::PhotonCore>              "photonCore"                ""                "RECO"    
#vector<reco::PreshowerCluster>        "multi5x5SuperClustersWithPreshower"   "preshowerXClusters"   "RECO"    
#vector<reco::PreshowerCluster>        "multi5x5SuperClustersWithPreshower"   "preshowerYClusters"   "RECO"    
#vector<reco::PreshowerClusterShape>    "multi5x5PreshowerClusterShape"   "multi5x5PreshowerXClustersShape"   "RECO"    
#vector<reco::PreshowerClusterShape>    "multi5x5PreshowerClusterShape"   "multi5x5PreshowerYClustersShape"   "RECO"    
#vector<reco::RecoChargedRefCandidate>    "trackRefsForJets"          ""                "RECO"    
#vector<reco::RecoEcalCandidate>       "hfRecoEcalCandidate"       ""                "RECO"    
#vector<reco::RecoTauPiZero>           "hpsPFTauProducer"          "pizeros"         "RECO"    
#vector<reco::SuperCluster>            "correctedHybridSuperClusters"   ""                "RECO"    
#vector<reco::SuperCluster>            "correctedMulti5x5SuperClustersWithPreshower"   ""                "RECO"    
#vector<reco::SuperCluster>            "hfEMClusters"              ""                "RECO"    
#vector<reco::SuperCluster>            "particleFlowEGamma"        ""                "RECO"    
#vector<reco::SuperCluster>            "particleFlowSuperClusterECAL"   "particleFlowSuperClusterECALBarrel"   "RECO"    
#vector<reco::SuperCluster>            "particleFlowSuperClusterECAL"   "particleFlowSuperClusterECALEndcapWithPreshower"   "RECO"    
#vector<reco::SuperCluster>            "hybridSuperClusters"       "uncleanOnlyHybridSuperClusters"   "RECO"    
#vector<reco::Track>                   "ckfInOutTracksFromConversions"   ""                "RECO"    
#vector<reco::Track>                   "ckfOutInTracksFromConversions"   ""                "RECO"    
#vector<reco::Track>                   "conversionStepTracks"      ""                "RECO"    
#vector<reco::Track>                   "cosmicMuons"               ""                "RECO"    
#vector<reco::Track>                   "cosmicMuons1Leg"           ""                "RECO"    
#vector<reco::Track>                   "displacedGlobalMuons"      ""                "RECO"    
#vector<reco::Track>                   "displacedStandAloneMuons"   ""                "RECO"    
#vector<reco::Track>                   "displacedTracks"           ""                "RECO"    
#vector<reco::Track>                   "generalTracks"             ""                "RECO"    
#vector<reco::Track>                   "globalMuons"               ""                "RECO"    
#vector<reco::Track>                   "refittedStandAloneMuons"   ""                "RECO"    
#vector<reco::Track>                   "standAloneMuons"           ""                "RECO"    
#vector<reco::Track>                   "uncleanedOnlyCkfInOutTracksFromConversions"   ""                "RECO"    
#vector<reco::Track>                   "uncleanedOnlyCkfOutInTracksFromConversions"   ""                "RECO"    
#vector<reco::Track>                   "refittedStandAloneMuons"   "UpdatedAtVtx"    "RECO"    
#vector<reco::Track>                   "standAloneMuons"           "UpdatedAtVtx"    "RECO"    
#vector<reco::Track>                   "tevMuons"                  "default"         "RECO"    
#vector<reco::Track>                   "tevMuons"                  "dyt"             "RECO"    
#vector<reco::Track>                   "tevMuons"                  "firstHit"        "RECO"    
#vector<reco::Track>                   "tevMuons"                  "picky"           "RECO"    
#vector<reco::TrackExtra>              "displacedGlobalMuons"      ""                "RECO"    
#vector<reco::TrackExtra>              "displacedStandAloneMuons"   ""                "RECO"    
#vector<reco::TrackExtra>              "globalMuons"               ""                "RECO"    
#vector<reco::TrackExtra>              "refittedStandAloneMuons"   ""                "RECO"    
#vector<reco::TrackExtra>              "standAloneMuons"           ""                "RECO"    
#vector<reco::TrackExtra>              "tevMuons"                  "default"         "RECO"    
#vector<reco::TrackExtra>              "tevMuons"                  "dyt"             "RECO"    
#vector<reco::TrackExtra>              "tevMuons"                  "firstHit"        "RECO"    
#vector<reco::TrackExtra>              "tevMuons"                  "picky"           "RECO"    
#vector<reco::TrackExtrapolation>      "trackExtrapolator"         ""                "RECO"    
#vector<reco::TrackJet>                "ak4TrackJets"              ""                "RECO"    
#vector<reco::Vertex>                  "inclusiveSecondaryVertices"   ""                "RECO"    
#vector<reco::Vertex>                  "offlinePrimaryVertices"    ""                "RECO"    
#vector<reco::Vertex>                  "offlinePrimaryVerticesWithBS"   ""                "RECO"    
#vector<reco::VertexCompositeCandidate>    "generalV0Candidates"       "Kshort"          "RECO"    
#vector<reco::VertexCompositeCandidate>    "generalV0Candidates"       "Lambda"          "RECO"    
#vector<reco::VertexCompositePtrCandidate>    "inclusiveCandidateSecondaryVertices"   ""                "RECO"    
#vector<reco::VertexCompositePtrCandidate>    "inclusiveCandidateSecondaryVerticesCvsL"   ""                "RECO"    

#Type                                  Module                      Label             Process   
#----------------------------------------------------------------------------------------------
#LHEEventProduct                       "externalLHEProducer"       ""                "LHE"     
#GenEventInfoProduct                   "generator"                 ""                "SIM"     
#edm::TriggerResults                   "TriggerResults"            ""                "SIM"     
#edm::TriggerResults                   "TriggerResults"            ""                "HLT"     
#HcalNoiseSummary                      "hcalnoise"                 ""                "RECO"    
#L1GlobalTriggerReadoutRecord          "gtDigis"                   ""                "RECO"    
#double                                "fixedGridRhoAll"           ""                "RECO"    
#double                                "fixedGridRhoFastjetAll"    ""                "RECO"    
#double                                "fixedGridRhoFastjetAllCalo"   ""                "RECO"    
#double                                "fixedGridRhoFastjetCentral"   ""                "RECO"    
#double                                "fixedGridRhoFastjetCentralCalo"   ""                "RECO"    
#double                                "fixedGridRhoFastjetCentralChargedPileUp"   ""                "RECO"    
#double                                "fixedGridRhoFastjetCentralNeutral"   ""                "RECO"    
#edm::TriggerResults                   "TriggerResults"            ""                "RECO"    
#reco::BeamHaloSummary                 "BeamHaloSummary"           ""                "RECO"    
#reco::BeamSpot                        "offlineBeamSpot"           ""                "RECO"    
#reco::CSCHaloData                     "CSCHaloData"               ""                "RECO"    
#vector<l1extra::L1EmParticle>         "l1extraParticles"          "Isolated"        "RECO"    
#vector<l1extra::L1EmParticle>         "l1extraParticles"          "NonIsolated"     "RECO"    
#vector<l1extra::L1EtMissParticle>     "l1extraParticles"          "MET"             "RECO"    
#vector<l1extra::L1EtMissParticle>     "l1extraParticles"          "MHT"             "RECO"    
#vector<l1extra::L1HFRings>            "l1extraParticles"          ""                "RECO"    
#vector<l1extra::L1JetParticle>        "l1extraParticles"          "Central"         "RECO"    
#vector<l1extra::L1JetParticle>        "l1extraParticles"          "Forward"         "RECO"    
#vector<l1extra::L1JetParticle>        "l1extraParticles"          "IsoTau"          "RECO"    
#vector<l1extra::L1JetParticle>        "l1extraParticles"          "Tau"             "RECO"    
#vector<l1extra::L1MuonParticle>       "l1extraParticles"          ""                "RECO"    
#edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "reducedEgamma"             "reducedEBRecHits"   "PAT"     
#edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "reducedEgamma"             "reducedEERecHits"   "PAT"     
#edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "reducedEgamma"             "reducedESRecHits"   "PAT"     
#edm::TriggerResults                   "TriggerResults"            ""                "PAT"     
#edm::ValueMap<float>                  "offlineSlimmedPrimaryVertices"   ""                "PAT"     
#pat::PackedTriggerPrescales           "patTrigger"                ""                "PAT"     
#pat::PackedTriggerPrescales           "patTrigger"                "l1max"           "PAT"     
#pat::PackedTriggerPrescales           "patTrigger"                "l1min"           "PAT"     
#vector<PileupSummaryInfo>             "slimmedAddPileupInfo"      ""                "PAT"     
#vector<pat::Electron>                 "slimmedElectrons"          ""                "PAT"     
#vector<pat::Jet>                      "slimmedJets"               ""                "PAT"     
#vector<pat::Jet>                      "slimmedJetsAK8"            ""                "PAT"     
#vector<pat::Jet>                      "slimmedJetsPuppi"          ""                "PAT"     
#vector<pat::Jet>                      "slimmedJetsAK8PFCHSSoftDropPacked"   "SubJets"         "PAT"     
#vector<pat::Jet>                      "slimmedJetsCMSTopTagCHSPacked"   "SubJets"         "PAT"     
#vector<pat::MET>                      "slimmedMETs"               ""                "PAT"     
#vector<pat::MET>                      "slimmedMETsNoHF"           ""                "PAT"     
#vector<pat::MET>                      "slimmedMETsPuppi"          ""                "PAT"     
#vector<pat::Muon>                     "slimmedMuons"              ""                "PAT"     
#vector<pat::PackedCandidate>          "lostTracks"                ""                "PAT"     
#vector<pat::PackedCandidate>          "packedPFCandidates"        ""                "PAT"     
#vector<pat::PackedGenParticle>        "packedGenParticles"        ""                "PAT"     
#vector<pat::Photon>                   "slimmedPhotons"            ""                "PAT"     
#vector<pat::Tau>                      "slimmedTaus"               ""                "PAT"     
#vector<pat::Tau>                      "slimmedTausBoosted"        ""                "PAT"     
#vector<pat::TriggerObjectStandAlone>    "selectedPatTrigger"        ""                "PAT"     
#vector<reco::CATopJetTagInfo>         "caTopTagInfosPAT"          ""                "PAT"     
#vector<reco::CaloCluster>             "reducedEgamma"             "reducedEBEEClusters"   "PAT"     
#vector<reco::CaloCluster>             "reducedEgamma"             "reducedESClusters"   "PAT"     
#vector<reco::Conversion>              "reducedEgamma"             "reducedConversions"   "PAT"     
#vector<reco::Conversion>              "reducedEgamma"             "reducedSingleLegConversions"   "PAT"     
#vector<reco::GenJet>                  "slimmedGenJets"            ""                "PAT"     
#vector<reco::GenJet>                  "slimmedGenJetsAK8"         ""                "PAT"     
#vector<reco::GenParticle>             "prunedGenParticles"        ""                "PAT"     
#vector<reco::GsfElectronCore>         "reducedEgamma"             "reducedGedGsfElectronCores"   "PAT"     
#vector<reco::PhotonCore>              "reducedEgamma"             "reducedGedPhotonCores"   "PAT"     
#vector<reco::SuperCluster>            "reducedEgamma"             "reducedSuperClusters"   "PAT"     
#vector<reco::Vertex>                  "offlineSlimmedPrimaryVertices"   ""                "PAT"     
#vector<reco::VertexCompositePtrCandidate>    "slimmedSecondaryVertices"   ""                "PA

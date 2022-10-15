//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 29 09:47:48 2013 by ROOT version 5.32/00
// from TTree EventTree/Event data
// found on file: /data1/ggNtuples/V05-03-11-00/job_hiele_2013pA_PRv1.root
//////////////////////////////////////////////////////////

#ifndef xAna_h
#define xAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TChainElement.h>

#include <iostream>
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxeleESEffSigmaRR = 1;
const Int_t kMaxphoESEffSigmaRR = 1;
const Int_t kMaxnPFPho = 36;
const Int_t kMaxPFPhoEt = 1;
const Int_t kMaxPFPhoEta = 1;
const Int_t kMaxPFPhoPhi = 1;
const Int_t kMaxPFPhoType = 1;
const Int_t kMaxPFPhoIso = 1;

class xAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nHLT;
   Int_t           HLT[164];   //[nHLT]
   Int_t           HLTIndex[70];
   Float_t         bspotPos[3];
   Int_t           nVtx;
   Float_t         vtx[7][3];   //[nVtx]
   Int_t           IsVtxGood;
   Int_t           nGoodVtx;
   Float_t         centrality[5];
   Int_t           nVtxBS;
   Float_t         vtxbs[6][3];   //[nVtxBS]
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         recoPfMET;
   Float_t         recoPfMETPhi;
   Float_t         recoPfMETsumEt;
   Float_t         recoPfMETmEtSig;
   Float_t         recoPfMETSig;
   Float_t         trkMETxPV;
   Float_t         trkMETyPV;
   Float_t         trkMETPhiPV;
   Float_t         trkMETPV;
   Float_t         trkMETx[6];   //[nVtxBS]
   Float_t         trkMETy[6];   //[nVtxBS]
   Float_t         trkMETPhi[6];   //[nVtxBS]
   Float_t         trkMET[6];   //[nVtxBS]
   Int_t           metFilters[10];
   Int_t           nEle;
   Int_t           eleTrg[9][16];   //[nEle]
   Int_t           eleClass[9];   //[nEle]
   Int_t           eleIsEcalDriven[9];   //[nEle]
   Int_t           eleCharge[9];   //[nEle]
   Float_t         eleEn[9];   //[nEle]
   Float_t         eleEcalEn[9];   //[nEle]
   Float_t         eleSCRawEn[9];   //[nEle]
   Float_t         eleSCEn[9];   //[nEle]
   Float_t         eleESEn[9];   //[nEle]
   Float_t         elePt[9];   //[nEle]
   Float_t         eleEta[9];   //[nEle]
   Float_t         elePhi[9];   //[nEle]
   Float_t         eleEtaVtx[9][100];   //[nEle]
   Float_t         elePhiVtx[9][100];   //[nEle]
   Float_t         eleEtVtx[9][100];   //[nEle]
   Float_t         eleSCEta[9];   //[nEle]
   Float_t         eleSCPhi[9];   //[nEle]
   Float_t         eleSCEtaWidth[9];   //[nEle]
   Float_t         eleSCPhiWidth[9];   //[nEle]
   Float_t         eleVtx[9][3];   //[nEle]
   Float_t         eleD0[9];   //[nEle]
   Float_t         eleDz[9];   //[nEle]
   Float_t         eleD0GV[9];   //[nEle]
   Float_t         eleDzGV[9];   //[nEle]
   Float_t         eleD0Vtx[9][100];   //[nEle]
   Float_t         eleDzVtx[9][100];   //[nEle]
   Float_t         eleHoverE[9];   //[nEle]
   Float_t         eleHoverE12[9];   //[nEle]
   Float_t         eleEoverP[9];   //[nEle]
   Float_t         elePin[9];   //[nEle]
   Float_t         elePout[9];   //[nEle]
   Float_t         eleTrkMomErr[9];   //[nEle]
   Float_t         eleBrem[9];   //[nEle]
   Float_t         eledEtaAtVtx[9];   //[nEle]
   Float_t         eledPhiAtVtx[9];   //[nEle]
   Float_t         eleSigmaIEtaIEta[9];   //[nEle]
   Float_t         eleSigmaIEtaIPhi[9];   //[nEle]
   Float_t         eleSigmaIPhiIPhi[9];   //[nEle]
   Float_t         eleEmax[9];   //[nEle]
   Float_t         eleE1x5[9];   //[nEle]
   Float_t         eleE3x3[9];   //[nEle]
   Float_t         eleE5x5[9];   //[nEle]
   Float_t         eleE2x5Max[9];   //[nEle]
   Float_t         eleRegrE[9];   //[nEle]
   Float_t         eleRegrEerr[9];   //[nEle]
   Float_t         elePhoRegrE[9];   //[nEle]
   Float_t         elePhoRegrEerr[9];   //[nEle]
   Float_t         eleSeedTime[9];   //[nEle]
   Int_t           eleRecoFlag[9];   //[nEle]
   Int_t           elePos[9];   //[nEle]
   Float_t         eleIsoTrkDR03[9];   //[nEle]
   Float_t         eleIsoEcalDR03[9];   //[nEle]
   Float_t         eleIsoHcalDR03[9];   //[nEle]
   Float_t         eleIsoHcalDR0312[9];   //[nEle]
   Float_t         eleIsoTrkDR04[9];   //[nEle]
   Float_t         eleIsoEcalDR04[9];   //[nEle]
   Float_t         eleIsoHcalDR04[9];   //[nEle]
   Float_t         eleIsoHcalDR0412[9];   //[nEle]
   Float_t         eleModIsoTrk[9];   //[nEle]
   Float_t         eleModIsoEcal[9];   //[nEle]
   Float_t         eleModIsoHcal[9];   //[nEle]
   Int_t           eleMissHits[9];   //[nEle]
   Float_t         eleConvDist[9];   //[nEle]
   Float_t         eleConvDcot[9];   //[nEle]
   Int_t           eleConvVtxFit[9];   //[nEle]
   Float_t         eleIP3D[9];   //[nEle]
   Float_t         eleIP3DErr[9];   //[nEle]
   Float_t         eleIDMVANonTrig[9];   //[nEle]
   Float_t         eleIDMVATrig[9];   //[nEle]
   Float_t         elePFChIso03[9];   //[nEle]
   Float_t         elePFPhoIso03[9];   //[nEle]
   Float_t         elePFNeuIso03[9];   //[nEle]
   Float_t         elePFChIso04[9];   //[nEle]
   Float_t         elePFPhoIso04[9];   //[nEle]
   Float_t         elePFNeuIso04[9];   //[nEle]
   Float_t         eleESEffSigmaRR[9][3];   //[nEle]
   Int_t           nPho;
   Int_t           phoTrg[10][8];   //[nPho]
   Int_t           phoTrgFilter[10][50];   //[nPho]
   Bool_t          phoIsPhoton[10];   //[nPho]
   Float_t         phoSCPos[10][3];   //[nPho]
   Float_t         phoCaloPos[10][3];   //[nPho]
   Float_t         phoE[10];   //[nPho]
   Float_t         phoEt[10];   //[nPho]
   Float_t         phoEta[10];   //[nPho]
   Float_t         phoVtx[10][3];   //[nPho]
   Float_t         phoPhi[10];   //[nPho]
   Float_t         phoEtVtx[10][100];   //[nPho]
   Float_t         phoEtaVtx[10][100];   //[nPho]
   Float_t         phoPhiVtx[10][100];   //[nPho]
   Float_t         phoR9[10];   //[nPho]
   Float_t         phoTrkIsoHollowDR03[10];   //[nPho]
   Float_t         phoEcalIsoDR03[10];   //[nPho]
   Float_t         phoHcalIsoDR03[10];   //[nPho]
   Float_t         phoHcalIsoDR0312[10];   //[nPho]
   Float_t         phoTrkIsoHollowDR04[10];   //[nPho]
   Float_t         phoCiCTrkIsoDR03[10][100];   //[nPho]
   Float_t         phoCiCTrkIsoDR04[10][100];   //[nPho]
   Float_t         phoCiCdRtoTrk[10];   //[nPho]
   Float_t         phoEcalIsoDR04[10];   //[nPho]
   Float_t         phoHcalIsoDR04[10];   //[nPho]
   Float_t         phoHcalIsoDR0412[10];   //[nPho]
   Float_t         phoHoverE[10];   //[nPho]
   Float_t         phoHoverE12[10];   //[nPho]
   Int_t           phoEleVeto[10];   //[nPho]
   Float_t         phoSigmaIEtaIEta[10];   //[nPho]
   Float_t         phoSigmaIEtaIPhi[10];   //[nPho]
   Float_t         phoSigmaIPhiIPhi[10];   //[nPho]
   Float_t         phoCiCPF4phopfIso03[10];   //[nPho]
   Float_t         phoCiCPF4phopfIso04[10];   //[nPho]
   Float_t         phoCiCPF4chgpfIso02[10][100];   //[nPho]
   Float_t         phoCiCPF4chgpfIso03[10][100];   //[nPho]
   Float_t         phoCiCPF4chgpfIso04[10][100];   //[nPho]
   Float_t         phoEmax[10];   //[nPho]
   Float_t         phoE3x3[10];   //[nPho]
   Float_t         phoE3x1[10];   //[nPho]
   Float_t         phoE1x3[10];   //[nPho]
   Float_t         phoE5x5[10];   //[nPho]
   Float_t         phoE1x5[10];   //[nPho]
   Float_t         phoE2x2[10];   //[nPho]
   Float_t         phoE2x5Max[10];   //[nPho]
   Float_t         phoPFChIso[10];   //[nPho]
   Float_t         phoPFPhoIso[10];   //[nPho]
   Float_t         phoPFNeuIso[10];   //[nPho]
   Float_t         phoSCRChIso[10];   //[nPho]
   Float_t         phoSCRPhoIso[10];   //[nPho]
   Float_t         phoSCRNeuIso[10];   //[nPho]
   Float_t         phoRegrE[10];   //[nPho]
   Float_t         phoRegrEerr[10];   //[nPho]
   Float_t         phoSeedTime[10];   //[nPho]
   Int_t           phoSeedDetId1[10];   //[nPho]
   Int_t           phoSeedDetId2[10];   //[nPho]
   Float_t         phoLICTD[10];   //[nPho]
   Int_t           phoRecoFlag[10];   //[nPho]
   Int_t           phoPos[10];   //[nPho]
   Float_t         phoSCE[10];   //[nPho]
   Float_t         phoSCRawE[10];   //[nPho]
   Float_t         phoESEn[10];   //[nPho]
   Float_t         phoSCEt[10];   //[nPho]
   Float_t         phoSCEta[10];   //[nPho]
   Float_t         phoSCPhi[10];   //[nPho]
   Float_t         phoSCEtaWidth[10];   //[nPho]
   Float_t         phoSCPhiWidth[10];   //[nPho]
   Float_t         phoSCBrem[10];   //[nPho]
   Int_t           phoOverlap[10];   //[nPho]
   Int_t           phohasPixelSeed[10];   //[nPho]
   Int_t           pho_hasConvPf[10];   //[nPho]
   Int_t           pho_hasSLConvPf[10];   //[nPho]
   Float_t         pho_pfconvVtxZ[10];   //[nPho]
   Float_t         pho_pfconvVtxZErr[10];   //[nPho]
   Int_t           pho_nSLConv[10];   //[nPho]
   Float_t         pho_pfSLConvPos[10][3];   //[nPho]
   Float_t         pho_pfSLConvVtxZ[10][20];   //[nPho]
   Int_t           phoIsConv[10];   //[nPho]
   Int_t           phoNConv[10];   //[nPho]
   Float_t         phoConvInvMass[10];   //[nPho]
   Float_t         phoConvCotTheta[10];   //[nPho]
   Float_t         phoConvEoverP[10];   //[nPho]
   Float_t         phoConvZofPVfromTrks[10];   //[nPho]
   Float_t         phoConvMinDist[10];   //[nPho]
   Float_t         phoConvdPhiAtVtx[10];   //[nPho]
   Float_t         phoConvdPhiAtCalo[10];   //[nPho]
   Float_t         phoConvdEtaAtCalo[10];   //[nPho]
   Float_t         phoConvTrkd0[10][2];   //[nPho]
   Float_t         phoConvTrkPin[10][2];   //[nPho]
   Float_t         phoConvTrkPout[10][2];   //[nPho]
   Float_t         phoConvTrkdz[10][2];   //[nPho]
   Float_t         phoConvTrkdzErr[10][2];   //[nPho]
   Float_t         phoConvChi2[10];   //[nPho]
   Float_t         phoConvChi2Prob[10];   //[nPho]
   Int_t           phoConvNTrks[10];   //[nPho]
   Float_t         phoConvCharge[10][2];   //[nPho]
   Float_t         phoConvValidVtx[10];   //[nPho]
   Float_t         phoConvLikeLihood[10];   //[nPho]
   Float_t         phoConvP4[10][4];   //[nPho]
   Float_t         phoConvVtx[10][3];   //[nPho]
   Float_t         phoConvVtxErr[10][3];   //[nPho]
   Float_t         phoConvPairMomentum[10][3];   //[nPho]
   Float_t         phoConvRefittedMomentum[10][3];   //[nPho]
   Int_t           SingleLegConv[10];   //[nPho]
   Float_t         phoPFConvVtx[10][3];   //[nPho]
   Float_t         phoPFConvMom[10][3];   //[nPho]
   Float_t         phoESEffSigmaRR[10][3];   //[nPho]
   Int_t           nMu;
   Int_t           muTrg[15][10];   //[nMu]
   Float_t         muEta[15];   //[nMu]
   Float_t         muPhi[15];   //[nMu]
   Int_t           muCharge[15];   //[nMu]
   Float_t         muPt[15];   //[nMu]
   Float_t         muPz[15];   //[nMu]
   Float_t         muVtx[15][3];   //[nMu]
   Float_t         muVtxGlb[15][3];   //[nMu]
   Float_t         mucktPt[15];   //[nMu]
   Float_t         mucktPtErr[15];   //[nMu]
   Float_t         mucktEta[15];   //[nMu]
   Float_t         mucktPhi[15];   //[nMu]
   Float_t         mucktdxy[15];   //[nMu]
   Float_t         mucktdz[15];   //[nMu]
   Float_t         muIsoTrk[15];   //[nMu]
   Float_t         muIsoCalo[15];   //[nMu]
   Float_t         muIsoEcal[15];   //[nMu]
   Float_t         muIsoHcal[15];   //[nMu]
   Float_t         muChi2NDF[15];   //[nMu]
   Float_t         muInnerChi2NDF[15];   //[nMu]
   Float_t         muPFIsoR04_CH[15];   //[nMu]
   Float_t         muPFIsoR04_NH[15];   //[nMu]
   Float_t         muPFIsoR04_Pho[15];   //[nMu]
   Float_t         muPFIsoR04_PU[15];   //[nMu]
   Float_t         muPFIsoR04_CPart[15];   //[nMu]
   Float_t         muPFIsoR04_NHHT[15];   //[nMu]
   Float_t         muPFIsoR04_PhoHT[15];   //[nMu]
   Float_t         muPFIsoR03_CH[15];   //[nMu]
   Float_t         muPFIsoR03_NH[15];   //[nMu]
   Float_t         muPFIsoR03_Pho[15];   //[nMu]
   Float_t         muPFIsoR03_PU[15];   //[nMu]
   Float_t         muPFIsoR03_CPart[15];   //[nMu]
   Float_t         muPFIsoR03_NHHT[15];   //[nMu]
   Float_t         muPFIsoR03_PhoHT[15];   //[nMu]
   Int_t           muType[15];   //[nMu]
   Float_t         muD0[15];   //[nMu]
   Float_t         muDz[15];   //[nMu]
   Float_t         muD0GV[15];   //[nMu]
   Float_t         muDzGV[15];   //[nMu]
   Float_t         muD0Vtx[15][100];   //[nMu]
   Float_t         muDzVtx[15][100];   //[nMu]
   Float_t         muInnerD0[15];   //[nMu]
   Float_t         muInnerDz[15];   //[nMu]
   Float_t         muInnerD0GV[15];   //[nMu]
   Float_t         muInnerDzGV[15];   //[nMu]
   Float_t         muInnerPt[15];   //[nMu]
   Float_t         muInnerPtErr[15];   //[nMu]
   Int_t           muNumberOfValidTrkLayers[15];   //[nMu]
   Int_t           muNumberOfValidTrkHits[15];   //[nMu]
   Int_t           muNumberOfValidPixelLayers[15];   //[nMu]
   Int_t           muNumberOfValidPixelHits[15];   //[nMu]
   Int_t           muNumberOfValidMuonHits[15];   //[nMu]
   Int_t           muStations[15];   //[nMu]
   Int_t           muChambers[15];   //[nMu]
   Float_t         muIP3D[15];   //[nMu]
   Float_t         muIP3DErr[15];   //[nMu]
   Int_t           nPFPho;
   Float_t         PFPhoEt[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFPhoEta[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFPhoPhi[kMaxnPFPho];   //[nPFPho_]
   Int_t           PFPhoType[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFPhoIso[kMaxnPFPho];   //[nPFPho_]
   Float_t         rho25;
   Float_t         rho25_neu;
   Float_t         rho25_muPFiso;
   Float_t         rho25_elePFiso;
   Float_t         rho2011;
   Float_t         rho2012;
   Int_t           nJet;
   Int_t           jetTrg[36][14];   //[nJet]
   Float_t         jetEn[36];   //[nJet]
   Float_t         jetPt[36];   //[nJet]
   Float_t         jetEta[36];   //[nJet]
   Float_t         jetPhi[36];   //[nJet]
   Float_t         jetCharge[36];   //[nJet]
   Float_t         jetEt[36];   //[nJet]
   Float_t         jetRawPt[36];   //[nJet]
   Float_t         jetRawEn[36];   //[nJet]
   Float_t         jetArea[36];   //[nJet]
   Float_t         jetCHF[36];   //[nJet]
   Float_t         jetNHF[36];   //[nJet]
   Float_t         jetCEF[36];   //[nJet]
   Float_t         jetNEF[36];   //[nJet]
   Int_t           jetNCH[36];   //[nJet]
   Float_t         jetHFHAE[36];   //[nJet]
   Float_t         jetHFEME[36];   //[nJet]
   Int_t           jetNConstituents[36];   //[nJet]
   Float_t         jetCombinedSecondaryVtxBJetTags[36];   //[nJet]
   Float_t         jetCombinedSecondaryVtxMVABJetTags[36];   //[nJet]
   Float_t         jetJetProbabilityBJetTags[36];   //[nJet]
   Float_t         jetJetBProbabilityBJetTags[36];   //[nJet]
   Float_t         jetTrackCountingHighPurBJetTags[36];   //[nJet]
   Float_t         jetBetaStar[36][100];   //[nJet]
   Bool_t          jetPFLooseId[36];   //[nJet]
   Float_t         jetDRMean[36];   //[nJet]
   Float_t         jetDR2Mean[36];   //[nJet]
   Float_t         jetDZ[36];   //[nJet]
   Float_t         jetFrac01[36];   //[nJet]
   Float_t         jetFrac02[36];   //[nJet]
   Float_t         jetFrac03[36];   //[nJet]
   Float_t         jetFrac04[36];   //[nJet]
   Float_t         jetFrac05[36];   //[nJet]
   Float_t         jetFrac06[36];   //[nJet]
   Float_t         jetFrac07[36];   //[nJet]
   Float_t         jetBeta[36];   //[nJet]
   Float_t         jetBetaStarCMG[36];   //[nJet]
   Float_t         jetBetaStarClassic[36];   //[nJet]
   Float_t         jetBetaExt[36][100];   //[nJet]
   Float_t         jetBetaStarCMGExt[36][100];   //[nJet]
   Float_t         jetBetaStarClassicExt[36][100];   //[nJet]
   Float_t         jetNNeutrals[36];   //[nJet]
   Float_t         jetNCharged[36];   //[nJet]
   Float_t         jetMVAs[36][4];   //[nJet]
   Int_t           jetWPLevels[36][4];   //[nJet]
   Float_t         jetMVAsExt[36][4][100];   //[nJet]
   Int_t           jetWPLevelsExt[36][4][100];   //[nJet]
   Float_t         jetMt[36];   //[nJet]
   Float_t         jetJECUnc[36];   //[nJet]
   Float_t         jetLeadTrackPt[36];   //[nJet]
   Float_t         jetVtxPt[36];   //[nJet]
   Float_t         jetVtxMass[36];   //[nJet]
   Float_t         jetVtx3dL[36];   //[nJet]
   Float_t         jetVtx3deL[36];   //[nJet]
   Float_t         jetSoftLeptPt[36];   //[nJet]
   Float_t         jetSoftLeptPtRel[36];   //[nJet]
   Float_t         jetSoftLeptdR[36];   //[nJet]
   Float_t         jetSoftLeptIdlooseMu[36];   //[nJet]
   Float_t         jetSoftLeptIdEle95[36];   //[nJet]
   Float_t         jetDPhiMETJet[36];   //[nJet]
   Float_t         jetPuJetIdL[36];   //[nJet]
   Float_t         jetPuJetIdM[36];   //[nJet]
   Float_t         jetPuJetIdT[36];   //[nJet]
   Int_t           nLowPtJet;
   Float_t         jetLowPtEn[49];   //[nLowPtJet]
   Float_t         jetLowPtPt[49];   //[nLowPtJet]
   Float_t         jetLowPtEta[49];   //[nLowPtJet]
   Float_t         jetLowPtPhi[49];   //[nLowPtJet]
   Float_t         jetLowPtCharge[49];   //[nLowPtJet]
   Float_t         jetLowPtEt[49];   //[nLowPtJet]
   Float_t         jetLowPtRawPt[49];   //[nLowPtJet]
   Float_t         jetLowPtRawEn[49];   //[nLowPtJet]
   Float_t         jetLowPtArea[49];   //[nLowPtJet]
   Int_t           nConv;
   Float_t         convP4[422][4];   //[nConv]
   Float_t         convVtx[422][3];   //[nConv]
   Float_t         convVtxErr[422][3];   //[nConv]
   Float_t         convPairMomentum[422][3];   //[nConv]
   Float_t         convRefittedMomentum[422][3];   //[nConv]
   Int_t           convNTracks[422];   //[nConv]
   Float_t         convPairInvMass[422];   //[nConv]
   Float_t         convPairCotThetaSep[422];   //[nConv]
   Float_t         convEoverP[422];   //[nConv]
   Float_t         convDistOfMinApproach[422];   //[nConv]
   Float_t         convDPhiTrksAtVtx[422];   //[nConv]
   Float_t         convDPhiTrksAtEcal[422];   //[nConv]
   Float_t         convDEtaTrksAtEcal[422];   //[nConv]
   Float_t         convDxy[422];   //[nConv]
   Float_t         convDz[422];   //[nConv]
   Float_t         convLxy[422];   //[nConv]
   Float_t         convLz[422];   //[nConv]
   Float_t         convZofPrimVtxFromTrks[422];   //[nConv]
   Int_t           convNHitsBeforeVtx[422][2];   //[nConv]
   Int_t           convNSharedHits[422];   //[nConv]
   Int_t           convValidVtx[422];   //[nConv]
   Float_t         convMVALikelihood[422];   //[nConv]
   Float_t         convChi2[422];   //[nConv]
   Float_t         convChi2Probability[422];   //[nConv]
   Float_t         convTk1Dz[422];   //[nConv]
   Float_t         convTk2Dz[422];   //[nConv]
   Float_t         convTk1DzErr[422];   //[nConv]
   Float_t         convTk2DzErr[422];   //[nConv]
   Int_t           convCh1Ch2[422];   //[nConv]
   Float_t         convTk1D0[422];   //[nConv]
   Float_t         convTk1Pout[422];   //[nConv]
   Float_t         convTk1Pin[422];   //[nConv]
   Float_t         convTk2D0[422];   //[nConv]
   Float_t         convTk2Pout[422];   //[nConv]
   Float_t         convTk2Pin[422];   //[nConv]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_HLT;   //!
   TBranch        *b_HLTIndex;   //!
   TBranch        *b_bspotPos;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_IsVtxGood;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_centrality;   //!
   TBranch        *b_nVtxBS;   //!
   TBranch        *b_vtxbs;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_recoPfMET;   //!
   TBranch        *b_recoPfMETPhi;   //!
   TBranch        *b_recoPfMETsumEt;   //!
   TBranch        *b_recoPfMETmEtSig;   //!
   TBranch        *b_recoPfMETSig;   //!
   TBranch        *b_trkMETxPV;   //!
   TBranch        *b_trkMETyPV;   //!
   TBranch        *b_trkMETPhiPV;   //!
   TBranch        *b_trkMETPV;   //!
   TBranch        *b_trkMETx;   //!
   TBranch        *b_trkMETy;   //!
   TBranch        *b_trkMETPhi;   //!
   TBranch        *b_trkMET;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleTrg;   //!
   TBranch        *b_eleClass;   //!
   TBranch        *b_eleIsEcalDriven;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleEcalEn;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleEtaVtx;   //!
   TBranch        *b_elePhiVtx;   //!
   TBranch        *b_eleEtVtx;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleVtx;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleD0GV;   //!
   TBranch        *b_eleDzGV;   //!
   TBranch        *b_eleD0Vtx;   //!
   TBranch        *b_eleDzVtx;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleHoverE12;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_elePin;   //!
   TBranch        *b_elePout;   //!
   TBranch        *b_eleTrkMomErr;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleEmax;   //!
   TBranch        *b_eleE1x5;   //!
   TBranch        *b_eleE3x3;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE2x5Max;   //!
   TBranch        *b_eleRegrE;   //!
   TBranch        *b_eleRegrEerr;   //!
   TBranch        *b_elePhoRegrE;   //!
   TBranch        *b_elePhoRegrEerr;   //!
   TBranch        *b_eleSeedTime;   //!
   TBranch        *b_eleRecoFlag;   //!
   TBranch        *b_elePos;   //!
   TBranch        *b_eleIsoTrkDR03;   //!
   TBranch        *b_eleIsoEcalDR03;   //!
   TBranch        *b_eleIsoHcalDR03;   //!
   TBranch        *b_eleIsoHcalDR0312;   //!
   TBranch        *b_eleIsoTrkDR04;   //!
   TBranch        *b_eleIsoEcalDR04;   //!
   TBranch        *b_eleIsoHcalDR04;   //!
   TBranch        *b_eleIsoHcalDR0412;   //!
   TBranch        *b_eleModIsoTrk;   //!
   TBranch        *b_eleModIsoEcal;   //!
   TBranch        *b_eleModIsoHcal;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleConvDist;   //!
   TBranch        *b_eleConvDcot;   //!
   TBranch        *b_eleConvVtxFit;   //!
   TBranch        *b_eleIP3D;   //!
   TBranch        *b_eleIP3DErr;   //!
   TBranch        *b_eleIDMVANonTrig;   //!
   TBranch        *b_eleIDMVATrig;   //!
   TBranch        *b_elePFChIso03;   //!
   TBranch        *b_elePFPhoIso03;   //!
   TBranch        *b_elePFNeuIso03;   //!
   TBranch        *b_elePFChIso04;   //!
   TBranch        *b_elePFPhoIso04;   //!
   TBranch        *b_elePFNeuIso04;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoTrg;   //!
   TBranch        *b_phoTrgFilter;   //!
   TBranch        *b_phoIsPhoton;   //!
   TBranch        *b_phoSCPos;   //!
   TBranch        *b_phoCaloPos;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoVtx;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEtVtx;   //!
   TBranch        *b_phoEtaVtx;   //!
   TBranch        *b_phoPhiVtx;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoTrkIsoHollowDR03;   //!
   TBranch        *b_phoEcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR0312;   //!
   TBranch        *b_phoTrkIsoHollowDR04;   //!
   TBranch        *b_phoCiCTrkIsoDR03;   //!
   TBranch        *b_phoCiCTrkIsoDR04;   //!
   TBranch        *b_phoCiCdRtoTrk;   //!
   TBranch        *b_phoEcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR0412;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoHoverE12;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoCiCPF4phopfIso03;   //!
   TBranch        *b_phoCiCPF4phopfIso04;   //!
   TBranch        *b_phoCiCPF4chgpfIso02;   //!
   TBranch        *b_phoCiCPF4chgpfIso03;   //!
   TBranch        *b_phoCiCPF4chgpfIso04;   //!
   TBranch        *b_phoEmax;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE3x1;   //!
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoSCRChIso;   //!
   TBranch        *b_phoSCRPhoIso;   //!
   TBranch        *b_phoSCRNeuIso;   //!
   TBranch        *b_phoRegrE;   //!
   TBranch        *b_phoRegrEerr;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedDetId1;   //!
   TBranch        *b_phoSeedDetId2;   //!
   TBranch        *b_phoLICTD;   //!
   TBranch        *b_phoRecoFlag;   //!
   TBranch        *b_phoPos;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoSCEt;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phoOverlap;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_pho_hasConvPf;   //!
   TBranch        *b_pho_hasSLConvPf;   //!
   TBranch        *b_pho_pfconvVtxZ;   //!
   TBranch        *b_pho_pfconvVtxZErr;   //!
   TBranch        *b_pho_nSLConv;   //!
   TBranch        *b_pho_pfSLConvPos;   //!
   TBranch        *b_pho_pfSLConvVtxZ;   //!
   TBranch        *b_phoIsConv;   //!
   TBranch        *b_phoNConv;   //!
   TBranch        *b_phoConvInvMass;   //!
   TBranch        *b_phoConvCotTheta;   //!
   TBranch        *b_phoConvEoverP;   //!
   TBranch        *b_phoConvZofPVfromTrks;   //!
   TBranch        *b_phoConvMinDist;   //!
   TBranch        *b_phoConvdPhiAtVtx;   //!
   TBranch        *b_phoConvdPhiAtCalo;   //!
   TBranch        *b_phoConvdEtaAtCalo;   //!
   TBranch        *b_phoConvTrkd0;   //!
   TBranch        *b_phoConvTrkPin;   //!
   TBranch        *b_phoConvTrkPout;   //!
   TBranch        *b_phoConvTrkdz;   //!
   TBranch        *b_phoConvTrkdzErr;   //!
   TBranch        *b_phoConvChi2;   //!
   TBranch        *b_phoConvChi2Prob;   //!
   TBranch        *b_phoConvNTrks;   //!
   TBranch        *b_phoConvCharge;   //!
   TBranch        *b_phoConvValidVtx;   //!
   TBranch        *b_phoConvLikeLihood;   //!
   TBranch        *b_phoConvP4;   //!
   TBranch        *b_phoConvVtx;   //!
   TBranch        *b_phoConvVtxErr;   //!
   TBranch        *b_phoConvPairMomentum;   //!
   TBranch        *b_phoConvRefittedMomentum;   //!
   TBranch        *b_SingleLegConv;   //!
   TBranch        *b_phoPFConvVtx;   //!
   TBranch        *b_phoPFConvMom;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muTrg;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muPz;   //!
   TBranch        *b_muVtx;   //!
   TBranch        *b_muVtxGlb;   //!
   TBranch        *b_mucktPt;   //!
   TBranch        *b_mucktPtErr;   //!
   TBranch        *b_mucktEta;   //!
   TBranch        *b_mucktPhi;   //!
   TBranch        *b_mucktdxy;   //!
   TBranch        *b_mucktdz;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muIsoCalo;   //!
   TBranch        *b_muIsoEcal;   //!
   TBranch        *b_muIsoHcal;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerChi2NDF;   //!
   TBranch        *b_muPFIsoR04_CH;   //!
   TBranch        *b_muPFIsoR04_NH;   //!
   TBranch        *b_muPFIsoR04_Pho;   //!
   TBranch        *b_muPFIsoR04_PU;   //!
   TBranch        *b_muPFIsoR04_CPart;   //!
   TBranch        *b_muPFIsoR04_NHHT;   //!
   TBranch        *b_muPFIsoR04_PhoHT;   //!
   TBranch        *b_muPFIsoR03_CH;   //!
   TBranch        *b_muPFIsoR03_NH;   //!
   TBranch        *b_muPFIsoR03_Pho;   //!
   TBranch        *b_muPFIsoR03_PU;   //!
   TBranch        *b_muPFIsoR03_CPart;   //!
   TBranch        *b_muPFIsoR03_NHHT;   //!
   TBranch        *b_muPFIsoR03_PhoHT;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muD0GV;   //!
   TBranch        *b_muDzGV;   //!
   TBranch        *b_muD0Vtx;   //!
   TBranch        *b_muDzVtx;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muInnerD0GV;   //!
   TBranch        *b_muInnerDzGV;   //!
   TBranch        *b_muInnerPt;   //!
   TBranch        *b_muInnerPtErr;   //!
   TBranch        *b_muNumberOfValidTrkLayers;   //!
   TBranch        *b_muNumberOfValidTrkHits;   //!
   TBranch        *b_muNumberOfValidPixelLayers;   //!
   TBranch        *b_muNumberOfValidPixelHits;   //!
   TBranch        *b_muNumberOfValidMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muChambers;   //!
   TBranch        *b_muIP3D;   //!
   TBranch        *b_muIP3DErr;   //!
   TBranch        *b_nPFPho_;   //!
   TBranch        *b_PFPhoEt;   //!
   TBranch        *b_PFPhoEta;   //!
   TBranch        *b_PFPhoPhi;   //!
   TBranch        *b_PFPhoType;   //!
   TBranch        *b_PFPhoIso;   //!
   TBranch        *b_rho25;   //!
   TBranch        *b_rho25_neu;   //!
   TBranch        *b_rho25_muPFiso;   //!
   TBranch        *b_rho25_elePFiso;   //!
   TBranch        *b_rho2011;   //!
   TBranch        *b_rho2012;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetTrg;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCharge;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetCombinedSecondaryVtxBJetTags;   //!
   TBranch        *b_jetCombinedSecondaryVtxMVABJetTags;   //!
   TBranch        *b_jetJetProbabilityBJetTags;   //!
   TBranch        *b_jetJetBProbabilityBJetTags;   //!
   TBranch        *b_jetTrackCountingHighPurBJetTags;   //!
   TBranch        *b_jetBetaStar;   //!
   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetDRMean;   //!
   TBranch        *b_jetDR2Mean;   //!
   TBranch        *b_jetDZ;   //!
   TBranch        *b_jetFrac01;   //!
   TBranch        *b_jetFrac02;   //!
   TBranch        *b_jetFrac03;   //!
   TBranch        *b_jetFrac04;   //!
   TBranch        *b_jetFrac05;   //!
   TBranch        *b_jetFrac06;   //!
   TBranch        *b_jetFrac07;   //!
   TBranch        *b_jetBeta;   //!
   TBranch        *b_jetBetaStarCMG;   //!
   TBranch        *b_jetBetaStarClassic;   //!
   TBranch        *b_jetBetaExt;   //!
   TBranch        *b_jetBetaStarCMGExt;   //!
   TBranch        *b_jetBetaStarClassicExt;   //!
   TBranch        *b_jetNNeutrals;   //!
   TBranch        *b_jetNCharged;   //!
   TBranch        *b_jetMVAs;   //!
   TBranch        *b_jetWPLevels;   //!
   TBranch        *b_jetMVAsExt;   //!
   TBranch        *b_jetWPLevelsExt;   //!
   TBranch        *b_jetMt;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetLeadTrackPt;   //!
   TBranch        *b_jetVtxPt;   //!
   TBranch        *b_jetVtxMass;   //!
   TBranch        *b_jetVtx3dL;   //!
   TBranch        *b_jetVtx3deL;   //!
   TBranch        *b_jetSoftLeptPt;   //!
   TBranch        *b_jetSoftLeptPtRel;   //!
   TBranch        *b_jetSoftLeptdR;   //!
   TBranch        *b_jetSoftLeptIdlooseMu;   //!
   TBranch        *b_jetSoftLeptIdEle95;   //!
   TBranch        *b_jetDPhiMETJet;   //!
   TBranch        *b_jetPuJetIdL;   //!
   TBranch        *b_jetPuJetIdM;   //!
   TBranch        *b_jetPuJetIdT;   //!
   TBranch        *b_nLowPtJet;   //!
   TBranch        *b_jetLowPtEn;   //!
   TBranch        *b_jetLowPtPt;   //!
   TBranch        *b_jetLowPtEta;   //!
   TBranch        *b_jetLowPtPhi;   //!
   TBranch        *b_jetLowPtCharge;   //!
   TBranch        *b_jetLowPtEt;   //!
   TBranch        *b_jetLowPtRawPt;   //!
   TBranch        *b_jetLowPtRawEn;   //!
   TBranch        *b_jetLowPtArea;   //!
   TBranch        *b_nConv;   //!
   TBranch        *b_convP4;   //!
   TBranch        *b_convVtx;   //!
   TBranch        *b_convVtxErr;   //!
   TBranch        *b_convPairMomentum;   //!
   TBranch        *b_convRefittedMomentum;   //!
   TBranch        *b_convNTracks;   //!
   TBranch        *b_convPairInvMass;   //!
   TBranch        *b_convPairCotThetaSep;   //!
   TBranch        *b_convEoverP;   //!
   TBranch        *b_convDistOfMinApproach;   //!
   TBranch        *b_convDPhiTrksAtVtx;   //!
   TBranch        *b_convDPhiTrksAtEcal;   //!
   TBranch        *b_convDEtaTrksAtEcal;   //!
   TBranch        *b_convDxy;   //!
   TBranch        *b_convDz;   //!
   TBranch        *b_convLxy;   //!
   TBranch        *b_convLz;   //!
   TBranch        *b_convZofPrimVtxFromTrks;   //!
   TBranch        *b_convNHitsBeforeVtx;   //!
   TBranch        *b_convNSharedHits;   //!
   TBranch        *b_convValidVtx;   //!
   TBranch        *b_convMVALikelihood;   //!
   TBranch        *b_convChi2;   //!
   TBranch        *b_convChi2Probability;   //!
   TBranch        *b_convTk1Dz;   //!
   TBranch        *b_convTk2Dz;   //!
   TBranch        *b_convTk1DzErr;   //!
   TBranch        *b_convTk2DzErr;   //!
   TBranch        *b_convCh1Ch2;   //!
   TBranch        *b_convTk1D0;   //!
   TBranch        *b_convTk1Pout;   //!
   TBranch        *b_convTk1Pin;   //!
   TBranch        *b_convTk2D0;   //!
   TBranch        *b_convTk2Pout;   //!
   TBranch        *b_convTk2Pin;   //!

   xAna(TTree *tree=0);
   virtual ~xAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     InitTreePass();
   virtual void     InitTreeFail();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   // Lepton selection
   virtual Int_t    eleKin(Int_t i = 0);
   virtual Float_t  geteIDEA(Int_t i = 0);
   virtual Float_t  getePFIso(Int_t i = 0);
   virtual Int_t    eID2012(Int_t i = 0, Int_t idWP = 0);
   virtual Int_t    muID2012(Int_t i = 0);
   virtual Int_t    jetID(Int_t i = 0);

   TFile *fout_;
   TTree *PassingWP85Tree;
   TTree *FailingWP85Tree;
   TString filename_;
   
   Double_t Weight;
   Double_t Tag_Pt;
   Double_t Tag_Eta;
   Double_t Tag_Phi;
   Double_t Tag_SCEta;
   Double_t Tag_isMCEle;
   Int_t    Tag_charge;
   Int_t    Tag_status;
   Double_t Probe_Pt;
   Double_t Probe_Eta;
   Double_t Probe_Phi;
   Double_t Probe_SCEta;
   Double_t Probe_isMCEle;
   Int_t    Probe_charge;
   Int_t    Probe_status;
   Double_t Zmass;
   Double_t cent;
   Int_t elepasstrg;
   Int_t elefailtrg;
};

#endif

#ifdef xAna_cxx
xAna::xAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data1/ggNtuples/V05-03-11-00/job_hiele_2013pA_PRv1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data1/ggNtuples/V05-03-11-00/job_hiele_2013pA_PRv1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/data1/ggNtuples/V05-03-11-00/job_hiele_2013pA_PRv1.root:/ggNtuplizer");
      dir->GetObject("EventTree",tree);

   }
   Init(tree);
}

xAna::~xAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t xAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t xAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void xAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("HLT", HLT, &b_HLT);
   fChain->SetBranchAddress("HLTIndex", HLTIndex, &b_HLTIndex);
   fChain->SetBranchAddress("bspotPos", bspotPos, &b_bspotPos);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("IsVtxGood", &IsVtxGood, &b_IsVtxGood);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("centrality", centrality, &b_centrality);
   fChain->SetBranchAddress("nVtxBS", &nVtxBS, &b_nVtxBS);
   fChain->SetBranchAddress("vtxbs", vtxbs, &b_vtxbs);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("recoPfMET", &recoPfMET, &b_recoPfMET);
   fChain->SetBranchAddress("recoPfMETPhi", &recoPfMETPhi, &b_recoPfMETPhi);
   fChain->SetBranchAddress("recoPfMETsumEt", &recoPfMETsumEt, &b_recoPfMETsumEt);
   fChain->SetBranchAddress("recoPfMETmEtSig", &recoPfMETmEtSig, &b_recoPfMETmEtSig);
   fChain->SetBranchAddress("recoPfMETSig", &recoPfMETSig, &b_recoPfMETSig);
   fChain->SetBranchAddress("trkMETxPV", &trkMETxPV, &b_trkMETxPV);
   fChain->SetBranchAddress("trkMETyPV", &trkMETyPV, &b_trkMETyPV);
   fChain->SetBranchAddress("trkMETPhiPV", &trkMETPhiPV, &b_trkMETPhiPV);
   fChain->SetBranchAddress("trkMETPV", &trkMETPV, &b_trkMETPV);
   fChain->SetBranchAddress("trkMETx", trkMETx, &b_trkMETx);
   fChain->SetBranchAddress("trkMETy", trkMETy, &b_trkMETy);
   fChain->SetBranchAddress("trkMETPhi", trkMETPhi, &b_trkMETPhi);
   fChain->SetBranchAddress("trkMET", trkMET, &b_trkMET);
   fChain->SetBranchAddress("metFilters", metFilters, &b_metFilters);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleTrg", eleTrg, &b_eleTrg);
   fChain->SetBranchAddress("eleClass", eleClass, &b_eleClass);
   fChain->SetBranchAddress("eleIsEcalDriven", eleIsEcalDriven, &b_eleIsEcalDriven);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleEn", eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleEcalEn", eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("eleSCRawEn", eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEn", eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleESEn", eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleEtaVtx", eleEtaVtx, &b_eleEtaVtx);
   fChain->SetBranchAddress("elePhiVtx", elePhiVtx, &b_elePhiVtx);
   fChain->SetBranchAddress("eleEtVtx", eleEtVtx, &b_eleEtVtx);
   fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCEtaWidth", eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleVtx", eleVtx, &b_eleVtx);
   fChain->SetBranchAddress("eleD0", eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleD0GV", eleD0GV, &b_eleD0GV);
   fChain->SetBranchAddress("eleDzGV", eleDzGV, &b_eleDzGV);
   fChain->SetBranchAddress("eleD0Vtx", eleD0Vtx, &b_eleD0Vtx);
   fChain->SetBranchAddress("eleDzVtx", eleDzVtx, &b_eleDzVtx);
   fChain->SetBranchAddress("eleHoverE", eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleHoverE12", eleHoverE12, &b_eleHoverE12);
   fChain->SetBranchAddress("eleEoverP", eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("elePin", elePin, &b_elePin);
   fChain->SetBranchAddress("elePout", elePout, &b_elePout);
   fChain->SetBranchAddress("eleTrkMomErr", eleTrkMomErr, &b_eleTrkMomErr);
   fChain->SetBranchAddress("eleBrem", eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleEmax", eleEmax, &b_eleEmax);
   fChain->SetBranchAddress("eleE1x5", eleE1x5, &b_eleE1x5);
   fChain->SetBranchAddress("eleE3x3", eleE3x3, &b_eleE3x3);
   fChain->SetBranchAddress("eleE5x5", eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE2x5Max", eleE2x5Max, &b_eleE2x5Max);
   fChain->SetBranchAddress("eleRegrE", eleRegrE, &b_eleRegrE);
   fChain->SetBranchAddress("eleRegrEerr", eleRegrEerr, &b_eleRegrEerr);
   fChain->SetBranchAddress("elePhoRegrE", elePhoRegrE, &b_elePhoRegrE);
   fChain->SetBranchAddress("elePhoRegrEerr", elePhoRegrEerr, &b_elePhoRegrEerr);
   fChain->SetBranchAddress("eleSeedTime", eleSeedTime, &b_eleSeedTime);
   fChain->SetBranchAddress("eleRecoFlag", eleRecoFlag, &b_eleRecoFlag);
   fChain->SetBranchAddress("elePos", elePos, &b_elePos);
   fChain->SetBranchAddress("eleIsoTrkDR03", eleIsoTrkDR03, &b_eleIsoTrkDR03);
   fChain->SetBranchAddress("eleIsoEcalDR03", eleIsoEcalDR03, &b_eleIsoEcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR03", eleIsoHcalDR03, &b_eleIsoHcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR0312", eleIsoHcalDR0312, &b_eleIsoHcalDR0312);
   fChain->SetBranchAddress("eleIsoTrkDR04", eleIsoTrkDR04, &b_eleIsoTrkDR04);
   fChain->SetBranchAddress("eleIsoEcalDR04", eleIsoEcalDR04, &b_eleIsoEcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR04", eleIsoHcalDR04, &b_eleIsoHcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR0412", eleIsoHcalDR0412, &b_eleIsoHcalDR0412);
   fChain->SetBranchAddress("eleModIsoTrk", eleModIsoTrk, &b_eleModIsoTrk);
   fChain->SetBranchAddress("eleModIsoEcal", eleModIsoEcal, &b_eleModIsoEcal);
   fChain->SetBranchAddress("eleModIsoHcal", eleModIsoHcal, &b_eleModIsoHcal);
   fChain->SetBranchAddress("eleMissHits", eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleConvDist", eleConvDist, &b_eleConvDist);
   fChain->SetBranchAddress("eleConvDcot", eleConvDcot, &b_eleConvDcot);
   fChain->SetBranchAddress("eleConvVtxFit", eleConvVtxFit, &b_eleConvVtxFit);
   fChain->SetBranchAddress("eleIP3D", eleIP3D, &b_eleIP3D);
   fChain->SetBranchAddress("eleIP3DErr", eleIP3DErr, &b_eleIP3DErr);
   fChain->SetBranchAddress("eleIDMVANonTrig", eleIDMVANonTrig, &b_eleIDMVANonTrig);
   fChain->SetBranchAddress("eleIDMVATrig", eleIDMVATrig, &b_eleIDMVATrig);
   fChain->SetBranchAddress("elePFChIso03", elePFChIso03, &b_elePFChIso03);
   fChain->SetBranchAddress("elePFPhoIso03", elePFPhoIso03, &b_elePFPhoIso03);
   fChain->SetBranchAddress("elePFNeuIso03", elePFNeuIso03, &b_elePFNeuIso03);
   fChain->SetBranchAddress("elePFChIso04", elePFChIso04, &b_elePFChIso04);
   fChain->SetBranchAddress("elePFPhoIso04", elePFPhoIso04, &b_elePFPhoIso04);
   fChain->SetBranchAddress("elePFNeuIso04", elePFNeuIso04, &b_elePFNeuIso04);
   fChain->SetBranchAddress("eleESEffSigmaRR", eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoTrg", phoTrg, &b_phoTrg);
   fChain->SetBranchAddress("phoTrgFilter", phoTrgFilter, &b_phoTrgFilter);
   fChain->SetBranchAddress("phoIsPhoton", phoIsPhoton, &b_phoIsPhoton);
   fChain->SetBranchAddress("phoSCPos", phoSCPos, &b_phoSCPos);
   fChain->SetBranchAddress("phoCaloPos", phoCaloPos, &b_phoCaloPos);
   fChain->SetBranchAddress("phoE", phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoVtx", phoVtx, &b_phoVtx);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEtVtx", phoEtVtx, &b_phoEtVtx);
   fChain->SetBranchAddress("phoEtaVtx", phoEtaVtx, &b_phoEtaVtx);
   fChain->SetBranchAddress("phoPhiVtx", phoPhiVtx, &b_phoPhiVtx);
   fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoTrkIsoHollowDR03", phoTrkIsoHollowDR03, &b_phoTrkIsoHollowDR03);
   fChain->SetBranchAddress("phoEcalIsoDR03", phoEcalIsoDR03, &b_phoEcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR03", phoHcalIsoDR03, &b_phoHcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR0312", phoHcalIsoDR0312, &b_phoHcalIsoDR0312);
   fChain->SetBranchAddress("phoTrkIsoHollowDR04", phoTrkIsoHollowDR04, &b_phoTrkIsoHollowDR04);
   fChain->SetBranchAddress("phoCiCTrkIsoDR03", phoCiCTrkIsoDR03, &b_phoCiCTrkIsoDR03);
   fChain->SetBranchAddress("phoCiCTrkIsoDR04", phoCiCTrkIsoDR04, &b_phoCiCTrkIsoDR04);
   fChain->SetBranchAddress("phoCiCdRtoTrk", phoCiCdRtoTrk, &b_phoCiCdRtoTrk);
   fChain->SetBranchAddress("phoEcalIsoDR04", phoEcalIsoDR04, &b_phoEcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR04", phoHcalIsoDR04, &b_phoHcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR0412", phoHcalIsoDR0412, &b_phoHcalIsoDR0412);
   fChain->SetBranchAddress("phoHoverE", phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoHoverE12", phoHoverE12, &b_phoHoverE12);
   fChain->SetBranchAddress("phoEleVeto", phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoCiCPF4phopfIso03", phoCiCPF4phopfIso03, &b_phoCiCPF4phopfIso03);
   fChain->SetBranchAddress("phoCiCPF4phopfIso04", phoCiCPF4phopfIso04, &b_phoCiCPF4phopfIso04);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso02", phoCiCPF4chgpfIso02, &b_phoCiCPF4chgpfIso02);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso03", phoCiCPF4chgpfIso03, &b_phoCiCPF4chgpfIso03);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso04", phoCiCPF4chgpfIso04, &b_phoCiCPF4chgpfIso04);
   fChain->SetBranchAddress("phoEmax", phoEmax, &b_phoEmax);
   fChain->SetBranchAddress("phoE3x3", phoE3x3, &b_phoE3x3);
   fChain->SetBranchAddress("phoE3x1", phoE3x1, &b_phoE3x1);
   fChain->SetBranchAddress("phoE1x3", phoE1x3, &b_phoE1x3);
   fChain->SetBranchAddress("phoE5x5", phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoE1x5", phoE1x5, &b_phoE1x5);
   fChain->SetBranchAddress("phoE2x2", phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE2x5Max", phoE2x5Max, &b_phoE2x5Max);
   fChain->SetBranchAddress("phoPFChIso", phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoSCRChIso", phoSCRChIso, &b_phoSCRChIso);
   fChain->SetBranchAddress("phoSCRPhoIso", phoSCRPhoIso, &b_phoSCRPhoIso);
   fChain->SetBranchAddress("phoSCRNeuIso", phoSCRNeuIso, &b_phoSCRNeuIso);
   fChain->SetBranchAddress("phoRegrE", phoRegrE, &b_phoRegrE);
   fChain->SetBranchAddress("phoRegrEerr", phoRegrEerr, &b_phoRegrEerr);
   fChain->SetBranchAddress("phoSeedTime", phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedDetId1", phoSeedDetId1, &b_phoSeedDetId1);
   fChain->SetBranchAddress("phoSeedDetId2", phoSeedDetId2, &b_phoSeedDetId2);
   fChain->SetBranchAddress("phoLICTD", phoLICTD, &b_phoLICTD);
   fChain->SetBranchAddress("phoRecoFlag", phoRecoFlag, &b_phoRecoFlag);
   fChain->SetBranchAddress("phoPos", phoPos, &b_phoPos);
   fChain->SetBranchAddress("phoSCE", phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEn", phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoSCEt", phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCEta", phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phoOverlap", phoOverlap, &b_phoOverlap);
   fChain->SetBranchAddress("phohasPixelSeed", phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("pho_hasConvPf", pho_hasConvPf, &b_pho_hasConvPf);
   fChain->SetBranchAddress("pho_hasSLConvPf", pho_hasSLConvPf, &b_pho_hasSLConvPf);
   fChain->SetBranchAddress("pho_pfconvVtxZ", pho_pfconvVtxZ, &b_pho_pfconvVtxZ);
   fChain->SetBranchAddress("pho_pfconvVtxZErr", pho_pfconvVtxZErr, &b_pho_pfconvVtxZErr);
   fChain->SetBranchAddress("pho_nSLConv", pho_nSLConv, &b_pho_nSLConv);
   fChain->SetBranchAddress("pho_pfSLConvPos", pho_pfSLConvPos, &b_pho_pfSLConvPos);
   fChain->SetBranchAddress("pho_pfSLConvVtxZ", pho_pfSLConvVtxZ, &b_pho_pfSLConvVtxZ);
   fChain->SetBranchAddress("phoIsConv", phoIsConv, &b_phoIsConv);
   fChain->SetBranchAddress("phoNConv", phoNConv, &b_phoNConv);
   fChain->SetBranchAddress("phoConvInvMass", phoConvInvMass, &b_phoConvInvMass);
   fChain->SetBranchAddress("phoConvCotTheta", phoConvCotTheta, &b_phoConvCotTheta);
   fChain->SetBranchAddress("phoConvEoverP", phoConvEoverP, &b_phoConvEoverP);
   fChain->SetBranchAddress("phoConvZofPVfromTrks", phoConvZofPVfromTrks, &b_phoConvZofPVfromTrks);
   fChain->SetBranchAddress("phoConvMinDist", phoConvMinDist, &b_phoConvMinDist);
   fChain->SetBranchAddress("phoConvdPhiAtVtx", phoConvdPhiAtVtx, &b_phoConvdPhiAtVtx);
   fChain->SetBranchAddress("phoConvdPhiAtCalo", phoConvdPhiAtCalo, &b_phoConvdPhiAtCalo);
   fChain->SetBranchAddress("phoConvdEtaAtCalo", phoConvdEtaAtCalo, &b_phoConvdEtaAtCalo);
   fChain->SetBranchAddress("phoConvTrkd0", phoConvTrkd0, &b_phoConvTrkd0);
   fChain->SetBranchAddress("phoConvTrkPin", phoConvTrkPin, &b_phoConvTrkPin);
   fChain->SetBranchAddress("phoConvTrkPout", phoConvTrkPout, &b_phoConvTrkPout);
   fChain->SetBranchAddress("phoConvTrkdz", phoConvTrkdz, &b_phoConvTrkdz);
   fChain->SetBranchAddress("phoConvTrkdzErr", phoConvTrkdzErr, &b_phoConvTrkdzErr);
   fChain->SetBranchAddress("phoConvChi2", phoConvChi2, &b_phoConvChi2);
   fChain->SetBranchAddress("phoConvChi2Prob", phoConvChi2Prob, &b_phoConvChi2Prob);
   fChain->SetBranchAddress("phoConvNTrks", phoConvNTrks, &b_phoConvNTrks);
   fChain->SetBranchAddress("phoConvCharge", phoConvCharge, &b_phoConvCharge);
   fChain->SetBranchAddress("phoConvValidVtx", phoConvValidVtx, &b_phoConvValidVtx);
   fChain->SetBranchAddress("phoConvLikeLihood", phoConvLikeLihood, &b_phoConvLikeLihood);
   fChain->SetBranchAddress("phoConvP4", phoConvP4, &b_phoConvP4);
   fChain->SetBranchAddress("phoConvVtx", phoConvVtx, &b_phoConvVtx);
   fChain->SetBranchAddress("phoConvVtxErr", phoConvVtxErr, &b_phoConvVtxErr);
   fChain->SetBranchAddress("phoConvPairMomentum", phoConvPairMomentum, &b_phoConvPairMomentum);
   fChain->SetBranchAddress("phoConvRefittedMomentum", phoConvRefittedMomentum, &b_phoConvRefittedMomentum);
   fChain->SetBranchAddress("SingleLegConv", SingleLegConv, &b_SingleLegConv);
   fChain->SetBranchAddress("phoPFConvVtx", phoPFConvVtx, &b_phoPFConvVtx);
   fChain->SetBranchAddress("phoPFConvMom", phoPFConvMom, &b_phoPFConvMom);
   fChain->SetBranchAddress("phoESEffSigmaRR", phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muTrg", muTrg, &b_muTrg);
   fChain->SetBranchAddress("muEta", muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", muCharge, &b_muCharge);
   fChain->SetBranchAddress("muPt", muPt, &b_muPt);
   fChain->SetBranchAddress("muPz", muPz, &b_muPz);
   fChain->SetBranchAddress("muVtx", muVtx, &b_muVtx);
   fChain->SetBranchAddress("muVtxGlb", muVtxGlb, &b_muVtxGlb);
   fChain->SetBranchAddress("mucktPt", mucktPt, &b_mucktPt);
   fChain->SetBranchAddress("mucktPtErr", mucktPtErr, &b_mucktPtErr);
   fChain->SetBranchAddress("mucktEta", mucktEta, &b_mucktEta);
   fChain->SetBranchAddress("mucktPhi", mucktPhi, &b_mucktPhi);
   fChain->SetBranchAddress("mucktdxy", mucktdxy, &b_mucktdxy);
   fChain->SetBranchAddress("mucktdz", mucktdz, &b_mucktdz);
   fChain->SetBranchAddress("muIsoTrk", muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muIsoCalo", muIsoCalo, &b_muIsoCalo);
   fChain->SetBranchAddress("muIsoEcal", muIsoEcal, &b_muIsoEcal);
   fChain->SetBranchAddress("muIsoHcal", muIsoHcal, &b_muIsoHcal);
   fChain->SetBranchAddress("muChi2NDF", muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerChi2NDF", muInnerChi2NDF, &b_muInnerChi2NDF);
   fChain->SetBranchAddress("muPFIsoR04_CH", muPFIsoR04_CH, &b_muPFIsoR04_CH);
   fChain->SetBranchAddress("muPFIsoR04_NH", muPFIsoR04_NH, &b_muPFIsoR04_NH);
   fChain->SetBranchAddress("muPFIsoR04_Pho", muPFIsoR04_Pho, &b_muPFIsoR04_Pho);
   fChain->SetBranchAddress("muPFIsoR04_PU", muPFIsoR04_PU, &b_muPFIsoR04_PU);
   fChain->SetBranchAddress("muPFIsoR04_CPart", muPFIsoR04_CPart, &b_muPFIsoR04_CPart);
   fChain->SetBranchAddress("muPFIsoR04_NHHT", muPFIsoR04_NHHT, &b_muPFIsoR04_NHHT);
   fChain->SetBranchAddress("muPFIsoR04_PhoHT", muPFIsoR04_PhoHT, &b_muPFIsoR04_PhoHT);
   fChain->SetBranchAddress("muPFIsoR03_CH", muPFIsoR03_CH, &b_muPFIsoR03_CH);
   fChain->SetBranchAddress("muPFIsoR03_NH", muPFIsoR03_NH, &b_muPFIsoR03_NH);
   fChain->SetBranchAddress("muPFIsoR03_Pho", muPFIsoR03_Pho, &b_muPFIsoR03_Pho);
   fChain->SetBranchAddress("muPFIsoR03_PU", muPFIsoR03_PU, &b_muPFIsoR03_PU);
   fChain->SetBranchAddress("muPFIsoR03_CPart", muPFIsoR03_CPart, &b_muPFIsoR03_CPart);
   fChain->SetBranchAddress("muPFIsoR03_NHHT", muPFIsoR03_NHHT, &b_muPFIsoR03_NHHT);
   fChain->SetBranchAddress("muPFIsoR03_PhoHT", muPFIsoR03_PhoHT, &b_muPFIsoR03_PhoHT);
   fChain->SetBranchAddress("muType", muType, &b_muType);
   fChain->SetBranchAddress("muD0", muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", muDz, &b_muDz);
   fChain->SetBranchAddress("muD0GV", muD0GV, &b_muD0GV);
   fChain->SetBranchAddress("muDzGV", muDzGV, &b_muDzGV);
   fChain->SetBranchAddress("muD0Vtx", muD0Vtx, &b_muD0Vtx);
   fChain->SetBranchAddress("muDzVtx", muDzVtx, &b_muDzVtx);
   fChain->SetBranchAddress("muInnerD0", muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muInnerD0GV", muInnerD0GV, &b_muInnerD0GV);
   fChain->SetBranchAddress("muInnerDzGV", muInnerDzGV, &b_muInnerDzGV);
   fChain->SetBranchAddress("muInnerPt", muInnerPt, &b_muInnerPt);
   fChain->SetBranchAddress("muInnerPtErr", muInnerPtErr, &b_muInnerPtErr);
   fChain->SetBranchAddress("muNumberOfValidTrkLayers", muNumberOfValidTrkLayers, &b_muNumberOfValidTrkLayers);
   fChain->SetBranchAddress("muNumberOfValidTrkHits", muNumberOfValidTrkHits, &b_muNumberOfValidTrkHits);
   fChain->SetBranchAddress("muNumberOfValidPixelLayers", muNumberOfValidPixelLayers, &b_muNumberOfValidPixelLayers);
   fChain->SetBranchAddress("muNumberOfValidPixelHits", muNumberOfValidPixelHits, &b_muNumberOfValidPixelHits);
   fChain->SetBranchAddress("muNumberOfValidMuonHits", muNumberOfValidMuonHits, &b_muNumberOfValidMuonHits);
   fChain->SetBranchAddress("muStations", muStations, &b_muStations);
   fChain->SetBranchAddress("muChambers", muChambers, &b_muChambers);
   fChain->SetBranchAddress("muIP3D", muIP3D, &b_muIP3D);
   fChain->SetBranchAddress("muIP3DErr", muIP3DErr, &b_muIP3DErr);
   fChain->SetBranchAddress("nPFPho", &nPFPho, &b_nPFPho_);
   fChain->SetBranchAddress("PFPhoEt", PFPhoEt, &b_PFPhoEt);
   fChain->SetBranchAddress("PFPhoEta", PFPhoEta, &b_PFPhoEta);
   fChain->SetBranchAddress("PFPhoPhi", PFPhoPhi, &b_PFPhoPhi);
   fChain->SetBranchAddress("PFPhoType", PFPhoType, &b_PFPhoType);
   fChain->SetBranchAddress("PFPhoIso", PFPhoIso, &b_PFPhoIso);
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("rho25_neu", &rho25_neu, &b_rho25_neu);
   fChain->SetBranchAddress("rho25_muPFiso", &rho25_muPFiso, &b_rho25_muPFiso);
   fChain->SetBranchAddress("rho25_elePFiso", &rho25_elePFiso, &b_rho25_elePFiso);
   fChain->SetBranchAddress("rho2011", &rho2011, &b_rho2011);
   fChain->SetBranchAddress("rho2012", &rho2012, &b_rho2012);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetTrg", jetTrg, &b_jetTrg);
   fChain->SetBranchAddress("jetEn", jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetCharge", jetCharge, &b_jetCharge);
   fChain->SetBranchAddress("jetEt", jetEt, &b_jetEt);
   fChain->SetBranchAddress("jetRawPt", jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetArea", jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetCHF", jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetHFHAE", jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxBJetTags", jetCombinedSecondaryVtxBJetTags, &b_jetCombinedSecondaryVtxBJetTags);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxMVABJetTags", jetCombinedSecondaryVtxMVABJetTags, &b_jetCombinedSecondaryVtxMVABJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetJetBProbabilityBJetTags", jetJetBProbabilityBJetTags, &b_jetJetBProbabilityBJetTags);
   fChain->SetBranchAddress("jetTrackCountingHighPurBJetTags", jetTrackCountingHighPurBJetTags, &b_jetTrackCountingHighPurBJetTags);
   fChain->SetBranchAddress("jetBetaStar", jetBetaStar, &b_jetBetaStar);
   fChain->SetBranchAddress("jetPFLooseId", jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetDRMean", jetDRMean, &b_jetDRMean);
   fChain->SetBranchAddress("jetDR2Mean", jetDR2Mean, &b_jetDR2Mean);
   fChain->SetBranchAddress("jetDZ", jetDZ, &b_jetDZ);
   fChain->SetBranchAddress("jetFrac01", jetFrac01, &b_jetFrac01);
   fChain->SetBranchAddress("jetFrac02", jetFrac02, &b_jetFrac02);
   fChain->SetBranchAddress("jetFrac03", jetFrac03, &b_jetFrac03);
   fChain->SetBranchAddress("jetFrac04", jetFrac04, &b_jetFrac04);
   fChain->SetBranchAddress("jetFrac05", jetFrac05, &b_jetFrac05);
   fChain->SetBranchAddress("jetFrac06", jetFrac06, &b_jetFrac06);
   fChain->SetBranchAddress("jetFrac07", jetFrac07, &b_jetFrac07);
   fChain->SetBranchAddress("jetBeta", jetBeta, &b_jetBeta);
   fChain->SetBranchAddress("jetBetaStarCMG", jetBetaStarCMG, &b_jetBetaStarCMG);
   fChain->SetBranchAddress("jetBetaStarClassic", jetBetaStarClassic, &b_jetBetaStarClassic);
   fChain->SetBranchAddress("jetBetaExt", jetBetaExt, &b_jetBetaExt);
   fChain->SetBranchAddress("jetBetaStarCMGExt", jetBetaStarCMGExt, &b_jetBetaStarCMGExt);
   fChain->SetBranchAddress("jetBetaStarClassicExt", jetBetaStarClassicExt, &b_jetBetaStarClassicExt);
   fChain->SetBranchAddress("jetNNeutrals", jetNNeutrals, &b_jetNNeutrals);
   fChain->SetBranchAddress("jetNCharged", jetNCharged, &b_jetNCharged);
   fChain->SetBranchAddress("jetMVAs", jetMVAs, &b_jetMVAs);
   fChain->SetBranchAddress("jetWPLevels", jetWPLevels, &b_jetWPLevels);
   fChain->SetBranchAddress("jetMVAsExt", jetMVAsExt, &b_jetMVAsExt);
   fChain->SetBranchAddress("jetWPLevelsExt", jetWPLevelsExt, &b_jetWPLevelsExt);
   fChain->SetBranchAddress("jetMt", jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetJECUnc", jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetLeadTrackPt", jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetVtxPt", jetVtxPt, &b_jetVtxPt);
   fChain->SetBranchAddress("jetVtxMass", jetVtxMass, &b_jetVtxMass);
   fChain->SetBranchAddress("jetVtx3dL", jetVtx3dL, &b_jetVtx3dL);
   fChain->SetBranchAddress("jetVtx3deL", jetVtx3deL, &b_jetVtx3deL);
   fChain->SetBranchAddress("jetSoftLeptPt", jetSoftLeptPt, &b_jetSoftLeptPt);
   fChain->SetBranchAddress("jetSoftLeptPtRel", jetSoftLeptPtRel, &b_jetSoftLeptPtRel);
   fChain->SetBranchAddress("jetSoftLeptdR", jetSoftLeptdR, &b_jetSoftLeptdR);
   fChain->SetBranchAddress("jetSoftLeptIdlooseMu", jetSoftLeptIdlooseMu, &b_jetSoftLeptIdlooseMu);
   fChain->SetBranchAddress("jetSoftLeptIdEle95", jetSoftLeptIdEle95, &b_jetSoftLeptIdEle95);
   fChain->SetBranchAddress("jetDPhiMETJet", jetDPhiMETJet, &b_jetDPhiMETJet);
   fChain->SetBranchAddress("jetPuJetIdL", jetPuJetIdL, &b_jetPuJetIdL);
   fChain->SetBranchAddress("jetPuJetIdM", jetPuJetIdM, &b_jetPuJetIdM);
   fChain->SetBranchAddress("jetPuJetIdT", jetPuJetIdT, &b_jetPuJetIdT);
   fChain->SetBranchAddress("nLowPtJet", &nLowPtJet, &b_nLowPtJet);
   fChain->SetBranchAddress("jetLowPtEn", jetLowPtEn, &b_jetLowPtEn);
   fChain->SetBranchAddress("jetLowPtPt", jetLowPtPt, &b_jetLowPtPt);
   fChain->SetBranchAddress("jetLowPtEta", jetLowPtEta, &b_jetLowPtEta);
   fChain->SetBranchAddress("jetLowPtPhi", jetLowPtPhi, &b_jetLowPtPhi);
   fChain->SetBranchAddress("jetLowPtCharge", jetLowPtCharge, &b_jetLowPtCharge);
   fChain->SetBranchAddress("jetLowPtEt", jetLowPtEt, &b_jetLowPtEt);
   fChain->SetBranchAddress("jetLowPtRawPt", jetLowPtRawPt, &b_jetLowPtRawPt);
   fChain->SetBranchAddress("jetLowPtRawEn", jetLowPtRawEn, &b_jetLowPtRawEn);
   fChain->SetBranchAddress("jetLowPtArea", jetLowPtArea, &b_jetLowPtArea);
   fChain->SetBranchAddress("nConv", &nConv, &b_nConv);
   fChain->SetBranchAddress("convP4", convP4, &b_convP4);
   fChain->SetBranchAddress("convVtx", convVtx, &b_convVtx);
   fChain->SetBranchAddress("convVtxErr", convVtxErr, &b_convVtxErr);
   fChain->SetBranchAddress("convPairMomentum", convPairMomentum, &b_convPairMomentum);
   fChain->SetBranchAddress("convRefittedMomentum", convRefittedMomentum, &b_convRefittedMomentum);
   fChain->SetBranchAddress("convNTracks", convNTracks, &b_convNTracks);
   fChain->SetBranchAddress("convPairInvMass", convPairInvMass, &b_convPairInvMass);
   fChain->SetBranchAddress("convPairCotThetaSep", convPairCotThetaSep, &b_convPairCotThetaSep);
   fChain->SetBranchAddress("convEoverP", convEoverP, &b_convEoverP);
   fChain->SetBranchAddress("convDistOfMinApproach", convDistOfMinApproach, &b_convDistOfMinApproach);
   fChain->SetBranchAddress("convDPhiTrksAtVtx", convDPhiTrksAtVtx, &b_convDPhiTrksAtVtx);
   fChain->SetBranchAddress("convDPhiTrksAtEcal", convDPhiTrksAtEcal, &b_convDPhiTrksAtEcal);
   fChain->SetBranchAddress("convDEtaTrksAtEcal", convDEtaTrksAtEcal, &b_convDEtaTrksAtEcal);
   fChain->SetBranchAddress("convDxy", convDxy, &b_convDxy);
   fChain->SetBranchAddress("convDz", convDz, &b_convDz);
   fChain->SetBranchAddress("convLxy", convLxy, &b_convLxy);
   fChain->SetBranchAddress("convLz", convLz, &b_convLz);
   fChain->SetBranchAddress("convZofPrimVtxFromTrks", convZofPrimVtxFromTrks, &b_convZofPrimVtxFromTrks);
   fChain->SetBranchAddress("convNHitsBeforeVtx", convNHitsBeforeVtx, &b_convNHitsBeforeVtx);
   fChain->SetBranchAddress("convNSharedHits", convNSharedHits, &b_convNSharedHits);
   fChain->SetBranchAddress("convValidVtx", convValidVtx, &b_convValidVtx);
   fChain->SetBranchAddress("convMVALikelihood", convMVALikelihood, &b_convMVALikelihood);
   fChain->SetBranchAddress("convChi2", convChi2, &b_convChi2);
   fChain->SetBranchAddress("convChi2Probability", convChi2Probability, &b_convChi2Probability);
   fChain->SetBranchAddress("convTk1Dz", convTk1Dz, &b_convTk1Dz);
   fChain->SetBranchAddress("convTk2Dz", convTk2Dz, &b_convTk2Dz);
   fChain->SetBranchAddress("convTk1DzErr", convTk1DzErr, &b_convTk1DzErr);
   fChain->SetBranchAddress("convTk2DzErr", convTk2DzErr, &b_convTk2DzErr);
   fChain->SetBranchAddress("convCh1Ch2", convCh1Ch2, &b_convCh1Ch2);
   fChain->SetBranchAddress("convTk1D0", convTk1D0, &b_convTk1D0);
   fChain->SetBranchAddress("convTk1Pout", convTk1Pout, &b_convTk1Pout);
   fChain->SetBranchAddress("convTk1Pin", convTk1Pin, &b_convTk1Pin);
   fChain->SetBranchAddress("convTk2D0", convTk2D0, &b_convTk2D0);
   fChain->SetBranchAddress("convTk2Pout", convTk2Pout, &b_convTk2Pout);
   fChain->SetBranchAddress("convTk2Pin", convTk2Pin, &b_convTk2Pin);
   Notify();
}

Bool_t xAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void xAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t xAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void xAna::InitTreePass() {

  PassingWP85Tree = new TTree("PassingWP85Tree", "Passing WP85 probes");

  Weight = 1;
  nGoodVtx = 1;
  Tag_isMCEle = false;
  Probe_isMCEle = false;

  PassingWP85Tree->Branch("run", &run, "run/I");
  PassingWP85Tree->Branch("event", &event, "event/I");
  PassingWP85Tree->Branch("Weight", &Weight, "Weight/D");
  PassingWP85Tree->Branch("nGoodVtx", &nGoodVtx, "nGoodVtx/I");
  PassingWP85Tree->Branch("Tag_Pt", &Tag_Pt, "Tag_Pt/D");
  PassingWP85Tree->Branch("Tag_Eta", &Tag_Eta, "Tag_Eta/D");
  PassingWP85Tree->Branch("Tag_SCEta", &Tag_SCEta, "Tag_SCEta/D");
  PassingWP85Tree->Branch("Tag_isMCEle",&Tag_isMCEle,"Tag_isMCEle/O");
  PassingWP85Tree->Branch("Tag_Phi", &Tag_Phi, "Tag_Pt/D");
  PassingWP85Tree->Branch("Tag_charge", &Tag_charge, "Tag_charge/I");
  PassingWP85Tree->Branch("Tag_status", &Tag_status, "Tag_status/I");
  PassingWP85Tree->Branch("Probe_Pt", &Probe_Pt, "Probe_Pt/D");
  PassingWP85Tree->Branch("Probe_Eta", &Probe_Eta, "Probe_Eta/D");
  PassingWP85Tree->Branch("Probe_SCEta", &Probe_SCEta, "Probe_SCEta/D");
  PassingWP85Tree->Branch("Probe_isMCEle",&Probe_isMCEle,"Probe_isMCEle/O");
  PassingWP85Tree->Branch("Probe_Phi", &Probe_Phi, "Probe_Pt/D");
  PassingWP85Tree->Branch("Probe_charge", &Probe_charge, "Probe_charge/I");
  PassingWP85Tree->Branch("Probe_status", &Probe_status, "Probe_status/I");
  PassingWP85Tree->Branch("Zmass", &Zmass, "Zmass/D");
  PassingWP85Tree->Branch("cent", &cent, "cent/D");
  PassingWP85Tree->Branch("elepasstrg", &elepasstrg, "elepasstrg/I");
  PassingWP85Tree->Branch("elefailtrg", &elefailtrg, "elefailtrg/I");
}

void xAna::InitTreeFail() {

  FailingWP85Tree = new TTree("FailingWP85Tree", "Failing WP85 probes");

  Weight = 1;
  nGoodVtx = 1;
  Tag_isMCEle = false;
  Probe_isMCEle= false;

  FailingWP85Tree->Branch("run", &run, "run/I");
  FailingWP85Tree->Branch("event", &event, "event/I");
  FailingWP85Tree->Branch("Weight", &Weight, "Weight/D");
  FailingWP85Tree->Branch("nGoodVtx", &nGoodVtx, "nGoodVtx/I");
  FailingWP85Tree->Branch("Tag_Pt", &Tag_Pt, "Tag_Pt/D");
  FailingWP85Tree->Branch("Tag_Eta", &Tag_Eta, "Tag_Eta/D");
  FailingWP85Tree->Branch("Tag_SCEta", &Tag_SCEta, "Tag_SCEta/D");
  FailingWP85Tree->Branch("Tag_isMCEle",&Tag_isMCEle,"Tag_isMCEle/O");
  FailingWP85Tree->Branch("Tag_Phi", &Tag_Phi, "Tag_Pt/D");
  FailingWP85Tree->Branch("Tag_charge", &Tag_charge, "Tag_charge/I");
  FailingWP85Tree->Branch("Tag_status", &Tag_status, "Tag_status/I");
  FailingWP85Tree->Branch("Probe_Pt", &Probe_Pt, "Probe_Pt/D");
  FailingWP85Tree->Branch("Probe_Eta", &Probe_Eta, "Probe_Eta/D");
  FailingWP85Tree->Branch("Probe_SCEta", &Probe_SCEta, "Probe_SCEta/D");
  FailingWP85Tree->Branch("Probe_isMCEle",&Probe_isMCEle,"Probe_isMCEle/O");
  FailingWP85Tree->Branch("Probe_Phi", &Probe_Phi, "Probe_Pt/D");
  FailingWP85Tree->Branch("Probe_charge", &Probe_charge, "Probe_charge/I");
  FailingWP85Tree->Branch("Probe_status", &Probe_status, "Probe_status/I");
  FailingWP85Tree->Branch("Zmass", &Zmass, "Zmass/D");
  FailingWP85Tree->Branch("cent", &cent, "cent/D");
  FailingWP85Tree->Branch("elepasstrg", &elepasstrg, "elepasstrg/I");
  FailingWP85Tree->Branch("elefailtrg", &elefailtrg, "elefailtrg/I");
}

#endif // #ifdef xAna_cxx

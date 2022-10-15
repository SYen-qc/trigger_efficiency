#ifndef xAna_tag_pro_h
#define xAna_tag_pro_h

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TChainElement.h>

#include <iostream>
using namespace std;

const Int_t kMaxeleESEffSigmaRR = 1;
const Int_t kMaxphoESEffSigmaRR = 1;
const Int_t kMaxnPFEle = 4;
const Int_t kMaxPFElePt = 1;
const Int_t kMaxPFEleEta = 1;
const Int_t kMaxPFElePhi = 1;
const Int_t kMaxPFEleEn = 1;
const Int_t kMaxPFEleCharge = 1;

const Int_t maxP = 500;

class xAna_tag_pro {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nHLT;
   Int_t           HLT[maxP];   //[nHLT]
   Int_t           HLTIndex[maxP];
   Float_t         bspotPos[maxP];
   Int_t           nVtx;
   Float_t         vtx[maxP][3];   //[nVtx]
   Int_t           IsVtxGood;
   Int_t           nGoodVtx;
   Float_t         centrality[5];
   Int_t           nVtxBS;
   Float_t         vtxbs[maxP][3];   //[nVtxBS]
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
   Float_t         trkMETx[maxP];   //[nVtxBS]
   Float_t         trkMETy[maxP];   //[nVtxBS]
   Float_t         trkMETPhi[maxP];   //[nVtxBS]
   Float_t         trkMET[maxP];   //[nVtxBS]
   Int_t           metFilters[maxP];
   Int_t           nEle;
   Int_t           eleTrg[9][maxP];   //[nEle]
   Int_t           eleClass[maxP];   //[nEle]
   Int_t           eleIsEcalDriven[maxP];   //[nEle]
   Int_t           eleCharge[maxP];   //[nEle]
   Float_t         eleEn[maxP];   //[nEle]
   Float_t         eleEcalEn[maxP];   //[nEle]
   Float_t         eleSCRawEn[maxP];   //[nEle]
   Float_t         eleSCEn[maxP];   //[nEle]
   Float_t         eleESEn[maxP9];   //[nEle]
   Float_t         elePt[maxP];   //[nEle]
   Float_t         eleEta[maxP];   //[nEle]
   Float_t         elePhi[maxP];   //[nEle]
   Float_t         eleEtaVtx[maxP][100];   //[nEle]
   Float_t         elePhiVtx[maxP][100];   //[nEle]
   Float_t         eleEtVtx[maxP][100];   //[nEle]
   Float_t         eleSCEta[maxP];   //[nEle]
   Float_t         eleSCPhi[maxP];   //[nEle]
   Float_t         eleSCEtaWidth[maxP];   //[nEle]
   Float_t         eleSCPhiWidth[maxP];   //[nEle]
   Float_t         eleVtx[maxP][3];   //[nEle]
   Float_t         eleD0[maxP];   //[nEle]
   Float_t         eleDz[maxP];   //[nEle]
   Float_t         eleD0GV[maxP];   //[nEle]
   Float_t         eleDzGV[maxP];   //[nEle]
   Float_t         eleD0Vtx[maxP][100];   //[nEle]
   Float_t         eleDzVtx[maxP][100];   //[nEle]
   Float_t         eleHoverE[maxP];   //[nEle]
   Float_t         eleHoverE12[maxP9];   //[nEle]
   Float_t         eleEoverP[maxP];   //[nEle]
   Float_t         elePin[maxP];   //[nEle]
   Float_t         elePout[maxP];   //[nEle]
   Float_t         eleTrkMomErr[maxP];   //[nEle]
   Float_t         eleBrem[maxP];   //[nEle]
   Float_t         eledEtaAtVtx[maxP];   //[nEle]
   Float_t         eledPhiAtVtx[maxP];   //[nEle]
   Float_t         eleSigmaIEtaIEta[maxP];   //[nEle]
   Float_t         eleSigmaIEtaIPhi[maxP];   //[nEle]
   Float_t         eleSigmaIPhiIPhi[maxP];   //[nEle]
   Float_t         eleEmax[maxP];   //[nEle]
   Float_t         eleE1x5[maxP];   //[nEle]
   Float_t         eleE3x3[maxP];   //[nEle]
   Float_t         eleE5x5[maxP];   //[nEle]
   Float_t         eleE2x5Max[maxP];   //[nEle]
   Float_t         eleRegrE[maxP];   //[nEle]
   Float_t         eleRegrEerr[maxP];   //[nEle]
   Float_t         elePhoRegrE[maxP];   //[nEle]
   Float_t         elePhoRegrEerr[maxP];   //[nEle]
   Float_t         eleSeedTime[maxP];   //[nEle]
   Int_t           eleRecoFlag[maxP];   //[nEle]
   Int_t           elePos[maxP];   //[nEle]
   Float_t         eleIsoTrkDR03[maxP];   //[nEle]
   Float_t         eleIsoEcalDR03[maxP];   //[nEle]
   Float_t         eleIsoHcalDR03[maxP];   //[nEle]
   Float_t         eleIsoHcalDR0312[maxP];   //[nEle]
   Float_t         eleIsoTrkDR04[maxP];   //[nEle]
   Float_t         eleIsoEcalDR04[maxP];   //[nEle]
   Float_t         eleIsoHcalDR04[maxP];   //[nEle]
   Float_t         eleIsoHcalDR0412[maxP];   //[nEle]
   Float_t         eleModIsoTrk[maxP];   //[nEle]
   Float_t         eleModIsoEcal[maxP];   //[nEle]
   Float_t         eleModIsoHcal[maxP];   //[nEle]
   Int_t           eleMissHits[maxP];   //[nEle]
   Float_t         eleConvDist[maxP];   //[nEle]
   Float_t         eleConvDcot[maxP];   //[nEle]
   Int_t           eleConvVtxFit[maxP];   //[nEle]
   Float_t         eleIP3D[maxP];   //[nEle]
   Float_t         eleIP3DErr[maxP];   //[nEle]
   Float_t         eleIDMVANonTrig[maxP];   //[nEle]
   Float_t         eleIDMVATrig[maxP];   //[nEle]
   Float_t         elePFChIso03[maxP];   //[nEle]
   Float_t         elePFPhoIso03[maxP];   //[nEle]
   Float_t         elePFNeuIso03[maxP];   //[nEle]
   Float_t         elePFChIso04[maxP];   //[nEle]
   Float_t         elePFPhoIso04[maxP];   //[nEle]
   Float_t         elePFNeuIso04[maxP];   //[nEle]
   Float_t         eleESEffSigmaRR[maxP][3];   //[nEle]
   Int_t           nPho;
   Int_t           phoTrg[maxP][8];   //[nPho]
   Int_t           phoTrgFilter[maxP][50];   //[nPho]
   Bool_t          phoIsPhoton[maxP];   //[nPho]
   Float_t         phoSCPos[maxP][3];   //[nPho]
   Float_t         phoCaloPos[maxP][3];   //[nPho]
   Float_t         phoE[maxP];   //[nPho]
   Float_t         phoEt[maxP];   //[nPho]
   Float_t         phoEta[maxP];   //[nPho]
   Float_t         phoVtx[maxP][3];   //[nPho]
   Float_t         phoPhi[maxP];   //[nPho]
   Float_t         phoEtVtx[maxP][100];   //[nPho]
   Float_t         phoEtaVtx[maxP][100];   //[nPho]
   Float_t         phoPhiVtx[maxP][100];   //[nPho]
   Float_t         phoR9[maxP];   //[nPho]
   Float_t         phoTrkIsoHollowDR03[maxP];   //[nPho]
   Float_t         phoEcalIsoDR03[maxP];   //[nPho]
   Float_t         phoHcalIsoDR03[maxP];   //[nPho]
   Float_t         phoHcalIsoDR0312[maxP];   //[nPho]
   Float_t         phoTrkIsoHollowDR04[maxP];   //[nPho]
   Float_t         phoCiCTrkIsoDR03[maxP][100];   //[nPho]
   Float_t         phoCiCTrkIsoDR04[maxP][100];   //[nPho]
   Float_t         phoCiCdRtoTrk[maxP];   //[nPho]
   Float_t         phoEcalIsoDR04[maxP];   //[nPho]
   Float_t         phoHcalIsoDR04[maxP];   //[nPho]
   Float_t         phoHcalIsoDR0412[maxP];   //[nPho]
   Float_t         phoHoverE[maxP];   //[nPho]
   Float_t         phoHoverE12[maxP];   //[nPho]
   Int_t           phoEleVeto[maxP];   //[nPho]
   Float_t         phoSigmaIEtaIEta[maxP];   //[nPho]
   Float_t         phoSigmaIEtaIPhi[maxP];   //[nPho]
   Float_t         phoSigmaIPhiIPhi[maxP];   //[nPho]
   Float_t         phoCiCPF4phopfIso03[maxP];   //[nPho]
   Float_t         phoCiCPF4phopfIso04[maxP];   //[nPho]
   Float_t         phoCiCPF4chgpfIso02[maxP][100];   //[nPho]
   Float_t         phoCiCPF4chgpfIso03[maxP][100];   //[nPho]
   Float_t         phoCiCPF4chgpfIso04[maxP][100];   //[nPho]
   Float_t         phoEmax[maxP];   //[nPho]
   Float_t         phoE3x3[maxP];   //[nPho]
   Float_t         phoE3x1[maxP];   //[nPho]
   Float_t         phoE1x3[maxP];   //[nPho]
   Float_t         phoE5x5[maxP];   //[nPho]
   Float_t         phoE1x5[maxP];   //[nPho]
   Float_t         phoE2x2[maxP];   //[nPho]
   Float_t         phoE2x5Max[maxP];   //[nPho]
   Float_t         phoPFChIso[maxP];   //[nPho]
   Float_t         phoPFPhoIso[maxP];   //[nPho]
   Float_t         phoPFNeuIso[maxP];   //[nPho]
   Float_t         phoSCRChIso[maxP];   //[nPho]
   Float_t         phoSCRPhoIso[maxP];   //[nPho]
   Float_t         phoSCRNeuIso[maxP];   //[nPho]
   Float_t         phoRegrE[maxP];   //[nPho]
   Float_t         phoRegrEerr[maxP];   //[nPho]
   Float_t         phoSeedTime[maxP];   //[nPho]
   Int_t           phoSeedDetId1[maxP];   //[nPho]
   Int_t           phoSeedDetId2[maxP];   //[nPho]
   Float_t         phoLICTD[maxP];   //[nPho]
   Int_t           phoRecoFlag[maxP];   //[nPho]
   Int_t           phoPos[maxP];   //[nPho]
   Float_t         phoSCE[maxP];   //[nPho]
   Float_t         phoSCRawE[maxP];   //[nPho]
   Float_t         phoESEn[maxP];   //[nPho]
   Float_t         phoSCEt[maxP];   //[nPho]
   Float_t         phoSCEta[maxP];   //[nPho]
   Float_t         phoSCPhi[maxP];   //[nPho]
   Float_t         phoSCEtaWidth[maxP];   //[nPho]
   Float_t         phoSCPhiWidth[maxP];   //[nPho]
   Float_t         phoSCBrem[maxP];   //[nPho]
   Int_t           phoOverlap[maxP];   //[nPho]
   Int_t           phohasPixelSeed[maxP];   //[nPho]
   Int_t           pho_hasConvPf[maxP];   //[nPho]
   Int_t           pho_hasSLConvPf[maxP];   //[nPho]
   Float_t         pho_pfconvVtxZ[maxP];   //[nPho]
   Float_t         pho_pfconvVtxZErr[maxP];   //[nPho]
   Int_t           pho_nSLConv[maxP];   //[nPho]
   Float_t         pho_pfSLConvPos[maxP][3];   //[nPho]
   Float_t         pho_pfSLConvVtxZ[maxP][20];   //[nPho]
   Int_t           phoIsConv[maxP];   //[nPho]
   Int_t           phoNConv[maxP];   //[nPho]
   Float_t         phoConvInvMass[maxP];   //[nPho]
   Float_t         phoConvCotTheta[maxP];   //[nPho]
   Float_t         phoConvEoverP[maxP];   //[nPho]
   Float_t         phoConvZofPVfromTrks[maxP];   //[nPho]
   Float_t         phoConvMinDist[maxP];   //[nPho]
   Float_t         phoConvdPhiAtVtx[maxP];   //[nPho]
   Float_t         phoConvdPhiAtCalo[maxP];   //[nPho]
   Float_t         phoConvdEtaAtCalo[maxP];   //[nPho]
   Float_t         phoConvTrkd0[maxP][2];   //[nPho]
   Float_t         phoConvTrkPin[maxP][2];   //[nPho]
   Float_t         phoConvTrkPout[maxP][2];   //[nPho]
   Float_t         phoConvTrkdz[maxP][2];   //[nPho]
   Float_t         phoConvTrkdzErr[maxP][2];   //[nPho]
   Float_t         phoConvChi2[maxP];   //[nPho]
   Float_t         phoConvChi2Prob[maxP];   //[nPho]
   Int_t           phoConvNTrks[maxP];   //[nPho]
   Float_t         phoConvCharge[maxP][2];   //[nPho]
   Float_t         phoConvValidVtx[maxP];   //[nPho]
   Float_t         phoConvLikeLihood[maxP];   //[nPho]
   Float_t         phoConvP4[10][4];   //[nPho]
   Float_t         phoConvVtx[10][3];   //[nPho]
   Float_t         phoConvVtxErr[10][3];   //[nPho]
   Float_t         phoConvPairMomentum[10][3];   //[nPho]
   Float_t         phoConvRefittedMomentum[10][3];   //[nPho]
   Int_t           SingleLegConv[maxP];   //[nPho]
   Float_t         phoPFConvVtx[maxP][3];   //[nPho]
   Float_t         phoPFConvMom[maxP][3];   //[nPho]
   Float_t         phoESEffSigmaRR[maxP][3];   //[nPho]
   Int_t           nMu;
   Int_t           muTrg[maxP][10];   //[nMu]
   Float_t         muEta[maxP];   //[nMu]
   Float_t         muPhi[maxP];   //[nMu]
   Int_t           muCharge[maxP];   //[nMu]
   Float_t         muPt[maxP];   //[nMu]
   Float_t         muPz[maxP];   //[nMu]
   Float_t         muVtx[maxP][3];   //[nMu]
   Float_t         muVtxGlb[maxP][3];   //[nMu]
   Float_t         mucktPt[maxP];   //[nMu]
   Float_t         mucktPtErr[maxP];   //[nMu]
   Float_t         mucktEta[maxP];   //[nMu]
   Float_t         mucktPhi[maxP];   //[nMu]
   Float_t         mucktdxy[maxP];   //[nMu]
   Float_t         mucktdz[maxP];   //[nMu]
   Float_t         muIsoTrk[maxP];   //[nMu]
   Float_t         muIsoCalo[maxP];   //[nMu]
   Float_t         muIsoEcal[maxP];   //[nMu]
   Float_t         muIsoHcal[maxP];   //[nMu]
   Float_t         muChi2NDF[maxP];   //[nMu]
   Float_t         muInnerChi2NDF[maxP];   //[nMu]
   Float_t         muPFIsoR04_CH[maxP];   //[nMu]
   Float_t         muPFIsoR04_NH[maxP];   //[nMu]
   Float_t         muPFIsoR04_Pho[maxP];   //[nMu]
   Float_t         muPFIsoR04_PU[maxP];   //[nMu]
   Float_t         muPFIsoR04_CPart[maxP];   //[nMu]
   Float_t         muPFIsoR04_NHHT[maxP];   //[nMu]
   Float_t         muPFIsoR04_PhoHT[maxP];   //[nMu]
   Float_t         muPFIsoR03_CH[maxP];   //[nMu]
   Float_t         muPFIsoR03_NH[maxP];   //[nMu]
   Float_t         muPFIsoR03_Pho[maxP];   //[nMu]
   Float_t         muPFIsoR03_PU[maxP];   //[nMu]
   Float_t         muPFIsoR03_CPart[maxP];   //[nMu]
   Float_t         muPFIsoR03_NHHT[maxP];   //[nMu]
   Float_t         muPFIsoR03_PhoHT[maxP];   //[nMu]
   Int_t           muType[maxP];   //[nMu]
   Float_t         muD0[maxP];   //[nMu]
   Float_t         muDz[maxP];   //[nMu]
   Float_t         muD0GV[maxP];   //[nMu]
   Float_t         muDzGV[maxP];   //[nMu]
   Float_t         muD0Vtx[maxP][100];   //[nMu]
   Float_t         muDzVtx[maxP][100];   //[nMu]
   Float_t         muInnerD0[maxP];   //[nMu]
   Float_t         muInnerDz[maxP];   //[nMu]
   Float_t         muInnerD0GV[maxP];   //[nMu]
   Float_t         muInnerDzGV[maxP];   //[nMu]
   Float_t         muInnerPt[maxP];   //[nMu]
   Float_t         muInnerPtErr[maxP];   //[nMu]
   Int_t           muNumberOfValidTrkLayers[maxP];   //[nMu]
   Int_t           muNumberOfValidTrkHits[maxP];   //[nMu]
   Int_t           muNumberOfValidPixelLayers[maxP];   //[nMu]
   Int_t           muNumberOfValidPixelHits[maxP];   //[nMu]
   Int_t           muNumberOfValidMuonHits[maxP];   //[nMu]
   Int_t           muStations[maxP];   //[nMu]
   Int_t           muChambers[maxP];   //[nMu]
   Float_t         muIP3D[maxP];   //[nMu]
   Float_t         muIP3DErr[maxP];   //[nMu]
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
   Int_t           jetTrg[maxP][14];   //[nJet]
   Float_t         jetEn[maxP];   //[nJet]
   Float_t         jetPt[maxP];   //[nJet]
   Float_t         jetEta[maxP];   //[nJet]
   Float_t         jetPhi[maxP];   //[nJet]
   Float_t         jetCharge[maxP];   //[nJet]
   Float_t         jetEt[maxP];   //[nJet]
   Float_t         jetRawPt[maxP];   //[nJet]
   Float_t         jetRawEn[maxP];   //[nJet]
   Float_t         jetArea[maxP];   //[nJet]
   Float_t         jetCHF[maxP];   //[nJet]
   Float_t         jetNHF[maxP];   //[nJet]
   Float_t         jetCEF[maxP];   //[nJet]
   Float_t         jetNEF[maxP];   //[nJet]
   Int_t           jetNCH[maxP];   //[nJet]
   Float_t         jetHFHAE[maxP];   //[nJet]
   Float_t         jetHFEME[maxP];   //[nJet]
   Int_t           jetNConstituents[maxP];   //[nJet]
   Float_t         jetCombinedSecondaryVtxBJetTags[maxP];   //[nJet]
   Float_t         jetCombinedSecondaryVtxMVABJetTags[maxP];   //[nJet]
   Float_t         jetJetProbabilityBJetTags[maxP];   //[nJet]
   Float_t         jetJetBProbabilityBJetTags[maxP];   //[nJet]
   Float_t         jetTrackCountingHighPurBJetTags[maxP];   //[nJet]
   Float_t         jetBetaStar[maxP][100];   //[nJet]
   Bool_t          jetPFLooseId[maxP];   //[nJet]
   Float_t         jetDRMean[maxP];   //[nJet]
   Float_t         jetDR2Mean[maxP];   //[nJet]
   Float_t         jetDZ[maxP];   //[nJet]
   Float_t         jetFrac01[maxP];   //[nJet]
   Float_t         jetFrac02[maxP];   //[nJet]
   Float_t         jetFrac03[maxP];   //[nJet]
   Float_t         jetFrac04[maxP];   //[nJet]
   Float_t         jetFrac05[maxP];   //[nJet]
   Float_t         jetFrac06[maxP];   //[nJet]
   Float_t         jetFrac07[maxP];   //[nJet]
   Float_t         jetBeta[maxP];   //[nJet]
   Float_t         jetBetaStarCMG[maxP];   //[nJet]
   Float_t         jetBetaStarClassic[maxP];   //[nJet]
   Float_t         jetBetaExt[maxP][100];   //[nJet]
   Float_t         jetBetaStarCMGExt[maxP][100];   //[nJet]
   Float_t         jetBetaStarClassicExt[maxP][100];   //[nJet]
   Float_t         jetNNeutrals[maxP];   //[nJet]
   Float_t         jetNCharged[maxP];   //[nJet]
   Float_t         jetMVAs[maxP][4];   //[nJet]
   Int_t           jetWPLevels[maxP][4];   //[nJet]
   Float_t         jetMVAsExt[maxP][4][100];   //[nJet]
   Int_t           jetWPLevelsExt[maxP][4][100];   //[nJet]
   Float_t         jetMt[maxP];   //[nJet]
   Float_t         jetJECUnc[maxP];   //[nJet]
   Float_t         jetLeadTrackPt[maxP];   //[nJet]
   Float_t         jetVtxPt[maxP];   //[nJet]
   Float_t         jetVtxMass[maxP];   //[nJet]
   Float_t         jetVtx3dL[maxP];   //[nJet]
   Float_t         jetVtx3deL[maxP];   //[nJet]
   Float_t         jetSoftLeptPt[maxP];   //[nJet]
   Float_t         jetSoftLeptPtRel[maxP];   //[nJet]
   Float_t         jetSoftLeptdR[maxP];   //[nJet]
   Float_t         jetSoftLeptIdlooseMu[maxP];   //[nJet]
   Float_t         jetSoftLeptIdEle95[maxP];   //[nJet]
   Float_t         jetDPhiMETJet[maxP];   //[nJet]
   Float_t         jetPuJetIdL[maxP];   //[nJet]
   Float_t         jetPuJetIdM[maxP];   //[nJet]
   Float_t         jetPuJetIdT[maxP];   //[nJet]
   Int_t           nLowPtJet;
   Float_t         jetLowPtEn[maxP];   //[nLowPtJet]
   Float_t         jetLowPtPt[maxP];   //[nLowPtJet]
   Float_t         jetLowPtEta[maxP];   //[nLowPtJet]
   Float_t         jetLowPtPhi[maxP];   //[nLowPtJet]
   Float_t         jetLowPtCharge[maxP];   //[nLowPtJet]
   Float_t         jetLowPtEt[maxP];   //[nLowPtJet]
   Float_t         jetLowPtRawPt[maxP];   //[nLowPtJet]
   Float_t         jetLowPtRawEn[maxP];   //[nLowPtJet]
   Float_t         jetLowPtArea[maxP];   //[nLowPtJet]
   Int_t           nConv;
   Float_t         convP4[maxP][4];   //[nConv]
   Float_t         convVtx[maxP][3];   //[nConv]
   Float_t         convVtxErr[maxP][3];   //[nConv]
   Float_t         convPairMomentum[maxP][3];   //[nConv]
   Float_t         convRefittedMomentum[maxP][3];   //[nConv]
   Int_t           convNTracks[maxP];   //[nConv]
   Float_t         convPairInvMass[maxP]   //[nConv]
   Float_t         convPairCotThetaSep[maxP];   //[nConv]
   Float_t         convEoverP[maxP];   //[nConv]
   Float_t         convDistOfMinApproach[maxP];   //[nConv]
   Float_t         convDPhiTrksAtVtx[maxP];   //[nConv]
   Float_t         convDPhiTrksAtEcal[maxP];   //[nConv]
   Float_t         convDEtaTrksAtEcal[maxP];   //[nConv]
   Float_t         convDxy[maxP];   //[nConv]
   Float_t         convDz[maxP];   //[nConv]
   Float_t         convLxy[maxP];   //[nConv]
   Float_t         convLz[maxP];   //[nConv]
   Float_t         convZofPrimVtxFromTrks[maxP];   //[nConv]
   Int_t           convNHitsBeforeVtx[maxP][2];   //[nConv]
   Int_t           convNSharedHits[maxP];   //[nConv]
   Int_t           convValidVtx[maxP];   //[nConv]
   Float_t         convMVALikelihood[maxP];   //[nConv]
   Float_t         convChi2[maxP];   //[nConv]
   Float_t         convChi2Probability[maxP];   //[nConv]
   Float_t         convTk1Dz[maxP];   //[nConv]
   Float_t         convTk2Dz[maxP];   //[nConv]
   Float_t         convTk1DzErr[maxP];   //[nConv]
   Float_t         convTk2DzErr[maxP];   //[nConv]
   Int_t           convCh1Ch2[maxP];   //[nConv]
   Float_t         convTk1D0[maxP];   //[nConv]
   Float_t         convTk1Pout[maxP];   //[nConv]
   Float_t         convTk1Pin[maxP];   //[nConv]
   Float_t         convTk2D0[maxP];   //[nConv]
   Float_t         convTk2Pout[maxP];   //[nConv]
   Float_t         convTk2Pin[maxP];   //[nConv]

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
   
   xAna_tag_pro(Int_t mode = 0, Int_t lepChannel = 0, TString dirname = "/scratch/cmkuo/", TString filename="job_photon_May23rereco_skim.root", Int_t skim=0);
   virtual ~xAna_tag_pro();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     InitHists();
   virtual void     InitTree();
   virtual void     InitTreePass();
   virtual void     InitTreeFail();
   virtual Double_t deltaPhi(Double_t phi1, Double_t phi2);
   virtual Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
   virtual vector<float> getPUTrueWeight(TH1* hDataPU, TH1* hMCPU);
   // Lepton selection
   virtual Int_t    eleKin(Int_t i = 0); 
   virtual Float_t  geteIDEA(Int_t i = 0);
   virtual Float_t  getePFIso(Int_t i = 0);
   virtual Int_t    eID2012(Int_t i = 0, Int_t idWP = 0);
   virtual Int_t    muID2012(Int_t i = 0);
   virtual Int_t    jetID(Int_t i = 0); 

   Int_t mode_;
   Int_t lepChannel_;
   Int_t skim_;

   TFile *fout_;
   TTree *outtree_;
   TTree *PassingWP85Tree;
   TTree *FailingWP85Tree;
   TString filename_;

   TH1F *hEvents_;
   TH1F *hMCPU_;
   TH1F *hDataPU_;

   TH1F *hnEle_;
   TH1F *hnMu_;
   TH1F *hnJets_;
   TH1F *hJetPt_;
   TH1F *hmZ_;
   TH1F *hJZPt_;
   TH1F *hJZPhi_;
     
   Int_t   trg1_;
   Int_t   trg2_;
   Float_t cent_;
   Float_t met_;
   Float_t metx_;
   Float_t mety_;
   Float_t mWT_;
   Float_t mZ_;
   Float_t ptZ_;
   Float_t etaZ_;
   Float_t phiZ_;
   Float_t yZ_;
   Int_t   nSelLep_;
   Int_t   ch1_;
   Float_t pt1_;
   Float_t eta1_;
   Float_t phi1_;
   Float_t px1_;
   Float_t py1_;
   Float_t pz1_;
   Float_t y1_;
   Float_t etacm1_;
   Int_t   ch2_;
   Float_t pt2_;
   Float_t eta2_;
   Float_t phi2_;
   Float_t px2_;
   Float_t py2_;
   Float_t pz2_;
   Float_t y2_;
   Float_t etacm2_;
   Int_t   nSelJets_;
   Float_t ptJ_[100];
   Float_t etaJ_[100];
   Float_t phiJ_[100];
   Float_t enJ_[100];
   Float_t ptCombJ_;
   Float_t etaCombJ_;
   Float_t phiCombJ_;
   Float_t enCombJ_;
   Float_t puwei_;
   Float_t iso1_[6];
   Float_t iso2_[6];
   Int_t elestatus_[100];
   Float_t TTmZ_;
   Float_t TTeta1_;
   Float_t TTphi1_;
   Float_t TTpt1_;
   Float_t TTcharge1_;
   Int_t TTstatus1_;
   Float_t TTeta2_;
   Float_t TTphi2_;
   Float_t TTpt2_;
   Float_t TTcharge2_;
   Int_t TTstatus2_;
   Float_t TPmZ_;
   Float_t TPeta1_;
   Float_t TPphi1_;
   Float_t TPpt1_;
   Float_t TPcharge1_;
   Int_t TPstatus1_;
   Float_t TPeta2_;
   Float_t TPphi2_;
   Float_t TPpt2_;
   Float_t TPcharge2_;
   Int_t TPstatus2_;
   Float_t TFmZ_;
   Float_t TFeta1_;
   Float_t TFphi1_;
   Float_t TFpt1_;
   Float_t TFcharge1_;
   Int_t TFstatus1_;
   Float_t TFeta2_;
   Float_t TFphi2_;
   Float_t TFpt2_;
   Float_t TFcharge2_;
   Int_t TFstatus2_;

   Double_t Weight;
   Double_t nGoodVtx;
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

};

//#endif

//#ifdef xAna_tag_pro_cxx
xAna_tag_pro::xAna_tag_pro(Int_t mode, Int_t lepChannel, TString dirname, TString filename, Int_t skim) {
  
  mode_ = mode; // 0 : data, 1 : MC    
  lepChannel_ = lepChannel; // 0 : electron, 1 : muon
                                                                                                                        
  skim_ = skim;
  filename_ = filename;
  
  TChain *chain = new TChain("ggNtuplizer/EventTree");
  cout<<"Opening ntuple : "<<dirname + filename<<endl;
  chain->Add(dirname + filename);
  
  Init(chain);

  // Get PU reference histogram
  if (mode != 0) {
    TFile *fpuData = new TFile("inputFiles/TruePU_69300_2012A.root", "read");
    hDataPU_ = (TH1F*) fpuData->Get("pileup");
    hMCPU_   = (TH1F*) hDataPU_->Clone();
    hMCPU_->Reset();
    hMCPU_->SetName("pileupMC"); 
    chain->Draw("puTrue[1]>>pileupMC","","goff");
  }

  // Get PU and # of events histogram  
  TIter nextfile(chain->GetListOfFiles());
  TChainElement *elem;
  while ((elem = (TChainElement*)nextfile())) {
    TFile *f = new TFile(elem->GetTitle());
    //hMCPU_ = (TH1F*) f->Get("ggNtuplizer/hPUTrue");
    //hMCPU_->SetDirectory(0);
    hEvents_ = (skim == 0) ? (TH1F*) f->Get("ggNtuplizer/hEvents") : (TH1F*) f->Get("ggNtuplizer/hskim");
    hEvents_->SetDirectory(0);
    delete f;
  }
  delete elem;

}

xAna_tag_pro::~xAna_tag_pro() {

   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t xAna_tag_pro::GetEntry(Long64_t entry) { 
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t xAna_tag_pro::LoadTree(Long64_t entry) {
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
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

Bool_t xAna_tag_pro::Notify() {
  return kTRUE;
}

void xAna_tag_pro::Show(Long64_t entry) {
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t xAna_tag_pro::Cut(Long64_t entry) {
  return 1;
}

vector<float> xAna_tag_pro::getPUTrueWeight(TH1* hDataPU, TH1* hMCPU) {

  vector<float> puWeight;

  hDataPU->Scale(1./hDataPU->Integral(-1, -1));
  hMCPU->Scale(1./hMCPU->Integral(-1, -1));

  Int_t nBinsX = hDataPU->GetXaxis()->GetNbins() ;
  for (Int_t iBin = 0; iBin < nBinsX; ++iBin) {
    if ( hMCPU->GetBinContent(iBin+1) > 0 )
      puWeight.push_back(hDataPU->GetBinContent(iBin+1)/hMCPU->GetBinContent(iBin+1));
    else puWeight.push_back(0.) ;
  }

  return puWeight;
}

Double_t xAna_tag_pro::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {

  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);

  return sqrt(dEta*dEta+dPhi*dPhi);
}

Double_t xAna_tag_pro::deltaPhi(Double_t phi1, Double_t phi2) {

  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();

  return dPhi;
}

void xAna_tag_pro::InitHists() {
  
  hnEle_   = new TH1F("hnEle",  "# of electrons",  10,  0,  10);
  hnMu_    = new TH1F("hnMu",   "# of muons",      10,  0,  10);
  hnJets_  = new TH1F("hnJets", "# of jets",       10,  0,  10);
  hJetPt_  = new TH1F("hJetPt", "jet pt",          40,  0, 200);
  hmZ_     = new TH1F("hmZ",    "Z mass",          55, 40, 150);
  hJZPt_   = new TH1F("hJZPt",  "Jet-Z Pt ratio",  50,  0,   5);
  hJZPhi_  = new TH1F("hJZPhi", "Jet-Z dPhi",      80,  4,  -4);
}

void xAna_tag_pro::InitTree() {
  
  outtree_ = new TTree("t", "mini tree");
  
  outtree_->Branch("run",        &run,         "run/I");
  outtree_->Branch("event",      &event,       "event/L");
  outtree_->Branch("trg1",       &trg1_,       "trg1/I");
  outtree_->Branch("trg2",       &trg2_,       "trg2/I");
  outtree_->Branch("cent",       &cent_,       "cent/F");
  outtree_->Branch("met",        &met_,        "met/F");
  outtree_->Branch("metx",       &metx_,       "metx/F");
  outtree_->Branch("mety",       &mety_,       "mety/F");
  outtree_->Branch("mWT",        &mWT_,        "mWT/F");
  outtree_->Branch("mZ",         &mZ_,         "mZ/F");
  outtree_->Branch("ptZ",        &ptZ_,        "ptZ/F");
  outtree_->Branch("etaZ",       &etaZ_,       "etaZ/F");
  outtree_->Branch("phiZ",       &phiZ_,       "phiZ/F");
  outtree_->Branch("yZ",         &yZ_,         "yZ/F");
  outtree_->Branch("nSelLep",    &nSelLep_,    "nSelLep/I");
  outtree_->Branch("ch1",        &ch1_,        "ch1/I");
  outtree_->Branch("pt1",        &pt1_,        "pt1/F");
  outtree_->Branch("eta1",       &eta1_,       "eta1/F");
  outtree_->Branch("phi1",       &phi1_,       "phi1/F");
  outtree_->Branch("px1",        &px1_,        "px1/F");
  outtree_->Branch("py1",        &py1_,        "py1/F");
  outtree_->Branch("pz1",        &pz1_,        "pz1/F");
  outtree_->Branch("y1",         &y1_,         "y1/F");
  outtree_->Branch("etacm1",     &etacm1_,     "etacm1/F");
  outtree_->Branch("ch2",        &ch2_,        "ch2/I");
  outtree_->Branch("pt2",        &pt2_,        "pt2/F");
  outtree_->Branch("eta2",       &eta2_,       "eta2/F");
  outtree_->Branch("phi2",       &phi2_,       "phi2/F");
  outtree_->Branch("px2",        &px2_,        "px2/F");
  outtree_->Branch("py2",        &py2_,        "py2/F");
  outtree_->Branch("pz2",        &pz2_,        "pz2/F");
  outtree_->Branch("y2",         &y2_,         "y2/F");
  outtree_->Branch("etacm2",     &etacm2_,     "etacm2/F");
  outtree_->Branch("nSelJets",   &nSelJets_,   "nSelJets/I");
  outtree_->Branch("ptJ",         ptJ_,        "ptJ[nSelJets]/F");
  outtree_->Branch("etaJ",        etaJ_,       "etaJ[nSelJets]/F");
  outtree_->Branch("phiJ",        phiJ_,       "phiJ[nSelJets]/F");
  outtree_->Branch("enJ",         enJ_,        "enJ[nSelJets]/F");
  outtree_->Branch("ptCombJ",    &ptCombJ_,    "ptCombJ/F");  
  outtree_->Branch("etaCombJ",   &etaCombJ_,   "etaCombJ/F");  
  outtree_->Branch("phiCombJ",   &phiCombJ_,   "phiCombJ/F");  
  outtree_->Branch("enCombJ",    &enCombJ_,    "enCombJ/F");  
  outtree_->Branch("puwei",      &puwei_,      "puwei/F");  
  outtree_->Branch("iso1",        iso1_,       "iso1[6]/F");
  outtree_->Branch("iso2",        iso2_,       "iso2[6]/F");  
  outtree_->Branch("elestatus",   elestatus_,  "elestatus[nEle]/F");
  outtree_->Branch("TTmz", &TTmZ_, "TTmZ/F");    
  outtree_->Branch("TTeta1", &TTeta1_, "TTeta1/F");  
  outtree_->Branch("TTphi1", &TTphi1_, "TTphi1/F");  
  outtree_->Branch("TTpt1", &TTpt1_, "TTpt1/F");
  outtree_->Branch("TTcharge1", &TTcharge1_, "TTcharge1/F");
  outtree_->Branch("TTeta2", &TTeta2_, "TTeta2/F");
  outtree_->Branch("TTphi2", &TTphi2_, "TTphi2/F");
  outtree_->Branch("TTpt2", &TTpt2_, "TTpt2/F");
  outtree_->Branch("TTcharge2", &TTcharge2_, "TTcharge2/F");
  outtree_->Branch("TPmz", &TPmZ_, "TPmZ/F");
  outtree_->Branch("TPeta1", &TPeta1_, "TPeta1/F");
  outtree_->Branch("TPphi1", &TPphi1_, "TPphi1/F");
  outtree_->Branch("TPpt1", &TPpt1_, "TPpt1/F");
  outtree_->Branch("TPcharge1", &TPcharge1_, "TPcharge1/F");
  outtree_->Branch("TPeta2", &TPeta2_, "TPeta2/F");
  outtree_->Branch("TPphi2", &TPphi2_, "TPphi2/F");
  outtree_->Branch("TPpt2", &TPpt2_, "TPpt2/F");
  outtree_->Branch("TPcharge2", &TPcharge2_, "TPcharge2/F");
  outtree_->Branch("TFmz", &TFmZ_, "TFmZ/F");
  outtree_->Branch("TFeta1", &TFeta1_, "TFeta1/F");
  outtree_->Branch("TFphi1", &TFphi1_, "TFphi1/F");
  outtree_->Branch("TFpt1", &TFpt1_, "TFpt1/F");
  outtree_->Branch("TFcharge1", &TFcharge1_, "TFcharge1/F");
  outtree_->Branch("TFeta2", &TFeta2_, "TFeta2/F");
  outtree_->Branch("TFphi2", &TFphi2_, "TFphi2/F");
  outtree_->Branch("TFpt2", &TFpt2_, "TFpt2/F");
  outtree_->Branch("TFcharge2", &TFcharge2_, "TFcharge2/F");
}

void xAna_tag_pro::InitTreePass() {

  PassingWP85Tree = new TTree("PassingWP85Tree", "Passing WP85 probes");

  /*t_ele1Pt    = 0 ;
  t_ele1Eta   = 0 ;
  t_ele1SCEta = 0 ;
  t_ele2Pt    = 0 ;
  t_ele2Eta   = 0 ;
  t_ele2SCEta = 0 ;
  EvtWeight   = 1 ;
  t_massZ     = 0 ;
  t_Ele1IsMCEle = false ;
  t_Ele2IsMCEle = false ;*/

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
  
}

void xAna_tag_pro::InitTreeFail() {

  FailingWP85Tree = new TTree("FailingWP85Tree", "Failing WP85 probes");

  /*t_ele1Pt    = 0 ;
  t_ele1Eta   = 0 ;
  t_ele1SCEta = 0 ;
  t_ele2Pt    = 0 ;
  t_ele2Eta   = 0 ;
  t_ele2SCEta = 0 ;
  EvtWeight   = 1 ;
  t_massZ     = 0 ;
  t_Ele1IsMCEle = false ;
  t_Ele2IsMCEle = false ;*/

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

}

#endif // #ifdef xAna_tag_pro_cxx

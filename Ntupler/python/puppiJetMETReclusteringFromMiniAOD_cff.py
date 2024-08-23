##################################################
# Adapted from https://github.com/nurfikri89/cmssw/blob/from14011_puppiJetsReclusteringFromMini/PhysicsTools/PatAlgos/python/tools/puppiJetMETReclusteringFromMiniAOD_cff.py
# Author: nurfikri89
##################################################

import FWCore.ParameterSet.Config as cms
from DeepNTuples.Ntupler.puppiJetMETReclusteringTools import puppiAK4METReclusterFromMiniAOD
from DeepNTuples.Ntupler.puppiJetMETReclusteringTools import puppiAK8ReclusterFromMiniAOD

def puppiJetMETReclusterFromMiniAOD(process, runOnMC, useExistingWeights=False, custom_dict=dict()):

  #
  # AK4 and MET
  #
  from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import _pfParticleNetFromMiniAODAK4PuppiCentralJetTagsAll as pfParticleNetFromMiniAODAK4PuppiCentralJetTagsAll
  from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import _pfParticleNetFromMiniAODAK4PuppiForwardJetTagsAll as pfParticleNetFromMiniAODAK4PuppiForwardJetTagsAll
  from RecoBTag.ONNXRuntime.pfUnifiedParticleTransformerAK4_cff import _pfUnifiedParticleTransformerAK4JetTagsAll as pfUnifiedParticleTransformerAK4JetTagsAll

  btagDiscriminatorsAK4 = cms.PSet(
   names=cms.vstring(
    'pfDeepFlavourJetTags:probb',
    'pfDeepFlavourJetTags:probbb',
    'pfDeepFlavourJetTags:problepb',
    'pfDeepFlavourJetTags:probc',
    'pfDeepFlavourJetTags:probuds',
    'pfDeepFlavourJetTags:probg')
    + pfParticleNetFromMiniAODAK4PuppiCentralJetTagsAll
    + pfParticleNetFromMiniAODAK4PuppiForwardJetTagsAll
    + pfUnifiedParticleTransformerAK4JetTagsAll
  )

  process = puppiAK4METReclusterFromMiniAOD(process, runOnMC,
    useExistingWeights=useExistingWeights,
    btagDiscriminatorsAK4=btagDiscriminatorsAK4
  )

  #
  # AK8
  #
  from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
  from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetMassRegressionOutputs as pfParticleNetMassRegressionOutputs
  from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetMassCorrelatedJetTagsAll as pfParticleNetMassCorrelatedJetTagsAll
  from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK8_cff import _pfParticleNetFromMiniAODAK8JetTagsAll as pfParticleNetFromMiniAODAK8JetTagsAll

  btagDiscriminatorsAK8 = pfParticleNetMassCorrelatedJetTagsAll+ \
      pfParticleNetFromMiniAODAK8JetTagsAll+ \
      pfParticleNetJetTagsAll+ \
      pfParticleNetMassRegressionOutputs+ \
      custom_dict.get('btagDiscriminatorsAK8', []) # add custom btag discriminators
  btagDiscriminatorsAK8 = cms.PSet(names = cms.vstring(btagDiscriminatorsAK8))

  btagDiscriminatorsAK8Subjets = cms.PSet(names = cms.vstring(
      'pfDeepCSVJetTags:probb',
      'pfDeepCSVJetTags:probbb',
    )
  )
  process = puppiAK8ReclusterFromMiniAOD(process, runOnMC,
    useExistingWeights=useExistingWeights,
    btagDiscriminatorsAK8=btagDiscriminatorsAK8,
    btagDiscriminatorsAK8Subjets=btagDiscriminatorsAK8Subjets,
    custom_dict=custom_dict,
  )

  return process

def puppiJetMETReclusterFromMiniAOD_MC(process):
  process = puppiJetMETReclusterFromMiniAOD(process, runOnMC=True)
  return process

def puppiJetMETReclusterFromMiniAOD_Data(process):
  process = puppiJetMETReclusterFromMiniAOD(process, runOnMC=False)
  return process
import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output.root'
# options.inputFiles = '/store/mc/RunIISummer20UL18MiniAODv2/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/001CA706-28B4-B24B-A2F7-63045C2BEB7F.root'
options.inputFiles = 'file:/eos/cms/store/cmst3/group/vhcc/sfTuples/GluGluToBulkGravitonToHHTo4W_JHUGen_M-1000_narrow/MiniAOD20UL/miniaod_20ul_35914-0.root'
options.maxEvents = -1

options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skip N events")
options.register('inputDataset',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input dataset")
options.register('isTrainSample', True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "if the sample is used for training")

options.parseArguments()

globalTagMap = {
    'auto': 'auto:phase1_2018_realistic',
    'UL18': '106X_upgrade2018_realistic_v16_L1v1',
    'UL17': '106X_mc2017_realistic_v9',
}

era = 'auto'
for k in globalTagMap:
    if k in options.inputDataset:
        era = k
if options.inputDataset == '':
    _inputfile = options.inputFiles[0]
    if '2017/' in _inputfile or 'UL17' in _inputfile:
        era = 'UL17'
    elif '2018/' in _inputfile or 'UL18' in _inputfile:
        era = 'UL18'
print('++era: ', era)
# ---------------------------------------------------------
process = cms.Process("DNNFiller")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(False)
)

print('Using output file ' + options.outputFile)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source('PoolSource',
                            fileNames=cms.untracked.vstring(options.inputFiles),
                            skipEvents=cms.untracked.uint32(options.skipEvents)
                            )
# ---------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTagMap[era], '')
print('Using global tag', process.GlobalTag.globaltag)
# ---------------------------------------------------------
# read JEC from sqlite
# if era == 'Summer19UL17':
#     import os
#     jecTag = 'Summer19UL17_V5_MC'
#     jecFile = '%s.db' % jecTag
#     if not os.path.exists(jecFile):
#         os.symlink('../data/' + jecFile, jecFile)
#     from CondCore.CondDB.CondDB_cfi import CondDB
#     CondDBJECFile = CondDB.clone(connect=cms.string('sqlite:%s' % jecFile))
#     process.jec = cms.ESSource('PoolDBESSource',
#                                CondDBJECFile,
#                                toGet=cms.VPSet(
#                                    cms.PSet(
#                                        record=cms.string('JetCorrectionsRecord'),
#                                        tag=cms.string('JetCorrectorParametersCollection_%s_AK4PFchs' % jecTag),
#                                        label=cms.untracked.string('AK4PFchs')
#                                    ),
#                                    cms.PSet(
#                                        record=cms.string('JetCorrectionsRecord'),
#                                        tag=cms.string('JetCorrectorParametersCollection_%s_AK4PFPuppi' % jecTag),
#                                        label=cms.untracked.string('AK4PFPuppi')
#                                    ),
#                                    # ...and so on for all jet types you need
#                                )
#                                )
#     print(jecTag, process.jec.toGet)
#     # Add an ESPrefer to override JEC that might be available from the global tag
#     process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')
# ---------------------------------------------------------
# Update to PuppiV14
# from CommonTools.PileupAlgos.customizePuppiTune_cff import UpdatePuppiTuneV14_MC
# UpdatePuppiTuneV14_MC(process)
# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from RecoBTag.ONNXRuntime.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsAll as pfDeepBoostedJetTagsAll
from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetMassRegressionOutputs as pfParticleNetMassRegressionOutputs

# !!! set `useReclusteredJets = True ` if you need to recluster jets (e.g., to adopt a new Puppi tune) !!!
useReclusteredJets = False
doCustomTaggerInference = True
jetR = 0.8

bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfBoostedDoubleSecondaryVertexAK8BJetTags',
    # 'pfDeepDoubleBvLJetTags:probHbb',
    # 'pfDeepDoubleCvLJetTags:probHcc',
    # 'pfDeepDoubleCvBJetTags:probHcc',
    # 'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
    # 'pfMassIndependentDeepDoubleCvLJetTags:probHcc',
    # 'pfMassIndependentDeepDoubleCvBJetTags:probHcc',
]

subjetBTagDiscriminators = ['None']

btagDiscriminatorsCustom = []
if doCustomTaggerInference:
    from DeepNTuples.Ntupler.jetTools import updateJetCollection # use custom updataJetCollection
    from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedDeepHWWV1_cff import _pfMassDecorrelatedDeepHWWV1JetTagsAll
    from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV1_cff import _pfMassDecorrelatedInclParticleTransformerV1JetTagsAll
    btagDiscriminatorsCustom = _pfMassDecorrelatedDeepHWWV1JetTagsAll + _pfMassDecorrelatedInclParticleTransformerV1JetTagsAll

if useReclusteredJets:
    JETCorrLevels = ['L2Relative', 'L3Absolute']

    from DeepNTuples.Ntupler.jetToolbox_cff import jetToolbox
    jetToolbox(process, 'ak8', 'dummySeq', 'noOutput', PUMethod='Puppi', JETCorrPayload='AK8PFPuppi',
               JETCorrLevels=JETCorrLevels, Cut='pt > 170.0 && abs(rapidity()) < 2.4', runOnMC=True, addNsub=True,
               maxTau=3, addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi',
               subJETCorrLevels=JETCorrLevels, bTagDiscriminators=['pfCombinedInclusiveSecondaryVertexV2BJetTags'],
               subjetBTagDiscriminators=subjetBTagDiscriminators)

    updateJetCollection(
        process,
        jetSource=cms.InputTag('packedPatJetsAK8PFPuppiSoftDrop'),
        rParam=jetR,
        jetCorrections=('AK8PFPuppi', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators=bTagDiscriminators + pfDeepBoostedJetTagsAll + pfParticleNetJetTagsAll + pfParticleNetMassRegressionOutputs + btagDiscriminatorsCustom,
        postfix='AK8WithPuppiDaughters',  # needed to tell the producers that the daughters are puppi-weighted
    )
    srcJets = cms.InputTag('selectedUpdatedPatJetsAK8WithPuppiDaughters')
else:
    updateJetCollection(
        process,
        jetSource=cms.InputTag('slimmedJetsAK8'),
        rParam=jetR,
        jetCorrections=('AK8PFPuppi', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators=pfParticleNetMassRegressionOutputs + btagDiscriminatorsCustom,
    )
    srcJets = cms.InputTag('selectedUpdatedPatJets')
# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask, addToProcessAndTask
patTask = getPatAlgosToolsTask(process)

from RecoJets.JetProducers.ak8GenJets_cfi import ak8GenJets
process.ak8GenJetsWithNu = ak8GenJets.clone(
    src='packedGenParticles',
    rParam=cms.double(jetR),
    jetPtMin=100.0
)
process.ak8GenJetsWithNuSoftDrop = process.ak8GenJetsWithNu.clone(
    useSoftDrop=cms.bool(True),
    zcut=cms.double(0.1),
    beta=cms.double(0.0),
    R0=cms.double(jetR),
    useExplicitGhosts=cms.bool(True)
)
process.ak8GenJetsWithNuMatch = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
                                               src=srcJets,  # RECO jets (any View<Jet> is ok)
                                               # GEN jets  (must be GenJetCollection)
                                               matched=cms.InputTag("ak8GenJetsWithNu"),
                                               mcPdgId=cms.vint32(),  # n/a
                                               mcStatus=cms.vint32(),  # n/a
                                               checkCharge=cms.bool(False),  # n/a
                                               maxDeltaR=cms.double(jetR),  # Minimum deltaR for the match
                                               # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
                                               # Forbid two RECO objects to match to the same GEN object
                                               resolveAmbiguities=cms.bool(True),
                                               # False = just match input in order; True = pick lowest deltaR pair first
                                               resolveByMatchQuality=cms.bool(False),
                                               )
process.ak8GenJetsWithNuSoftDropMatch = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
                                                       src=srcJets,  # RECO jets (any View<Jet> is ok)
                                                       # GEN jets  (must be GenJetCollection)
                                                       matched=cms.InputTag("ak8GenJetsWithNuSoftDrop"),
                                                       mcPdgId=cms.vint32(),  # n/a
                                                       mcStatus=cms.vint32(),  # n/a
                                                       checkCharge=cms.bool(False),  # n/a
                                                       maxDeltaR=cms.double(jetR),  # Minimum deltaR for the match
                                                       # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
                                                       # Forbid two RECO objects to match to the same GEN object
                                                       resolveAmbiguities=cms.bool(True),
                                                       # False = just match input in order; True = pick lowest deltaR pair first
                                                       resolveByMatchQuality=cms.bool(False),
                                                       )
process.genJetTask = cms.Task(
    process.ak8GenJetsWithNu,
    process.ak8GenJetsWithNuMatch,
    process.ak8GenJetsWithNuSoftDrop,
    process.ak8GenJetsWithNuSoftDropMatch,
)

# DeepNtuplizer
process.load("DeepNTuples.Ntupler.DeepNtuplizer_cfi")
process.deepntuplizer.jets = srcJets
process.deepntuplizer.useReclusteredJets = useReclusteredJets
process.deepntuplizer.bDiscriminators = bTagDiscriminators + pfDeepBoostedJetTagsAll + pfParticleNetJetTagsAll + pfParticleNetMassRegressionOutputs + btagDiscriminatorsCustom

process.deepntuplizer.genJetsMatch = 'ak8GenJetsWithNuMatch'
process.deepntuplizer.genJetsSoftDropMatch = 'ak8GenJetsWithNuSoftDropMatch'

# determine sample type with inputFiles name
_inputfile = options.inputFiles[0]
process.deepntuplizer.isQCDSample = '/QCD_' in _inputfile
process.deepntuplizer.isTTBarSample = 'tott' in _inputfile.lower() or 'ttbar' in _inputfile.lower()
process.deepntuplizer.isPythia = 'pythia' in _inputfile.lower()
process.deepntuplizer.isHerwig = 'herwig' in _inputfile.lower()
# note: MG can be interfaced w/ either pythia or herwig
process.deepntuplizer.isMadGraph = 'madgraph' in _inputfile.lower()

process.deepntuplizer.isTrainSample = options.isTrainSample
# adding MET
process.deepntuplizer.addMET = True
if not options.inputDataset:
    # interactive running
    process.deepntuplizer.isTrainSample = False
#==============================================================================================================================#
process.p = cms.Path(process.deepntuplizer)
process.p.associate(patTask)
process.p.associate(process.genJetTask)

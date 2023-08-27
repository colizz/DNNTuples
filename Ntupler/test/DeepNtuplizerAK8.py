import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output.root'
# options.inputFiles = '/store/mc/RunIISummer20UL17MiniAODv2/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/240000/82A9346F-78D8-AE45-A0A9-F3BC0A33480B.root' ## QCD test
# options.inputFiles = '/store/mc/RunIISummer20UL17MiniAODv2/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/230000/99B7EC88-541E-3143-8892-6DDD179CBC4B.root' ## QCD test 2
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/BulkGravitonToHHTo4QTau_MX-600to6000_MH-15to250/20UL17MiniAODv2/part2/miniv2_14998562-1403.root' ## H->qq/tautau
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/BulkGravitonToHHTo4Glu_MX-600to6000_MH-15to250/20UL17MiniAODv2/part3/miniv2_15687298-876.root' ## H->gg
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/ChargedHiggs_HplusToBC_HminusToBC_HT-600to6000_MH-15to250/20UL17MiniAODv2/miniv2_15164446-2044.root' ## H->bc
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/ChargedHiggs_HplusToCS_HminusToCS_HT-600to6000_MH-15to250/20UL17MiniAODv2/miniv2_15687299-2655.root' ## H->cs
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/H1H1_2HDM_H1ToBS_HT-600to6000_MH-15to250/20UL17MiniAODv2/miniv2_15723624-285.root' ## H->bs
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/BulkGravitonToHHTo4L_MX-600to6000_MH-15to250/20UL17MiniAODv2/part3/miniv2_15744104-989.root' ## H->ll
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/Spin0ToTT_VariableMass_WhadWhad_MX-600to6000_MH-15to250/20UL17MiniAODv2/miniv2_15094466-3566.root' ## top hadronic
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/Spin0ToTT_VariableMass_WlepWlep_MX-600to6000_MH-15to250_ext1/20UL17MiniAODv2/miniv2_15094468-7131.root' ## top leptonic
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/20220523_HWW_JHUVariableWMass/2017/mc/BulkGravitonToHHTo4W_MX-600to6000_MH-15to250_JHUVariableWMass_part2/DNNTuples_PrivateMC/220523_095639/0002/miniv2_2419.root' ## HWW 1
# options.inputFiles = '/store/group/cmst3/group/vhcc/sfTuples/20220523_HWW_JHUVariableWMass/2017/mc/BulkGravitonToHHTo4W_MX-600to6000_MH-15to250_JHUVariableWMass_part2/DNNTuples_PrivateMC/220523_095639/0002/miniv2_2420.root' ## HWW 2
# options.inputFiles = '/store/cmst3/group/vhcc/sfTuples/BulkGravitonToHHTo4W_MX-600to6000_MH-15to250_JHUVariableWMass2DMesh/20UL17MiniAODv2/part2/miniv2_15683018-4322.root' ## HWW 2dmesh
# options.inputFiles = '/store/cmst3/group/vhcc/sfTuples/BulkGravitonToHHTo4Z_MX-600to6000_MH-15to250_JHUVariableZMass/20UL17MiniAODv2/part3/miniv2_15820048-548.root' ## HZZ 1
# options.inputFiles = '/store/cmst3/group/vhcc/sfTuples/BulkGravitonToHHTo4Z_MX-600to6000_MH-15to250_JHUVariableZMass/20UL17MiniAODv2/part3/miniv2_15820048-549.root' ## HZZ 2
# options.inputFiles = '/store/cmst3/group/vhcc/sfTuples/BulkGravitonToHHTo4Z_MX-600to6000_MH-15to250_JHUVariableZMass2DMesh/20UL17MiniAODv2/part2/miniv2_15822201-2585.root' ## HZZ 2dmesh
# options.inputFiles = '/store/cmst3/group/vhcc/sfTuples/PairVectorLQ_LQToBTau_HT-600to6000_M-15to250/20UL17MiniAODv2/miniv2_16076647-1.root' ## customized btau
# options.inputFiles = '/store/mc/RunIISummer20UL17MiniAODv2/BulkGravToZZToZhadZhad_narrow_M-1000_TuneCP5_13TeV-madgraph-pythia/MINIAODSIM/106X_mc2017_realistic_v9-v2/110000/DABA0ABE-8F97-9747-9A7A-4E31435442E1.root' ## Zqq inference
# options.inputFiles = '/store/mc/RunIISummer20UL17MiniAODv2/BulkGravToWWToWhadWhad_narrow_M-1000_TuneCP5_13TeV-madgraph-pythia/MINIAODSIM/106X_mc2017_realistic_v9-v2/260000/F1F668E3-CB4A-ED4E-9493-3485628D5059.root'  ## Wqq inference
options.inputFiles = '/store/mc/RunIISummer20UL17MiniAODv2/GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/50000/352A1D3C-77B8-1A42-A964-87E486EDBB7E.root' # AD ggHbb sample

options.maxEvents = -1

options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skip N events")
options.register('inputDataset',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input dataset")
options.register('isTrainSample', True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "if the sample is used for training")
# special output configs
options.register('addMET', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "add MET vars to output file")
options.register('addLowLevel', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "add low-level vars to output file")

options.parseArguments()

# test command: cmsRun DeepNtuplizerAK8.py maxEvents=100 isTrainSample=1

globalTagMap = {
    'auto': 'auto:phase1_2018_realistic',
    'UL18': '106X_upgrade2018_realistic_v16_L1v1',
    'UL17': '106X_mc2017_realistic_v9',
}

## current workflow: only run interactively
assert options.inputDataset == '' and len(options.inputFiles) == 1, 'Only run interactively with file len=1'
print("Running", options.inputFiles)
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
jetR = 0.8

# only inference the new taggers when not producing the training sample
doCustomTaggerInference = True if not options.isTrainSample else False

bTagDiscriminators = [
]

subjetBTagDiscriminators = ['None']

# from dnntuple v9: infer the new tagger so as to store the hidden layer scores in a special branch jet_custom_discs
from DeepNTuples.Ntupler.jetTools import updateJetCollection # use custom updataJetCollection
from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV2_cff import _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbsHidNeurons
from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV2_cff import _pfMassDecorrelatedInclParticleTransformerV2M125HidLayerJetTagsProbsHidNeurons
btagDiscriminatorsCustomSaveAsCompact1 = _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbsHidNeurons
btagDiscriminatorsCustomSaveAsCompact2 = _pfMassDecorrelatedInclParticleTransformerV2M125HidLayerJetTagsProbsHidNeurons
btagDiscriminatorsCustomSaveAsSeparate = []

if useReclusteredJets:
    raise NotImplementedError('Reclustering not implemented')

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
        btagDiscriminators=bTagDiscriminators + pfDeepBoostedJetTagsAll + pfParticleNetJetTagsAll + pfParticleNetMassRegressionOutputs + btagDiscriminatorsCustomSaveAsCompact + btagDiscriminatorsCustomSaveAsSeparate,
        postfix='AK8WithPuppiDaughters',  # needed to tell the producers that the daughters are puppi-weighted
    )
    srcJets = cms.InputTag('selectedUpdatedPatJetsAK8WithPuppiDaughters')
else:
    updateJetCollection(
        process,
        jetSource=cms.InputTag('slimmedJetsAK8'),
        rParam=jetR,
        jetCorrections=('AK8PFPuppi', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators=pfParticleNetMassRegressionOutputs + btagDiscriminatorsCustomSaveAsCompact1 + btagDiscriminatorsCustomSaveAsCompact2 + btagDiscriminatorsCustomSaveAsSeparate,
    )
    srcJets = cms.InputTag('selectedUpdatedPatJets')
# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask, addToProcessAndTask
patTask = getPatAlgosToolsTask(process)

from RecoJets.JetProducers.ak8GenJets_cfi import ak8GenJets
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
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
process.packedGenParticlesForJetsNoNu = genParticlesForJetsNoNu.clone(
    src='packedGenParticles'
)
process.ak8GenJetsNoNu = process.ak8GenJetsWithNu.clone(
    src='packedGenParticlesForJetsNoNu' # using packedGenParticles with no neutrinos as input
)
process.ak8GenJetsNoNuSoftDrop = process.ak8GenJetsWithNuSoftDrop.clone(
    src='packedGenParticlesForJetsNoNu'
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
process.ak8GenJetsNoNuMatch = process.ak8GenJetsWithNuMatch.clone(matched=cms.InputTag("ak8GenJetsNoNu"))
process.ak8GenJetsNoNuSoftDropMatch = process.ak8GenJetsWithNuSoftDropMatch.clone(matched=cms.InputTag("ak8GenJetsNoNuSoftDrop"))

process.genJetTask = cms.Task(
    process.ak8GenJetsWithNu,
    process.ak8GenJetsWithNuMatch,
    process.ak8GenJetsWithNuSoftDrop,
    process.ak8GenJetsWithNuSoftDropMatch,
    process.packedGenParticlesForJetsNoNu,
    process.ak8GenJetsNoNu,
    process.ak8GenJetsNoNuMatch,
    process.ak8GenJetsNoNuSoftDrop,
    process.ak8GenJetsNoNuSoftDropMatch,
)

# DeepNtuplizer
process.load("DeepNTuples.Ntupler.DeepNtuplizer_cfi")
process.deepntuplizer.jets = srcJets
process.deepntuplizer.useReclusteredJets = useReclusteredJets
process.deepntuplizer.bDiscriminators = bTagDiscriminators + pfDeepBoostedJetTagsAll + pfParticleNetJetTagsAll + pfParticleNetMassRegressionOutputs + btagDiscriminatorsCustomSaveAsSeparate
process.deepntuplizer.bDiscriminatorsCompactSave1 = btagDiscriminatorsCustomSaveAsCompact1
process.deepntuplizer.bDiscriminatorsCompactSave2 = btagDiscriminatorsCustomSaveAsCompact2

process.deepntuplizer.genJetsWithNuMatch = 'ak8GenJetsWithNuMatch'
process.deepntuplizer.genJetsWithNuSoftDropMatch = 'ak8GenJetsWithNuSoftDropMatch'
process.deepntuplizer.genJetsNoNuMatch = 'ak8GenJetsNoNuMatch'
process.deepntuplizer.genJetsNoNuSoftDropMatch = 'ak8GenJetsNoNuSoftDropMatch'

# determine sample type with inputFiles name
_inputfile = options.inputFiles[0]
process.deepntuplizer.isQCDSample = '/QCD_' in _inputfile
process.deepntuplizer.isTTBarSample = 'tott' in _inputfile.lower() or 'ttbar' in _inputfile.lower()
process.deepntuplizer.isHVV2DVarMassSample = '2DMesh' in _inputfile
process.deepntuplizer.isPythia = 'pythia' in _inputfile.lower()
process.deepntuplizer.isHerwig = 'herwig' in _inputfile.lower()
# note: MG can be interfaced w/ either pythia or herwig
process.deepntuplizer.isMadGraph = 'madgraph' in _inputfile.lower()

process.deepntuplizer.isTrainSample = options.isTrainSample

# special output configs
process.deepntuplizer.addMET = options.addMET
process.deepntuplizer.addLowLevel = options.addLowLevel

# AD configs
process.deepntuplizer.addHLT = True
process.deepntuplizer.HLTList = cms.untracked.vstring(     
    "HLT_AK8PFHT800_TrimMass50_v",
    "HLT_AK8PFJet400_TrimMass30_v",
    "HLT_AK8PFJet420_TrimMass30_v",
    "HLT_AK8PFJet500_v",
    "HLT_PFHT1050_v",
    "HLT_PFJet500_v",
)

#==============================================================================================================================#
process.p = cms.Path(process.deepntuplizer)
process.p.associate(patTask)
process.p.associate(process.genJetTask)
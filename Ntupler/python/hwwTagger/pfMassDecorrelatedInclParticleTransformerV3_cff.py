import FWCore.ParameterSet.Config as cms

# use ParticleTransformerV2JetTagInfosProducer (include lost track and pixel inputs)
from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer

## GloParT 3 models: ##
# downloaded from https://coli.web.cern.ch/coli/tmp/.230626-003937_partv2_model
# V03pre1: "ak8_MD_inclv8_part_splitreg_extmass_rmhbs_manual.useamp.large_fc2048.gm5.ddp-bs384-lr6e-4.nepoch100.v3fixsplit2", no H->bs class, split-class regression, no pixel inputs
# V03beta1: "ak8_MD_inclv8std_rmhbs_manual.useamp.large_fc2048.pemb64_block10.gm5.ddp-bs512-lr7e-4.nepoch100", no H->bs class, entire model (313+313+28 nodes), selected output (24) + hidden neurons

pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTags = boostedJetONNXJetTagsProducer.clone(
    src = 'pfMassDecorrelatedInclParticleTransformerV2TagInfos', # v3 still uses v2 tag infos
    preprocess_json = 'DeepNTuples/Ntupler/data/InclParticleTransformer-MD/ak8/V03beta1/preprocess.json',
    model_path = 'DeepNTuples/Ntupler/data/InclParticleTransformer-MD/ak8/V03beta1/model.onnx',
    flav_names = [
        'probXbb', 'probXcc', 'probXcs', 'probXqq', 'probXtauhtaue', 'probXtauhtaum', 'probXtauhtauh', 'probXWW4q', 'probXWW3q', 'probXWWqqev', 'probXWWqqmv', 'probTopbWqq', 'probTopbWq', 'probTopbWev', 'probTopbWmv', 'probTopbWtauhv', 'probQCD', 'massCorrX2p', 'massCorrXWW', 'massCorrTopbW', 'massCorrQCD', 'probWithMassTopvsQCD', 'probWithMassWvsQCD', 'probWithMassZvsQCD'
    ] + ['hidNeuron' + str(i).zfill(3) for i in range(256)],
    debugMode = False,
)

# declare all the discriminators
# probs
_pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTagsProbs = ['pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTags:' + flav_name
                                 for flav_name in pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTags.flav_names]
_pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTagsProbsRawScores = [disc for disc in _pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTagsProbs if 'hidNeuron' not in disc]
_pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTagsProbsHidNeurons = [disc for disc in _pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTagsProbs if 'hidNeuron' in disc]

# meta-taggers
_pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTagsMetaDiscrs = []
_pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTagsAll = _pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTagsProbs + _pfMassDecorrelatedInclParticleTransformerV3HidLayerJetTagsMetaDiscrs

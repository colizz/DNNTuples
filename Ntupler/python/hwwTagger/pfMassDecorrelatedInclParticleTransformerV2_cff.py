import FWCore.ParameterSet.Config as cms

# use ParticleTransformerV2JetTagInfosProducer (include lost track and pixel inputs)
from DeepNTuples.Ntupler.pfParticleTransformerAK8TagInfos_cfi import pfParticleTransformerAK8TagInfos as pfParticleTransformerV2JetTagInfos
from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer
from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV2DiscriminatorsJetTags_cfi import pfMassDecorrelatedInclParticleTransformerV2DiscriminatorsJetTags

pfMassDecorrelatedInclParticleTransformerV2TagInfos = pfParticleTransformerV2JetTagInfos.clone(
    use_puppiP4 = False
)

pfMassDecorrelatedInclParticleTransformerV2JetTags = boostedJetONNXJetTagsProducer.clone(
    src = 'pfMassDecorrelatedInclParticleTransformerV2TagInfos',
    preprocess_json = 'DeepNTuples/Ntupler/data/InclParticleTransformer-MD/ak8/V02/preprocess_corr.json',
    model_path = 'DeepNTuples/Ntupler/data/InclParticleTransformer-MD/ak8/V02/model.onnx',
    flav_names = [
        'probTopbWcs', 'probTopbWqq', 'probTopbWc', 'probTopbWs', 'probTopbWq', 'probTopbWev', 'probTopbWmv', 'probTopbWtauev', 'probTopbWtaumv', 'probTopbWtauhv', 
        'probTopWcs', 'probTopWqq', 'probTopWev', 'probTopWmv', 'probTopWtauev', 'probTopWtaumv', 'probTopWtauhv', 
        'probHbb', 'probHcc', 'probHss', 'probHqq', 'probHbc', 'probHbs', 'probHcs', 'probHgg', 'probHee', 'probHmm', 'probHtauhtaue', 'probHtauhtaum', 'probHtauhtauh', 
        'probHWWcscs', 'probHWWcsqq', 'probHWWqqqq', 'probHWWcsc', 'probHWWcss', 'probHWWcsq', 'probHWWqqc', 'probHWWqqs', 'probHWWqqq', 
        'probHWWcsev', 'probHWWqqev', 'probHWWcsmv', 'probHWWqqmv', 'probHWWcstauev', 'probHWWqqtauev', 'probHWWcstaumv', 'probHWWqqtaumv', 'probHWWcstauhv', 'probHWWqqtauhv', 
        'probHWxWxcscs', 'probHWxWxcsqq', 'probHWxWxqqqq', 'probHWxWxcsc', 'probHWxWxcss', 'probHWxWxcsq', 'probHWxWxqqc', 'probHWxWxqqs', 'probHWxWxqqq', 
        'probHWxWxcsev', 'probHWxWxqqev', 'probHWxWxcsmv', 'probHWxWxqqmv', 'probHWxWxcstauev', 'probHWxWxqqtauev', 'probHWxWxcstaumv', 'probHWxWxqqtaumv', 'probHWxWxcstauhv', 'probHWxWxqqtauhv', 
        'probHWxWxStarcscs', 'probHWxWxStarcsqq', 'probHWxWxStarqqqq', 'probHWxWxStarcsc', 'probHWxWxStarcss', 'probHWxWxStarcsq', 'probHWxWxStarqqc', 'probHWxWxStarqqs', 'probHWxWxStarqqq', 
        'probHWxWxStarcsev', 'probHWxWxStarqqev', 'probHWxWxStarcsmv', 'probHWxWxStarqqmv', 'probHWxWxStarcstauev', 'probHWxWxStarqqtauev', 'probHWxWxStarcstaumv', 'probHWxWxStarqqtaumv', 'probHWxWxStarcstauhv', 'probHWxWxStarqqtauhv', 
        'probHZZbbbb', 'probHZZbbcc', 'probHZZbbss', 'probHZZbbqq', 'probHZZcccc', 'probHZZccss', 'probHZZccqq', 'probHZZssss', 'probHZZssqq', 'probHZZqqqq', 'probHZZbbb', 'probHZZbbc', 'probHZZbbs', 'probHZZbbq', 'probHZZccb', 'probHZZccc', 'probHZZccs', 'probHZZccq', 'probHZZssb', 'probHZZssc', 'probHZZsss', 'probHZZssq', 'probHZZqqb', 'probHZZqqc', 'probHZZqqs', 'probHZZqqq', 
        'probHZZbbee', 'probHZZbbmm', 'probHZZbbe', 'probHZZbbm', 'probHZZbee', 'probHZZbmm', 'probHZZbbtauhtaue', 'probHZZbbtauhtaum', 'probHZZbbtauhtauh', 'probHZZbtauhtaue', 'probHZZbtauhtaum', 'probHZZbtauhtauh', 'probHZZccee', 'probHZZccmm', 'probHZZcce', 'probHZZccm', 'probHZZcee', 'probHZZcmm', 'probHZZcctauhtaue', 'probHZZcctauhtaum', 'probHZZcctauhtauh', 'probHZZctauhtaue', 'probHZZctauhtaum', 'probHZZctauhtauh', 'probHZZssee', 'probHZZssmm', 'probHZZsse', 'probHZZssm', 'probHZZsee', 'probHZZsmm', 'probHZZsstauhtaue', 'probHZZsstauhtaum', 'probHZZsstauhtauh', 'probHZZstauhtaue', 'probHZZstauhtaum', 'probHZZstauhtauh', 'probHZZqqee', 'probHZZqqmm', 'probHZZqqe', 'probHZZqqm', 'probHZZqee', 'probHZZqmm', 'probHZZqqtauhtaue', 'probHZZqqtauhtaum', 'probHZZqqtauhtauh', 'probHZZqtauhtaue', 'probHZZqtauhtaum', 'probHZZqtauhtauh', 
        'probHZxZxbbbb', 'probHZxZxbbcc', 'probHZxZxbbss', 'probHZxZxbbqq', 'probHZxZxcccc', 'probHZxZxccss', 'probHZxZxccqq', 'probHZxZxssss', 'probHZxZxssqq', 'probHZxZxqqqq', 'probHZxZxbbb', 'probHZxZxbbc', 'probHZxZxbbs', 'probHZxZxbbq', 'probHZxZxccb', 'probHZxZxccc', 'probHZxZxccs', 'probHZxZxccq', 'probHZxZxssb', 'probHZxZxssc', 'probHZxZxsss', 'probHZxZxssq', 'probHZxZxqqb', 'probHZxZxqqc', 'probHZxZxqqs', 'probHZxZxqqq', 
        'probHZxZxbbee', 'probHZxZxbbmm', 'probHZxZxbbe', 'probHZxZxbbm', 'probHZxZxbee', 'probHZxZxbmm', 'probHZxZxbbtauhtaue', 'probHZxZxbbtauhtaum', 'probHZxZxbbtauhtauh', 'probHZxZxbtauhtaue', 'probHZxZxbtauhtaum', 'probHZxZxbtauhtauh', 'probHZxZxccee', 'probHZxZxccmm', 'probHZxZxcce', 'probHZxZxccm', 'probHZxZxcee', 'probHZxZxcmm', 'probHZxZxcctauhtaue', 'probHZxZxcctauhtaum', 'probHZxZxcctauhtauh', 'probHZxZxctauhtaue', 'probHZxZxctauhtaum', 'probHZxZxctauhtauh', 'probHZxZxssee', 'probHZxZxssmm', 'probHZxZxsse', 'probHZxZxssm', 'probHZxZxsee', 'probHZxZxsmm', 'probHZxZxsstauhtaue', 'probHZxZxsstauhtaum', 'probHZxZxsstauhtauh', 'probHZxZxstauhtaue', 'probHZxZxstauhtaum', 'probHZxZxstauhtauh', 'probHZxZxqqee', 'probHZxZxqqmm', 'probHZxZxqqe', 'probHZxZxqqm', 'probHZxZxqee', 'probHZxZxqmm', 'probHZxZxqqtauhtaue', 'probHZxZxqqtauhtaum', 'probHZxZxqqtauhtauh', 'probHZxZxqtauhtaue', 'probHZxZxqtauhtaum', 'probHZxZxqtauhtauh', 
        'probHZxZxStarbbbb', 'probHZxZxStarbbcc', 'probHZxZxStarbbss', 'probHZxZxStarbbqq', 'probHZxZxStarcccc', 'probHZxZxStarccss', 'probHZxZxStarccqq', 'probHZxZxStarssss', 'probHZxZxStarssqq', 'probHZxZxStarqqqq', 'probHZxZxStarbbb', 'probHZxZxStarbbc', 'probHZxZxStarbbs', 'probHZxZxStarbbq', 'probHZxZxStarccb', 'probHZxZxStarccc', 'probHZxZxStarccs', 'probHZxZxStarccq', 'probHZxZxStarssb', 'probHZxZxStarssc', 'probHZxZxStarsss', 'probHZxZxStarssq', 'probHZxZxStarqqb', 'probHZxZxStarqqc', 'probHZxZxStarqqs', 'probHZxZxStarqqq', 
        'probHZxZxStarbbee', 'probHZxZxStarbbmm', 'probHZxZxStarbbe', 'probHZxZxStarbbm', 'probHZxZxStarbee', 'probHZxZxStarbmm', 'probHZxZxStarbbtauhtaue', 'probHZxZxStarbbtauhtaum', 'probHZxZxStarbbtauhtauh', 'probHZxZxStarbtauhtaue', 'probHZxZxStarbtauhtaum', 'probHZxZxStarbtauhtauh', 'probHZxZxStarccee', 'probHZxZxStarccmm', 'probHZxZxStarcce', 'probHZxZxStarccm', 'probHZxZxStarcee', 'probHZxZxStarcmm', 'probHZxZxStarcctauhtaue', 'probHZxZxStarcctauhtaum', 'probHZxZxStarcctauhtauh', 'probHZxZxStarctauhtaue', 'probHZxZxStarctauhtaum', 'probHZxZxStarctauhtauh', 'probHZxZxStarssee', 'probHZxZxStarssmm', 'probHZxZxStarsse', 'probHZxZxStarssm', 'probHZxZxStarsee', 'probHZxZxStarsmm', 'probHZxZxStarsstauhtaue', 'probHZxZxStarsstauhtaum', 'probHZxZxStarsstauhtauh', 'probHZxZxStarstauhtaue', 'probHZxZxStarstauhtaum', 'probHZxZxStarstauhtauh', 'probHZxZxStarqqee', 'probHZxZxStarqqmm', 'probHZxZxStarqqe', 'probHZxZxStarqqm', 'probHZxZxStarqee', 'probHZxZxStarqmm', 'probHZxZxStarqqtauhtaue', 'probHZxZxStarqqtauhtaum', 'probHZxZxStarqqtauhtauh', 'probHZxZxStarqtauhtaue', 'probHZxZxStarqtauhtaum', 'probHZxZxStarqtauhtauh', 
        'probQCDbb', 'probQCDcc', 'probQCDb', 'probQCDc', 'probQCDothers', 
        'resonanceMassCorr', 'visiableMassCorr',
    ],
    debugMode = False,
)

# declare all the discriminators
# probs
_pfMassDecorrelatedInclParticleTransformerV2JetTagsProbs = ['pfMassDecorrelatedInclParticleTransformerV2JetTags:' + flav_name
                                 for flav_name in pfMassDecorrelatedInclParticleTransformerV2JetTags.flav_names]
# meta-taggers
_pfMassDecorrelatedInclParticleTransformerV2JetTagsMetaDiscrs = ['pfMassDecorrelatedInclParticleTransformerV2DiscriminatorsJetTags:' + disc.name.value()
                                      for disc in pfMassDecorrelatedInclParticleTransformerV2DiscriminatorsJetTags.discriminators]

_pfMassDecorrelatedInclParticleTransformerV2JetTagsAll = _pfMassDecorrelatedInclParticleTransformerV2JetTagsProbs + _pfMassDecorrelatedInclParticleTransformerV2JetTagsMetaDiscrs

# ===========================================================================================
# For AK8 jets with hidden layer neurons

pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTags = boostedJetONNXJetTagsProducer.clone(
    src = 'pfMassDecorrelatedInclParticleTransformerV2TagInfos',
    preprocess_json = 'DeepNTuples/Ntupler/data/InclParticleTransformer-MD/ak8/V02-HidLayer/preprocess_corr.json',
    model_path = 'DeepNTuples/Ntupler/data/InclParticleTransformer-MD/ak8/V02-HidLayer/model_embed.onnx',
    flav_names = pfMassDecorrelatedInclParticleTransformerV2JetTags.flav_names + ['hidNeuron%s' % str(idx).zfill(3) for idx in range(256)],
    debugMode = False,
)

# declare all the discriminators
# probs
_pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbs = ['pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTags:' + flav_name
                                 for flav_name in pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTags.flav_names]
_pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbsRawScores = [disc for disc in _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbs if 'hidNeuron' not in disc]
_pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbsHidNeurons = [disc for disc in _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbs if 'hidNeuron' in disc]

# meta-taggers
_pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsMetaDiscrs = ['pfMassDecorrelatedInclParticleTransformerV2HidLayerDiscriminatorsJetTags:' + disc.name.value()
                                      for disc in pfMassDecorrelatedInclParticleTransformerV2DiscriminatorsJetTags.discriminators]

_pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsAll = _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbs + _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsMetaDiscrs

# ===========================================================================================
# For AK8 jets with hidden layer neurons, rescaled to targetd mass

mass_edges = [0 + k * 10 for k in range(100)]

## special tag info producer
pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdxM1TagInfos = pfParticleTransformerV2JetTagInfos.clone(
    use_puppiP4 = False,
    mass_edges = mass_edges,
    mass_edge_index = -1,
)
pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdx0TagInfos = pfParticleTransformerV2JetTagInfos.clone(
    use_puppiP4 = False,
    mass_edges = mass_edges,
    mass_edge_index = 0,
)
pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdx1TagInfos = pfParticleTransformerV2JetTagInfos.clone(
    use_puppiP4 = False,
    mass_edges = mass_edges,
    mass_edge_index = 1,
)
pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdx2TagInfos = pfParticleTransformerV2JetTagInfos.clone(
    use_puppiP4 = False,
    mass_edges = mass_edges,
    mass_edge_index = 2,
)

pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdxM1HidLayerJetTags = pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTags.clone(
    src = 'pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdxM1TagInfos',
)
pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdx0HidLayerJetTags = pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTags.clone(
    src = 'pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdx0TagInfos',
)
pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdx1HidLayerJetTags = pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTags.clone(
    src = 'pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdx1TagInfos',
)
pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdx2HidLayerJetTags = pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTags.clone(
    src = 'pfMassDecorrelatedInclParticleTransformerV2JetP4ScalingIdx2TagInfos',
)

# declare all the discriminators
# probs
def getJetTagsProbs(name, requireHidNeurons=False):
    jetTagsProbs = [name + 'JetTags:' + flav_name for flav_name in pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTags.flav_names]
    if requireHidNeurons:
        jetTagsProbs = [disc for disc in jetTagsProbs if 'hidNeuron' in disc]
    return jetTagsProbs

# meta-taggers
def getJetTagsMetaDiscrs(name):
    return [name + 'DiscriminatorsJetTags:' + disc.name.value() for disc in pfMassDecorrelatedInclParticleTransformerV2DiscriminatorsJetTags.discriminators]

# ===========================================================================================
# For AK15 jets:
# obtain tag infos with radius 1.5, then inference with the model trained with radius 1.5 jets

pfMassDecorrelatedInclParticleTransformerAK15V2TagInfos = pfParticleTransformerV2JetTagInfos.clone(
    use_puppiP4 = False,
    jet_radius = 1.5,
)

pfMassDecorrelatedInclParticleTransformerAK15V2JetTags = pfMassDecorrelatedInclParticleTransformerV2JetTags.clone(
    src = 'pfMassDecorrelatedInclParticleTransformerAK15V2TagInfos',
    preprocess_json = 'DeepNTuples/Ntupler/data/InclParticleTransformer-MD/ak15/V02/preprocess_corr.json',
    model_path = 'DeepNTuples/Ntupler/data/InclParticleTransformer-MD/ak15/V02/model.onnx',
)

# declare all the discriminators
# probs
_pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsProbs = ['pfMassDecorrelatedInclParticleTransformerAK15V2JetTags:' + flav_name
                                 for flav_name in pfMassDecorrelatedInclParticleTransformerAK15V2JetTags.flav_names]

_pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsProbsSelected = ['pfMassDecorrelatedInclParticleTransformerAK15V2JetTags:' + flav_name
                                 for flav_name in ['probHbb', 'probHcc', 'probHss', 'probHqq', 'probHbc', 'probHbs', 'probHcs', 'probHgg', 'probHee', 'probHmm', 'probHtauhtaue', 'probHtauhtaum', 'probHtauhtauh', 
                                                   'probQCDbb', 'probQCDcc', 'probQCDb', 'probQCDc', 'probQCDothers', 
                                                   'resonanceMassCorr', 'visiableMassCorr']]
# meta-taggers
_pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsMetaDiscrs = ['pfMassDecorrelatedInclParticleTransformerAK15V2DiscriminatorsJetTags:' + disc.name.value()
                                      for disc in pfMassDecorrelatedInclParticleTransformerV2DiscriminatorsJetTags.discriminators]

_pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsAll = _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsProbs + _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsMetaDiscrs
_pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsAllSelected = _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsProbsSelected + _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsMetaDiscrs

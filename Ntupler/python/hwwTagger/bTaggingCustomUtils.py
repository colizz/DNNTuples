def updateSupportedBtagDiscr(supportedBtagInfos, supportedBtagDiscr, supportedMetaDiscr):
    
    ## Update taggers in DeepHWWV1, InclParticleTransformerV1, InclParticleTransformerV2
    from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedDeepHWWV1_cff import _pfMassDecorrelatedDeepHWWV1JetTagsProbs, _pfMassDecorrelatedDeepHWWV1JetTagsMetaDiscrs
    from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV1_cff import _pfMassDecorrelatedInclParticleTransformerV1JetTagsProbs, _pfMassDecorrelatedInclParticleTransformerV1JetTagsMetaDiscrs
    from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV2_cff import _pfMassDecorrelatedInclParticleTransformerV2JetTagsProbs, _pfMassDecorrelatedInclParticleTransformerV2JetTagsMetaDiscrs
    from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV2_cff import _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsProbs, _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsMetaDiscrs # AK15 tagger
    from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV2_cff import _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbs, _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsMetaDiscrs
    from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV2_cff import getJetTagsProbs, getJetTagsMetaDiscrs
    
    # update supportedBtagDiscr
    supportedBtagInfos.extend(["pfMassDecorrelatedDeepHWWV1TagInfos"])
    supportedBtagInfos.extend(["pfMassDecorrelatedInclParticleTransformerV1TagInfos"])
    supportedBtagInfos.extend(["pfMassDecorrelatedInclParticleTransformerV2TagInfos"])
    supportedBtagInfos.extend(["pfMassDecorrelatedInclParticleTransformerV2RemoveLeadLeptonTagInfos"])
    supportedBtagInfos.extend(["pfMassDecorrelatedInclParticleTransformerAK15V2TagInfos"])
    for disc in _pfMassDecorrelatedDeepHWWV1JetTagsProbs + _pfMassDecorrelatedDeepHWWV1JetTagsMetaDiscrs:
        supportedBtagDiscr[disc] = [["pfMassDecorrelatedDeepHWWV1TagInfos"]]
    for disc in _pfMassDecorrelatedInclParticleTransformerV1JetTagsProbs + _pfMassDecorrelatedInclParticleTransformerV1JetTagsMetaDiscrs:
        supportedBtagDiscr[disc] = [["pfMassDecorrelatedInclParticleTransformerV1TagInfos"]]
    for disc in _pfMassDecorrelatedInclParticleTransformerV2JetTagsProbs + _pfMassDecorrelatedInclParticleTransformerV2JetTagsMetaDiscrs:
        supportedBtagDiscr[disc] = [["pfMassDecorrelatedInclParticleTransformerV2TagInfos"]]
    for disc in _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsProbs + _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsMetaDiscrs:
        supportedBtagDiscr[disc] = [["pfMassDecorrelatedInclParticleTransformerAK15V2TagInfos"]]
    for disc in _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbs + _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsMetaDiscrs:
        supportedBtagDiscr[disc] = [["pfMassDecorrelatedInclParticleTransformerV2TagInfos"]]
    for n in ['RemoveLeadLepton']:
        for disc in getJetTagsProbs('pfMassDecorrelatedInclParticleTransformerV2%sHidLayer' % n) + getJetTagsMetaDiscrs('pfMassDecorrelatedInclParticleTransformerV2%sHidLayer' % n):
            supportedBtagDiscr[disc] = [["pfMassDecorrelatedInclParticleTransformerV2%sTagInfos" % n]]
    # update supportedMetaDiscr
    for disc in _pfMassDecorrelatedDeepHWWV1JetTagsMetaDiscrs:
        supportedMetaDiscr[disc] = _pfMassDecorrelatedDeepHWWV1JetTagsProbs
    for disc in _pfMassDecorrelatedInclParticleTransformerV1JetTagsMetaDiscrs:
        supportedMetaDiscr[disc] = _pfMassDecorrelatedInclParticleTransformerV1JetTagsProbs
    for disc in _pfMassDecorrelatedInclParticleTransformerV2JetTagsMetaDiscrs:
        supportedMetaDiscr[disc] = _pfMassDecorrelatedInclParticleTransformerV2JetTagsProbs
    for disc in _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsMetaDiscrs:
        supportedMetaDiscr[disc] = _pfMassDecorrelatedInclParticleTransformerAK15V2JetTagsProbs
    for disc in _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsMetaDiscrs:
        supportedMetaDiscr[disc] = _pfMassDecorrelatedInclParticleTransformerV2HidLayerJetTagsProbs
    for n in ['RemoveLeadLepton']:
        for disc in getJetTagsMetaDiscrs('pfMassDecorrelatedInclParticleTransformerV2%sHidLayer' % n):
            supportedMetaDiscr[disc] = getJetTagsProbs('pfMassDecorrelatedInclParticleTransformerV2%sHidLayer' % n)

    return supportedBtagInfos, supportedBtagDiscr, supportedMetaDiscr

## Import TagInfos additional to RecoBTag_cff
from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedDeepHWWV1_cff import *
from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV1_cff import *
from DeepNTuples.Ntupler.hwwTagger.pfMassDecorrelatedInclParticleTransformerV2_cff import *
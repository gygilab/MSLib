#include "Match.hpp"

map<string,SpectralMatch::ptr2scoringFnc> SpectralMatch::scoringFncPtrMap;

// Match peptide and manage storage in peptideSet/peptideSequenceSet on the fly as needed
SpectralMatch::SpectralMatch(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, Proteome& proteome, tr1::unordered_map<Peptide,int>::iterator p_it, ModIterator &mIt) {
    int max_charge;
    double dta_lo, dta_hi;
    initScoringFunctions();
    getSpectralRangeParams(specItLo, specItHi, max_charge, dta_lo, dta_hi);
    ft = new PeptideFragmentTable(mIt, max_charge, dta_lo, dta_hi);
    for (multimap<double, Spectrum>::iterator specIt = specItLo; specIt != specItHi; specIt++) {
        dta = &((*specIt).second);
        computeScore();
        if (score1 > 0 && score1 > dta->_hits.getBottomHitScore1()) {
            populateFragmentDevItyVecs();
            computeSumFragIty();

            sumFilteredIty = dta->getFilteredSumIty(DigestOptions::maxMS2ItyCutoffForTIC);
            proteome.incrementPeptideHit(p_it);
            SearchResult * psr = getSearchResult(p_it->first);
            dta->_hits.storeHit(psr, proteome);
        }
    }
    delete ft;
    return;
}

//// Match peptide for a peptide-indexed search where peptideSet is already pre-populated and peptideSequenceSet is not needed
//SpectralMatch::SpectralMatch(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, ModIterator &mIt, const Peptide &pep) {
//    int max_charge;
//    double dta_lo, dta_hi;
//    initScoringFunctions();
//    getSpectralRangeParams(specItLo, specItHi, max_charge, dta_lo, dta_hi);
//    ft = new PeptideFragmentTable(mIt, max_charge, dta_lo, dta_hi);
//    for (multimap<double, Spectrum>::iterator specIt = specItLo; specIt != specItHi; specIt++) {
//        dta = &((*specIt).second);
//        computeScore();
//        if (score1 > 0.0 && score1 > dta->_hits.getBottomHitScore1()) {
//            SearchResult * psr = getSearchResult(pep);
//            dta->_hits.storeHit(psr);
//        }
//    }
//    delete ft;
//    return;
//}

// Match x-linked peptide to get matching info without storing hits; currently used only for noModDeltaScore calc
SpectralMatch::SpectralMatch(multimap<double, Spectrum>::iterator specIt, ModIterator &mIt, ModIterator &mIt2) {
    int max_charge = specIt->second._charge;
    double dta_lo = specIt->second.getMinMZ();
    double dta_hi = max(specIt->second.getMaxMZ(), specIt->second.getM_Z());
    initScoringFunctions();
    ft = new XLinkedFragmentTable(mIt, mIt2, -1, -1, max_charge, dta_lo, dta_hi);
    
    dta = &((*specIt).second);
    computeScore();

    delete ft;
    return;
}

// Match x-linked peptides (assumes peptides themselves are already pre-indexed in peptideSet)
SpectralMatch::SpectralMatch(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, Proteome& proteome, ModIterator &mIt, ModIterator &mIt2, int &pos1, int &pos2) {
    int max_charge;
    double dta_lo, dta_hi;
    initScoringFunctions();
    getSpectralRangeParams(specItLo, specItHi, max_charge, dta_lo, dta_hi);
    ft = new XLinkedFragmentTable(mIt, mIt2, pos1, pos2, max_charge, dta_lo, dta_hi);
    bool mirror_image_xlinks = false;
    if (mIt.modSequence.compare(mIt2.modSequence) == 0 && pos1 == pos2) {
        mirror_image_xlinks = true;
    }
    for (multimap<double, Spectrum>::iterator specIt = specItLo; specIt != specItHi; specIt++) {
        dta = &((*specIt).second);
        computeScore();
        if (score1 > 0.0 && score1 > dta->_hits.getBottomHitScore1()) {
            populateFragmentDevItyVecs();
            meanFragmentDeviation1 = JMUtil::vectorMean(fragmentDeviation1, fragmentIty1);
            meanFragmentDeviation2 = JMUtil::vectorMean(fragmentDeviation2, fragmentIty2);
            stdFragmentDeviation1 = JMUtil::vectorStd(fragmentDeviation1, fragmentIty1);
            stdFragmentDeviation2 = JMUtil::vectorStd(fragmentDeviation2, fragmentIty2);
            
            computeSumFragIty();

            sumFilteredIty = dta->getFilteredSumIty(DigestOptions::maxMS2ItyCutoffForTIC);

            double hi_score = max(OM1Score1, OM2Score1);
            maxOM_delta_score = (score1 - hi_score)/score1;

            if (maxOM_delta_score <= DigestOptions::xFilterOMDiff && !mirror_image_xlinks) {
                continue;
            }

            SpectralMatch xsm(specIt, mIt, mIt2);
            noModDeltaScoreX = (score1 - xsm.score1)/score1;
            noModnMatchedPeaks = xsm.nMatchedPeaks;

            // Peaks from 2nd peptide always match first, so make some (occasionally arbitrary) adjustments
            if (mirror_image_xlinks) {
                meanFragmentDeviation1 = meanFragmentDeviation2;
                stdFragmentDeviation1 = 1.0;
                OM1Score1 = OM2Score1/2.0;
                OM2Score1 = OM1Score1;
                nMatchedPeaks1 = floor(nMatchedPeaks2/2.0);
                nMatchedPeaks2 = nMatchedPeaks1;
                nIons1 = nIons2;
            }
            
            if (nMatchedPeaks1 < 1 || nMatchedPeaks2 < 1) {
                continue;
            }

            if (OM1Score1 >= OM2Score1) {
                frag_z = fabs( (meanFragmentDeviation1 - meanFragmentDeviation2) / stdFragmentDeviation1);
            } else {
                frag_z = fabs( (meanFragmentDeviation1 - meanFragmentDeviation2) / stdFragmentDeviation2);
            }
            
            // TODO: check here that it's a hi-res run??
            if (frag_z > DigestOptions::xFilterFragDevZ) {
                // This includes cases where the frag_z = inf because std = 0
                // These turn out to be cases where a single peak was matched AND this was the highest-scoring peptide of the two
                continue;
            }

            proteome.incrementPeptideHit(*mIt.pep);
            proteome.incrementPeptideHit(*mIt2.pep);
            SearchResult * xsr = getSearchResult(*(mIt.pep), *(mIt2.pep));
            dta->_hits.storeHit(xsr);
        }
    }
    delete ft;
    return;
}

// Match loop-linked peptides (assumes peptides themselves are already pre-indexed in peptideSet)
SpectralMatch::SpectralMatch(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, Proteome& proteome, tr1::unordered_map<Peptide,int>::iterator p_it, ModIterator &mIt, int &pos1, int &pos2) {
    int max_charge;
    double dta_lo, dta_hi;
    initScoringFunctions();
    getSpectralRangeParams(specItLo, specItHi, max_charge, dta_lo, dta_hi);
    ft = new LoopLinkedFragmentTable(mIt, pos1, pos2, max_charge, dta_lo, dta_hi);
    for (multimap<double, Spectrum>::iterator specIt = specItLo; specIt != specItHi; specIt++) {
        dta = &((*specIt).second);
        computeScore();
        if (score1 > 0.0 && score1 > dta->_hits.getBottomHitScore1()) {
            populateFragmentDevItyVecs();
            computeSumFragIty();

            sumFilteredIty = dta->getFilteredSumIty(DigestOptions::maxMS2ItyCutoffForTIC);
            proteome.incrementPeptideHit(p_it);
            SearchResult * xsr = getLoopLinkedSearchResult(*(mIt.pep));
            dta->_hits.storeHit(xsr, proteome);
        }
    }
    delete ft;
    return;
}

void SpectralMatch::computeScore() {
    ////cout << "... Initializing SpectralMatch ... " << endl; ////
    SpectralMatch::ptr2scoringFnc scoring_func = scoringFncPtrMap[ DigestOptions::scoringFunction ];
//    for (int i = 0; i < 1000000; i++) {
        nMatchedPeaks = 0;
        nIons = 0;
        nMatchedPeaks1 = 0;
        nIons1 = 0;
        nMatchedPeaks2 = 0;
        nIons2 = 0;
        nIonsContainingXlink = 0;
        meanFragmentDeviation1 = 0.0;
        meanFragmentDeviation2 = 0.0;
        stdFragmentDeviation1 = 0.0;
        stdFragmentDeviation2 = 0.0;
        fragmentDeviation1.clear();
        fragmentDeviation2.clear();
        fragmentIty1.clear();
        fragmentIty2.clear();
        frag_z = 0.0;
        sumFragIty1 = 0.0;
        sumFragIty2 = 0.0;
        sumFilteredIty = 0.0;
        OM1Score1 = 0.0;
        OM2Score1 = 0.0;
        maxOM_delta_score = 0.0;
        score1 = 0.0;
       (this->*scoring_func) ();
       dta->_numPSMs++;
//    }
}

void SpectralMatch::initScoringFunctions() {
    scoringFncPtrMap["Binomial"] = &SpectralMatch::scoreBinomial;
    scoringFncPtrMap["StepwiseMaxBinomial"] = &SpectralMatch::scoreStepwiseMaxBinomial;
    return;
}

SearchResult * SpectralMatch::getSearchResult(const Peptide &pep) {
    SearchResult * sr = new PeptideSearchResult(ft->peptide, pep);
    sr->score1 = score1;
//    sr->meanFragmentDeviation = getMeanFragmentDeviation();
    sr->dtaPath = dta->_dtaPath;
    sr->nMatches = nMatchedPeaks;
    sr->nIons = nIons;
    sr->theo_m_h = ft->m_h_mass;
    sr->ppm = getMatchedPrecursorDeviation(); ////
//    sr->debugInfo = getDebugInfo();
    sr->fracItyExplained = sumFragIty1/sumFilteredIty;
    return sr;
}

SearchResult * SpectralMatch::getLoopLinkedSearchResult(const Peptide &pep) {
    string pos_info = toStr(ft->pos1 + 1) + "-" + toStr(ft->pos2 + 1) + ":" + toStr(pep._start + ft->pos1) + "-" + toStr(pep._start + ft->pos2);
    SearchResult * sr = new LoopLinkedSearchResult(ft->peptide, pep, pos_info);

    sr->score1 = score1;
//    sr->meanFragmentDeviation = getMeanFragmentDeviation();
    sr->dtaPath = dta->_dtaPath;
    sr->nMatches = nMatchedPeaks;
    sr->nIons = nIons;
    sr->theo_m_h = ft->m_h_mass;
    sr->ppm = getMatchedPrecursorDeviation(); ////
    sr->fracItyExplained = sumFragIty1/sumFilteredIty;
    return sr;
}

SearchResult * SpectralMatch::getSearchResult(const Peptide &pep, const Peptide &pep2) {
    SearchResult * xsr;
    string frac1 = toStr(nMatchedPeaks1) + "/" + toStr(nIons1);
    string frac2 = toStr(nMatchedPeaks2) + "/" + toStr(nIons2);
    string matched_info = toStr(frag_z, 2) + " " + toStr(nIonsContainingXlink) + " " + toStr(maxOM_delta_score,2) + " " + toStr(noModDeltaScoreX,2) + "_" + toStr(nMatchedPeaks1 + nMatchedPeaks2 - noModnMatchedPeaks);
    string matched_frags_12 = frac1 + " " + frac2 + " " + toStrE(JMUtil::vectorSum(fragmentIty1), 3) + "_" + toStrE(JMUtil::vectorSum(fragmentIty2), 3) + " " + toStr(meanFragmentDeviation1, 2) + "_" + toStr(meanFragmentDeviation2, 2) + " " + toStr(OM1Score1,2) + "_" + toStr(OM2Score1,2);
    string matched_frags_21 = frac2 + " " + frac1 + " " + toStrE(JMUtil::vectorSum(fragmentIty2), 3) + "_" + toStrE(JMUtil::vectorSum(fragmentIty1), 3) + " " + toStr(meanFragmentDeviation2, 2) + "_" + toStr(meanFragmentDeviation1, 2) + " " + toStr(OM2Score1,2) + "_" + toStr(OM1Score1,2);
    string pos_12_info = toStr(ft->pos1 + 1) + "-" + toStr(ft->pos2 + 1) + ":" + toStr(pep._start + ft->pos1) + "-" + toStr(pep2._start + ft->pos2);
    string pos_21_info = toStr(ft->pos2 + 1) + "-" + toStr(ft->pos1 + 1) + ":" + toStr(pep2._start + ft->pos2) + "-" + toStr(pep._start + ft->pos1);
    // Order by _reference, peptide._start, peptide.length()
    if (pep.getProteinReference() < pep2.getProteinReference()) {
        string pos_info = pos_12_info;
        xsr = new XLinkedSearchResult(ft->peptide, ft->peptide2, pep, pep2, pos_info);
        xsr->addMatchInfo = matched_frags_12 + " " + matched_info;
    } else if (pep.getProteinReference() > pep2.getProteinReference()) {
        string pos_info = pos_21_info;
        xsr = new XLinkedSearchResult(ft->peptide2, ft->peptide, pep2, pep, pos_info);
        xsr->addMatchInfo = matched_frags_21 + " " + matched_info;
    } else { // (*it).first._reference == (*it2).first._reference
        if (pep._start < pep2._start
                ||
                (pep._start == pep2._start && pep._sequence.length() < pep2._sequence.length())) {
            string pos_info = pos_12_info;
            xsr = new XLinkedSearchResult(ft->peptide, ft->peptide2, pep, pep2, pos_info);
            xsr->addMatchInfo = matched_frags_12 + " " + matched_info;
        } else {
            string pos_info = pos_21_info;
            xsr = new XLinkedSearchResult(ft->peptide2, ft->peptide, pep2, pep, pos_info);
            xsr->addMatchInfo = matched_frags_21 + " " + matched_info;
        }
    }
    xsr->score1 = score1;
//    xsr->meanFragmentDeviation = getMeanFragmentDeviation();
    xsr->dtaPath = dta->_dtaPath;
    xsr->nMatches = nMatchedPeaks;
    xsr->nIons = nIons;
    xsr->theo_m_h = ft->m_h_mass;
    xsr->ppm = getMatchedPrecursorDeviation(); ////
    
    xsr->fracItyExplained = (sumFragIty1 + sumFragIty2)/sumFilteredIty;
//    xsr->debugInfo = getDebugInfo();
    return xsr;
}

void SpectralMatch::computeSumFragIty() {
    sumFragIty1 = 0.0;
    for (vector<double>::iterator vecIt = fragmentIty1.begin(); vecIt != fragmentIty1.end(); vecIt++) {
        if (*vecIt > DigestOptions::maxMS2ItyCutoffForTIC)
            sumFragIty1 += *vecIt;
    }
    sumFragIty2 = 0.0;
    for (vector<double>::iterator vecIt = fragmentIty2.begin(); vecIt != fragmentIty2.end(); vecIt++) {
        if (*vecIt > DigestOptions::maxMS2ItyCutoffForTIC)
            sumFragIty2 += *vecIt;
    }
    return;
}

void SpectralMatch::scoreBinomial() {
    matchTheoreticalPeaks();
    double p = (double) dta->_topPeaks * (MSMass::fragmentMass2Diff(1000.0) * 2) / ((double) DigestOptions::windowSize);
    double binomial_cdf = JMUtil::calcBinomialCdfUpper(nMatchedPeaks, p, nIons);
    score1 = -10 * log10(binomial_cdf);

    if (nIons2 > 0) {
        double binomial_cdf1 = JMUtil::calcBinomialCdfUpper(nMatchedPeaks1, p, nIons1);
        OM1Score1 = -10 * log10(binomial_cdf1);
        double binomial_cdf2 = JMUtil::calcBinomialCdfUpper(nMatchedPeaks2, p, nIons2);
        OM2Score1 = -10 * log10(binomial_cdf2);
    }
}

void SpectralMatch::scoreStepwiseMaxBinomial() {
    //ft->generatePeakList(dta->getMinMZ(), dta->getMaxMZ(), dta->_charge);
    matchTheoreticalPeaks();

    vector<int> n_matches_by_peak_depth(DigestOptions::peakDepth + 1, 0);
    for (tr1::unordered_map<TheoreticalPeak*,SpectrumPeak*>::iterator teppmIt = theoExpPtrPairMap.begin(); teppmIt != theoExpPtrPairMap.end(); teppmIt++) {
        n_matches_by_peak_depth[(teppmIt->second)->rank]++;
    }
    score1 = 0;
    double max_binomial = 0.0;
    double binomial_score = 0.0;
    double p = 0.0;
    int n_cum_matches = 0;
    int n_trials = nIons; //
    int peak_depth = 1;
    while (peak_depth <= DigestOptions::peakDepth && n_cum_matches < nMatchedPeaks) {
        n_cum_matches += n_matches_by_peak_depth[peak_depth];
        p = (double) peak_depth * (MSMass::fragmentMass2Diff(1000.0) * 2) / ((double) DigestOptions::windowSize);
        binomial_score = -10 * log10(calcBinomialCdfUpper(n_cum_matches, p, n_trials));
        if (binomial_score > max_binomial || binomial_score == 0) {
            max_binomial = binomial_score;
            peak_depth++;
        } else {
            score1 += max_binomial;
            n_trials -= (n_cum_matches - n_matches_by_peak_depth[peak_depth]);
            n_cum_matches = 0;
            max_binomial = 0.0;
        }
    }
    score1 += max_binomial;
    
    if (nIons2 > 0) {
        scoreStepwiseMaxBinomialX1();
        scoreStepwiseMaxBinomialX2();
    }
    return;
}

void SpectralMatch::scoreStepwiseMaxBinomialX1() {
    vector<int> n_matches_by_peak_depth(DigestOptions::peakDepth + 1, 0);
    for (tr1::unordered_map<TheoreticalPeak*,SpectrumPeak*>::iterator teppmIt = theoExpPtrPairMap.begin(); teppmIt != theoExpPtrPairMap.end(); teppmIt++) {
        if (teppmIt->first->pep_1_2 == 1) { // make sure this behavior of counting linker fragments matches what's counted in matchTheoreticalPeaks
            n_matches_by_peak_depth[(teppmIt->second)->rank]++;
        }
    }
    OM1Score1 = 0;
    double max_binomial = 0.0;
    double binomial_score = 0.0;
    double p = 0.0;
    int n_cum_matches = 0;
    int n_trials = nIons1; //
    int peak_depth = 1;
    while (peak_depth <= DigestOptions::peakDepth && n_cum_matches < nMatchedPeaks1) {
        n_cum_matches += n_matches_by_peak_depth[peak_depth];
        p = (double) peak_depth * (MSMass::fragmentMass2Diff(1000.0) * 2) / ((double) DigestOptions::windowSize);
        binomial_score = -10 * log10(calcBinomialCdfUpper(n_cum_matches, p, n_trials));
        if (binomial_score > max_binomial || binomial_score == 0) {
            max_binomial = binomial_score;
            peak_depth++;
        } else {
            OM1Score1 += max_binomial;
            n_trials -= (n_cum_matches - n_matches_by_peak_depth[peak_depth]);
            n_cum_matches = 0;
            max_binomial = 0.0;
        }
    }
    OM1Score1 += max_binomial;
    return;
}

void SpectralMatch::scoreStepwiseMaxBinomialX2() {
    vector<int> n_matches_by_peak_depth(DigestOptions::peakDepth + 1, 0);
    for (tr1::unordered_map<TheoreticalPeak*,SpectrumPeak*>::iterator teppmIt = theoExpPtrPairMap.begin(); teppmIt != theoExpPtrPairMap.end(); teppmIt++) {
        if (teppmIt->first->pep_1_2 == 2) { // make sure this behavior of counting linker fragments matches what's counted in matchTheoreticalPeaks
            n_matches_by_peak_depth[(teppmIt->second)->rank]++;
        }
    }
    OM2Score1 = 0;
    double max_binomial = 0.0;
    double binomial_score = 0.0;
    double p = 0.0;
    int n_cum_matches = 0;
    int n_trials = nIons2; //
    int peak_depth = 1;
    while (peak_depth <= DigestOptions::peakDepth && n_cum_matches < nMatchedPeaks2) {
        n_cum_matches += n_matches_by_peak_depth[peak_depth];
        p = (double) peak_depth * (MSMass::fragmentMass2Diff(1000.0) * 2) / ((double) DigestOptions::windowSize);
        binomial_score = -10 * log10(calcBinomialCdfUpper(n_cum_matches, p, n_trials));
        if (binomial_score > max_binomial || binomial_score == 0) {
            max_binomial = binomial_score;
            peak_depth++;
        } else {
            OM2Score1 += max_binomial;
            n_trials -= (n_cum_matches - n_matches_by_peak_depth[peak_depth]);
            n_cum_matches = 0;
            max_binomial = 0.0;
        }
    }
    OM2Score1 += max_binomial;
    return;
}

/*
 * Find the closest matching experimental/observed peak pairs within a tolerance window
 */
void SpectralMatch::matchTheoreticalPeaks() {
    ////cout << "== SpectralMatch::matchTheoreticalPeaks ==" << endl; ////
    vector<SpectrumPeak>::iterator spec_peak_it;
    vector<SpectrumPeak>::iterator prev_spec_peak_it = dta->_peakVector.begin();
    vector<TheoreticalPeak>::iterator theo_peak_it;
    multimap<int,expThPtrPair> deviation_2_exp_th_ptr_pair;
    theoExpPtrPairMap.clear();
    tr1::unordered_multimap<SpectrumPeak*,TheoreticalPeak*> exp_th_ptr_pair_map;
    vector<TheoreticalPeak>::iterator vec_it_lo = lower_bound(ft->fullPeakList.begin(), ft->fullPeakList.end(), MSMass::double2int( dta->getMinMZ() ), TheoreticalPeakCompare());
    vector<TheoreticalPeak>::iterator vec_it_hi = upper_bound(ft->fullPeakList.begin(), ft->fullPeakList.end(), MSMass::double2int( max(dta->getMaxMZ(), dta->getM_Z()) ), TheoreticalPeakCompare());
    int max_fragment_charge = dta->_charge > 1 ? (dta->_charge - 1) : 1;
    for (theo_peak_it = vec_it_lo; theo_peak_it != vec_it_hi; theo_peak_it++) {
        /* Theoretical fragments were created for a spectral range which means we may have charge states here that are 
         * above what's appropriate for this particular spectrum */
        if (theo_peak_it->charge > max_fragment_charge)
            continue;
        
        nIons++;
        TheoreticalPeak tp = *theo_peak_it;
        if        (tp.pep_1_2 == 2 && tp.ion_type != 'L') { // make sure this behavior of counting linker fragments matches what's counted in scoreStepwiseMaxBinomialX1 and 2
            nIons2++;
        } else if (tp.pep_1_2 == 1 && tp.ion_type != 'L') { // make sure this behavior of counting linker fragments matches what's counted in scoreStepwiseMaxBinomialX1 and 2
            nIons1++;
        }
        int mass_diff = MSMass::fragmentMass2Diff(tp.m_z);
        
        // advance the iterator to start of match window
        spec_peak_it = prev_spec_peak_it;
        while (spec_peak_it != dta->_peakVector.end() && spec_peak_it->m_z < (tp.m_z - mass_diff)) {
            spec_peak_it++;
        }
        prev_spec_peak_it = spec_peak_it;

        while (spec_peak_it != dta->_peakVector.end() && spec_peak_it->m_z <= (tp.m_z + mass_diff)) {
            // if charge state was assigned during de-isotoping, make sure it matches theoretical
            if (spec_peak_it->charge == 0 || (spec_peak_it->charge > 0 && tp.charge == spec_peak_it->charge) ) {
                int obs_peak_deviation = abs(tp.getPeakDeviation((*spec_peak_it)));
                // store experimental->theoretical map pairings and the deviation between them
                deviation_2_exp_th_ptr_pair.insert(pair<int,expThPtrPair> (obs_peak_deviation, expThPtrPair(&(*spec_peak_it), &(*theo_peak_it))) );
                exp_th_ptr_pair_map.insert(expThPtrPair(&(*spec_peak_it), &(*theo_peak_it)));
            }
            spec_peak_it++;
        }
    }
    for (multimap<int,expThPtrPair>::iterator d2etppIt = deviation_2_exp_th_ptr_pair.begin(); d2etppIt != deviation_2_exp_th_ptr_pair.end(); d2etppIt++) {
        if (exp_th_ptr_pair_map.find( (d2etppIt->second).first ) != exp_th_ptr_pair_map.end()) {
            // store new theoretical-experimental pairs, thus keeping only the best match per theoretical peak
            theoExpPtrPairMap.insert( thExpPtrPair( (d2etppIt->second).second, (d2etppIt->second).first ) );
            exp_th_ptr_pair_map.erase( (d2etppIt->second).first );
        }
    }
    nMatchedPeaks = theoExpPtrPairMap.size();
    if (nIons2 > 0) {
        for (tr1::unordered_map<TheoreticalPeak*,SpectrumPeak*>::iterator teppMIt = theoExpPtrPairMap.begin(); teppMIt != theoExpPtrPairMap.end(); teppMIt++) {
            if        (teppMIt->first->pep_1_2 == 2) {
                if (teppMIt->first->ion_type != 'L') { // make sure this behavior of counting linker fragments matches what's counted in scoreStepwiseMaxBinomialX1 and 2
                    nMatchedPeaks2++;
                }
            } else if (teppMIt->first->pep_1_2 == 1) {
                if (teppMIt->first->ion_type != 'L') { // make sure this behavior of counting linker fragments matches what's counted in scoreStepwiseMaxBinomialX1 and 2
                    nMatchedPeaks1++;
                }
            }
            if (teppMIt->first->ion_type == 'B' || teppMIt->first->ion_type == 'Y')
                nIonsContainingXlink++;
        }
    }
    return;
}


// DEBUG
void SpectralMatch::printMatchedPeaks(string outfile) {
    cout << "   = SpectralMatch::printMatchedPeaks =  " << endl;
    std::ofstream outFile(outfile.c_str(), ios::out);
    cout << dta->_dtaPath << "\t\t" << ft->peptide << " - " << ft->peptide2 << "\t" << ft->pos1 << "-" << ft->pos2 << "\t" << ft->id << endl;
    vector<TheoreticalPeak>::iterator vec_it_lo = lower_bound(ft->fullPeakList.begin(), ft->fullPeakList.end(), MSMass::double2int( dta->getMinMZ() ), TheoreticalPeakCompare());
    vector<TheoreticalPeak>::iterator vec_it_hi = upper_bound(ft->fullPeakList.begin(), ft->fullPeakList.end(), MSMass::double2int( dta->getMaxMZ() ), TheoreticalPeakCompare());
    cout << fixed << setprecision(5) << "SpectralMatch: Obz: m/z Ity Rank\t Expected: m/z z fragDev ion_type pepNum" << endl;
    for (vector<TheoreticalPeak>::iterator vecIt = vec_it_lo; vecIt != vec_it_hi; vecIt++) {
        if ((*vecIt).charge > DigestOptions::maxTheoFragCharge)
            continue;
        tr1::unordered_map<TheoreticalPeak*,SpectrumPeak*>::iterator temp_it = theoExpPtrPairMap.find( &(*vecIt) );
        if (temp_it != theoExpPtrPairMap.end()) {
            outFile << fixed << setprecision(5) << "SpectralMatch: " 
                                             << MSMass::int2double( (temp_it->second)->m_z ) << " "
                                             << (temp_it->second)->intensity << " "
                                             << (temp_it->second)->rank << "\t\t"

                                             << MSMass::int2double( (temp_it->first)->m_z ) << " "
                                             <<                     (temp_it->first)->charge << " "
                                             << MSMass::int2double( (temp_it->first)->getPeakDeviation(*(temp_it->second)) ) << " "
                                             << (temp_it->first)->ion_type << " "
                                             << (temp_it->first)->pep_1_2 << endl;
        }
    }
    outFile << nMatchedPeaks << " / " << nIons << " / " << dta->_peakVector.size() << endl;
    cout << nMatchedPeaks << " / " << nIons << " / " << dta->_peakVector.size() << endl;
    outFile.close();
    return;
}

// DEBUG
string SpectralMatch::getDebugInfo() {
//    cout << "   = SpectralMatch::getDebugInfo =  " << endl;
    stringstream oss;
    oss << dta->_dtaPath << "\t\t" << ft->peptide << " - " << ft->peptide2 << "\t" << ft->pos1 << "-" << ft->pos2 << "\t"
        << nMatchedPeaks << " / " << nIons << " / " << dta->_peakVector.size() << endl;
    vector<TheoreticalPeak>::iterator vec_it_lo = lower_bound(ft->fullPeakList.begin(), ft->fullPeakList.end(), MSMass::double2int( dta->getMinMZ() ), TheoreticalPeakCompare());
    vector<TheoreticalPeak>::iterator vec_it_hi = upper_bound(ft->fullPeakList.begin(), ft->fullPeakList.end(), MSMass::double2int( dta->getMaxMZ() ), TheoreticalPeakCompare());
    for (vector<TheoreticalPeak>::iterator vecIt = vec_it_lo; vecIt != vec_it_hi; vecIt++) {
        if ((*vecIt).charge > dta->_charge - 1)
            continue;
        tr1::unordered_map<TheoreticalPeak*,SpectrumPeak*>::iterator temp_it = theoExpPtrPairMap.find( &(*vecIt) );
        if (temp_it != theoExpPtrPairMap.end()) {
            oss << fixed << setprecision(3) << MSMass::int2double( (temp_it->second)->m_z ) << "\t"
                                            << (temp_it->second)->intensity << "\t"
                                            << (temp_it->second)->rank << "\t"

                                            << MSMass::int2double( (temp_it->first)->m_z ) << "\t"
                                            <<                     (temp_it->first)->charge << "\t"
                                            << MSMass::int2double( (temp_it->first)->getPeakDeviation(*(temp_it->second)) ) << "\t"
                                            << (temp_it->first)->ion_type << "\t"
                                            << (temp_it->first)->pep_1_2 << endl;
        }
    }
    return oss.str();
}

/* Get charge and min/max m/z for a single spectrum or max charge, min m/z, max m/z for a continuous range of spectra. 
 * This helps in reducing calls to generateFullPeakList in Fragment by creating a somewhat generalized FragmentTable 
 * that can be re-used to match against multiple spectra.  MatchTheoreticalPeaks will only consider peaks appropriate for
 * the spectrum at hand if higher charged fragments are given. 
 */
void SpectralMatch::getSpectralRangeParams(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, int &c, double &lo, double &hi) {
    c = specItLo->second._charge;
    lo = specItLo->second.getMinMZ();
    hi = max(specItLo->second.getMaxMZ(), specItLo->second.getM_Z());
    specItLo++; //// ??
    while (specItLo != specItHi) {
        double cur_charge = specItLo->second._charge;
        if (cur_charge > c)
            c = cur_charge;

        double cur_lo = specItLo->second.getMinMZ();
        if (cur_lo < lo)
            lo = cur_lo;

        double cur_hi = max(specItLo->second.getMaxMZ(), specItLo->second.getM_Z());
        if (cur_hi > hi)
            hi = cur_hi;

        specItLo++;
    }
    return;
}

double SpectralMatch::getMatchedPrecursorDeviation() {
    double diff = dta->_M_H - ft->m_h_mass;
//    if (DigestOptions::precursorTolRelative)
        return diff / ft->m_h_mass * 1000000.0;
//    else
//        return diff;
}


void SpectralMatch::populateFragmentDevItyVecs() {
    for (tr1::unordered_map<TheoreticalPeak*,SpectrumPeak*>::iterator teppmIt = theoExpPtrPairMap.begin(); teppmIt != theoExpPtrPairMap.end(); teppmIt++) {
        if (teppmIt->first->pep_1_2 == 2) {
            fragmentDeviation2.push_back(MSMass::int2double((teppmIt->first)->getPeakDeviation(*(teppmIt->second))));
            fragmentIty2.push_back( teppmIt->second->intensity );
        } else if (teppmIt->first->pep_1_2 == 1) {
            fragmentDeviation1.push_back(MSMass::int2double((teppmIt->first)->getPeakDeviation(*(teppmIt->second))));
            fragmentIty1.push_back( teppmIt->second->intensity );
       }
    }
    return;
}

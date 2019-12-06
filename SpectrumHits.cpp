#include "SpectrumHits.hpp"

// Peptide-indexed search.  Spectra may or may not be indexed
void SpectrumHits::storeHit (SearchResult * xsr) {
    _write_mutex->lock();
    hitMap.insert(xsr);
    if (hitMap.size() > (DigestOptions::maxHitRank + 1)) {
        delete *(hitMap.begin());
        hitMap.erase(hitMap.begin());
    }
    _write_mutex->unlock();
    return;
}

// "standard" spectrum-indexed peptide search, where only top-N scoring peptides get stored/removed on the fly as needed
void SpectrumHits::storeHit (SearchResult * xsr, Proteome &proteome) {
    _write_mutex->lock();
    hitMap.insert(xsr);
    if (hitMap.size() > (DigestOptions::maxHitRank + 1)) {
        // Only storing top N matches so remove the peptide from peptideSet when if falls off the priority queue
        // the extra map.find() here doesn't appear to cost much
        proteome.decrementPeptideHit(*(*hitMap.begin())->pepPtr);
        delete *(hitMap.begin());
        hitMap.erase(hitMap.begin());
    }
    _write_mutex->unlock();
    return;
}

void SpectrumHits::output(ofstream &outFile) {
    if (hitMap.empty())
        return;
    calcDeltaScore();
    multiset<SearchResult *, SearchResultPtrCompare>::reverse_iterator XSRmapIt;
    int top_hit_count = 0;
    for (XSRmapIt = hitMap.rbegin(); XSRmapIt != hitMap.rend() && top_hit_count < DigestOptions::maxHitOutputRank; XSRmapIt++) {
        outFile << (*XSRmapIt)->getOutputString() << endl;
//        (*XSRmapIt)->printDebugInfo(); ////
        top_hit_count++;
    }
    return;
}

// Debugging
void SpectrumHits::output() {
    if (hitMap.empty())
        return;
    multiset<SearchResult *, SearchResultPtrCompare>::reverse_iterator XSRmapIt;
    for (XSRmapIt = hitMap.rbegin(); XSRmapIt != hitMap.rend(); XSRmapIt++) {
        cout << (*XSRmapIt)->getOutputString() << endl;
    }
    return;
}

double SpectrumHits::getTopHitScore1() {
    return (*(hitMap.rbegin()))->score1;
}

double SpectrumHits::getBottomHitScore1() {
    double score = -1.0;
    _write_mutex->lock();
    if (hitMap.size() < (DigestOptions::maxHitRank + 1)) {
        _write_mutex->unlock();
    } else {
        score = (*hitMap.begin())->score1; // data race here?
        _write_mutex->unlock();
    }
    return score;
}

void SpectrumHits::calcDeltaScore() {
    double diffSeqDeltaScore = -1.0;

    multiset<SearchResult *, SearchResultPtrCompare>::reverse_iterator XSRmapItTop = hitMap.rbegin();
    double topScore = (*XSRmapItTop)->score1;

    /* Deal with edge cases where there is only one hit for a spectrum.  This is more likely for small db searches
     * Theoretically, in such cases, the default delta/dsDelta values should be 1.0 but that might affect the delta/dsDelta
     * score distributions for the LDA (-1.0 is not even a possible delta value), so it's safer to default them to 0.0
     */
    multiset<SearchResult *, SearchResultPtrCompare>::iterator XSRmapItNplusOne = hitMap.begin();
    (*XSRmapItNplusOne)->deltaScore1 = 0.0;
    (*XSRmapItNplusOne)->diffSeqDeltaScore1 = 0.0;
    if (hitMap.size() == 1) {
        return;
    }

    // Calculate deltas from the bottom up, carrying the diffSeqDelta when sequences are the same
    multiset<SearchResult *, SearchResultPtrCompare>::iterator XSRmapItN = XSRmapItNplusOne;
    XSRmapItN++;
    for (; XSRmapItN != hitMap.end(); XSRmapItN++) {
        (*XSRmapItN)->deltaScore1 = (topScore - (*XSRmapItNplusOne)->score1) / topScore;

        string sequence_n_plus_one = (*XSRmapItNplusOne)->getSequence();
        sequence_n_plus_one = MSMass::makeIsobaric(sequence_n_plus_one);

        string sequence_n = (*XSRmapItN)->getSequence();
        sequence_n = MSMass::makeIsobaric(sequence_n);

        if (sequence_n_plus_one.compare(sequence_n) != 0 || diffSeqDeltaScore == -1.0) {
            (*XSRmapItN)->diffSeqDeltaScore1 = (*XSRmapItN)->deltaScore1;
            diffSeqDeltaScore = (*XSRmapItN)->deltaScore1;
        } else {
            (*XSRmapItN)->diffSeqDeltaScore1 = diffSeqDeltaScore;
        }

        XSRmapItNplusOne++;
    }
    return;
}

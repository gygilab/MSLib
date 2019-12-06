#include "SearchResult.hpp"

string SearchResult::getHeaderString() {
    string header = string("#dtaPath") + "\t"
            + "Score" + "\t"
            + "fracItyExplained" + "\t"
            + "diffSeqDeltaScore" + "\t"
            + "nMatches/nIons" +"\t" //  + " " + "additionalMatchInfo"
            + "theo_m_h" + "\t"
            + "ppm" + "\t"
            + "peptide" + "\t"
            + "protein" + "\t"
            + "redun" + "\t"
            + "posInfo" + "\t"
            + "searchType";
    return header;
}

PeptideSearchResult::PeptideSearchResult(string &p, const Peptide &pep) {
    deltaScore1 = -1;
    diffSeqDeltaScore1 = -1;
    searchType = "peptide";
    pepPtr = &pep;
    peptide = p;
}

string PeptideSearchResult::getSequence() {
    return pepPtr->_sequence;
}

string PeptideSearchResult::getOutputString() const {
    string output = (dtaPath + "\t"
            + toStr(score1) + "\t"
            + toStr(fracItyExplained) + "\t"
//            + toStr(deltaScore1) + "\t"
//            + toStr(meanFragmentDeviation) + "\t"
            + toStr(diffSeqDeltaScore1) + "\t"
            + toStr(nMatches) + "/" + toStr(nIons) + "\t"
            + toStr(theo_m_h) + "\t"
            + toStr(ppm) + "\t"
            + pepPtr->_flankResNterm + "." + peptide + "." + pepPtr->_flankResCterm + " - \t"
            + pepPtr->getProteinReference() + "\t"
            + toStr(pepPtr->getNumRedundantRefs()) + "\t"
            + toStr(pepPtr->_start) + "\t"
            + searchType);
//    output = output + "\n" + peakDeviationInfo; ////
    return output;
}


XLinkedSearchResult::XLinkedSearchResult(string &p, string &p2, const Peptide &pep, const Peptide &pep2, string pos_info) {
    deltaScore1 = -1;
    diffSeqDeltaScore1 = -1;
    searchType = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkType;
    pepPtr = &pep;
    pepPtr2 = &pep2;
    peptide = p;
    peptide2 = p2;
    posInfo = pos_info;
}

string XLinkedSearchResult::getSequence() {
    return (pepPtr->_sequence + "-" + pepPtr2->_sequence);
}

string XLinkedSearchResult::getOutputString() const {
    string output = (dtaPath + "\t"
            + toStr(score1) + "\t"
            + toStr(fracItyExplained) + "\t"
//            + toStr(deltaScore1) + "\t"
//            + toStr(meanFragmentDeviation) + "\t"
            + toStr(diffSeqDeltaScore1) + "\t"
            + toStr(nMatches) + "/" + toStr(nIons) + " " + addMatchInfo +"\t"
            + toStr(theo_m_h) + "\t"
            + toStr(ppm) + "\t"
            + pepPtr->_flankResNterm + "." + peptide + "." + pepPtr->_flankResCterm + " - "
            + pepPtr2->_flankResNterm + "." + peptide2 + "." + pepPtr2->_flankResCterm + "\t"
            + pepPtr->getProteinReference() + " - " + pepPtr2->getProteinReference() + "\t"
            + toStr(pepPtr->getNumRedundantRefs(DigestOptions::xFilterTopProts)) + "-" + toStr(pepPtr2->getNumRedundantRefs(DigestOptions::xFilterTopProts)) + "\t"
            + posInfo + "\t"
            + searchType);
//    output = output + "\n" + peakDeviationInfo; ////
    return output;
}

LoopLinkedSearchResult::LoopLinkedSearchResult(string &p, const Peptide &pep, string pos_info) {
    deltaScore1 = -1;
    diffSeqDeltaScore1 = -1;
    searchType = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkType + "-loop";
    pepPtr = &pep;
    peptide = p;
    posInfo = pos_info;
}

string LoopLinkedSearchResult::getSequence() {
    return pepPtr->_sequence;
}

string LoopLinkedSearchResult::getOutputString() const {
    string output = (dtaPath + "\t"
            + toStr(score1) + "\t"
            + toStr(fracItyExplained) + "\t"
//            + toStr(deltaScore1) + "\t"
//            + toStr(meanFragmentDeviation) + "\t"
            + toStr(diffSeqDeltaScore1) + "\t"
            + toStr(nMatches) + "/" + toStr(nIons) + "\t"
            + toStr(theo_m_h) + "\t"
            + toStr(ppm) + "\t"
            + pepPtr->_flankResNterm + "." + peptide + "." + pepPtr->_flankResCterm + " - \t"
            + pepPtr->getProteinReference() + "\t"
            + toStr(pepPtr->getNumRedundantRefs()) + "\t"
            + posInfo + "\t"
            + searchType);
//    output = output + "\n" + peakDeviationInfo; ////
    return output;
}

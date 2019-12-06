#ifndef _MATCH_HPP
#define	_MATCH_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <tr1/unordered_map>
#include <algorithm>

#include "../JMUtil/JMUtil.hpp"
#include "Proteome.hpp"
#include "Peptide.hpp"
#include "Peak.hpp"
#include "MSMass.hpp"
#include "Spectrum.hpp"
#include "Fragment.hpp"
#include "DigestOptions.hpp"
#include "SearchResult.hpp"

using namespace std;

class Spectrum;
// can be called with two arguments thus ensuring that we store specific PeakMatch info
// OR it can be called with no arguments and thus we can store just the scores and means without
// wasting memory on vector<PeakMatch> matches
class SpectralMatch {
public:
    SpectralMatch(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, Proteome& proteome, tr1::unordered_map<Peptide,int>::iterator p_it, ModIterator &mIt);
//    SpectralMatch(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, ModIterator &mIt, const Peptide &pep);
    SpectralMatch(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, Proteome& proteome, ModIterator &mIt, ModIterator &mIt2, int &pos1, int &pos2);
    SpectralMatch(multimap<double, Spectrum>::iterator specIt, ModIterator &mIt, ModIterator &mIt2);
    SpectralMatch(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, Proteome& proteome, tr1::unordered_map<Peptide,int>::iterator p_it, ModIterator &mIt, int &pos1, int &pos2);

    void getSpectralRangeParams(multimap<double, Spectrum>::iterator specItLo, multimap<double, Spectrum>::iterator specItHi, int &c, double &lo, double &hi);

    SearchResult * getSearchResult(const Peptide &pep);
    SearchResult * getSearchResult(const Peptide &pep, const Peptide &pep2);
    SearchResult * getLoopLinkedSearchResult(const Peptide &pep);

private:
    typedef void (SpectralMatch::*ptr2scoringFnc)();
    typedef pair<SpectrumPeak*,TheoreticalPeak*> expThPtrPair;
    typedef pair<TheoreticalPeak*,SpectrumPeak*> thExpPtrPair;

    Spectrum * dta;
    FragmentTable * ft;
    tr1::unordered_map<TheoreticalPeak*,SpectrumPeak*> theoExpPtrPairMap;

    int nMatchedPeaks;
    int nIons;
    int nMatchedPeaks1;
    int nIons1;
    int nMatchedPeaks2;
    int nIons2;
    int nIonsContainingXlink;
    double meanFragmentDeviation1;
    double meanFragmentDeviation2;
    double stdFragmentDeviation1;
    double stdFragmentDeviation2;
    vector<double> fragmentDeviation1;
    vector<double> fragmentDeviation2;
    vector<double> fragmentIty1;
    vector<double> fragmentIty2;
    double sumFragIty1;
    double sumFragIty2;
    double sumFilteredIty;
    double frag_z;
    double OM1Score1;
    double OM2Score1;
    double noModDeltaScoreX;
    int noModnMatchedPeaks;
    double maxOM_delta_score;
    double score1;

    static map<string, ptr2scoringFnc> scoringFncPtrMap;
    static void initScoringFunctions();

    void computeScore();
    void matchTheoreticalPeaks();
    void scoreBinomial();
    void scoreStepwiseMaxBinomial();
    void scoreStepwiseMaxBinomialX1();
    void scoreStepwiseMaxBinomialX2();

    void printMatchedPeaks(string outfile);
    string getDebugInfo();
    double getMatchedPrecursorDeviation();
    void computeSumFragIty();
    void populateFragmentDevItyVecs();
};

#endif	/* _MATCH_HPP */


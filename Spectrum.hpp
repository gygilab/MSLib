#ifndef _SPECTRUM_HPP
#define	_SPECTRUM_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <math.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "MSMass.hpp"
#include "Peak.hpp"
#include "DigestOptions.hpp"
#include "Peptide.hpp"
#include "SpectrumHits.hpp"
#include "Isotope.hpp"

using namespace std;

class Spectrum {
public:
    string _dtaPath;
    unsigned int _firstScanNum;
    int _charge;
    double _M_H;
    int _topPeaks; // store the current number of top filtered peaks
    unsigned int _numPSMs;
    double _meanLogIntensities;
    double _stdLogIntensities;
    int _fullPeakCount;
//    vector<SpectrumPeak*> _mzLookup;
    vector<SpectrumPeak> _peakVector;
    SpectrumHits _hits;

    Spectrum(string file) {
        readDta(file);
    }

    void readDta(string dtafile);
    int deIsotope(Isotope &iso);
    void filterDta(int window, int num_peaks, bool flg_deisotope, bool flg_nl, bool flg_nl_water, bool flg_nl_phospho);
    int filterRanges(set<Range,RangeCompare> fdr_passed_bin_ranges);

    int enforcePeakFilter();
    void flipPeakFilterFlagsTrue();
    void flipPeakFilterFlagsFalse();
    double getMinMZ();
    double getMaxMZ();
    double getM_Z();
    void adjustM_H_at(int z);
    void setFirstScan();
    int getFirstScan();
    double getFilteredSumIty(double min_ity = 0.0, double min_m_z = 0.0);

    void printDta(string outfile);
    void printFilteredRankedDta(string outfile);
};

#endif	/* _SPECTRUM_HPP */


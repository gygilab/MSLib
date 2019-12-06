#ifndef _RUN_HPP
#define	_RUN_HPP

#include <iostream>
#include <map>
#include <set>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/detail/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>

#include "Spectrum.hpp"
#include "DigestOptions.hpp"
#include "Isotope.hpp"
#include "../JMUtil/JMUtil.hpp"

namespace bfs = boost::filesystem;

class SpectrumHits;
class SearchResult;

class Run {
public:
    vector<string> msnSpectraPaths;
    multimap<double, Spectrum> msnSpectra;
    static set<Range,RangeCompare> fdr_passed_bin_ranges;
    bool _spectraLoaded;

    Run ();

    void loadAllSpectra();
    void loadAllSpectraMT(int cur_thread_num = 0);
    void adjustChargeStates(int z);
    void loadSpectrum(string path);
    void outputTab (string out_file);
    double getMinM_Z();
    double getMaxM_Z();
    void deIsotope();
    void filterCommonPeaks();
    void filterWindowPeaks();
    void outputFilteredDtas();
    double getNumTotalPSMs ();
    static bool isInFilteredRange(unsigned int m_z);

private:
    boost::shared_ptr<boost::mutex> _write_mutex;
    void findSpectra();
};

#endif	/* _RUN_HPP */


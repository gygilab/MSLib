#include "Run.hpp"

set<Range,RangeCompare> Run::fdr_passed_bin_ranges;

Run::Run () : _write_mutex(new boost::mutex()) {
    _spectraLoaded = false;
    findSpectra();
}

void Run::findSpectra() {
    bfs::path p(DigestOptions::searchDir); // p reads clearer than argv[1] in the following code
    const boost::regex re("^\\w+\\.\\d+\\.\\d+\\.\\d+\\.dta$");
    try {
        if (bfs::exists(p)) {
            if (bfs::is_directory(p)) {
                cout << "In search directory " << p << " \n";
                bfs::directory_iterator dIt(p), dIt_end;
                for (; dIt != dIt_end; dIt++) {
                    if (boost::regex_match(dIt->path().filename().string(), re)) {
                        if (bfs::is_regular_file(*dIt)) {
                            msnSpectraPaths.push_back( dIt->path().string() );
                        } else {
                            cerr << endl << *dIt << " doesn't appear to be a regular file" << endl;
                        }
                    }
                }
            } else {
                cout << p << " exists, but is not a directory\n";
            }
        } else {
            cout << p << " does not exist\n";
        }
    } catch (const bfs::filesystem_error& ex) {
        cout << ex.what() << '\n';
    }

    if (msnSpectraPaths.size() > 0) {
        cout << "...  found " << msnSpectraPaths.size() << " spectra" << endl;
    } else {
        cerr << endl << "No spectra found" << endl;
        exit(0);
    }
    return;
}

// external wrapper around loadAllSpectraMT
void Run::loadAllSpectra() {
    // Check if already loaded
    if (_spectraLoaded)
        return;
    
    cout << "...Loading spectra" << endl;
    if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
        loadAllSpectraMT();
    } else {
        boost::thread_group btg;
        for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
            if (DigestOptions::verbose) {
                cout << "... starting thread # " << i << endl; ////
            }
            btg.add_thread(new boost::thread(boost::bind(&Run::loadAllSpectraMT, this, i)));
        }
        btg.join_all();
    }
    msnSpectraPaths.clear();

    if (msnSpectra.size() > 0) {
        cout << "... loaded " << msnSpectra.size() << " spectra" << endl;
        // make sure the data is appropriate for de-isotoping - for now use a somewhat arbitrary cutoff
        if (DigestOptions::deIsotope && DigestOptions::fragmentTolInt < MSMass::double2int(0.05)) {
            deIsotope();
        }
        if (DigestOptions::filterCommonCorrel) {
            filterCommonPeaks();
        }
        filterWindowPeaks();
 
        multimap<double,Spectrum>::iterator msMmIt;

        // remove spectra with few peaks remaining
        cout << "... filtering sparse spectra" << endl;
        for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end();) {
            vector<SpectrumPeak>::iterator pv_it_lo = lower_bound(msMmIt->second._peakVector.begin(), msMmIt->second._peakVector.end(), MSMass::double2int(200.0));
            int n_peaks_over_cutoff = distance(pv_it_lo, msMmIt->second._peakVector.end());
//            cout << "JMDEBUG: " << msMmIt->second._firstScanNum << "\t" << msMmIt->second._peakVector.size() << "\t" << n_peaks_over_cutoff << endl; ////
            if (n_peaks_over_cutoff < 10) {
                msnSpectra.erase(msMmIt++);
            } else {
                ++msMmIt;
            }
        }
  
    if (DigestOptions::printDeisoSpectra) {
            outputFilteredDtas();
        }

    } else {
        cerr << endl << "No spectra loaded" << endl;
        exit(0);
    }

    cout << "... " << msnSpectra.size() << " spectra remaining" << endl;
    return;
}

// Load spectra
void Run::loadAllSpectraMT(int cur_thread_num) {
    int local_cnt = 0;
    for (vector<string>::iterator vecIt = msnSpectraPaths.begin(); vecIt != msnSpectraPaths.end(); vecIt++) {
        local_cnt++;
        if ((local_cnt - cur_thread_num - 1) % DigestOptions::maxNumThreads != 0)
            continue;
        
        Spectrum * dta = new Spectrum( *vecIt );
        if (dta->_peakVector.size() < DigestOptions::minRequiredPeaks) {
            delete dta;
            continue;
        }
        _write_mutex->lock();
        msnSpectra.insert(pair<double,Spectrum>(dta->_M_H, *dta));
//        cout << "... loaded " << msnSpectra.size() << " spectra\r" << flush; // outputting progress here somehow sometimes causes CPU inefficiency
        _write_mutex->unlock();
        delete dta;

    }
    return;
}

// Use this to find cases of incorrect charge assignment by Thermo
void Run::adjustChargeStates(int z) {
    multimap<double, Spectrum> temp_spectra;
    multimap<double,Spectrum>::iterator msMmIt;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        Spectrum s = msMmIt->second;
        s.adjustM_H_at(z);
        temp_spectra.insert(pair<double,Spectrum>(s._M_H, s));
    }
    msnSpectra = temp_spectra;
    return;
}

// Load specific spectrum
void Run::loadSpectrum(string path) {
    Spectrum * dta = new Spectrum(path);
    dta->filterDta(DigestOptions::windowSize, DigestOptions::peakDepth, DigestOptions::deIsotope, DigestOptions::removePrecursorNL, DigestOptions::removeWaterNLoss, DigestOptions::removePhosphoNLoss);
    if (dta->_peakVector.size() < DigestOptions::minRequiredPeaks) {
        delete dta;
        return;
    }
    msnSpectra.insert(pair<double, Spectrum > (dta->_M_H, *dta));
    //        cout << "... loaded " << msnSpectra.size() << " spectra\r" << flush; // outputting progress here somehow sometimes causes CPU inefficiency
    delete dta;
    return;
}

/* Output search results in tab-delimited format
 */
void Run::outputTab(string out_file) {
    cout << "... Writing to file " << out_file << endl; ////
    multimap<double, Spectrum *> resultMap;
    // could just sort the vector with a simple custom function
    multimap<double,Spectrum>::iterator msMmIt;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        if (((*msMmIt).second)._hits.hitMap.empty()) {
            continue;
        }
        //cout << ((*msMmIt).second)._dtaPath << endl; ////
        resultMap.insert(pair<double, Spectrum *>( ((*msMmIt).second)._hits.getTopHitScore1(), &(*msMmIt).second));
    }

    std::ofstream outFile(out_file.c_str(), ios::out);
//    multimap<double,Spectrum>::iterator msMmIt;
//    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
//        ((*msMmIt).second)._hits.output(outFile);
//    }
    outFile << SearchResult::getHeaderString() << endl;
    multimap<double, Spectrum *>::reverse_iterator RmapIt;
    for (RmapIt = resultMap.rbegin(); RmapIt != resultMap.rend(); RmapIt++) {
        RmapIt->second->_hits.output(outFile);
    }

    outFile.close();
    return;
}

void Run::deIsotope() {
    cout << "... de-isotoping spectra" << endl; ////
    int num_ms2_peaks = 0;
    int n_peaks_removed = 0;
    Isotope iso;
    multimap<double,Spectrum>::iterator msMmIt;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        num_ms2_peaks += ((*msMmIt).second)._peakVector.size();
        n_peaks_removed += ((*msMmIt).second).deIsotope(iso);
    }
    cout << n_peaks_removed << " isotopic peaks removed (" << (double) n_peaks_removed / (double) num_ms2_peaks * 100.0 << "%)" << endl; ////
    //exit(0); ////
    return;
}

void Run::filterWindowPeaks() {
    cout << "...Filtering top N peaks in m/z windows" << endl; ////
    multimap<double,Spectrum>::iterator msMmIt;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        ((*msMmIt).second).filterDta(DigestOptions::windowSize, DigestOptions::peakDepth, DigestOptions::deIsotope, DigestOptions::removePrecursorNL, DigestOptions::removeWaterNLoss, DigestOptions::removePhosphoNLoss);
    }
}

void Run::outputFilteredDtas() {
    cout << "... Writing filtered dtas" << endl; ////
    multimap<double,Spectrum>::iterator msMmIt;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        ((*msMmIt).second).printDta((*msMmIt).second._dtaPath + ".dei");
    }
    return;
}


// Correlation-based common peak filter.
void Run::filterCommonPeaks() {
    cout << "...Filtering common peaks" << endl; ////
    unsigned int min_mz = MSMass::double2int( getMinM_Z() );
    unsigned int max_mz = MSMass::double2int( getMaxM_Z() );
    vector<double> mz_bins;
    gsl_histogram * mz_hist;

    // build the histogram with appropriate bins
    mz_bins.push_back(min_mz);
    int bin_width = DigestOptions::fragmentTolInt;
    cout << "... binning between " << min_mz/MSMass::intDoubleConvertFactor << " and " << max_mz/MSMass::intDoubleConvertFactor; ////
    while (*mz_bins.rbegin() <= max_mz) {
        if (DigestOptions::fragmentTolRelative) {
            bin_width = DigestOptions::fragmentTol/1000000.0 * (*mz_bins.rbegin());
        }
        mz_bins.push_back(*mz_bins.rbegin() + bin_width);
    }
    if (!DigestOptions::fragmentTolRelative) {
        cout << " at bin width " << bin_width/MSMass::intDoubleConvertFactor << endl; ////
    } else {
        cout << " at variable bin width "
        << DigestOptions::fragmentTol/1000000.0 * (*mz_bins.begin()/MSMass::intDoubleConvertFactor) << " - "
        << DigestOptions::fragmentTol/1000000.0 * (*mz_bins.rbegin()/MSMass::intDoubleConvertFactor) << endl; ////
    }
    mz_hist = gsl_histogram_alloc (mz_bins.size() - 1);
    gsl_histogram_set_ranges (mz_hist, &mz_bins[0], mz_bins.size());
    // build the histogram with appropriate bins

    // populate mz_hist
    int total_ms2_peaks = 0;
    map<size_t,set<int> > mz_bin_scans;
    multimap<double,Spectrum>::iterator msMmIt;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        int scan_num = ((*msMmIt).second).getFirstScan();
        int num_ms2_peaks = ((*msMmIt).second)._peakVector.size();
        total_ms2_peaks += num_ms2_peaks;
        for (vector<SpectrumPeak>::iterator pvIt = ((*msMmIt).second)._peakVector.begin(); pvIt != ((*msMmIt).second)._peakVector.end(); pvIt++) {
            gsl_histogram_increment (mz_hist, (*pvIt).m_z);
            size_t bin_num;
            gsl_histogram_find (mz_hist, (*pvIt).m_z, &bin_num);
            mz_bin_scans[bin_num].insert(scan_num);
        }
    }
    cout << gsl_histogram_bins(mz_hist) << " m/z bins in histogram" << endl << endl; ////
//    gsl_histogram_fprintf (stdout, mz_hist, "%f", "%g"); ////


    //// First find bins containing more than expected counts
    double expected_prob = 1.0 / (double) mz_bin_scans.size();
    map<size_t,double> mz_bin_counts_p_val;
    for (map<size_t,set<int> >::iterator mbIt = mz_bin_scans.begin(); mbIt != mz_bin_scans.end(); mbIt++) {
        double p_val = calcBinomialCdfUpper(mz_bin_scans[mbIt->first].size(), expected_prob, total_ms2_peaks);
        mz_bin_counts_p_val[mbIt->first] = p_val;
    }

    double mz_bin_counts_q_value_cutoff = 1.0 / (double) gsl_histogram_bins(mz_hist);
    vector<size_t> fdr_passed_mz_bin_counts = JMUtil::benjaminiHochbergFDR(mz_bin_counts_p_val, mz_bin_counts_q_value_cutoff);
    cout << mz_bin_counts_p_val.size() << " non-empty m/z bins" << endl;
    cout << fdr_passed_mz_bin_counts.size() << " m/z bins passing fdr cutoff at " << mz_bin_counts_q_value_cutoff << endl;
    //// First find bins containing more than expected counts

    /* Should ideally fit mz_hist (non-zero values) to weibull or exp and then get p-values from that.
     * For now just go with a heuristic with this setup, if the fragment ion (and mz bin) tolerance
     * is set low, might not work
     */
    int min_bin_cnt = msnSpectra.size() * 0.25; ////

//    cout << "binNum\tbinCnt\tbinLo\tbinHi\tpValueBinom" << endl; ////
    size_t prev_bin = -1; ////
    map<size_t,size_t> merge_bins;
    for (vector<size_t>::iterator vecIt = fdr_passed_mz_bin_counts.begin(); vecIt != fdr_passed_mz_bin_counts.end(); vecIt++) {
       double lower, upper; ////
        gsl_histogram_get_range (mz_hist, *vecIt, &lower, &upper); ////
        int cnt = gsl_histogram_get(mz_hist, *vecIt); ////
//        cout << *vecIt << "\t" << cnt << "\t" << lower/MSMass::intDoubleConvertFactor << "\t" << upper/MSMass::intDoubleConvertFactor << "\t" << mz_bin_counts_p_val[*vecIt] << endl; ////
        if (*vecIt == prev_bin + 1 && cnt >= min_bin_cnt ||
            *vecIt == prev_bin + 1 && mz_bin_scans[prev_bin].size() >= min_bin_cnt) {
//            cout << "merge" << endl; ////
            merge_bins[*vecIt] = prev_bin;
        }
        prev_bin = *vecIt;
    }
////    exit(0); ////

    //// Merge adjacent bins
    vector<double> mz_bins_merged;
    for (int i = mz_bins.size() - 1; i >= 0; i--) {
        if (merge_bins.count(i) == 0) {
            mz_bins_merged.push_back( mz_bins[i] );
        }
    }
    reverse(mz_bins_merged.begin(), mz_bins_merged.end());
    cout << mz_bins_merged.size() << " bins after merging" << endl;
    //// Merge adjacent bins

    //// populate mz_hist again with merged bins
    gsl_histogram_free(mz_hist);
    mz_hist = gsl_histogram_alloc (mz_bins_merged.size() - 1);
    gsl_histogram_set_ranges (mz_hist, &mz_bins_merged[0], mz_bins_merged.size());
    vector<int> scan_nums;
    map<int,int> scan_2_ms2_peak_nums;
    mz_bin_scans.erase(mz_bin_scans.begin(), mz_bin_scans.end());
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        int scan_num = ((*msMmIt).second).getFirstScan();
        scan_nums.push_back(scan_num);
        int num_ms2_peaks = ((*msMmIt).second)._peakVector.size();
        scan_2_ms2_peak_nums[scan_num] = num_ms2_peaks;
        for (vector<SpectrumPeak>::iterator pvIt = ((*msMmIt).second)._peakVector.begin(); pvIt != ((*msMmIt).second)._peakVector.end(); pvIt++) {
            gsl_histogram_increment (mz_hist, (*pvIt).m_z);
            size_t bin_num;
            gsl_histogram_find (mz_hist, (*pvIt).m_z, &bin_num);
            mz_bin_scans[bin_num].insert(scan_num);
        }
    }
    cout << gsl_histogram_bins(mz_hist) << " m/z bins in merged histogram" << endl << endl; ////
//    gsl_histogram_fprintf (stdout, mz_hist, "%f", "%g"); ////
//    exit(0);


    sort(scan_nums.begin(), scan_nums.end());
    pair<double,double> bins_width = JMUtil::numBinsFreedmanDiaconis(scan_nums);
    int num_scan_bins = bins_width.first;
    double scan_binwidth = bins_width.second;

    if (scan_binwidth < 1.0 || scan_nums.size() < num_scan_bins*5) {
        cerr << "Not enough scans for correlation-based filtering" << endl;
        return;
    }
    
    cout << "... binning ms2 scan numbers between " << *scan_nums.begin() << " and " << *scan_nums.rbegin() << " at bin width " << scan_binwidth << endl; ////
    gsl_histogram * all_scan_hist = gsl_histogram_alloc (num_scan_bins);
    gsl_histogram_set_ranges_uniform (all_scan_hist, *scan_nums.begin(), *scan_nums.rbegin() + 1);
    for (vector<int>::iterator vecIt = scan_nums.begin(); vecIt != scan_nums.end(); vecIt++) {
        gsl_histogram_increment(all_scan_hist, *vecIt);
    }
    cout << gsl_histogram_bins(all_scan_hist) << " ms2 scan bins in histogram" << endl; ////
    double all_scans_bin_count[num_scan_bins]; ////
    double all_scans_bin_peak_count[num_scan_bins]; // all scan ms2 peak count histogram to run correlations against ////
    fill_n(all_scans_bin_peak_count, num_scan_bins, 0);
    for (size_t i = 0; i < gsl_histogram_bins(all_scan_hist); i++) {
        all_scans_bin_count[i] = gsl_histogram_get(all_scan_hist, i);
    }
    for (size_t i = 0; i < scan_nums.size(); i++) {
        size_t bin_num;
        gsl_histogram_find (all_scan_hist, scan_nums[i], &bin_num);
        all_scans_bin_peak_count[bin_num] += scan_2_ms2_peak_nums[ scan_nums[i] ];
    }
//    copy(all_scan_bin_count.begin(), all_scan_bin_count.end(), ostream_iterator<int>(cout, " ")); cout << endl; ////

    gsl_histogram * scan_hist = gsl_histogram_clone (all_scan_hist);

    map<size_t,double> mz_bin_p_val;
    map<size_t,double> mz_bin_corr_coef;
    for (map<size_t,set<int> >::iterator mbIt = mz_bin_scans.begin(); mbIt != mz_bin_scans.end(); mbIt++) {
        if (mbIt->second.size() < min_bin_cnt) ////
            continue; ////

        gsl_histogram_reset(scan_hist);

        for (set<int>::iterator setIt = mbIt->second.begin(); setIt != mbIt->second.end(); setIt++) {
            gsl_histogram_increment(scan_hist, *setIt);
        }
        double scans_bin_count[num_scan_bins]; // scan histogram for each bin to correlate against all_scans_bin_count
        for (size_t i = 0; i < gsl_histogram_bins(scan_hist); i++) {
            scans_bin_count[i] = gsl_histogram_get(scan_hist, i);
        }
//        copy(scan_bin_count.begin(), scan_bin_count.end(), ostream_iterator<int>(cout, " ")); cout << endl; ////

        double scans_bin_peak_count[num_scan_bins]; // all scan ms2 peak count histogram to run correlations against
        fill_n(scans_bin_peak_count, num_scan_bins, 0);
        for (set<int>::iterator setIt = mbIt->second.begin(); setIt != mbIt->second.end(); setIt++) {
            size_t bin_num;
            gsl_histogram_find(scan_hist, *setIt, &bin_num);
            scans_bin_peak_count[bin_num] += scan_2_ms2_peak_nums[ *setIt ];
        }

        double r = gsl_stats_correlation (all_scans_bin_peak_count, 1, scans_bin_peak_count, 1, num_scan_bins);
        double t_stat = r * sqrt((num_scan_bins - 2) / (1 - pow(r,2)));
        //double p_val = 2 * gsl_cdf_tdist_P(-fabs(t_stat), num_scan_bins - 2); // 2-tailed
        double p_val = 2 * calcTdistCdfLower(-fabs(t_stat), num_scan_bins - 2); // 2-tailed       
        mz_bin_corr_coef[mbIt->first] = r;
        mz_bin_p_val[mbIt->first] = p_val;
    }

    double q_value_cutoff = 1.0 / (double) gsl_histogram_bins(mz_hist);
    vector<size_t> fdr_passed = JMUtil::benjaminiHochbergFDR(mz_bin_p_val, q_value_cutoff);
    cout << mz_bin_p_val.size() << " non-empty m/z bins" << endl;
    cout << fdr_passed.size() << " m/z bins passing fdr cutoff at " << q_value_cutoff << endl;

//    vector< pair<unsigned int,unsigned int> > fdr_passed_bin_ranges;
    fdr_passed_bin_ranges.clear();
    
    if (DigestOptions::verbose) {
        cout << "binNum\tbinCnt\tbinLo\tbinHi\tCorrCoef\tpValue" << endl; ////
    }
    for (vector<size_t>::iterator vecIt = fdr_passed.begin(); vecIt != fdr_passed.end(); vecIt++) {
        double lower, upper; ////
        gsl_histogram_get_range (mz_hist, *vecIt, &lower, &upper); ////
        fdr_passed_bin_ranges.insert(Range(lower,upper));
        int cnt = gsl_histogram_get(mz_hist, *vecIt); ////
        if (DigestOptions::verbose) {
            cout << *vecIt << "\t" << cnt << "\t" << lower/MSMass::intDoubleConvertFactor << "\t" << upper/MSMass::intDoubleConvertFactor << "\t" << mz_bin_corr_coef[*vecIt] << "\t" << mz_bin_p_val[*vecIt] << endl; ////
        }
    }
    //exit(0);
    
    // do the actual filtering
    int n_peaks_removed = 0;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
//        cout << ((*msMmIt).second)._M_H << "\t" << ((*msMmIt).second)._peakVector.size() << "\t" << ((*msMmIt).second)._dtaPath << endl; ////
        n_peaks_removed += ((*msMmIt).second).filterRanges(fdr_passed_bin_ranges);
    }
    cout << n_peaks_removed << " peaks filtered from " << fdr_passed_bin_ranges.size() << " m/z bins (" << (double) n_peaks_removed / (double) total_ms2_peaks * 100.0 << "%)" << endl; ////
//    gsl_histogram_fprintf (stdout, mz_hist, "%g", "%g");
    gsl_histogram_free (mz_hist);
    gsl_histogram_free (all_scan_hist);
    gsl_histogram_free (scan_hist);
    return;
}

bool Run::isInFilteredRange(unsigned int m_z) {
    if (fdr_passed_bin_ranges.empty()) {
        return false;
    }
    return inRange(fdr_passed_bin_ranges, m_z);
}

double Run::getMinM_Z() {
    double min_mz = MSMass::numericMaxLimit;
    double temp;
    multimap<double,Spectrum>::iterator msMmIt;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        temp = ((*msMmIt).second).getMinMZ();
        if (temp < min_mz) {
            min_mz = temp;
        }
    }
    return min_mz;
}

double Run::getMaxM_Z() {
    double max_mz = 0.0;
    double temp;
    multimap<double,Spectrum>::iterator msMmIt;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        temp = ((*msMmIt).second).getMaxMZ();
        if (temp > max_mz) {
            max_mz = temp;
        }
    }
    return max_mz;
}

double Run::getNumTotalPSMs() {
    double c = 0;
    multimap<double,Spectrum>::iterator msMmIt;
    for (msMmIt = msnSpectra.begin(); msMmIt != msnSpectra.end(); msMmIt++) {
        c += msMmIt->second._numPSMs;
    }
    return c;
}

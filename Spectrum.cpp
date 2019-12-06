#include "Spectrum.hpp"

// Read dta file
void Spectrum::readDta (string file) {
  ifstream dtafile (file.c_str());
  if (! dtafile.is_open()) {
    cout << "Error opening dta file " << file << endl;
    exit (1);
  }
  string line;
  int line_count = 0;
  while (! dtafile.eof() ) {
    getline (dtafile, line);
    if (line == "")
      break;
    string buf; // Have a buffer string
    stringstream ss(line); // Insert the string into a stream
    vector<string> tokens; // vector to hold tokens
    while (ss >> buf) {
      tokens.push_back(buf);
    }
    if (tokens.size() > 3) {
      cout << "Error: " << file << " does not appear to be a dta file" << endl;
      exit(0);
    }
    if (line_count == 0) {
      _M_H = atof(tokens[0].c_str());
      _charge = atoi(tokens[1].c_str());
    } else {
        double m_z = atof(tokens[0].c_str());
        if (m_z > MSMass::numericMaxLimit) {
            //cerr << "WARNING: skipping " << m_z << " m/z peak from " << file << endl;
        } else if (m_z < DigestOptions::filterLowMzMass) {
            //cerr << "WARNING: skipping " << m_z << " m/z peak from " << file << endl; //
        } else {
            double ity = atof(tokens[1].c_str());
            SpectrumPeak sp(MSMass::double2int( m_z ), ity, 0, false);
            _peakVector.push_back(sp);
        }
    }
    line_count ++;
  }
  dtafile.close();

  _topPeaks = -1;
//  enforcePeakFilter();

  _dtaPath = file;
  _fullPeakCount = _peakVector.size();
  setFirstScan();
  _numPSMs = 0;
  return;
}

// Filter the dta to top num_peaks.  This is pretty much a translation of sean's trim_dta_opt php version.
void Spectrum::filterDta(int window, int num_peaks, bool flg_deisotope, bool flg_nl, bool flg_nl_water, bool flg_nl_phospho) {
    flipPeakFilterFlagsFalse();
    int bin_number = 0;
    map<int, int> window_count;

    unsigned int m_z = MSMass::double2int( getM_Z() );
    int water = MSMass::double2int(18.01056); //mass of water
    int phosphate = MSMass::double2int(97.97632); // mass of phosphate
    int nl_water = m_z - ( water / _charge );
    int nl_2waters = m_z - (water * 2 / _charge );
    int phospho_nl = m_z - (phosphate / _charge);
    bool flg_nl_phospho_water = flg_nl_phospho && flg_nl_water;
    int nl_phospho_water = phospho_nl - nl_water;
    int nl_phospho_2waters = phospho_nl - nl_2waters;

    set<unsigned int> filtered10;
    filtered10.clear();

    sort(_peakVector.begin(), _peakVector.end(), sortByIntensity());

    for (vector<SpectrumPeak>::iterator pvIt = _peakVector.begin(); pvIt != _peakVector.end(); pvIt++) {
        unsigned int ion_mz = (*pvIt).m_z;
//        double ion_intensity = (*pvIt).intensity;
        bin_number = ion_mz / window / MSMass::intDoubleConvertFactor;
        if (window_count[bin_number] >= num_peaks) {
            continue;
        } else {
            int mass_diff = MSMass::fragmentMass2Diff(ion_mz);
            // check for isotopes
            if (flg_deisotope) {
                set<unsigned int>::iterator left = filtered10.lower_bound(ion_mz - mass_diff);
                set<unsigned int>::iterator right = filtered10.upper_bound(ion_mz + mass_diff);
                if (left != right) {
                    // cout << "... removing Isotope " << ion_mz << " " << ion_intensity << endl;
                    continue;
                }
            }

            // check for precursor
            if (flg_nl) {
                if (ion_mz >= (m_z - mass_diff) && ion_mz <= (m_z + mass_diff)) {
                    // cout << "... removing precursor " << ion_mz << " " << ion_intensity << endl;
                    continue;
                }
            }

            // check for waters
            if (flg_nl_water) {
                if (ion_mz >= (nl_water - mass_diff) && ion_mz <= (nl_water + mass_diff)) {
                    // cout << "... removing Water " << ion_mz << " " << ion_intensity << endl;
                    continue;
                }
                if (ion_mz >= (nl_2waters - mass_diff) && ion_mz <= (nl_2waters + mass_diff)) {
                    // cout << "... removing Water2 " << ion_mz << " " << ion_intensity << endl;
                    continue;
                }
            }

            if (flg_nl_phospho) {
                if (ion_mz >= (phospho_nl - mass_diff) && ion_mz <= (phospho_nl + mass_diff)) {
                    // cout << "... removing Phos NL " << ion_mz << " " << ion_intensity << endl;
                    continue;
                }
            }
            if (flg_nl_phospho_water) {
                if (ion_mz >= (nl_phospho_water - mass_diff) && ion_mz <= (nl_phospho_water + mass_diff)) {
                    // cout << "... removing water_nl " << ion_mz << " " << ion_intensity << endl;
                    continue;
                }
                if (ion_mz >= (nl_phospho_2waters - mass_diff) && ion_mz <= (nl_phospho_2waters + mass_diff)) {
                    // cout << "... removing water2_nl " << ion_mz << " " << ion_intensity << endl;
                    continue;
                }
            }
            
        }

        (*pvIt).rank = window_count[bin_number] + 1;
        (*pvIt).passedFilter = true;
        filtered10.insert(ion_mz); ////
        window_count[bin_number]++;
    }

    _topPeaks = num_peaks;
    sort(_peakVector.begin(), _peakVector.end()); // re-sort by m_z again
    enforcePeakFilter();
//    printDta(_dtaPath + ".out"); ////
//    printFilteredRankedDta(_dtaPath + ".filtered"); ////
    return;
}

// remember to flip passedFilter bool flags on/off depending on usage
int Spectrum::enforcePeakFilter() {
    int n_peaks_start = _peakVector.size();
    _peakVector.erase( std::remove_if(_peakVector.begin(), _peakVector.end(), [](const SpectrumPeak & s) { return !s.passedFilter; }), _peakVector.end() );
    int n_peaks_end = _peakVector.size();
    return (n_peaks_start - n_peaks_end);
}

void Spectrum::flipPeakFilterFlagsTrue() {
    for (vector<SpectrumPeak>::iterator vec_it = _peakVector.begin(); vec_it != _peakVector.end(); vec_it++) {
        vec_it->passedFilter = true;
    }
}

void Spectrum::flipPeakFilterFlagsFalse() {
    for (vector<SpectrumPeak>::iterator vec_it = _peakVector.begin(); vec_it != _peakVector.end(); vec_it++) {
        vec_it->passedFilter = false;
    }
}

int Spectrum::filterRanges(set<Range,RangeCompare> fdr_passed_bin_ranges) {
    flipPeakFilterFlagsTrue();
    vector<SpectrumPeak>::iterator pv_it = _peakVector.begin();
    int n_peaks_to_remove = 0;
    for (set< pair<unsigned int,unsigned int> >::iterator set_it = fdr_passed_bin_ranges.begin(); set_it != fdr_passed_bin_ranges.end(); set_it++) {
//        cout << "... removing from range\t" << (double) set_it->first/MSMass::intDoubleConvertFactor << "\t" << (double) set_it->second/MSMass::intDoubleConvertFactor << endl; ////
        vector<SpectrumPeak>::iterator pv_it_lo = lower_bound(pv_it, _peakVector.end(), set_it->first);
        if (pv_it_lo == _peakVector.end())
            break;
        while (pv_it_lo != _peakVector.end() && pv_it_lo->m_z < set_it->second) {
            pv_it_lo->passedFilter = false;
            n_peaks_to_remove++;
//            cout << "\t" << (double) pv_it_lo->m_z/MSMass::intDoubleConvertFactor << endl; ////
            pv_it_lo++;
        }
        pv_it = pv_it_lo;
    }
    enforcePeakFilter();
    flipPeakFilterFlagsFalse();
    return n_peaks_to_remove;
}

/**
 * Get (filtered) sum of intensities, i.e. peaks that passed filters, de-isotoping etc
 * @return sum_ity
 */
double Spectrum::getFilteredSumIty(double min_ity, double min_m_z) {
    double sum_ity = 0.0;
    for (vector<SpectrumPeak>::iterator vec_it = _peakVector.begin(); vec_it != _peakVector.end(); vec_it++) {
        if (vec_it->passedFilter && vec_it->intensity > min_ity && vec_it->m_z > min_m_z) {
            sum_ity += (vec_it->intensity);
        }
    }
    return sum_ity;
}

double Spectrum::getMinMZ() {
    return MSMass::int2double( (_peakVector.begin()->m_z - MSMass::fragmentMass2Diff(_peakVector.begin()->m_z)) );
}

double Spectrum::getMaxMZ() {
    return MSMass::int2double( (_peakVector.rbegin()->m_z + MSMass::fragmentMass2Diff(_peakVector.rbegin()->m_z)) );
}

// calculate precursor m/z
double Spectrum::getM_Z() {
    return MSMass::m_h2m_z(_M_H, _charge);
}

// Use this to find cases of incorrect charge assignment by Thermo
void Spectrum::adjustM_H_at(int z) {
    double m_z = getM_Z();
    double m_h = m_z * z - z * MSMass::protonMass + MSMass::protonMass;
    _M_H = m_h;
    _charge = z;
    return;
}

void Spectrum::setFirstScan() {
    boost::filesystem::path pathname(_dtaPath);
    string basename = pathname.filename().string();
    vector<string> tokens;
    boost::split(tokens, basename, boost::is_any_of("."));
    _firstScanNum = atoi( tokens[1].c_str() );
}

int Spectrum::getFirstScan() {
    return _firstScanNum;
}

// Print dta file
void Spectrum::printDta(string outfile) {
    std::ofstream outFile(outfile.c_str(), ios::out);
    outFile << fixed << setprecision(5) << _M_H << "\t" << _charge << endl;
    // outFile << m_z_filtered10.size() << " " << Int_filtered10.size() << " " << rank_filtered10.size() << endl;
    for (vector<SpectrumPeak>::iterator pvIt = _peakVector.begin(); pvIt != _peakVector.end(); pvIt++) {
        outFile << fixed << setprecision(5) << MSMass::int2double( (*pvIt).m_z ) << "\t";
        outFile << fixed << setprecision(2) << (*pvIt).intensity << "\t";
        outFile << fixed << (*pvIt).charge << endl;
    }
    outFile.close();
    return;
}

// Print dta file
void Spectrum::printFilteredRankedDta(string outfile) {
    std::ofstream outFile(outfile.c_str(), ios::out);
    outFile << fixed << setprecision(5) << _M_H << "\t" << _charge << endl;
    for (vector<SpectrumPeak>::iterator pvIt = _peakVector.begin(); pvIt != _peakVector.end(); pvIt++) {
        if ((*pvIt).passedFilter) {
            outFile << fixed << setprecision(5) << MSMass::int2double( (*pvIt).m_z ) << "\t";
            outFile << fixed << setprecision(2) << (*pvIt).intensity << "\t";
            outFile << (*pvIt).rank << endl;
        }
    }
    outFile.close();
    return;
}

int Spectrum::deIsotope(Isotope &iso) {
    vector<IsotopicEnvelopeMatch> envelope_vec;
    vector<SpectrumPeak>::iterator pvIt;
    for (pvIt = _peakVector.begin(); pvIt != _peakVector.end(); pvIt++) {
        int start_m_z = (*pvIt).m_z;
        for (int z = 1; z <= _charge; z++) {
            //TheoreticalEnvelope t_e = iso.getTheoreticalEnvelope(MSMass::int2double(start_m_z), z);
           TheoreticalEnvelope t_e = iso.getHashedTheoreticalEnvelope(MSMass::int2double(start_m_z), z);
            //t_e.print();
            IsotopicPeakMatch ipm(*pvIt, t_e._peakVector[0]);
            IsotopicEnvelopeMatch env_match;
            env_match._pmVec.push_back(ipm);
            iso.buildEnvelopeRecursive(t_e, env_match, *this, envelope_vec);
        }
    }
//    cout << iso.envelope_vec.size() << " envelopes" << endl;

    // retain the highest-scoring charges state for each envelope (mono-iso peak)
    set<IsotopicEnvelopeMatch> e_m_set;
    for (vector<IsotopicEnvelopeMatch>::iterator em_vec_it = envelope_vec.begin(); em_vec_it != envelope_vec.end(); em_vec_it++) {
        //em_vec_it->print(); ////
        pair< set<IsotopicEnvelopeMatch>::iterator, bool> ret_pair = e_m_set.insert(*em_vec_it);
        if (ret_pair.second == false) {
            if (ret_pair.first->_pevznerScore < em_vec_it->_pevznerScore) {
                e_m_set.erase(ret_pair.first);
                e_m_set.insert(*em_vec_it);
            }
        }
    }
//    cout << e_m_set.size() << " envelopes" << endl;

    // fold sub-envelopes into higher-scoring ones to the left
    for (set<IsotopicEnvelopeMatch>::iterator set_it = e_m_set.begin(); set_it != e_m_set.end(); set_it++) {
        //set_it->print(); ////
        vector<IsotopicPeakMatch>::const_iterator vec_it = set_it->_pmVec.begin();
        vec_it++;
        for (; vec_it != set_it->_pmVec.end(); vec_it++) {
            set<IsotopicEnvelopeMatch>::const_iterator temp_set_it = find_if(e_m_set.begin(), e_m_set.end(), findIEM(vec_it->_sp->m_z));
            if (temp_set_it != e_m_set.end() && temp_set_it->_pmVec[0]._tp.charge == set_it->_pmVec[0]._tp.charge) {
                if (temp_set_it->_pevznerScore < set_it->_pevznerScore) {
                    e_m_set.erase(temp_set_it);
                }
            }
        }
    }
//    cout << e_m_set.size() << " envelopes" << endl;

//    // pick the best scoring option of what remains - ensure one envelope per peak
    map<unsigned int, set<IsotopicEnvelopeMatch,compareScore> > peak_2_envelope_map;
    for (set<IsotopicEnvelopeMatch>::iterator set_it = e_m_set.begin(); set_it != e_m_set.end(); set_it++) {
        //set_it->print(); ////
        for (vector<IsotopicPeakMatch>::const_iterator vec_it = set_it->_pmVec.begin(); vec_it != set_it->_pmVec.end(); vec_it++) {
            if (peak_2_envelope_map.count(vec_it->_sp->m_z) > 0) {
                peak_2_envelope_map[ vec_it->_sp->m_z ].insert(*set_it);
            } else {
                set<IsotopicEnvelopeMatch,compareScore> overlap_set;
                overlap_set.insert(*set_it);
                peak_2_envelope_map.insert( pair<unsigned int, set<IsotopicEnvelopeMatch,compareScore> >(vec_it->_sp->m_z, overlap_set) );
            }
        }
    }
    for (map<unsigned int, set<IsotopicEnvelopeMatch,compareScore> >::iterator map_it = peak_2_envelope_map.begin(); map_it != peak_2_envelope_map.end(); map_it++) {
        if (map_it->second.size() > 1) {
            set<IsotopicEnvelopeMatch,compareScore>::iterator temp_set_it = map_it->second.begin();
            set<IsotopicEnvelopeMatch>::const_iterator set_it_to_erase;
            do { // go past envelopes already erased to top-scoring existing one
                set_it_to_erase = find_if(e_m_set.begin(), e_m_set.end(), findIEM(temp_set_it->_pmVec[0]._sp->m_z));
                temp_set_it++;
            } while (set_it_to_erase == e_m_set.end() && temp_set_it != map_it->second.end());
            for (; temp_set_it != map_it->second.end(); temp_set_it++) {
                set_it_to_erase = find_if(e_m_set.begin(), e_m_set.end(), findIEM(temp_set_it->_pmVec[0]._sp->m_z));
                if (set_it_to_erase != e_m_set.end()) {
                    e_m_set.erase(set_it_to_erase);
                }
            }
        }
    }

    // final (arbitrary) filter
    for (set<IsotopicEnvelopeMatch>::iterator set_it = e_m_set.begin(); set_it != e_m_set.end();) {
        if (set_it->_myNormalizedPpmScore < 0.6 || set_it->_myNormalizedItyScore < 0.6) {
            e_m_set.erase(set_it++);
        } else {
            ++set_it;
        }
    }
//    cout << e_m_set.size() << " envelopes" << endl;

    // remove/combine envelope peaks
    flipPeakFilterFlagsTrue();
    for (set<IsotopicEnvelopeMatch>::iterator set_it = e_m_set.begin(); set_it != e_m_set.end(); set_it++) {
        //set_it->print(); ////
        vector<SpectrumPeak>::iterator mono_iso_pv_it = find_if(_peakVector.begin(), _peakVector.end(), findSpectrumPeak(set_it->_pmVec[0]._sp->m_z));
        mono_iso_pv_it->charge = set_it->_pmVec[0]._tp.charge;
        for (int i = 1; i < set_it->_pmVec.size(); i++) {
            vector<SpectrumPeak>::iterator pv_it = find_if(_peakVector.begin(), _peakVector.end(), findSpectrumPeak(set_it->_pmVec[i]._sp->m_z));
            mono_iso_pv_it->intensity += pv_it->intensity;
            //cout << "... removing " << fixed << setprecision(6) << MSMass::int2double(pv_it->m_z) << "\t";
            //cout << fixed << setprecision(2) << pv_it->intensity << endl; ////
            pv_it->passedFilter = false;
        }
        //cout << endl; ////
    }
    int n_peaks_removed = enforcePeakFilter();
    return n_peaks_removed;
}


#include <vector>

#include "Isotope.hpp"
#include "Spectrum.hpp"

const string Isotope::isotopefn = "../MSLib/Isotope/ISOTOPE.DAT";

Isotope::Isotope() {
    if (!init_data(isotopefn, sad, em)) {
        string tmp(isotopefn);
        wcerr << "Could not open isotope file: "
                << wstring(tmp.begin(), tmp.end()) << endl;
        exit(0);
    }
    
    limit = 0;
    digits = 2;

    // default empty theoretical envelope
    TheoreticalEnvelope t_e;

    //Parser p_tmp(em);
    //p = &p_tmp;
    //Parser p(em);
}

IsotopicPeakMatch::IsotopicPeakMatch(const SpectrumPeak &sp, TheoreticalPeak &tp) {
    _sp = &sp;
    _tp = tp;
    calcPeakDeviation();
}

void IsotopicPeakMatch::calcPeakDeviation() {
    _peakDeviation = _tp.getPeakDeviation(*_sp);
}

Pattern Isotope::getIsotopePattern(string mol_formula, long charge) {
    wstring formula(mol_formula.length(), ' ');
    copy(mol_formula.begin(), mol_formula.end(), formula.begin());

    Parser p(em);

    // read a formula
    fm.clear();
    try {
        p.write_compound(formula, fm);
    } catch (Parser::Error e) {
        wcerr << e._message;
        exit(0); // skip to next compound
    }

    // initialize the result
    result.resize(1);
    result.front().mass = 0.0;
    result.front().rel_area = 1.0;

    calculate(tmp, result, fm, limit, charge, sad);
//    print_pattern(result, digits);
    return result;
}

void Isotope::buildEnvelopeRecursive(TheoreticalEnvelope &theo_e, IsotopicEnvelopeMatch &env_match, const Spectrum &dta, vector<IsotopicEnvelopeMatch> &envelope_vec) {
    int i = env_match._pmVec.size();
    //int mass_diff = MSMass::double2int(0.02);
    int mass_diff = MSMass::fragmentMass2Diff( theo_e._peakVector[0].m_z );
    set<const SpectrumPeak*,SpectrumPeakPtrCompare> idx;
    if (i < theo_e._peakVector.size()) {
        vector<SpectrumPeak>::const_iterator vec_it_lo = upper_bound(dta._peakVector.begin(), dta._peakVector.end(), theo_e._peakVector[i].m_z - mass_diff, SpectrumPeakCompare());
        vector<SpectrumPeak>::const_iterator vec_it_hi = lower_bound(dta._peakVector.begin(), dta._peakVector.end(), theo_e._peakVector[i].m_z + mass_diff, SpectrumPeakCompare());
        for (vector<SpectrumPeak>::const_iterator vec_it = vec_it_lo; vec_it != vec_it_hi; vec_it++) {
            idx.insert(&(*vec_it));
        }
    }

    // make sure the envelope doesn't have any really weird minima
    if (env_match._pmVec.size() >= 2) {
        double max_ity = 0.0;
        double max_rise_ity = 0.0;
        max_ity = env_match._pmVec[0]._sp->intensity;
        for (i = 1; i < env_match._pmVec.size(); i++) {
            if (env_match._pmVec[i]._sp->intensity > max_ity) {
                max_ity = env_match._pmVec[i]._sp->intensity;
            }
            double diff = env_match._pmVec[i]._sp->intensity - env_match._pmVec[i-1]._sp->intensity;
            if (diff > 0 && diff > max_rise_ity) {
                max_rise_ity = diff;
            }
        }

        size_t last_idx = env_match._pmVec.size() - 1;
        for (set<const SpectrumPeak*,SpectrumPeakPtrCompare>::iterator set_it = idx.begin(); set_it != idx.end();) {
            double diff = (*set_it)->intensity - env_match._pmVec[last_idx]._sp->intensity;
            if (env_match._pmVec[last_idx]._sp->intensity < max_ity && diff > 0) { // allow no minima
                idx.erase(set_it++);
            } else {
                ++set_it;
            }
        }
    }

    
    if (idx.empty()) {
        if (env_match._pmVec.size() > 1) {
            env_match.calcScores();
            // store envelope
            envelope_vec.push_back(env_match);
            return;
        }
    } else {
        for (set<const SpectrumPeak*,SpectrumPeakPtrCompare>::iterator set_it = idx.begin(); set_it != idx.end(); set_it++) {
            const SpectrumPeak * spPtr = *set_it;
            IsotopicPeakMatch ipm(*spPtr, theo_e._peakVector[i]);
            IsotopicEnvelopeMatch next_env_match = env_match;
            next_env_match._pmVec.push_back(ipm);
            buildEnvelopeRecursive(theo_e, next_env_match, dta, envelope_vec);
        }
    }
    return;
}


// retrieve pre-computed TheoreticalEnvelope or calculate it and add to hash
// store envelopes at whatever charge state comes first and then recalc to work around
// the UINT_MAX limitations on m_z
TheoreticalEnvelope Isotope::getHashedTheoreticalEnvelope(double start_m_z, long charge) {
    TheoreticalEnvelope t_e;
    double m_h = MSMass::m_z2m_h(start_m_z, charge);
    double m_h_rounded = JMUtil::roundTo(m_h, 0);
    unsigned long long int hash_key = (unsigned long long int) m_h_rounded;
    if (theoreticEnvelopeHash.count(hash_key)) {
        t_e = theoreticEnvelopeHash[hash_key];
    } else {
        t_e = getTheoreticalEnvelope(MSMass::m_h2m_z(m_h_rounded, charge), charge);
        theoreticEnvelopeHash[hash_key] = t_e;
    }
    // adjust charge state
    for (vector<TheoreticalPeak>::iterator vec_it = t_e._peakVector.begin(); vec_it != t_e._peakVector.end(); vec_it++) {
        vec_it->m_z = MSMass::m_h2m_z( MSMass::m_z2m_h(vec_it->m_z, vec_it->charge), charge);
        vec_it->charge = charge;
    }
    // adjust m/z
    unsigned int mass_diff = MSMass::double2int(start_m_z) - t_e._peakVector[0].m_z;
    for (vector<TheoreticalPeak>::iterator vec_it = t_e._peakVector.begin(); vec_it != t_e._peakVector.end(); vec_it++) {
        vec_it->m_z += mass_diff;
    }
    return t_e;
}

TheoreticalEnvelope Isotope::getTheoreticalEnvelope(double start_m_z, long charge) {
    TheoreticalEnvelope t_e;
    string mol_formula = MSMass::averagineModel( MSMass::m_z2m_h(start_m_z, charge) );
    Pattern result = getIsotopePattern(mol_formula, charge);
    // find the maximum
    double max_area = 0;
    double sum_area = 0;
    for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    if(max_area < (*i).rel_area)
        max_area = (*i).rel_area;
        sum_area += (*i).rel_area;
    }
    if(max_area == 0)
        return t_e; // empty pattern

    double print_limit = pow(10.0, -digits); //// / 2;

    Pattern::iterator i = result.begin();
    double mz_diff = start_m_z - (*i).mass;
    for (; i != result.end(); ++i) {
        double mass = (*i).mass + mz_diff;
        double rel_area = (*i).rel_area;
        double val_perc = rel_area / max_area * 100.0;
        //double val_norm = rel_area / sum_area;
        if (mass > DBL_MIN && val_perc >= print_limit) {
            TheoreticalPeak p(MSMass::double2int( mass ), val_perc, charge, 'i', 0);
            t_e._peakVector.push_back(p);
        }
    }

    return t_e;
}

void TheoreticalEnvelope::print() {
    cout << (*_peakVector.begin()).charge << "+" << endl;
    for (vector<TheoreticalPeak>::iterator pvIt = _peakVector.begin(); pvIt != _peakVector.end(); pvIt++) {
        cout << fixed << setprecision(6) << MSMass::int2double( (*pvIt).m_z ) << "\t";
        cout << fixed << setprecision(6) << (*pvIt).intensity << endl;
    }
}

void IsotopicEnvelopeMatch::calcScores() {
    adjustIntensityScale();
    calcIntensityScores();
    calcPpmScores();
    _pevznerScore = 0.0;
    _myNormalizedPpmScore = 0.0;
    _myNormalizedItyScore = 0.0;
    for (vector<IsotopicPeakMatch>::iterator vec_it = _pmVec.begin(); vec_it != _pmVec.end(); vec_it++) {
        vec_it->_pevznerScore = sqrt(vec_it->_scaledTheoreticIntensity) * vec_it->_ppmScore * vec_it->_ityScore;
        _pevznerScore += vec_it->_pevznerScore;
        _myNormalizedPpmScore += vec_it->_ppmScore;
        _myNormalizedItyScore += vec_it->_ityScore;
    }
    _myNormalizedPpmScore /= ((double) _pmVec.size());
    _myNormalizedItyScore /= ((double) _pmVec.size());
    _myNormalizedPevznerScore = (_myNormalizedPpmScore + _myNormalizedItyScore)/2.0;
    return;
}

// like Pevzner's, but shift the observed m/z by mean_diff to prevent skews (1st peak always a perfect match)
void IsotopicEnvelopeMatch::calcPpmScores() {
    //double mz_tol = 0.02; ////
    int mass_diff = MSMass::fragmentMass2Diff( _pmVec[0]._tp.m_z );
    double mz_tol = MSMass::int2double(mass_diff);
    double mean_diff = 0.0;
    for (vector<IsotopicPeakMatch>::iterator vec_it = _pmVec.begin(); vec_it != _pmVec.end(); vec_it++) {
        mean_diff += ( MSMass::int2double(vec_it->_sp->m_z) - MSMass::int2double(vec_it->_tp.m_z) );
    }
    mean_diff /= (double) _pmVec.size();

    for (vector<IsotopicPeakMatch>::iterator vec_it = _pmVec.begin(); vec_it != _pmVec.end(); vec_it++) {
        vec_it->_ppmScore = 1.0 - fabs(MSMass::int2double(vec_it->_sp->m_z) - mean_diff - MSMass::int2double(vec_it->_tp.m_z))/mz_tol;
    }
    return;
}

// scale theoretic intensities for subsequent intensity scoring
// This differs a bit from the Pevzner approach, but seems to be much more straight-forward
void IsotopicEnvelopeMatch::adjustIntensityScale() {
    double sum_observed_ity = 0.0;
    double sum_theoretic_ity = 0.0;
    for (vector<IsotopicPeakMatch>::iterator vec_it = _pmVec.begin(); vec_it != _pmVec.end(); vec_it++) {
        sum_observed_ity += vec_it->_sp->intensity;
        sum_theoretic_ity += vec_it->_tp.intensity;
    }
    double r = sum_theoretic_ity / sum_observed_ity;
    for (vector<IsotopicPeakMatch>::iterator vec_it = _pmVec.begin(); vec_it != _pmVec.end(); vec_it++) {
        vec_it->_scaledTheoreticIntensity = vec_it->_tp.intensity / r;
    }
    return;
}

void IsotopicEnvelopeMatch::calcIntensityScores() {
    for (vector<IsotopicPeakMatch>::iterator vec_it = _pmVec.begin(); vec_it != _pmVec.end(); vec_it++) {
        double observed_ity = vec_it->_sp->intensity;
        double expected_ity = vec_it->_scaledTheoreticIntensity;
        if        (observed_ity <  expected_ity && (expected_ity - observed_ity)/observed_ity <= 1) {
            vec_it->_ityScore =       1.0 - (expected_ity - observed_ity)/observed_ity;
        } else if (observed_ity >= expected_ity && (observed_ity - expected_ity)/observed_ity <= 1) {
            vec_it->_ityScore = sqrt( 1.0 - (observed_ity - expected_ity)/observed_ity );
        } else {
            vec_it->_ityScore = 0.0;
        }
    }
    return;
}

void IsotopicEnvelopeMatch::print() const {
    cout << _pmVec[0]._tp.charge << "+\t" << _pevznerScore << "\t" << _myNormalizedPevznerScore << endl;
    for (vector<IsotopicPeakMatch>::const_iterator ipmVecIt = _pmVec.begin(); ipmVecIt != _pmVec.end(); ipmVecIt++) {
        cout << fixed << setprecision(6) << MSMass::int2double( ipmVecIt->_sp->m_z ) << "\t";
        cout << fixed << setprecision(2) << ipmVecIt->_sp->intensity << "\t";
        cout << fixed << setprecision(6) << MSMass::int2double( ipmVecIt->_tp.m_z ) << "\t";
        ////cout << fixed << setprecision(2) << ipmVecIt->_tp->intensity << "\t";
        cout << fixed << setprecision(2) << ipmVecIt->_scaledTheoreticIntensity << "\t";
        cout << fixed << setprecision(6) << MSMass::int2double(ipmVecIt->_sp->m_z) - MSMass::int2double(ipmVecIt->_tp.m_z) << "\t"; ////
//        cout << fixed << setprecision(2) << MSMass::int2double( ipmVecIt->_peakDeviation ) << "\t";
        cout << fixed << setprecision(3) << ipmVecIt->_ppmScore << "\t";
        cout << fixed << setprecision(3) << ipmVecIt->_ityScore << "\t";
        cout << fixed << setprecision(3) << ipmVecIt->_pevznerScore << endl;
    }
}
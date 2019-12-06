#include "Peak.hpp"
#include "DigestOptions.hpp"

SpectrumPeak::SpectrumPeak(unsigned int mz, float ity, short int z, bool pass) {
    m_z = mz;
    intensity = ity;
    charge = z;
    rank = 0;
    passedFilter = pass;
}

TheoreticalPeak::TheoreticalPeak(unsigned int mz, float ity, short int z, char it, int short p12) {
    m_z = mz;
    intensity = ity;
    charge = z;
    ion_type = it;
    pep_1_2 = p12;
}

PeakMatch::PeakMatch(const SpectrumPeak &sp, const TheoreticalPeak &tp) {
    _sp = &sp;
    _tp = &tp;
    calcPeakDeviation();
}

void PeakMatch::calcPeakDeviation() {
    _peakDeviation = _tp->getPeakDeviation(*_sp);
}

// When called from TheoreticalPeak, this will give (observed - expected) / expected sign
int Peak::getPeakDeviation(const Peak &p) const {
    long long diff = (long long) p.m_z - (long long) m_z;
    if (DigestOptions::fragmentTolRelative)
        return (1000000 * diff * (long long) MSMass::intDoubleConvertFactor) / (long long) m_z;
    else
        return diff;
}

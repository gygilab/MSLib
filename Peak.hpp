#ifndef PEAK_HPP_
#define PEAK_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <climits>
#include "MSMass.hpp"

using namespace std;

class Peak {
public:
    unsigned int m_z;
    float intensity;
    short int charge;

    bool operator<(const Peak& s) const {
        return (m_z < s.m_z); // Ascending, default
    }
    bool operator<(unsigned int d) const {
        return (m_z < d); // Ascending, default
    }
    bool operator>(unsigned int d) const {
        return (m_z > d); // Ascending, default
    }
    bool operator==(const Peak& s) const {
        return (m_z == s.m_z); // Ascending, default
    }

    int getPeakDeviation(const Peak &obs) const;
};

class SpectrumPeak : public Peak {
public:
    short int rank;
    bool passedFilter;

    SpectrumPeak() {};
    SpectrumPeak(unsigned int mz, float ity, short int z, bool pass);

    bool operator<(const SpectrumPeak& s) const {
        return (m_z < s.m_z); // Ascending, default
    }
    bool operator<(unsigned int d) const {
        return (m_z < d); // Ascending, default
    }
    bool operator>(unsigned int d) const {
        return (m_z > d); // Ascending, default
    }
    bool operator==(const SpectrumPeak& s) const {
        return (m_z == s.m_z); // Ascending, default
    }
};

struct SpectrumPeakCompare {
    bool operator() (const SpectrumPeak lhs, const SpectrumPeak rhs) const {
        return lhs.m_z < rhs.m_z;
    }
    bool operator() (SpectrumPeak lhs, const int &d) const {
        return lhs.m_z < d;
    }
    bool operator() (const int &d, const SpectrumPeak rhs) const {
        return d < rhs.m_z;
    }
};

struct SpectrumPeakPtrCompare {
    bool operator() (const SpectrumPeak * lhs, const SpectrumPeak * rhs) const {
        return lhs->m_z < rhs->m_z;
    }
};

struct sortByIntensity { // Descending
    inline bool operator() (const SpectrumPeak& s1, const SpectrumPeak & s2) const {
        return s1.intensity > s2.intensity;
    }
};

struct sortByRank {
    inline bool operator() (const SpectrumPeak& s1, const SpectrumPeak & s2) const {
        return s1.rank < s2.rank; // Ascending
    }
};

struct findSpectrumPeak {
    unsigned int d;
    findSpectrumPeak(unsigned int d) : d(d) {}
    bool operator () ( const SpectrumPeak& sp ) const {
        return sp.m_z == d;
    }
};

class TheoreticalPeak : public Peak {
public:
    char ion_type; // for now, b = 0, y = 1
    short int pep_1_2; // peptide 1 or 2

    TheoreticalPeak() {};
    TheoreticalPeak(unsigned int mz, float ity, short int z, char it, short int p12);

};

// Still need this because compiler throws an error in Match.cpp upper_bound when using the operator< type comparison
// The reason for this apparently is I sort fullPeakList using sortByM_ZandCharge in Fragment and therefore the call
// to upper_bound is confused. Could add the sort by charge to default operator function but need to test to make sure
// it's not slowing things down. 
struct TheoreticalPeakCompare {
    bool operator() (TheoreticalPeak lhs, TheoreticalPeak rhs) {
        return lhs.m_z < rhs.m_z;
    }
    bool operator() (TheoreticalPeak lhs, const int &d) {
        return lhs.m_z < d;
    }
    bool operator() (const int &d, TheoreticalPeak rhs) {
        return d < rhs.m_z;
    }
};

struct sortByM_ZandCharge {
    inline bool operator() (const TheoreticalPeak &p1, const TheoreticalPeak &p2) const {
        return p1.m_z < p2.m_z || (p1.m_z == p2.m_z && p1.charge < p2.charge); // Ascending
    }
};

// resolve ties by assigning ambiguous peaks to the longer of xlinked peptides.  Help avoid high scores for spurious short peptides
struct sortByM_ZandCharge_Longer1 {
    inline bool operator() (const TheoreticalPeak &p1, const TheoreticalPeak &p2) const {
        return p1.m_z < p2.m_z || (p1.m_z == p2.m_z && p1.charge < p2.charge) || (p1.m_z == p2.m_z && p1.charge == p2.charge && p1.pep_1_2 < p2.pep_1_2); // Ascending
    }
};

// resolve ties by assigning ambiguous peaks to the longer of xlinked peptides.  Help avoid high scores for spurious short peptides
struct sortByM_ZandCharge_Longer2 {
    inline bool operator() (const TheoreticalPeak &p1, const TheoreticalPeak &p2) const {
        return p1.m_z < p2.m_z || (p1.m_z == p2.m_z && p1.charge < p2.charge) || (p1.m_z == p2.m_z && p1.charge == p2.charge && p2.pep_1_2 < p1.pep_1_2); // Ascending
    }
};

class PeakMatch {
public:
    const SpectrumPeak * _sp;
    const TheoreticalPeak * _tp;
    int _peakDeviation;

    PeakMatch(const SpectrumPeak &sp, const TheoreticalPeak &tp);

    // Make sure that both experimental and theoretical m/z values are different
    bool operator<(const PeakMatch& pm) const {
        if (_sp->m_z != pm._sp->m_z && _tp->m_z != pm._tp->m_z)
            return _sp->m_z < pm._sp->m_z;
        else
            return false;
    }

    void calcPeakDeviation();
};

struct compareAbsPMDeviation {
    inline bool operator() (const PeakMatch &pm1, const PeakMatch &pm2) const {
        return abs(pm1._peakDeviation) < abs(pm2._peakDeviation); // Ascending
    }
};

#endif /*PEAK_HPP_*/

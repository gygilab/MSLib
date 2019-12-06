#ifndef _ISOTOPE_HPP
#define	_ISOTOPE_HPP

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cwchar>
#include <cfloat>
#include <algorithm>

#include "Peak.hpp"
#include "Isotope/emass.h"
#include "Isotope/parser.h"

class TheoreticalEnvelope {
public:
    vector<TheoreticalPeak> _peakVector;

    void print();
//
//    bool operator<(const TheoreticalEnvelope& te) const {
//        return (_peakVector[0].m_z < te._peakVector[0].m_z);
//    }
//    bool operator<(unsigned int& d) const {
//        return (_peakVector[0].m_z < d);
//    }
//    bool operator>(unsigned int& d) const {
//        return (_peakVector[0].m_z > d);
//    }
//    bool operator==(const TheoreticalEnvelope& te) const {
//        return (_peakVector[0].m_z == te._peakVector[0].m_z);
//    }
};

class IsotopicPeakMatch {
public:
    const SpectrumPeak * _sp;
    TheoreticalPeak _tp;
    int _peakDeviation;
    int _obsPeakIdx; //
    double _scaledTheoreticIntensity;
    double _ppmScore;
    double _ityScore;
    double _pevznerScore;

    IsotopicPeakMatch(const SpectrumPeak &sp, TheoreticalPeak &tp);
    void calcPeakDeviation();
};

class IsotopicEnvelopeMatch {
public:
    string _obsPeakIdx;
    vector<IsotopicPeakMatch> _pmVec;
    double _pevznerScore;
    double _myNormalizedPpmScore;
    double _myNormalizedItyScore;
    double _myNormalizedPevznerScore;

    /** The scoring system to match theoretical/observed isotopic envelopes follows with some modifications
     *  the one described in
     *  Deconvolution and database search of complex tandem mass spectra of intact proteins: a combinatorial approach.
     *  Liu X, Inbar Y, Dorrestein PC, Wynne C, Edwards N, Souda P, Whitelegge JP, Bafna V, Pevzner PA
     *  Mol Cell Proteomics. 2010 Dec;9(12):2772-82
     */
    void adjustIntensityScale();
    void calcPpmScores();
    void calcIntensityScores();
    void calcScores();
    void print() const;

    bool operator<(const IsotopicEnvelopeMatch& em) const {
        return (_pmVec[0]._sp->m_z < em._pmVec[0]._sp->m_z);
    }
    bool operator<(unsigned int& d) const {
        return (_pmVec[0]._sp->m_z < d);
    }
    bool operator>(unsigned int& d) const {
        return (_pmVec[0]._sp->m_z > d);
    }
    bool operator==(const IsotopicEnvelopeMatch& em) const {
        return (_pmVec[0]._sp->m_z == em._pmVec[0]._sp->m_z);
    }
};

struct findIEM {
    unsigned int d;
    findIEM(unsigned int d) : d(d) {}
    bool operator () ( const IsotopicEnvelopeMatch& iem ) const
    {
        return iem._pmVec[0]._sp->m_z == d;
    }
};

// pevznerScore descending
struct compareScore {
    inline bool operator() (const IsotopicEnvelopeMatch &em1, const IsotopicEnvelopeMatch &em2) const {
        return (em1._pevznerScore > em2._pevznerScore);
        //return (em1._myNormalizedPevznerScore > em2._myNormalizedPevznerScore);
    }
};

class Spectrum;

class Isotope {
public:
    map<unsigned long int,TheoreticalEnvelope> theoreticEnvelopeHash;

    Isotope ();
    TheoreticalEnvelope getTheoreticalEnvelope(double start_m_z, long charge);
    TheoreticalEnvelope getHashedTheoreticalEnvelope(double start_m_z, long charge);
    void addToEnvelopeHash (TheoreticalEnvelope &theo_e);
    void buildEnvelopeRecursive(TheoreticalEnvelope &theo_e, IsotopicEnvelopeMatch &env_match, const Spectrum &dta, vector<IsotopicEnvelopeMatch> &envelope_vec);
    
private:
    SuperAtomData sad;
    ElemMap em;

    const static string isotopefn;
    double limit;
    int digits;

    //Parser p;

    FormMap fm;

    Pattern result;
    Pattern tmp;

    Pattern getIsotopePattern(string mol_formula, long charge);
};


#endif	/* _ISOTOPE_HPP */


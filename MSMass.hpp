#ifndef _MSMASS_HPP
#define	_MSMASS_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>
#include <string>
#include <climits>
#include "DigestOptions.hpp"
#include "Modification.hpp"
using namespace std;

#include <tr1/unordered_map>
#include <boost/regex.hpp>

class MSMass {
public:
    static double numericMaxLimit;
    static double protonMass;
    static double intDoubleConvertFactor;
    static tr1::unordered_map<char, double> elementMass;
    static tr1::unordered_map<char, string> aaComp;
    static tr1::unordered_map<char, double> aa2mass;
    static tr1::unordered_map<string, double> commonMass;
    static tr1::unordered_map<char, double> averagineAtomicComp;
    
//    MSMass ();
    MSMass (DigestOptions opt);
    
    static double elComp2Mass (string el_comp);
    static double seq2MHMass (string seq);
    static double seq2MHMass (char flank_res_nterm, string seq, char flank_res_cterm);
//    static void seq2MassVec (string seq, vector<double> &massVec);
    static void seq2MassVec (char flank_res_nterm, string seq, char flank_res_cterm, vector<double> &massVec);

    static string makeIsobaric (const string &peptide);

    template<class T> static inline double calcPpm(T m1, T m2) {
        double mass1 = (double) m1;
        double mass2 = (double) m2;
        double ppm = (mass1 - mass2) / mass2 * 1000000.0;
        return ppm;
    }

    static inline double precursorMass2Diff (double query_mass) {
        if (DigestOptions::precursorTolRelative) {
            return query_mass * DigestOptions::precursorTol / 1000000.0;
        } else {
            return DigestOptions::precursorTol;
        }
    }

    static inline double fragmentMass2Diff (double query_mass) {
        if (DigestOptions::fragmentTolRelative) {
            return query_mass * DigestOptions::fragmentTol / 1000000.0;
        } else {
            return DigestOptions::fragmentTol;
        }
    }

    static inline int fragmentMass2Diff (unsigned int query_mass) {
        if (DigestOptions::fragmentTolRelative) {
            return (long long int) query_mass * (long long) DigestOptions::fragmentTol / 1000000;
        } else {
            return DigestOptions::fragmentTolInt;
        }
    }

    static inline unsigned int double2int (double d) {
        if (d >= int2double(UINT_MAX)) {
            cerr << "WARNING: " << d << " is too large to be converted to integer" << endl;
//            exit(0);
        }
        return (unsigned int) round(d*intDoubleConvertFactor);
    }

    static inline double int2double (double n) {
        return n / intDoubleConvertFactor;
    }

    static inline double m_h2m_z (double m_h, int z) {
        return ((m_h - protonMass) + ((double) z) * protonMass) / (double) z;
    }

    static inline double m_z2m_h (double m_z, int z) {
        return m_z * ((double) z) - ((double) z) * protonMass + protonMass;
    }

    static string averagineModel (double m);

private:
    void initElementMass ();
    void initCommonMass ();
    void addElementMassMods (DigestOptions opt);

    void initAAComp();

    void initAAMass ();
    void addStaticAAMods (DigestOptions opt);
    void checkVariableMods (DigestOptions opt);

    static char i2l (char c);
};

#endif	/* _MSMASS_HPP */


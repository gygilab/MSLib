#ifndef _FRAGMENT_HPP
#define	_FRAGMENT_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

#include "Peak.hpp"
#include "MSMass.hpp"
#include "DigestOptions.hpp"
#include "../JMUtil/JMUtil.hpp"
#include "Peptide.hpp"
#include "Run.hpp"

using namespace std;
using namespace JMUtil;

class FragmentTable {
public:
    string peptide;
    string peptide2;
    string id;
    double m_h_mass;
    int max_charge;
    int pos1;
    int pos2;
    int pos1_y; //// ??
    int pos2_y; //// ??
    vector<double> massVec;
    vector<double> massVec2;
    vector< vector<double> > b;
    vector< vector<double> > y;
    vector< vector<double> > b2;
    vector< vector<double> > y2;
    vector< vector<double> > lf; // linker fragments
    vector<TheoreticalPeak> fullPeakList; // keyed on int to prevent rounding issues
    int n1Peaks; ////
    int n2Peaks; ////
    vector< pair<unsigned int,unsigned int> > fdr_passed_bin_ranges;

    virtual void generateFullPeakList(double lo, double hi) = 0;

    virtual void prettyPrint(string outfile) = 0;
    void printPeakList(string outfile); // should be able to just get away with base implementation??
    int filterRanges(vector< pair<unsigned int,unsigned int> > fdr_passed_bin_ranges);

    virtual void fragment() = 0;
};


// fragmented peptide table
class PeptideFragmentTable : public FragmentTable {
public:

    PeptideFragmentTable(ModIterator &mIt, int c, double lo, double hi) {
        peptide = mIt.modSequence;
        id = mIt.pep->getProteinReference();
        m_h_mass = mIt.pep->_mPlusH + mIt.modMass;
        MSMass::seq2MassVec(mIt.pep->_flankResNterm, mIt.pep->_sequence, mIt.pep->_flankResCterm, massVec);

        for (int i = (mIt._modIndexVector.size() - 1); i >= 0; i--) {
            pair<int,Modification*> position_mod_ptr_pair = mIt.pep->_modSites[ mIt._modIndexVector[i] - 1 ];
            // check against multiple mods per residue was already done by MSMass
            massVec[ position_mod_ptr_pair.first ] += position_mod_ptr_pair.second->_mass;
        }

        max_charge = c > 1 ? (c - 1) : 1;
        max_charge = (max_charge > DigestOptions::maxTheoFragCharge) ? DigestOptions::maxTheoFragCharge : max_charge;

        b.resize(max_charge, vector<double> (massVec.size() - 1));
        y.resize(max_charge, vector<double> (massVec.size() - 1));
        fragment();
        generateFullPeakList(lo, hi);
    }


    void generateFullPeakList(double lo, double hi) {
        /* THEORETICALLY this could become a problem on some platforms where INT_MAX is lower than 2e9
         * causing loss of theoretical fragments.  Made this uint just in case and put in this check
         * to at least make sure no rogue peaks get in.  Should really be fine in vast majority of situations
         */
        if (hi > MSMass::numericMaxLimit) {
//            cerr << "WARNING: " << hi << " m/z value appears to exceed numeric limits " << hi << endl;
            hi = MSMass::numericMaxLimit;
        }
        fullPeakList.clear();
        fullPeakList.reserve((massVec.size() - 1) * b.size() * 2);
        for (int c = 0; c < b.size(); c++) {
            for (int i = 0; i < b[0].size(); i++) {
                if (b[c][i] >= lo  && b[c][i] <= hi) {
                    TheoreticalPeak p(MSMass::double2int( b[c][i] ), 0, (c + 1), 'b', 1);
                    fullPeakList.push_back(p);
                }
                if (y[c][i] >= lo && y[c][i] <= hi) {
                    TheoreticalPeak p(MSMass::double2int( y[c][i] ), 0, (c + 1), 'y', 1);
                    fullPeakList.push_back(p);
                }
            }
        }

        vector<TheoreticalPeak>::iterator thVecIt;
        // sort by m/z AND THEN charge. Occasionally see identical m/z for multiple charge states
        // and want to keep the lower one
        sort(fullPeakList.begin(), fullPeakList.end(), sortByM_ZandCharge());
        thVecIt = unique(fullPeakList.begin(), fullPeakList.end());
        fullPeakList.resize(thVecIt - fullPeakList.begin());
        return;
    }

    void prettyPrint(string outfile);

    void fragment();
};


// fragmented peptide table
class XLinkedFragmentTable : public FragmentTable {
private:
    double m_h_mass1;
    double m_h_mass2;
    char res1;
    char res2;
    bool pos1_term;
    bool pos2_term;
    
public:
    XLinkedFragmentTable(ModIterator &mIt, ModIterator &mIt2, int cloc1, int cloc2, int c, double lo, double hi) {
        id = mIt.pep->getProteinReference() + " - " + mIt2.pep->getProteinReference();;
        peptide = mIt.modSequence;
        peptide2 = mIt2.modSequence;
        pos1 = cloc1;
        pos2 = cloc2;
        res1 = mIt.pep->_sequence[pos1];
        res2 = mIt2.pep->_sequence[pos2];
        m_h_mass1 = mIt.pep->_mPlusH + mIt.modMass;
        m_h_mass2 = mIt2.pep->_mPlusH + mIt2.modMass;
        m_h_mass = m_h_mass1 + m_h_mass2 + DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass - MSMass::protonMass;
        MSMass::seq2MassVec(mIt.pep->_flankResNterm, mIt.pep->_sequence, mIt.pep->_flankResCterm, massVec);
        MSMass::seq2MassVec(mIt2.pep->_flankResNterm, mIt2.pep->_sequence, mIt2.pep->_flankResCterm, massVec2);

        for (int i = (mIt._modIndexVector.size() - 1); i >= 0; i--) {
            pair<int,Modification*> position_mod_ptr_pair = mIt.pep->_modSites[ mIt._modIndexVector[i] - 1 ];
            // check against multiple mods per residue was already done by MSMass
            massVec[ position_mod_ptr_pair.first ] += position_mod_ptr_pair.second->_mass;
        }
        for (int i = (mIt2._modIndexVector.size() - 1); i >= 0; i--) {
            pair<int,Modification*> position_mod_ptr_pair = mIt2.pep->_modSites[ mIt2._modIndexVector[i] - 1 ];
            // check against multiple mods per residue was already done by MSMass
            massVec2[ position_mod_ptr_pair.first ] += position_mod_ptr_pair.second->_mass;
        }
        pos1_y = massVec.size() - pos1 - 1;
        pos2_y = massVec2.size() - pos2 - 1;

        pos1_term = 0;
        if ((pos1 == 0                    && mIt.pep->_flankResNterm == '+') ||
            (pos1 == (massVec.size() - 1) && mIt.pep->_flankResCterm == '-')) {
            pos1_term = 1;
        }
        pos2_term = 0;
        if ((pos2 == 0                     && mIt2.pep->_flankResNterm == '+') ||
            (pos2 == (massVec2.size() - 1) && mIt2.pep->_flankResCterm == '-')) {
            pos2_term = 1;
        }

        max_charge = c > 1 ? (c - 1) : 1;
        max_charge = (max_charge > DigestOptions::maxTheoFragCharge) ? DigestOptions::maxTheoFragCharge : max_charge;

        b.resize(max_charge, vector<double> (massVec.size() - 1));
        y.resize(max_charge, vector<double> (massVec.size() - 1));
        b2.resize(max_charge, vector<double> (massVec2.size() - 1));
        y2.resize(max_charge, vector<double> (massVec2.size() - 1));
        fragment();
        generateFullPeakList(lo, hi);
    }


    void generateFullPeakList(double lo, double hi) {
        /* THEORETICALLY this could become a problem on some platforms where INT_MAX is lower than 2e9
         * causing loss of theoretical fragments.  Made this uint just in case and put in this check
         * to at least make sure no rogue peaks get in.  Should really be fine in vast majority of situations
         */
        if (hi > MSMass::numericMaxLimit) {
//            cerr << "WARNING: " << hi << " m/z value appears to exceed numeric limits " << hi << endl;
            hi = MSMass::numericMaxLimit;
        }
        fullPeakList.clear();
        fullPeakList.reserve((massVec.size() - 1) * (b.size() + b2.size()) * 2);
        for (int c = 0; c < b.size(); c++) {
            for (int i = 0; i < b[0].size(); i++) {
                if (b[c][i] >= lo  && b[c][i] <= hi) {
                    char ion_type = (i >= pos1) ? 'B' : 'b';
                    TheoreticalPeak p(MSMass::double2int( b[c][i] ), 0, (c + 1), ion_type, 1);
                    if (! Run::isInFilteredRange( p.m_z )) {
                        fullPeakList.push_back(p);
                    }
                }
                if (y[c][i] >= lo && y[c][i] <= hi) {
                    char ion_type = (i >= pos1_y) ? 'Y' : 'y';
                    TheoreticalPeak p(MSMass::double2int( y[c][i] ), 0, (c + 1), ion_type, 1);
                    if (! Run::isInFilteredRange( p.m_z )) {
                        fullPeakList.push_back(p);
                    }
                }
            }
            for (int i = 0; i < b2[0].size(); i++) {
                if (b2[c][i] >= lo  && b2[c][i] <= hi) {
                    char ion_type = (i >= pos2) ? 'B' : 'b';
                    TheoreticalPeak p(MSMass::double2int( b2[c][i] ), 0, (c + 1), ion_type, 2);
                    if (! Run::isInFilteredRange( p.m_z )) {
                        fullPeakList.push_back(p);
                    }
                }
                if (y2[c][i] >= lo && y2[c][i] <= hi) {
                    char ion_type = (i >= pos2_y) ? 'Y' : 'y';
                    TheoreticalPeak p(MSMass::double2int( y2[c][i] ), 0, (c + 1), ion_type, 2);
                    if (! Run::isInFilteredRange( p.m_z )) {
                        fullPeakList.push_back(p);
                    }
                }
            }
        }
        for (int c = 0; c < lf.size(); c++) {
            for (int i = 0; i < lf[0].size(); i++) {
                if (lf[c][i] >= lo  && lf[c][i] <= hi) {
                    int pep_12 = (i+1) % 2 == 0 ? 2 : 1;
                    TheoreticalPeak p(MSMass::double2int( lf[c][i] ), 0, (c + 1), 'L', pep_12);
                    if (! Run::isInFilteredRange( p.m_z )) {
                        fullPeakList.push_back(p);
                    }
                }
            }
        }

        vector<TheoreticalPeak>::iterator thVecIt;
        // sort by m/z AND THEN charge. Occasionally see identical m/z for multiple charge states
        // and want to keep the lower one. If still tied, assign the peak to the longer of the 2 peptides
        if (peptide.length() > peptide2.length())
            sort(fullPeakList.begin(), fullPeakList.end(), sortByM_ZandCharge_Longer1());
        else
            sort(fullPeakList.begin(), fullPeakList.end(), sortByM_ZandCharge_Longer2());
        thVecIt = unique(fullPeakList.begin(), fullPeakList.end());
        fullPeakList.resize(thVecIt - fullPeakList.begin());
        return;
    }

    void prettyPrint(string outfile);

    void fragment();
};


// fragmented peptide table
class LoopLinkedFragmentTable : public FragmentTable {
public:
    LoopLinkedFragmentTable(ModIterator &mIt, int cloc1, int cloc2, int c, double lo, double hi) {
        id = mIt.pep->getProteinReference();
        peptide = mIt.modSequence;
        pos1 = cloc1;
        pos2 = cloc2;
        m_h_mass = mIt.pep->_mPlusH + mIt.modMass + DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
        MSMass::seq2MassVec(mIt.pep->_flankResNterm, mIt.pep->_sequence, mIt.pep->_flankResCterm, massVec);

        for (int i = (mIt._modIndexVector.size() - 1); i >= 0; i--) {
            pair<int,Modification*> position_mod_ptr_pair = mIt.pep->_modSites[ mIt._modIndexVector[i] - 1 ];
            // check against multiple mods per residue was already done by MSMass
            massVec[ position_mod_ptr_pair.first ] += position_mod_ptr_pair.second->_mass;
        }

        pos1_y = massVec.size() - pos1 - 1;
        pos2_y = massVec.size() - pos2 - 1;
        max_charge = c > 1 ? (c - 1) : 1;
        max_charge = (max_charge > DigestOptions::maxTheoFragCharge) ? DigestOptions::maxTheoFragCharge : max_charge;

        b.resize(max_charge, vector<double> (massVec.size() - 1));
        y.resize(max_charge, vector<double> (massVec.size() - 1));
        fragment();
        generateFullPeakList(lo, hi);
    }


    void generateFullPeakList(double lo, double hi) {
        /* THEORETICALLY this could become a problem on some platforms where INT_MAX is lower than 2e9
         * causing loss of theoretical fragments.  Made this uint just in case and put in this check
         * to at least make sure no rogue peaks get in.  Should really be fine in vast majority of situations
         */
        if (hi > MSMass::numericMaxLimit) {
//            cerr << "WARNING: " << hi << " m/z value appears to exceed numeric limits " << hi << endl;
            hi = MSMass::numericMaxLimit;
        }
        fullPeakList.clear();
        fullPeakList.reserve((massVec.size() - 1) * b.size() * 2);
        for (int c = 0; c < b.size(); c++) {
            for (int i = 0; i < b[0].size(); i++) {
                if (b[c][i] >= lo  && b[c][i] <= hi) {
                    TheoreticalPeak p(MSMass::double2int( b[c][i] ), 0, (c + 1), 'b', 1);
                    fullPeakList.push_back(p);
                }
                if (y[c][i] >= lo && y[c][i] <= hi) {
                    TheoreticalPeak p(MSMass::double2int( y[c][i] ), 0, (c + 1), 'y', 1);
                    fullPeakList.push_back(p);
                }
            }
        }

        vector<TheoreticalPeak>::iterator thVecIt;
        // sort by m/z AND THEN charge. Occasionally see identical m/z for multiple charge states
        // and want to keep the lower one
        sort(fullPeakList.begin(), fullPeakList.end(), sortByM_ZandCharge());
        thVecIt = unique(fullPeakList.begin(), fullPeakList.end());
        fullPeakList.resize(thVecIt - fullPeakList.begin());
        return;
    }

    void prettyPrint(string outfile);

    void fragment();
};

#endif	/* _FRAGMENT_HPP */


#ifndef _PEPTIDE_HPP
#define	_PEPTIDE_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include "MSMass.hpp"
#include "DigestOptions.hpp"
#include "Modification.hpp"
#include "../JMUtil/Combinatorics.hpp"

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;

class Fasta;

class Peptide {
public:
    friend class ModIterator;
    friend class FragmentTable;
    double _mPlusH; // can this be float?? Should I not store this at all and compute on the fly??
    mutable string _sequence;
    // liveIndex helps handle multi-search cases, distinguishing Peptides stored in SearchResults but not to be considered
    // It's easier to turn the flag off prior to each re-indexing and then turn it on for relevant peptides than check different properties every time
    mutable bool _liveIndex; // determines if the peptide is to be considered in the current indexed search
    mutable int _start;
    mutable char _flankResNterm;
    mutable char _flankResCterm;
    static vector<string> proteinRefVec;
    mutable set<int> _redundantRefs; // this is cheating a bit. Should be a pointer instead
    mutable vector< pair<int,Modification*> > _modSites; // store position-mod pairs
//    unsigned short int m_c;
    vector<short int> _missedCleavPosArray;

    Peptide(string sequence, string protein, int startPos);
    Peptide(const Fasta& prot, int site1, int site2);
    Peptide(string modified_pep_seq);
    Peptide();

    // This is more memory efficient than keeping a variable around and doesn't appear to cost any appreciable cpu time
    bool isReverse () const {
        if (getProteinReference().substr(0, DigestOptions::reverseTag.size()).compare(DigestOptions::reverseTag) == 0) {
            return true;
        }
        return false;
    }

    string getProteinReference () const {
        if (proteinRefVec.size() > 0 && _redundantRefs.size() > 0)
            return proteinRefVec[ *(_redundantRefs.begin()) ];
        else
            return "UNKNOWN_PROTEIN";
    }
    static void storeRefVec (string ref, int cur_cnt);
    void updateRedundant (const Peptide &p) const;
    int getNumRedundantRefs (int top_n = 0) const;

    bool operator<(const Peptide& b) const {
        return MSMass::makeIsobaric(_sequence) < MSMass::makeIsobaric(b._sequence);
    }

    bool operator==(const Peptide& b) const {
        return MSMass::makeIsobaric(_sequence).compare(MSMass::makeIsobaric(b._sequence)) == 0;
    }

    void output();

private:
    static boost::mutex _write_mutex;
    void getMplusH();
};

namespace std { namespace tr1 {
    template <>
        class hash<Peptide>{
        public :
            size_t operator()(const Peptide &p ) const {
                return hash<string>()( MSMass::makeIsobaric(p._sequence) );
            }
    };
}}

class ModIterator : public std::iterator<std::input_iterator_tag, std::vector<int> > {
private:
    bool _last;
    bool _lastCombo;
    int _numPossibleMods,_numSimultaneousMods;

    void findNextMod() {
        Combinatorics::next_combination(_numPossibleMods, _modIndexVector, _lastCombo);
        if (_lastCombo) {
            if (_numSimultaneousMods < DigestOptions::nMaxVariableMods && _numSimultaneousMods < _numPossibleMods) {
                _numSimultaneousMods++;
                _modIndexVector.resize(_numSimultaneousMods);
                buildComboRestrictions();
                Combinatorics::next_combination(_numPossibleMods, _modIndexVector, _lastCombo);
            } else {
                _last = true;
                return;
            }
        }
    }
    
    /*
     * Build 2D vector for looking up disallowed combinations
     * These are represented by 1's in the upper right triangle only
     */
    void buildComboRestrictions() {
        _comboRestriction.resize(_numPossibleMods, vector<int>(_numPossibleMods));
        for (int i = 0; i < _numPossibleMods - 1; i++) {
            for (int j = i + 1; j < _numPossibleMods; j++) {
                if (pep->_modSites[i].first == pep->_modSites[j].first) {
                    _comboRestriction[i][j] = 1;
                }
            }
        }
        return;
    }
    
    /*
     * Check if the current _modIndexVector passes combo restrictions
     */
    bool passComboRestriction() {
        if (_modIndexVector.size() <= 1) {
            return true;
        }
        for (int i = 0; i < _modIndexVector.size() - 1; i++) {
            if (_comboRestriction[ _modIndexVector[i]-1 ][ _modIndexVector[i] ] > 0) {
                return false;
            }
        }
        return true;
    }

public:
    const Peptide * pep;
    vector<int> _modIndexVector; // 1-based
    string modSequence;
    double modMass;
    vector<double> massVec;
    short int numUnmodifiedCleavageSites;
    vector< vector<int> > _comboRestriction;

    ModIterator() : _last(1) {
    };

    ModIterator(const Peptide &p) : pep(&p) {
        _last = false;
        _lastCombo = true;
        _numSimultaneousMods = 0;
        if (DigestOptions::variableAAMod.size() == 0) {
            _numPossibleMods = 0;
            return;
        }

        // Should clear modSites here in case this peptide had already been used to initiate a ModIterator
        pep->_modSites.clear();

        /*
         * Currently protein terminal peptides are not going to get peptide terminal mods.
         * This can be fairly easily changed but not bothering atm
         */
        
        // look for N-terminal variable mods and assign to _modSites
        char vAAmIt_key = pep->_flankResNterm == '+' ? '+' : ']';
        multimap<char, Modification>::iterator vAAmIt; // = DigestOptions::variableAAMod.find(vAAmIt_key);
        pair <multimap<char,Modification>::iterator, multimap<char,Modification>::iterator> range_it = DigestOptions::variableAAMod.equal_range(vAAmIt_key);
        for (vAAmIt = range_it.first; vAAmIt != range_it.second; ++vAAmIt) {
            if (vAAmIt != DigestOptions::variableAAMod.end()) {
                if (vAAmIt->second._proteinIDs.size() == 0 ||
                   (vAAmIt->second._proteinIDs.size() > 0 && binary_search(vAAmIt->second._proteinIDs.begin(), vAAmIt->second._proteinIDs.end(), pep->getProteinReference()))) 
                {
                    pep->_modSites.push_back( pair<int,Modification*>(0, &(vAAmIt->second)) );
                }
            }
        }
        
        string::const_iterator start_pos_it = pep->_sequence.begin();
        string::const_iterator stop_pos_it = pep->_sequence.end();

        // prevent mod on terminal cleavage sites if desired.  Usually is, except isotopic labels maybe
        if (DigestOptions::enzymeAfterSites.second && pep->_flankResCterm != '-' && DigestOptions::noModSiteCleavage) {
            stop_pos_it--;
        }
        if (DigestOptions::enzymeBeforeSites.second && pep->_flankResNterm != '+' && DigestOptions::noModSiteCleavage) {
            start_pos_it++;
        }

        // look for regular residue mods and assign to _modSites
        const boost::regex re("[" + DigestOptions::variablyModAas + "]");
        boost::sregex_iterator reIt(start_pos_it, stop_pos_it, re);
        boost::sregex_iterator end;
        while (reIt != end) {
            //vAAmIt = DigestOptions::variableAAMod.find( reIt->str()[0] );
            range_it = DigestOptions::variableAAMod.equal_range( reIt->str()[0] );
            for (vAAmIt = range_it.first; vAAmIt != range_it.second; ++vAAmIt) {
                if (vAAmIt->second._proteinIDs.size() == 0 ||
                   (vAAmIt->second._proteinIDs.size() > 0 && binary_search(vAAmIt->second._proteinIDs.begin(), vAAmIt->second._proteinIDs.end(), pep->getProteinReference()))) 
                {
                    pep->_modSites.push_back( pair<int,Modification*>( reIt->position(), &(vAAmIt->second ) ) );
                }
            }
            reIt++;
        }

        // look for C-terminal variable mods and assign to _modSites
        vAAmIt_key = pep->_flankResCterm == '-' ? '-' : '[';
        range_it = DigestOptions::variableAAMod.equal_range(vAAmIt_key);
        //vAAmIt = DigestOptions::variableAAMod.find(vAAmIt_key);
        for (vAAmIt = range_it.first; vAAmIt != range_it.second; ++vAAmIt) {
            if (vAAmIt != DigestOptions::variableAAMod.end()) {
                if (vAAmIt->second._proteinIDs.size() == 0 ||
                   (vAAmIt->second._proteinIDs.size() > 0 && binary_search(vAAmIt->second._proteinIDs.begin(), vAAmIt->second._proteinIDs.end(), pep->getProteinReference()))) 
                {
                    pep->_modSites.push_back( pair<int,Modification*>(pep->_sequence.length() - 1, &(vAAmIt->second ) ) );
                }
            }
        }

        _numPossibleMods = pep->_modSites.size();
        buildComboRestrictions();
    }

    ModIterator & operator=(const Peptide &p) {
        pep = &p;
        findNextMod();
        while (!passComboRestriction() && !_last) {
            findNextMod();
        }
        return (*this);
    }

    vector<int> operator*() const {
        return _modIndexVector;
    }

    bool operator==(ModIterator rhs) const {
        bool retv = false;
        if ((_last == 1 && rhs._last == 1) || (_modIndexVector == rhs._modIndexVector && _numSimultaneousMods > 0 && rhs._numSimultaneousMods > 0)) {
            retv = true;
        }
        return (retv);
    }

    bool operator!=(ModIterator rhs) const {
        return !operator==(rhs);
    }

    ModIterator & operator++() {
        findNextMod();
        while (!passComboRestriction() && !_last) {
            findNextMod();
        }
        return (*this);
    }

    ModIterator operator++(int) {
        ModIterator ret(*this);
        findNextMod();
        while (!passComboRestriction() && !_last) {
            findNextMod();
        }
        return (ret);
    }

    void prepOutput();
    void output();
    string getFlankedSeq();
};


#endif	/* _PEPTIDE_HPP */


#ifndef _PROTEIN_HPP
#define	_PROTEIN_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <cctype>

#include "DigestOptions.hpp"
#include "Peptide.hpp"

#include <boost/regex.hpp>

using namespace std;

class Fasta {
public:
    friend class DigestIterator;

    Fasta (string id, string full_header, string seq, int cnt = 1);

    void setId(string id);
    void setFullHeader(string fullHeader);
    void setSequence(string sequence);
    void setCnt(int cnt);

    string getId() const;
    string getFullHeader() const;
    string getSequence() const;
    int getCnt() const;
    string getFormattedFasta() const;
    void outputDigest() const;

    bool isReverse() const;
    void swapCutSites();
    void swapRevInitMet();

private:
    string _id;
    string _fullHeader;
    string _sequence;
    int _sequentialCnt;
};

class DigestIterator : public std::iterator<std::input_iterator_tag, Peptide> { // change start
private:
    bool last;
    bool _cleavedInitMet;
    const Fasta * prot;
    int peptideCTermSite;
    int peptideNTermSite;
    int nTermCleavageSiteIdx;
    int cTermCleavageSiteIdx;
    int peptideNTermPartialSite;
    int peptideCTermPartialSite;
    int max_missed_cleavages;
    vector<int> sites;
    Peptide _peptide;

    bool evaluatePeptide() {
//        cout << "... evaluatePeptide()\t" << prot->getSequence().substr(peptideNTermSite, (peptideCTermSite - peptideNTermSite)) << endl; ////
//        cout << nTermCleavageSiteIdx << "-" << cTermCleavageSiteIdx << "\t" << peptideNTermSite << "-" << peptideCTermSite << endl; ////
//        cout << peptideNTermPartialSite << "-" << peptideCTermPartialSite << endl; ////
        if ((peptideCTermSite - peptideNTermSite) >= DigestOptions::minLength &&
                (peptideCTermSite - peptideNTermSite) <= DigestOptions::maxLength) {
            // Check mass constraints  and whether the peptide sequence contains any unrecognized residues/chars (M_H < 0)
            double m_h = MSMass::seq2MHMass( prot->getSequence().substr(peptideNTermSite, (peptideCTermSite - peptideNTermSite)) );
            if (m_h < DigestOptions::minMass || m_h > DigestOptions::maxMass) {
                return false;
            }
            _peptide = Peptide(*prot, peptideNTermSite, peptideCTermSite);
            for (int i = nTermCleavageSiteIdx + 1; i < cTermCleavageSiteIdx; i++) {
                _peptide._missedCleavPosArray.push_back(sites[i] - peptideNTermSite - 1);
            }
            _cleavedInitMet = false;
            return true;
        } else {
            return false;
        }
    }

    void findNextPeptide() {
        if (DigestOptions::enzymeSpecificity == 1 && peptideCTermPartialSite > (sites[cTermCleavageSiteIdx - 1] + 1)) {
            peptideCTermPartialSite--;
            peptideCTermSite = peptideCTermPartialSite;
        } else if (DigestOptions::enzymeSpecificity == 1 && peptideNTermPartialSite < (sites[nTermCleavageSiteIdx + 1] - 1)) {
            peptideCTermSite = sites[cTermCleavageSiteIdx]; // return if slipped back in the above clause
            peptideNTermPartialSite++;
            peptideNTermSite = peptideNTermPartialSite;
        } else if (cTermCleavageSiteIdx <= (nTermCleavageSiteIdx + max_missed_cleavages) && cTermCleavageSiteIdx < (sites.size() - 1)) {
            cTermCleavageSiteIdx++;
            peptideNTermSite = sites[nTermCleavageSiteIdx];
            peptideCTermSite = sites[cTermCleavageSiteIdx];
            peptideCTermPartialSite = peptideCTermSite;
            peptideNTermPartialSite = peptideNTermSite;
        } else if (nTermCleavageSiteIdx < (sites.size() - 2)) {
            nTermCleavageSiteIdx++;
            cTermCleavageSiteIdx = nTermCleavageSiteIdx + 1;
            peptideNTermSite = sites[nTermCleavageSiteIdx];
            peptideCTermSite = sites[cTermCleavageSiteIdx];
            peptideCTermPartialSite = peptideCTermSite;
            peptideNTermPartialSite = peptideNTermSite;
        } else {
            last = 1;
        }
        return;
    }

public:

    DigestIterator() : last(1) {};

    DigestIterator(const Fasta &S) : prot(&S) {
        last = false;
        _cleavedInitMet = false;
        sites.reserve(prot->_sequence.size());
        if (DigestOptions::enzymeSpecificity == 0) {
            for (int i = 0; i <= prot->_sequence.size(); i++) {
                sites.push_back(i);
            }
            max_missed_cleavages = DigestOptions::maxLength;
            /* TODO: Need to rethink this.  Maybe what we want is 2 definitions of non-specific search:
             * 1) A search with an enzyme definition and max_missed_cleavages defined where when specificity is set to 0
             *    the peptide termini are allowed to be non-specific but missed cleavages are still capped
             * 2) A true "no enzyme" search which will be triggered by a non-specific regex and user will have to remember
             *    to set max_missed_cleavages to maxLength on their own
             */
        } else {
            const boost::regex re(DigestOptions::enzymeRegEx);
            boost::sregex_iterator reIt((prot->_sequence).begin(), (prot->_sequence).end(), re);
            boost::sregex_iterator end;

            // deal with the start position
            if (reIt != end) {
                int pos = (*reIt++).position();
                if (pos > 0) {
                    sites.push_back(0);
                    sites.push_back(pos);
                }
            } else {
                sites.push_back(0);
            }
            while (reIt != end) {
                sites.push_back((*reIt++).position());
            }
            // Deal with end position
            if (sites.back() != prot->_sequence.size()) {
                sites.push_back(prot->_sequence.size());
            }
            max_missed_cleavages = DigestOptions::maxMissedCleavages;
        }
        nTermCleavageSiteIdx = 0;
        cTermCleavageSiteIdx = 1;
        peptideCTermSite = sites[cTermCleavageSiteIdx];
        peptideNTermSite = sites[nTermCleavageSiteIdx];
        peptideCTermPartialSite = peptideCTermSite;
        peptideNTermPartialSite = peptideNTermSite;
        // define the first peptide
        while (!evaluatePeptide() && !last) {
            findNextPeptide();
        }
        return;
    }

//    DigestIterator & operator=(const Fasta &S) {
//        prot = &S;
//        findNextPeptide();
//        return (*this);
//    }

    Peptide operator*() const {
        return _peptide;
    }

    bool operator==(DigestIterator rhs) const {
        bool retv = false;
        if ((nTermCleavageSiteIdx == rhs.nTermCleavageSiteIdx && cTermCleavageSiteIdx == rhs.cTermCleavageSiteIdx) || (last == 1 && rhs.last == 1)) {
            retv = true;
        }
        return (retv);
    }

    bool operator!=(DigestIterator rhs) const {
        return !operator==(rhs);
    }

//    DigestIterator & operator++() {
//        findNextPeptide();
//        return (*this);
//    }

    DigestIterator operator++(int) {
        // only produce initiator Met cleaved product when enzymeSpecificity == 2, otherwise will wind up with 2 versions, with different flanking Nterm residue
        if (DigestOptions::cleaveInitMet && DigestOptions::enzymeSpecificity == 2 && !_cleavedInitMet && _peptide._start == 1 && _peptide._sequence[0] == 'M') {
            peptideNTermSite++; // skip the Met
            if (evaluatePeptide()) {
                _peptide._flankResNterm = '+';
                _cleavedInitMet = true;
                DigestIterator ret(*this);
                return (ret);
            }
        }
        
        do {
            findNextPeptide();
        } while (!evaluatePeptide() && !last);
        DigestIterator ret(*this);
        return (ret);
    }
};

#endif	/* _PROTEIN_HPP */


#ifndef _PROTEOME_HPP
#define	_PROTEOME_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <tr1/unordered_set>
#include <tr1/unordered_map>

#include "DigestOptions.hpp"
#include "Protein.hpp"
#include "Peptide.hpp"

#include <boost/regex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/detail/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

using namespace std;

class Proteome {
public:
    tr1::unordered_map<string,int> proteinSet; // keep track of unique protein ids seen
    set<string> xProteinSet; // FWD proteins to be x-linked
    tr1::unordered_map<Peptide,int> peptideSet; // save matched Peptide records
    /* this is a temporary multi-threading hack. Keeps peptide in limbo while all mods
     * are being evaluated.  Prevents from being deleted by concurrent threads
     */
    tr1::unordered_set<Peptide> threadSafeLimbo;
    tr1::unordered_set<string> peptideSequenceSet; // keep track of unique sequences already seen
    multimap<double,ModIterator> modMap; // index all mod masses

    void indexPeptides();
    void indexPeptidesMT(int cur_thread_num = 0);
    void reIndexPeptides();
    void resetLiveIndex();
    void keepOnlyMatched();

    void storeProteinSet (string ref, int cur_cnt);
    void storePeptideSequenceSet (string seq);
    pair< tr1::unordered_map<Peptide,int>::iterator, bool > storePeptide(const Peptide& p);
    void incrementPeptideHit (tr1::unordered_map<Peptide,int>::iterator pIt);
    void incrementPeptideHit (Peptide p);
    void decrementPeptideHit (tr1::unordered_map<Peptide,int>::iterator);
    void decrementPeptideHit (Peptide p);
    void checkRemovePeptide (tr1::unordered_map<Peptide,int>::iterator pIt);
    void outputDigest();

/*    Fasta getFastaFromId (string id);
    Fasta getFastaFromfullHeader (string header);

    string getSequenceFromId (string id);
    string getSequenceFromfullHeader (string header);*/
    Proteome() : _write_mutex(new boost::mutex()) {}

private:
    boost::shared_ptr<boost::mutex> _write_mutex;
};


class FastaStreamIterator : public std::iterator<std::input_iterator_tag, Fasta> {
private:
    Proteome * _proteome;
    ifstream _infile;
    string _fullHeader;
    string _id;
    string _sequence;
    int _cur_cnt;

    void readNextFasta() {
        string line;

        boost::regex re(DigestOptions::fastaHeaderRegEx);
        boost::match_results<string::const_iterator> regex_result;

        while (getline(_infile, line, '\n')) {
            if (line[0] == '>') {
                if (boost::regex_match(line, regex_result, re)) {
                    _fullHeader = regex_result[0];
                    _id = regex_result[1];
                    _sequence = "";
                    _cur_cnt ++;
                    _proteome->storeProteinSet(_id, _cur_cnt);
                    Peptide::storeRefVec(_id, _cur_cnt);
                } else {
                    cerr << "Error parsing " << DigestOptions::fastaFile << ", entry: " << endl;
                    cerr << line << endl;
                    cerr << "|" << regex_result[0] << "|" << DigestOptions::fastaHeaderRegEx << "|" << endl; ////
                    exit(1);
                }
            } else if (isalnum(line[0])) {
                _sequence += line;
            } else {
                cerr << "WARNING: Fasta appears to contain empty lines or is otherwise improperly formatted" << endl;
                cerr << "         Attempting to read through ..." << endl;
            }
            if (_infile.peek() == EOF || _infile.peek() == '>') {
                transform(_sequence.begin(), _sequence.end(), _sequence.begin(), ::toupper);
                return;
            }
        }

        return;
    }

public:
    FastaStreamIterator(const FastaStreamIterator &b); //not defined, so no copy

    FastaStreamIterator & operator=(const FastaStreamIterator &b); //not defined, so no copy
    // don't allow copies, because streams aren't copiable. make sure to always pass by reference
    // unfortunately, this means no stl algorithms either

    FastaStreamIterator(Proteome &p) : _proteome(&p) {
        _infile.open(DigestOptions::fastaFile.c_str());
        if (!_infile.is_open()) {
            cout << "Error opening fasta file " << DigestOptions::fastaFile << endl;
            exit(0);
        }
        _cur_cnt = 0;
        // initial iterator is empty, so need to increment before doing anything. Needed for elegant processing of last fasta entry
    }

    FastaStreamIterator() {
        _infile.setstate(ios::eofbit);
    }

    FastaStreamIterator & operator++() {
        //cout << "FastaStreamIterator &operator++()" << endl; ////
        readNextFasta();
        return (*this);
    }

    // normally you'd want operator++(int) as well, but that requires making a copy
    // and std::ifstream isn't copiable.

    bool operator==(const FastaStreamIterator &b) const {
        if (_infile.bad() == b._infile.bad() &&
                _infile.eof() == b._infile.eof() &&
                _infile.good() == b._infile.good() &&
                _infile.fail() == b._infile.fail())
            return true; //all end() and uninitialized iterators equal
        // so we can use FastaStreamIterator() as end()
        return false; //since they all are seperate strings, never equal
    }

    bool operator!=(const FastaStreamIterator &b) const {
        return !operator==(b);
    }

    Fasta operator*() const {
        Fasta f(_id, _fullHeader, _sequence, _cur_cnt);
        if (DigestOptions::swapReverseInitMet && f.isReverse()) {
            f.swapRevInitMet();
        }
        if (DigestOptions::swapReverseCutSites && f.isReverse()) {
            f.swapCutSites();
        }
        return f;
    }
};

#endif	/* _PROTEOME_HPP */


#ifndef _SEARCHRESULT_HPP
#define	_SEARCHRESULT_HPP

#include <tr1/unordered_set>
#include "Peptide.hpp"

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "../JMUtil/JMUtil.hpp"

using namespace JMUtil;

class SearchResult {
public:
    string dtaPath;
    double score1;
    double deltaScore1;
    double diffSeqDeltaScore1;
    int nMatches;
    int nIons;
    string addMatchInfo;
    double fracItyExplained;
//    string debugInfo;
    double theo_m_h;
    double ppm;
    string peptide;
    const Peptide * pepPtr;
    string peptide2;
    const Peptide * pepPtr2;
    string posInfo;
    string searchType; ////

    virtual ~SearchResult() {};
    virtual string getSequence() = 0;
    virtual string getOutputString() const = 0;
    static string getHeaderString();
//    void printDebugInfo();
};

class PeptideSearchResult : public SearchResult {
public:

    PeptideSearchResult(string &p, const Peptide &pep);
    string getSequence();
    string getOutputString() const;

    //bool operator< (const SearchResult* b) const;
};

class XLinkedSearchResult : public SearchResult {
public:
    XLinkedSearchResult(string &p, string &p2, const Peptide &pep, const Peptide &pep2, string pos_info);
    string getSequence();
    string getOutputString() const;
};

class LoopLinkedSearchResult : public SearchResult {
public:
    LoopLinkedSearchResult(string &p, const Peptide &pep, string pos_info);
    string getSequence();
    string getOutputString() const;
};

#endif	/* _SEARCHRESULT_HPP */


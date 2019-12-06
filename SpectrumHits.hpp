#ifndef _SPECTRUMHITS_HPP
#define	_SPECTRUMHITS_HPP

#include <iostream>
#include <map>
#include <cmath>
#include <tr1/unordered_set>
#include "Proteome.hpp"
#include "SearchResult.hpp"
#include "DigestOptions.hpp"

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

struct SearchResultPtrCompare {
  bool operator() (const SearchResult * lhs, const SearchResult * rhs) const {
      return lhs->score1 < rhs->score1 ||
      (lhs->score1 == rhs->score1 && abs(lhs->ppm) > abs(rhs->ppm)) ||
      (lhs->score1 == rhs->score1 && abs(lhs->ppm) == abs(rhs->ppm) && *(lhs->pepPtr->_redundantRefs.begin()) > *(rhs->pepPtr->_redundantRefs.begin())) ||
      (lhs->score1 == rhs->score1 && abs(lhs->ppm) == abs(rhs->ppm) && *(lhs->pepPtr->_redundantRefs.begin()) == *(rhs->pepPtr->_redundantRefs.begin()) && lhs->peptide2.length() > 0 && rhs->peptide2.length() > 0 && *(lhs->pepPtr2->_redundantRefs.begin()) > *(rhs->pepPtr2->_redundantRefs.begin()));
  }
};

class SpectrumHits {
public:
    multiset<SearchResult *, SearchResultPtrCompare> hitMap;
    void storeHit (SearchResult * xsr);
    void storeHit (SearchResult * xsr, Proteome &proteome);
    void calcDeltaScore();
    double getTopHitScore1();
    double getBottomHitScore1();
    void output(ofstream &outFile);
    void output();

    SpectrumHits () : _write_mutex(new boost::mutex()) {
    }

private:
    boost::shared_ptr<boost::mutex> _write_mutex;
};

#endif	/* _SPECTRUMHITS_HPP */


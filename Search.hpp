#ifndef _SEARCH_HPP
#define	_SEARCH_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <functional>

#include "../JMUtil/JMUtil.hpp"
#include "DigestOptions.hpp"
#include "Peptide.hpp"
#include "Run.hpp"
#include "Proteome.hpp"
#include "Protein.hpp"
#include "Match.hpp"

#include <boost/thread/thread.hpp>
#include <boost/thread/detail/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

class SearchResult;
class FragmentTable;

class xLinkPosition {
public:
    ModIterator mIt;
    set<int> pos;

    xLinkPosition (ModIterator m, set<int> p) {
        mIt = m;
        pos = p;
    }
};

struct sortMHxLinkPositionPair {
    inline bool operator() (pair<double,xLinkPosition> i, pair<double,xLinkPosition> j) {
        return i.first < j.first;
    }
    inline bool operator() (pair<double,xLinkPosition> i, const double &d) {
        return i.first < d;
    }
    inline bool operator() (const double &d, pair<double,xLinkPosition> j) {
        return d < j.first;
    }
};

/* Functor that can be used with stl containers that keeps an index counter
 * and returns true for even indices.  Note that to use this properly, need c++11 and std::reference_wrapper
 * Currently used without those when removing half of (rev)modMaps incorrectly.  gcc implementation lead to loss of 
 * the first even index which is fine given the size of that vector
 */
class isIndexEvenFunctor {
public:
  explicit isIndexEvenFunctor() : k(0) {}

  bool operator()(const pair<double,xLinkPosition> &i) {
    if(k++ % 2 == 0) {
      return true;
    } else {
      return false;
    }
  }
  
private:
  int k;
};

/* Sample uniformly from an container.  Keeps an index counter as well as a counter
 * of removed elements.  Requires c++11 to work fully although c++98 should only
 * skip the first occurrence and be produce very close results for large containers
 */
class uniformSamplingFunctor {
public:
  long int num_total;
  long int num_to_remove;
  
  explicit uniformSamplingFunctor(long int n_t, long int n_t_r) : k(0), r(1) {
    num_total = n_t;
    num_to_remove = n_t_r;
    pos_to_remove = JMUtil::roundTo(r * num_total / num_to_remove, 0);
    shift_for_pos_to_remove = JMUtil::roundTo(pos_to_remove / 2.0, 0);
    pos_to_remove -= shift_for_pos_to_remove;
  }

  bool operator()(const pair<double, xLinkPosition> &i) {
    //cout << "k: " << k << "; r: " << r << endl; ////
    if(k++ == pos_to_remove) {
      ++r;
      //cout << "removing " << pos_to_remove << endl; ////
      pos_to_remove = JMUtil::roundTo(r * num_total / num_to_remove, 0) - shift_for_pos_to_remove;
      return true;
    } else {
      return false;
    }
  }
private:
  long int k;
  long int r;
  long int pos_to_remove;
  long int shift_for_pos_to_remove;
};


class Search {
public:
    typedef void (Search::*ptr2searchFnc)(Run &run, Proteome &proteome);
    map<string, ptr2searchFnc> searchFncPtrMap;

    Search(Run &run, Proteome &proteome);
    void searchPeptidesIndexSpectra(Run &run, Proteome &proteome);
    void searchPeptidesIndexSpectraMT(Run &run, Proteome &proteome, int cur_thread_num = 0);
    void searchPeptidesIndexMods(Run &run, Proteome &proteome);
    void searchPeptidesIndexAll(Run &run, Proteome &proteome);
    void searchPeptidesIndexAllMT(Run &run, Proteome &proteome, int cur_thread_num = 0);

    void searchLoopLinkedPeptides(Run &run, Proteome &proteome);
    void searchHomoLoopLinkedMT(Run &run, Proteome &proteome, int cur_thread_num = 0);
    void searchHeteroLoopLinkedMT(Run &run, Proteome &proteome, int cur_thread_num = 0);

    void searchXLinkedPeptides(Run &run, Proteome &proteome);
    void searchHomoXLinked(Run &run, Proteome &proteome);
    void searchHomoXLinkedIndexedDB(Run &run, Proteome &proteome, vector<pair< double,xLinkPosition> > &massModMap, int cur_thread_num = 0);
    void searchHeteroXLinked(Run &run, Proteome &proteome);
    void searchHeteroXLinkedIndexedDB(Run &run, Proteome &proteome, vector<pair< double,xLinkPosition> > &massModMap1, vector<pair< double,xLinkPosition> > &massModMap2, int cur_thread_num = 0);
    void matchXlinkedIndexedDBRange(Run &run, Proteome &proteome, vector<pair< double,xLinkPosition> >::iterator mmMIt, vector<pair< double,xLinkPosition> >::iterator mmMIt2, int cur_thread_num = 0);
    void filterRevModToMatchFwd (vector<pair< double,xLinkPosition> > &revModMap, vector<pair< double,xLinkPosition> > &modMap);

    vector< pair<int,int> > getXLinkPairs(set<int> &pos1, set<int> &pos2);
    vector< pair<int,int> > getXLinkPairs(set<int> &pos);
    set<int> getXLinkPositions(ModIterator &mIt, string linked_res);

private:
    static vector<int> fastaCntVec;
    static vector<double> pepCntVec;
    static vector<double> psmCntVec;
    static vector<double> modCntVec;
//    static map<int,int> count_mc;
    double linkerMass;

    static void zeroOutCnts () {
        fill(fastaCntVec.begin(), fastaCntVec.end(), 0);
        fill(pepCntVec.begin(), pepCntVec.end(), 0);
        fill(psmCntVec.begin(), psmCntVec.end(), 0);
        fill(modCntVec.begin(), modCntVec.end(), 0);
    }

    static int getTotalFastaCnt () {
        return accumulate(fastaCntVec.begin(), fastaCntVec.end(), 0);
    }
    static double getTotalPepCnt () {
        return accumulate(pepCntVec.begin(), pepCntVec.end(), 0.0);
    }
    static double getTotalPsmCnt () {
        return accumulate(psmCntVec.begin(), psmCntVec.end(), 0.0);
    }
    static double getTotalModCnt () {
        return accumulate(modCntVec.begin(), modCntVec.end(), 0.0);
    }
};


#endif	/* _SEARCH_HPP */


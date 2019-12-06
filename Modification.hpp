/* 
 * File:   Modification.hpp
 * Author: jmintser
 *
 * Created on June 7, 2011, 1:35 AM
 */

#ifndef _MODIFICATION_HPP
#define	_MODIFICATION_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include <boost/algorithm/string.hpp>

#include "DigestOptions.hpp"

using namespace std;

class DigestOptions;

class Modification {
public:
    double _mass;
    char _symbol;
//    string _name;
    vector<string> _proteinIDs;

    Modification (double mass, char symbol, string prot_refs = "") {
        _mass = mass;
        _symbol = symbol;
        if (prot_refs.length() > 0) {
            boost::split(_proteinIDs, prot_refs, boost::is_any_of(DigestOptions::varModProteinSpecificDelimiter));
            int prot_cnt = _proteinIDs.size();
            for (int i = 0; i < prot_cnt; i++) {
                _proteinIDs.push_back(DigestOptions::reverseTag + _proteinIDs[i]);
            }
            sort(_proteinIDs.begin(), _proteinIDs.end()); // have to make sure it's sorted for binary search to work later in modIterator
        }
    }
};

#endif	/* _MODIFICATION_HPP */


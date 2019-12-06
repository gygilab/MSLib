#include "Peptide.hpp"
#include "Protein.hpp"

vector<string> Peptide::proteinRefVec;
boost::mutex Peptide::_write_mutex;

Peptide::Peptide (const Fasta& prot, int site1, int site2) {
    _sequence = prot.getSequence().substr(site1, (site2 - site1));
    _liveIndex = false;
    _start = site1 + 1;
    if (site1 == 0) {
        _flankResNterm = '+';
    } else {
        _flankResNterm = prot.getSequence()[site1 - 1];
    }
    if (site2 == prot.getSequence().length()) {
        _flankResCterm = '-';
    } else {
        _flankResCterm = prot.getSequence()[site2];
    }
    _redundantRefs.insert(prot.getCnt());
    getMplusH();
}

Peptide::Peptide(string sequence, string protein, int startPos) {
    _sequence = sequence;
    _liveIndex = false;
    _start = startPos;
    _redundantRefs.insert(1);
    getMplusH();
}

Peptide::Peptide(string modified_pep) {
    boost::regex re("^(.\\.){0,1}(.+?)(\\..){0,1}$");
    boost::match_results<string::const_iterator> regex_result;
    if (boost::regex_match(modified_pep, regex_result, re)) {
        string flank_res_n_term = regex_result[1];
        string modified_seq = regex_result[2];
        string flank_res_c_term = regex_result[3];
        boost::regex re_repl("[^A-Z]");
        _sequence = boost::regex_replace(modified_seq, re_repl, "");
        _liveIndex = false;
        _start = 1;
        _flankResNterm = flank_res_n_term[0];
        _flankResCterm = flank_res_c_term[1];
    } else {
        cerr << "Error parsing modified sequence: " << modified_pep << endl;
        exit(1);
    }
    getMplusH();
}

Peptide::Peptide() {
    _sequence = "";
    _liveIndex = false;
    _mPlusH = 0.0;
    _start = 1;
    _flankResNterm = '+';
    _flankResCterm = '-';
}

void Peptide::getMplusH() {
    _mPlusH = MSMass::seq2MHMass(_flankResNterm, _sequence, _flankResCterm);
}

int Peptide::getNumRedundantRefs(int top_n) const {
    if (top_n == 0) {
        return _redundantRefs.size();
    } else {
        return count_if(_redundantRefs.begin(), _redundantRefs.end(), bind2nd(less_equal<int>(), top_n));
    }
}

void Peptide::updateRedundant (const Peptide &p) const {
    _write_mutex.lock();
    if (*(p._redundantRefs.begin()) < *(_redundantRefs.begin()) ||
        (*(p._redundantRefs.begin()) == *(_redundantRefs.begin()) && p._start < _start)) {
        _start = p._start;
        _flankResNterm = p._flankResNterm;
        _flankResCterm = p._flankResCterm;
        _sequence = p._sequence;
    }
    _redundantRefs.insert(p._redundantRefs.begin(), p._redundantRefs.end());
    _write_mutex.unlock();
    return;
}

void Peptide::storeRefVec (string ref, int cur_cnt) {
    _write_mutex.lock();
    if (cur_cnt >= proteinRefVec.size()) { // grow the vector as needed
        proteinRefVec.resize(cur_cnt + 1);
    }
    proteinRefVec[cur_cnt] = ref;
    _write_mutex.unlock();
    return;
}

// construct the final sequence with all modifications
void ModIterator::prepOutput() {
    modSequence = pep->_sequence;
    modMass = 0.0;
    numUnmodifiedCleavageSites = pep->_missedCleavPosArray.size();
    for (int i = (_modIndexVector.size() - 1); i >= 0; i--) {
        int mod_position = pep->_modSites[ _modIndexVector[i] - 1 ].first;
        char aa = modSequence[ mod_position ];
        string aa_mod = "";
        if (DigestOptions::lowerCaseMods) {
            aa_mod.push_back(tolower(aa));
        } else {
            aa_mod.push_back(aa);
        }

        if (find(pep->_missedCleavPosArray.begin(), pep->_missedCleavPosArray.end(), mod_position) != pep->_missedCleavPosArray.end()) {
            numUnmodifiedCleavageSites --;
        }
        
        // check against multiple was done by MSMass
        Modification * mod_ptr = pep->_modSites[ _modIndexVector[i] - 1 ].second;
        if (!DigestOptions::lowerCaseMods) {
            aa_mod.push_back((mod_ptr->_symbol));
        }
        modMass += mod_ptr->_mass;

        modSequence = modSequence.replace(mod_position, 1, aa_mod);
    }
}

void ModIterator::output() {
    prepOutput();
    cout << fixed << setprecision(4)
            << (pep->_mPlusH + modMass) << "\t"
            << pep->_start << "\t"
            << (pep->_start + pep->_sequence.length() - 1) << "\t"
            //                 << modSequence << "\t"
            << pep->_flankResNterm << "." << modSequence << "." << pep->_flankResCterm << "\t"
            << pep->getProteinReference() << endl;
//    cout << "JMDEBUG: _n="  << _numPossibleMods << "; _k=" << _numSimultaneousMods << "\t"; ////
//    copy(_modIndexVector.begin(), _modIndexVector.end(), ostream_iterator<int>(cout, " ")); cout << endl; ////
}

void Peptide::output() {
    ModIterator mIt(*this), lastModIter;
    while (mIt != lastModIter) {
        //print(mIt.I);
        mIt.output();
        mIt++;
    }
    return;
}

string ModIterator::getFlankedSeq() {
    return pep->_flankResNterm + string(".") + modSequence + string(".") + pep->_flankResCterm;
}

//bool Peptide::operator<(const Peptide& b) const {
//    return _sequence < b._sequence;
//}
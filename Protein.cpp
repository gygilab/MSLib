#include "Protein.hpp"

Fasta::Fasta (string id, string full_header, string seq, int cnt) {
    _id = id;
    _fullHeader = full_header;
    _sequence = seq;
    _sequentialCnt = cnt;
    return;
}

void Fasta::setId (string id) {
    _id = id;
    return;
}
void Fasta::setFullHeader (string fullHeader) {
    _fullHeader = fullHeader;
    return;
}
void Fasta::setSequence (string sequence) {
    _sequence = sequence;
    return;
}

string Fasta::getId() const {
    return _id;
}
string Fasta::getFullHeader() const {
    return _fullHeader;
}
string Fasta::getSequence() const {
    return _sequence;
}
int Fasta::getCnt() const {
    return _sequentialCnt;
}

bool Fasta::isReverse() const {
    if (_id.substr(0, DigestOptions::reverseTag.size()).compare(DigestOptions::reverseTag) == 0) {
        return true;
    }
    return false;
}

/* Make sure REV proteins have fake init Mets
 */
void Fasta::swapRevInitMet() {
    if (*_sequence.rbegin() == 'M') {
        _sequence = "M" + _sequence.substr(0, _sequence.length() - 1);
    }
    return;
}

/* Swap cut sites for reversed proteins to avoid the issue of reverse peptides having the 
   same exact mass as the forward.  Following MaxQuant paper 2008
 * Turns out to be important to target all eg [KR] sites, including those before a P
 * If you don't, you introduce a small FWD bias for KP/PK cases
 */
void Fasta::swapCutSites() {
    string s;
    if (DigestOptions::enzymeAfterSites.second) {
        s = DigestOptions::enzymeAfterSites.first;
    } else if (DigestOptions::enzymeBeforeSites.second) {
        s = DigestOptions::enzymeBeforeSites.first;
    }

    for (string::iterator i = s.begin(); i != s.end(); i++) {
        if (!isalpha(s.at(i - s.begin()))) {
            s.erase(i);
            i--;
        }
    }

    if (DigestOptions::enzymeAfterSites.second) {
        s = "(?<=[" + s + "])";
    } else if (DigestOptions::enzymeBeforeSites.second) {
        s = "(?=[" + s + "])";
    }
    
    const boost::regex re(s);
    boost::sregex_iterator reIt((_sequence).begin(), (_sequence).end(), re);
    boost::sregex_iterator end;

    string old_seq = _sequence;
    vector<int> sites;
    sites.reserve(_sequence.size());
    while (reIt != end) {
        sites.push_back((*reIt++).position());
    }
    for (vector<int>::iterator vecIt = sites.begin(); vecIt != sites.end(); vecIt++) {
        if (*vecIt > 1) {
            swap(_sequence[*vecIt - 2], _sequence[*vecIt - 1]);
        }
    }
    return;
}

void Fasta::outputDigest() const {
    cout << _fullHeader << endl;
    DigestIterator dIt(*this);
    while (dIt != DigestIterator()) {
        Peptide pep = *dIt;
        pep.output();
        dIt++;
    }
    return;
}

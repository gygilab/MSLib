#include "Proteome.hpp"

void Proteome::indexPeptides() {
    cout << "...Indexing peptides" << endl; ////
    DigestOptions::indexAllPeptides = 1; 
    DigestOptions::indexSeenSequences = 0; /// TODO: this should be reset in DigestOptions
    if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
        indexPeptidesMT();
    } else {
        boost::thread_group btg;
        for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
            if (DigestOptions::verbose) {
                cout << "... starting thread # " << i << endl; ////
            }
            btg.add_thread(new boost::thread(boost::bind(&Proteome::indexPeptidesMT, this, i)));
        }
        btg.join_all();
    }
            
    // Some sanity checks on the fasta
    if (proteinSet.size() % 2 != 0) {
        cerr << "Error: Number of proteins " << proteinSet.size() << " appears to be not even suggesting that not all FWD proteins match a corresponding REV" << endl;
        exit(0);
    }
    for(tr1::unordered_map<string,int>::const_iterator fpsIt = proteinSet.begin(); fpsIt != proteinSet.end(); fpsIt++) {
        if (fpsIt->first.substr(0, DigestOptions::reverseTag.size()).compare(DigestOptions::reverseTag) != 0) {
            // make sure all forwards have a corresponding reverse
            if (proteinSet.find( DigestOptions::reverseTag + fpsIt->first ) == proteinSet.end()) {
                cerr << "Error: FWD protein " << fpsIt->first << " appears to not have a corresponding REV" << endl;
                exit(0);
            }
        }
    }
    
    /* make xProteinSet to contain only fwd proteins in the order they were first seen in the fasta
     * This will be used to enforce xFilterTopProts
     */
    if (xProteinSet.size() == 0 && DigestOptions::xFilterTopProts > 0) {
        vector<string>::const_iterator prvIt = Peptide::proteinRefVec.begin();
        while (prvIt != Peptide::proteinRefVec.end() && xProteinSet.size() < DigestOptions::xFilterTopProts) {
            if (prvIt->substr(0, DigestOptions::reverseTag.size()).compare(DigestOptions::reverseTag) != 0) {
                xProteinSet.insert(*prvIt);
            }
            prvIt++;
        }
    }

    cout << endl << peptideSet.size() << " peptides stored in peptideMap" << endl; ////
    return;
}

void Proteome::indexPeptidesMT(int cur_thread_num) {
    int local_cnt = 0;
    FastaStreamIterator fIt(*this), lastFasta;
    while (fIt != lastFasta) {
        ++fIt;

        local_cnt++;
        if ((local_cnt - cur_thread_num - 1) % DigestOptions::maxNumThreads != 0)
            continue;

        Fasta f = *fIt;
        DigestIterator dIt(f), lastPeptide;
        while (dIt != lastPeptide) {

            Peptide p = (*dIt);
            pair < tr1::unordered_map<Peptide,int>::iterator, bool> p_pair = storePeptide(p);
            (*p_pair.first).first._liveIndex = true;
            dIt++;

        }
    }
    return;
}

void Proteome::reIndexPeptides() {
    cout << "...re-indexing peptides" << endl; ////
    // clear unmatched peptides
    cout << "Trimming peptideMap from " << peptideSet.size() << " to ";
    for (tr1::unordered_map<Peptide,int>::iterator map_it = peptideSet.begin(); map_it != peptideSet.end();) {
        if (map_it->second <= 0) {
            peptideSet.erase(map_it++);
        } else {
            map_it->first._liveIndex = false; // make sure to temporarily remove from consideration all currently stored peptides
            ++map_it;
        }
    }
    cout << peptideSet.size() << endl;
    threadSafeLimbo.clear(); ////

    indexPeptides();
    return;
}

/**
 * Just re-set peptides considered for next search without deleting anything
 */
void Proteome::resetLiveIndex() {
    cout << "...re-setting liveIndexed peptides" << endl; ////
    // clear unmatched peptides
    for (tr1::unordered_map<Peptide,int>::iterator map_it = peptideSet.begin(); map_it != peptideSet.end(); map_it++) {
        map_it->first._liveIndex = false; // make sure to temporarily remove from consideration all currently stored peptides
    }
    threadSafeLimbo.clear(); ////

    indexPeptides();
    return;
}

void Proteome::keepOnlyMatched() {
    cout << "Trimming peptideMap from " << peptideSet.size() << " to ";
    for (tr1::unordered_map<Peptide,int>::iterator map_it = peptideSet.begin(); map_it != peptideSet.end(); map_it++) {
        if (map_it->second <= 0) {
            peptideSet.erase(map_it->first);
        }
    }
    cout << peptideSet.size() << endl;
    threadSafeLimbo.clear(); ////
}

// Digest each protein based on the Params
void Proteome::outputDigest() {
    FastaStreamIterator fIt(*this), lastFasta;
    while (fIt != lastFasta) {
        ++fIt;
        (*fIt).outputDigest();
    }
    return;
}

void Proteome::storePeptideSequenceSet (string seq) {
    _write_mutex->lock();
    peptideSequenceSet.insert(seq);
    _write_mutex->unlock();
    return;
}

void Proteome::storeProteinSet (string ref, int cur_cnt) {
    _write_mutex->lock();
    pair < tr1::unordered_map<string, int>::iterator, bool > ret_pair = proteinSet.insert(pair<string, int>(ref, cur_cnt));
    _write_mutex->unlock();
    if (!ret_pair.second && ret_pair.first->second != cur_cnt) { // This check is not necessary but also inexpensive
        cerr << endl << "ERROR: Multiple fasta sequences detected with same header " << ref << " " << cur_cnt << endl;
        exit(0);
    }
    return;
}

// return false if sequence has been seen and/or scored
pair< tr1::unordered_map<Peptide,int>::iterator, bool > Proteome::storePeptide(const Peptide& p) {
    _write_mutex->lock();
    if (DigestOptions::indexSeenSequences) {
        tr1::unordered_set<string>::iterator pssIt = peptideSequenceSet.find(p._sequence);
        if (pssIt != peptideSequenceSet.end()) { // sequence has been scored but did not score well enough to be stored in top N
            _write_mutex->unlock();
            return pair< tr1::unordered_map<Peptide,int>::iterator, bool >(peptideSet.end(), false);
        }
    }
    pair< tr1::unordered_map<Peptide,int>::iterator, bool> p_ret = peptideSet.insert(pair<Peptide,int>(p,0));
    if (p_ret.second == false) { // key exists
        (*p_ret.first).first.updateRedundant(p);
    }
    if (!DigestOptions::indexAllPeptides) {
        threadSafeLimbo.insert(p);
    }
    _write_mutex->unlock();
    return p_ret;
}

void Proteome::incrementPeptideHit (tr1::unordered_map<Peptide,int>::iterator pIt) {
    _write_mutex->lock();
    pIt->second++;
    _write_mutex->unlock();
    return;
}

void Proteome::incrementPeptideHit (Peptide p) {
    tr1::unordered_map<Peptide,int>::iterator pIt = peptideSet.find(p);
    incrementPeptideHit(pIt);
    return;
}

void Proteome::checkRemovePeptide (tr1::unordered_map<Peptide,int>::iterator pIt) {
    _write_mutex->lock();
    threadSafeLimbo.erase(pIt->first);
    if (pIt->second <= 0) {
        if (DigestOptions::indexSeenSequences) {
            peptideSequenceSet.insert(pIt->first._sequence);
        }
        peptideSet.erase(pIt->first);
    }
    _write_mutex->unlock();
    return;
}

void Proteome::decrementPeptideHit (tr1::unordered_map<Peptide,int>::iterator pIt) {
    _write_mutex->lock();
    pIt->second--;
    // boost::mutex isn't recursive so need to duplicate code from checkRemovePeptide for now
    tr1::unordered_set<Peptide>::iterator l_it = threadSafeLimbo.find(pIt->first);
    if (pIt->second <= 0 && l_it == threadSafeLimbo.end()) {
        if (DigestOptions::indexSeenSequences) {
            peptideSequenceSet.insert(pIt->first._sequence);
        }
        if (!DigestOptions::indexAllPeptides) {
            peptideSet.erase(pIt->first);
        }
    }
    _write_mutex->unlock();
    return;
}

void Proteome::decrementPeptideHit (Peptide p) {
    tr1::unordered_map<Peptide,int>::iterator pIt = peptideSet.find(p);
    decrementPeptideHit(pIt);
    return;
}

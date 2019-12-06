#include "Search.hpp"

vector<int> Search::fastaCntVec;
vector<double> Search::pepCntVec;
vector<double> Search::psmCntVec;
vector<double> Search::modCntVec;
//map<int,int> Search::count_mc;

/**
 * For regular peptide searches have the following options, ranked by increasing memory usage:
 *      - search while indexing spectra only (DigestOptions::indexAllPeptides and DigestOptions::indexSeenSequences should be set to 0)
 *        use this when searching huge fastas like NR
 *      - search while indexing spectra, not indexing all the peptides but indexing peptide sequences
 *      - search while indexing spectra and peptides
 * @param run
 * @param proteome
 */
Search::Search(Run& run, Proteome& proteome) {
    // Assume we always index spectra
    double vm, rss;
    run.loadAllSpectra(); // already multi-threaded in Run
    JMUtil::procMemUsage(vm, rss); cout << "VM: " << vm << "; RSS: " << rss << endl; ////
    
    searchFncPtrMap["PeptideIndexSpectra"] = &Search::searchPeptidesIndexSpectra;
    searchFncPtrMap["PeptideIndexAll"]     = &Search::searchPeptidesIndexAll;
    searchFncPtrMap["LoopLinkedPeptide"]   = &Search::searchLoopLinkedPeptides;
    searchFncPtrMap["XLink"]               = &Search::searchXLinkedPeptides;

    fastaCntVec.resize(DigestOptions::maxNumThreads);
    pepCntVec.resize(DigestOptions::maxNumThreads);
    psmCntVec.resize(DigestOptions::maxNumThreads);
    modCntVec.resize(DigestOptions::maxNumThreads);

    for (int s = 0; s < DigestOptions::multiSearchOptions.size(); s++) {
        DigestOptions::currentSearchOpt = s;
        if (searchFncPtrMap.find(DigestOptions::multiSearchOptions[s].searchType) != searchFncPtrMap.end()) {
            Search::ptr2searchFnc search_func = searchFncPtrMap[ DigestOptions::multiSearchOptions[s].searchType ];
            (this->*search_func) (run, proteome);
        } else {
            cerr << "ERROR: SearchType '" << DigestOptions::multiSearchOptions[s].searchType << "' not recognized!" << endl;
            exit(0);
        }
    }
    run.outputTab(DigestOptions::outputFile);
}

void Search::searchPeptidesIndexSpectra(Run& run, Proteome& proteome) {
    double vm, rss;
    cout << "Search::searchPeptidesIndexSpectra" << "\t"
         << DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].searchType << endl;
    DigestOptions::indexAllPeptides = 0;

    if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
        searchPeptidesIndexSpectraMT(run, proteome, 0);
    } else {
        boost::thread_group btg;
        for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
            if (DigestOptions::verbose) {
                cout << "... starting thread # " << i << endl; ////
            }
            btg.add_thread(new boost::thread(boost::bind(&Search::searchPeptidesIndexSpectraMT, this, boost::ref(run), boost::ref(proteome), i)));
        }
        btg.join_all();
    }
    cout << endl << proteome.peptideSet.size() << " peptides stored in peptideMap" << endl; ////
    cout << endl << proteome.peptideSequenceSet.size() << " sequences stored in peptideSequenceSet" << endl; ////
    cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalPepCnt() << " peptides (" << getTotalModCnt() << " mods)" << endl; ////
    JMUtil::procMemUsage(vm, rss); cout << "VM: " << vm << "; RSS: " << rss << endl; ////
    return;
}

void Search::searchPeptidesIndexSpectraMT(Run& run, Proteome& proteome, int cur_thread_num) {
    zeroOutCnts();
    int local_cnt = 0;
    FastaStreamIterator fIt(proteome), lastFasta;
    while (fIt != lastFasta) {
        ++fIt;

        local_cnt++;
        if ((local_cnt - cur_thread_num - 1) % DigestOptions::maxNumThreads != 0)
            continue;

        Fasta f = *fIt;
        if (getTotalFastaCnt() % 1000 == 0) ////
            cout << "." << flush; ////

        DigestIterator dIt(f), lastPeptide;
        while (dIt != lastPeptide) {

            pepCntVec[cur_thread_num] ++; ////
            Peptide p = (*dIt);

            pair < tr1::unordered_map<Peptide,int>::iterator, bool> p_pair = proteome.storePeptide(p);
            if (p_pair.second == false) { // peptide is already a scored and stored top N search result
                dIt++;
                continue;
            }
            
            ModIterator mIt(p), lastMod;
            while (mIt != lastMod) {
                modCntVec[cur_thread_num] ++; ////
                mIt.prepOutput();
                double query_mass = p._mPlusH + mIt.modMass;
                double mass_diff = MSMass::precursorMass2Diff(query_mass);
                for (int step = 0; step < DigestOptions::precursorSteps.size(); step++) {
                    multimap<double, Spectrum>::iterator specItHi = run.msnSpectra.lower_bound(query_mass + DigestOptions::precursorSteps[step] + mass_diff);
                    multimap<double, Spectrum>::iterator specItLo = run.msnSpectra.upper_bound(query_mass + DigestOptions::precursorSteps[step] - mass_diff);
                    if (specItLo != specItHi) {
                        SpectralMatch xsm(specItLo, specItHi, proteome, p_pair.first, mIt);
                        psmCntVec[cur_thread_num] += distance(specItLo, specItHi); ////
                    }
                }
                mIt++;
            }
            proteome.checkRemovePeptide(p_pair.first);
            dIt++;
        }
        fastaCntVec[cur_thread_num] ++; ////
    }
    return;
}

void Search::searchPeptidesIndexAll(Run& run, Proteome& proteome) {
    cout << "Search::searchPeptidesIndexAll" << "\t"
         << DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].searchType << endl;

    double vm, rss;
    proteome.indexPeptides();
    JMUtil::procMemUsage(vm, rss); cout << "VM: " << vm << "; RSS: " << rss << endl; ////
    
    if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
        searchPeptidesIndexAllMT(run, proteome, 0);
    } else {
        boost::thread_group btg;
        for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
            if (DigestOptions::verbose) {
                cout << "... starting thread # " << i << endl; ////
            }
            btg.add_thread(new boost::thread(boost::bind(&Search::searchPeptidesIndexAllMT, this, boost::ref(run), boost::ref(proteome), i)));
        }
        btg.join_all();
    }
    cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalPepCnt() << " peptides (" << getTotalModCnt() << " mods)" << endl; ////
    JMUtil::procMemUsage(vm, rss); cout << "VM: " << vm << "; RSS: " << rss << endl; ////
}

void Search::searchPeptidesIndexAllMT(Run& run, Proteome& proteome, int cur_thread_num) {
    zeroOutCnts();
    int local_cnt = 0;
    for (tr1::unordered_map<Peptide,int>::iterator pIt = proteome.peptideSet.begin(); pIt != proteome.peptideSet.end(); pIt++) {
        local_cnt++;
        if ((local_cnt - cur_thread_num - 1) % DigestOptions::maxNumThreads != 0)
            continue;

        if (fmod(getTotalPepCnt(), 100000.0) == 0.0) ////
            cout << "." << flush; ////
        Peptide p = (*pIt).first;

        ModIterator mIt(p), lastMod;
        while (mIt != lastMod) {
            modCntVec[cur_thread_num] ++; ////
            mIt.prepOutput();
            double query_mass = p._mPlusH + mIt.modMass;
            double mass_diff = MSMass::precursorMass2Diff(query_mass);
            for (int step = 0; step < DigestOptions::precursorSteps.size(); step++) {
                multimap<double, Spectrum>::iterator specItHi = run.msnSpectra.lower_bound(query_mass + DigestOptions::precursorSteps[step] + mass_diff);
                multimap<double, Spectrum>::iterator specItLo = run.msnSpectra.upper_bound(query_mass + DigestOptions::precursorSteps[step] - mass_diff);
                if (specItLo != specItHi) {
                    SpectralMatch psm(specItLo, specItHi, proteome, pIt, mIt);
                    psmCntVec[cur_thread_num] += distance(specItLo, specItHi); ////
                }
            }

            mIt++;
        }

        pepCntVec[cur_thread_num] ++; ////
    }
    return;
}

// More of a spectrum-indexed approach.  Iterate over all possible xlinkable mod pairs
void Search::searchXLinkedPeptides(Run& run, Proteome& proteome) {
    DigestOptions::filterVarModsForXlink();
            
    if (DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].isHomoXLink) {
        searchHomoXLinked(run, proteome);
    } else {
        searchHeteroXLinked(run, proteome);
    }
}

// More of a spectrum-indexed approach.  Iterate over all possible xlinkable mod pairs
void Search::searchHomoXLinked(Run& run, Proteome& proteome) {
    cout << "Search::searchXLinkedPeptides" << "\t"
         << DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].join2string() << endl;

    if (DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].enzymeSpecificity != DigestOptions::enzymeSpecificity ||
        DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].effectiveMaxMissedCleavages1 != DigestOptions::maxMissedCleavages) {
        
        // re-init options
        cout << "Re-initializing [Digest] options: " << endl;
        DigestOptions::enzymeSpecificity = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].enzymeSpecificity;
        DigestOptions::maxMissedCleavages = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].effectiveMaxMissedCleavages1;

        proteome.reIndexPeptides();
    }

    double vm, rss;
    zeroOutCnts();
    vector<pair< double,xLinkPosition> > modMap;
    vector<pair< double,xLinkPosition> > revModMap;

    // Approximate the upper bound on the indexed mods to save memory (for non-specific digest searches)
    double max_query_mass = run.msnSpectra.rbegin()->first + MSMass::precursorMass2Diff(run.msnSpectra.rbegin()->first);
    double min_pep_mass = DigestOptions::minMass;
    double max_mod_pep_mass = max_query_mass - min_pep_mass + MSMass::protonMass - linkerMass;

    for (tr1::unordered_map<Peptide, int>::iterator dIt = proteome.peptideSet.begin(); dIt != proteome.peptideSet.end(); dIt++) {
        if (!(*dIt).first._liveIndex)
            continue;
        
        string cur_ref = (*dIt).first.getProteinReference();
        if (cur_ref.find("_contaminant") != string::npos)
            continue;

        if ((*dIt).first.isReverse()) {
            cur_ref = cur_ref.substr(DigestOptions::reverseTag.size());
        }
        if (proteome.xProteinSet.find(cur_ref) == proteome.xProteinSet.end())
            continue;

        ModIterator lastMod;
        for (ModIterator mIt((*dIt).first); mIt != lastMod; mIt++) {
            mIt.prepOutput();
            if ((*dIt).first._mPlusH + mIt.modMass > max_mod_pep_mass)
                continue;

            set<int> pos = Search::getXLinkPositions(mIt, DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkedRes1);
            if (pos.size() > 0) {
                if (mIt.pep->isReverse()) {
                    revModMap.push_back(pair<double,xLinkPosition> ((*dIt).first._mPlusH + mIt.modMass, xLinkPosition(mIt,pos)));
                } else {
                    modMap.push_back(pair<double,xLinkPosition> ((*dIt).first._mPlusH + mIt.modMass, xLinkPosition(mIt,pos)));
                }
            }
        }
    }
    sort(modMap.begin(), modMap.end(), sortMHxLinkPositionPair()); cout << modMap.size() << " mods stored in modMap" << endl;
//    for (vector<pair< double,xLinkPosition> >::iterator vecIt = modMap.begin(); vecIt != modMap.end(); vecIt++) {
//        cout << fixed << setprecision(4) << "JMDEBUG: modMap H:\t" << vecIt->first << "\t" << vecIt->second.mIt.modSequence << endl;
//    }
    sort(revModMap.begin(), revModMap.end(), sortMHxLinkPositionPair()); cout << revModMap.size() << " mods stored in revModMap" << endl;
//    for (vector<pair< double,xLinkPosition> >::iterator vecIt = revModMap.begin(); vecIt != revModMap.end(); vecIt++) {
//        cout << fixed << setprecision(4) << "JMDEBUG: revModMap H:\t" << vecIt->first << "\t" << vecIt->second.mIt.modSequence << endl;
//    }
//    return; ////
    if (!DigestOptions::skipREVREVsearch) {
        if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
            searchHomoXLinkedIndexedDB(run, proteome, revModMap, 0);
        } else {
            boost::thread_group btg;
            for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
                if (DigestOptions::verbose) {
                    cout << "... starting thread # " << i << endl; ////
                }
                btg.add_thread(new boost::thread(boost::bind(&Search::searchHomoXLinkedIndexedDB, this, boost::ref(run), boost::ref(proteome), boost::ref(revModMap), i)));
            }
            btg.join_all();
        }
        cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalModCnt() << " mod pairs" << endl << endl; ////
    }
    
    if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
        searchHomoXLinkedIndexedDB(run, proteome, modMap, 0);
    } else {
        boost::thread_group btg;
        for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
            if (DigestOptions::verbose) {
                cout << "... starting thread # " << i << endl; ////
            }
            btg.add_thread(new boost::thread(boost::bind(&Search::searchHomoXLinkedIndexedDB, this, boost::ref(run), boost::ref(proteome), boost::ref(modMap), i)));
        }
        btg.join_all();
    }
    cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalModCnt() << " mod pairs" << endl << endl; ////

    if (!DigestOptions::skipREVFWDsearch) {
        if (DigestOptions::skipHalfREVFWDsearch) {
            filterRevModToMatchFwd(revModMap, modMap);
            //filterRevModToMatchFwd(revModMap);            
            cout << revModMap.size() << " mods stored in revModMap" << endl;
        }

        cout << "... searching REV-FWD/FWD-REV pairs" << endl; ////
        if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
            searchHeteroXLinkedIndexedDB(run, proteome, revModMap, modMap, 0);
        } else {
            boost::thread_group btg;
            for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
                if (DigestOptions::verbose) {
                    cout << "... starting thread # " << i << endl; ////
                }
                btg.add_thread(new boost::thread(boost::bind(&Search::searchHeteroXLinkedIndexedDB, this, boost::ref(run), boost::ref(proteome), boost::ref(revModMap), boost::ref(modMap), i)));
            }
            btg.join_all();
        }
        cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalModCnt() << " mod pairs" << endl << endl; ////
    //    cout << "JMDEBUG: MC 0: " << count_mc[0] << "; MC 1: " << count_mc[1] << "; MC 2: " << count_mc[2] << "; MC 3: " << count_mc[3] << "; MC 4: " << count_mc[4] << "; MC 5: " << count_mc[5] << endl; ////
    }
    
    JMUtil::procMemUsage(vm, rss); cout << "VM: " << vm << "; RSS: " << rss << endl; ////
    return;
}


/* The basic idea is to remove half of the revModMap.  The simplest thing to do is to go through
 * the sorted map and remove every other entry.  This turns out leave distributions that are slightly 
 * biased when binned by M+H because short peptides ~< 7 are more likely to map to FWD than REV.
 * To deal with that, we bin modMaps into 100 M+H bins and remove enough revMods from each bin
 * to get to the target of 1/2 of however many mods are in corresponding modMap bin
 */
void Search::filterRevModToMatchFwd(vector<pair<double,xLinkPosition> >& revModMap, vector<pair<double,xLinkPosition> >& modMap) {
    // initialize histograms
    map<int, int> modMap_hist;
    for (int i = 0; i <= floor(modMap.rbegin()->first / 100.0) * 100; i += 100) {
        modMap_hist[i] = 0;
    }
    map<int, int> revModMap_hist = modMap_hist;
    // populate histograms
    for (vector<pair< double, xLinkPosition> >::iterator mm_it = modMap.begin(); mm_it != modMap.end(); mm_it++) {
        modMap_hist[ floor(mm_it->first / 100.0) * 100 ]++;
    }
    for (vector<pair< double, xLinkPosition> >::iterator mm_it = revModMap.begin(); mm_it != revModMap.end(); mm_it++) {
        revModMap_hist[ floor(mm_it->first / 100.0) * 100 ]++;
    }

    // merge bins if too small
    map<int, int>::iterator mm_hist_it = modMap_hist.begin();
    map<int, int>::iterator mm_hist_it_next;
    while (mm_hist_it->first < 700 && mm_hist_it != modMap_hist.end()) {
        while (mm_hist_it->second < 10 && mm_hist_it != modMap_hist.end()) {
            mm_hist_it_next = mm_hist_it;
            mm_hist_it_next++;
            modMap_hist[mm_hist_it_next->first] += modMap_hist[mm_hist_it->first];
            revModMap_hist[mm_hist_it_next->first] += revModMap_hist[mm_hist_it->first];
//            cout << "JMDEBUG: merging bins " << mm_hist_it->first << " into " << mm_hist_it_next->first << endl; ////
            revModMap_hist.erase(mm_hist_it->first);
            modMap_hist.erase(mm_hist_it++);
        }
        mm_hist_it++;
    }
    mm_hist_it = modMap_hist.lower_bound(800);
    map<int, int>::iterator mm_hist_it_penultimate = modMap_hist.end();
    --mm_hist_it_penultimate;
    while (mm_hist_it != mm_hist_it_penultimate && mm_hist_it != modMap_hist.end()) {
        while (mm_hist_it->second < 10 && mm_hist_it != modMap_hist.end()) {
            mm_hist_it_next = mm_hist_it;
            mm_hist_it_next++;
            modMap_hist[mm_hist_it_next->first] += modMap_hist[mm_hist_it->first];
            revModMap_hist[mm_hist_it_next->first] += revModMap_hist[mm_hist_it->first];
//            cout << "JMDEBUG: merging bins " << mm_hist_it->first << " into " << mm_hist_it_next->first << endl; ////
            revModMap_hist.erase(mm_hist_it->first);
            modMap_hist.erase(mm_hist_it++);
        }
        mm_hist_it++;
    }

    vector<pair< double, xLinkPosition> >::iterator mm_it_lo;
    vector<pair< double, xLinkPosition> >::iterator mm_it_hi;
    for (map<int, int>::iterator mm_hist_it = revModMap_hist.begin(); mm_hist_it != revModMap_hist.end(); mm_hist_it++) {
        mm_it_lo = lower_bound(revModMap.begin(), revModMap.end(), mm_hist_it->first, sortMHxLinkPositionPair());
        map<int, int>::iterator mm_hist_it_next = mm_hist_it;
        mm_hist_it_next++;
        if (mm_hist_it_next != revModMap_hist.end()) {
            mm_it_hi = upper_bound(revModMap.begin(), revModMap.end(), mm_hist_it_next->first, sortMHxLinkPositionPair());
        } else {
            mm_it_hi = revModMap.end();
        }
//        cout << "JMDEBUG: cutting revModMap bin " << mm_hist_it->first << " between " << mm_it_lo->first << " and " << mm_it_hi->first << endl; ////
        int i = 0;
        int num_revMods_to_remove = JMUtil::roundTo((double) revModMap_hist[ mm_hist_it->first ] - (double) modMap_hist[ mm_hist_it->first ] / 2.0, 0);
        
        if (DigestOptions::verbose) {
            cout << "JMDEBUG: modMap bin " << mm_hist_it->first << ": " << modMap_hist[ mm_hist_it->first ] << "\trevModMap bin: " << revModMap_hist[ mm_hist_it->first ] << "\tremoving: " << num_revMods_to_remove << endl; ////
        }    
        if (num_revMods_to_remove < 1) continue;
        
        int total_mods_in_bin = revModMap_hist[ mm_hist_it->first ];
        int pos_to_remove = JMUtil::roundTo(1 * total_mods_in_bin / num_revMods_to_remove, 0);

        uniformSamplingFunctor usf = uniformSamplingFunctor(total_mods_in_bin, num_revMods_to_remove);
        revModMap.erase(remove_if(mm_it_lo, mm_it_hi, std::ref(usf)), mm_it_hi);
    }
    return;
}

void Search::searchHeteroXLinked(Run& run, Proteome& proteome) {
    cout << "Search::searchXLinkedPeptides" << "\t"
         << DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].join2string() << endl;

    double vm, rss;
    zeroOutCnts();
    vector<pair< double,xLinkPosition> > modMap;
    vector<pair< double,xLinkPosition> > revModMap;
    vector<pair< double,xLinkPosition> > modMap2;
    vector<pair< double,xLinkPosition> > revModMap2;

    // Approximate the upper bound on the indexed mods to save memory (for non-specific digest searches)
    double max_query_mass = run.msnSpectra.rbegin()->first + MSMass::precursorMass2Diff(run.msnSpectra.rbegin()->first);
    double min_pep_mass = DigestOptions::minMass;
    double max_mod_pep_mass = max_query_mass - min_pep_mass + MSMass::protonMass - linkerMass;

    if (DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].enzymeSpecificity != DigestOptions::enzymeSpecificity ||
        DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].effectiveMaxMissedCleavages1 != DigestOptions::maxMissedCleavages) {
        
        // re-init options
        cout << "Re-initializing [Digest] options for peptide 1: " << endl;
        DigestOptions::enzymeSpecificity = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].enzymeSpecificity;
        DigestOptions::maxMissedCleavages = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].effectiveMaxMissedCleavages1;

        proteome.reIndexPeptides();
    }
    
    // Generate mods separately for peptide1 and peptide2 sets because they may have different max missed cleavages
    // peptide 1
    for (tr1::unordered_map<Peptide,int>::iterator dIt = proteome.peptideSet.begin(); dIt != proteome.peptideSet.end(); dIt++) {
        if (!(*dIt).first._liveIndex)
            continue;

        string cur_ref = (*dIt).first.getProteinReference();
        if (cur_ref.find("_contaminant") != string::npos)
            continue;
        if ((*dIt).first.isReverse()) {
            cur_ref = cur_ref.substr(DigestOptions::reverseTag.size());
        }
        if (proteome.xProteinSet.find(cur_ref) == proteome.xProteinSet.end())
            continue;
        
        ModIterator lastMod;
        for (ModIterator mIt((*dIt).first); mIt != lastMod; mIt++) {
            mIt.prepOutput();
            if ((*dIt).first._mPlusH + mIt.modMass > max_mod_pep_mass)
                continue;
            set<int> pos1 = Search::getXLinkPositions(mIt, DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkedRes1);
            if (pos1.size() > 0) {
                if (mIt.pep->isReverse()) {
                    revModMap.push_back(pair<double,xLinkPosition>((*dIt).first._mPlusH + mIt.modMass, xLinkPosition(mIt,pos1)));
                } else {
                    modMap.push_back(pair<double,xLinkPosition>((*dIt).first._mPlusH + mIt.modMass, xLinkPosition(mIt,pos1)));
                }
            }
        }
    }
    
    if (DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].enzymeSpecificity != DigestOptions::enzymeSpecificity ||
        DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].effectiveMaxMissedCleavages2 != DigestOptions::maxMissedCleavages) {
        
        // re-init options
        cout << "Re-initializing [Digest] options for peptide 2: " << endl;
        DigestOptions::enzymeSpecificity = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].enzymeSpecificity;
        DigestOptions::maxMissedCleavages = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].effectiveMaxMissedCleavages2;

        proteome.resetLiveIndex();
    }
    
    // Generate mods separately for peptide1 and peptide2 sets because they may have different max missed cleavages
    // peptide 2
    for (tr1::unordered_map<Peptide,int>::iterator dIt = proteome.peptideSet.begin(); dIt != proteome.peptideSet.end(); dIt++) {
        if (!(*dIt).first._liveIndex)
            continue;

        string cur_ref = (*dIt).first.getProteinReference();
        if (cur_ref.find("_contaminant") != string::npos)
            continue;
        if ((*dIt).first.isReverse()) {
            cur_ref = cur_ref.substr(DigestOptions::reverseTag.size());
        }
        if (proteome.xProteinSet.find(cur_ref) == proteome.xProteinSet.end())
            continue;
        
        ModIterator lastMod;
        for (ModIterator mIt((*dIt).first); mIt != lastMod; mIt++) {
            mIt.prepOutput();
            if ((*dIt).first._mPlusH + mIt.modMass > max_mod_pep_mass)
                continue;
            set<int> pos2 = Search::getXLinkPositions(mIt, DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkedRes2);
            if (pos2.size() > 0) {
                if (mIt.pep->isReverse()) {
                    revModMap2.push_back(pair<double,xLinkPosition>((*dIt).first._mPlusH + mIt.modMass, xLinkPosition(mIt,pos2)));
                } else {
                    modMap2.push_back(pair<double,xLinkPosition>((*dIt).first._mPlusH + mIt.modMass, xLinkPosition(mIt,pos2)));
                }
            }
        }
    }
    
    sort(modMap.begin(), modMap.end(), sortMHxLinkPositionPair()); cout << modMap.size() << " mods stored in modMap" << endl;
    sort(revModMap.begin(), revModMap.end(), sortMHxLinkPositionPair()); cout << revModMap.size() << " mods stored in revModMap" << endl;
    sort(modMap2.begin(), modMap2.end(), sortMHxLinkPositionPair()); cout << modMap2.size() << " mods stored in modMap2" << endl;
    sort(revModMap2.begin(), revModMap2.end(), sortMHxLinkPositionPair()); cout << revModMap2.size() << " mods stored in revModMap2" << endl;

    if (!DigestOptions::skipREVREVsearch) {
        if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
            searchHeteroXLinkedIndexedDB(run, proteome, revModMap, revModMap2, 0);
        } else {
            boost::thread_group btg;
            for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
                if (DigestOptions::verbose) {
                    cout << "... starting thread # " << i << endl; ////
                }
                btg.add_thread(new boost::thread(boost::bind(&Search::searchHeteroXLinkedIndexedDB, this, boost::ref(run), boost::ref(proteome), boost::ref(revModMap), boost::ref(revModMap2), i)));
            }
            btg.join_all();
        }
        cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalModCnt() << " mod pairs" << endl << endl; ////
    }

    if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
        searchHeteroXLinkedIndexedDB(run, proteome, modMap, modMap2, 0);
    } else {
        boost::thread_group btg;
        for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
            if (DigestOptions::verbose) {
                cout << "... starting thread # " << i << endl; ////
            }
            btg.add_thread(new boost::thread(boost::bind(&Search::searchHeteroXLinkedIndexedDB, this, boost::ref(run), boost::ref(proteome), boost::ref(modMap), boost::ref(modMap2), i)));
        }
        btg.join_all();
    }
    cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalModCnt() << " mod pairs" << endl << endl; ////

    if (!DigestOptions::skipREVFWDsearch) {
        if (DigestOptions::skipHalfREVFWDsearch) {
            
            cout << "... cutting revModMap in half" << endl;
            filterRevModToMatchFwd(revModMap, modMap);
            //filterRevModToMatchFwd(revModMap);
            cout << endl;
            cout << "... cutting revModMap2 in half" << endl;
            filterRevModToMatchFwd(revModMap2, modMap2);
            //filterRevModToMatchFwd(revModMap);
            
            cout << revModMap.size() << " mods stored in revModMap" << endl;
            cout << revModMap2.size() << " mods stored in revModMap2" << endl;
        }

        cout << "... searching REV-FWD/FWD-REV pairs" << endl; ////
        if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
            searchHeteroXLinkedIndexedDB(run, proteome, revModMap, modMap2, 0);
        } else {
            boost::thread_group btg;
            for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
                if (DigestOptions::verbose) {
                    cout << "... starting thread # " << i << endl; ////
                }
                btg.add_thread(new boost::thread(boost::bind(&Search::searchHeteroXLinkedIndexedDB, this, boost::ref(run), boost::ref(proteome), boost::ref(revModMap), boost::ref(modMap2), i)));
            }
            btg.join_all();
        }
        cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalModCnt() << " mod pairs" << endl << endl; ////

        cout << "... searching REV-FWD/FWD-REV pairs" << endl; ////
        if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
            searchHeteroXLinkedIndexedDB(run, proteome, modMap, revModMap2, 0);
        } else {
            boost::thread_group btg;
            for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
                if (DigestOptions::verbose) {
                    cout << "... starting thread # " << i << endl; ////
                }
                btg.add_thread(new boost::thread(boost::bind(&Search::searchHeteroXLinkedIndexedDB, this, boost::ref(run), boost::ref(proteome), boost::ref(modMap), boost::ref(revModMap2), i)));
            }
            btg.join_all();
        }
        cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalModCnt() << " mod pairs" << endl << endl; ////
    //    cout << "JMDEBUG: MC 0: " << count_mc[0] << "; MC 1: " << count_mc[1] << "; MC 2: " << count_mc[2] << "; MC 3: " << count_mc[3] << "; MC 4: " << count_mc[4] << "; MC 5: " << count_mc[5] << endl; ////
    }
        
    JMUtil::procMemUsage(vm, rss); cout << "VM: " << vm << "; RSS: " << rss << endl; ////
    return;
}

// Parallelize by stepping every n_threads over the outer loop starting with cur_thread
// Could adjust the upperIt to make sure it's divisible by n_threads.  This turns out to be a bit complicated:
// >> upper_it_parallel = -n_th*( [0:(n_th -1)]' >= modulus ) - modulus*( [0:(n_th -1)]' < modulus ) + mod([0:(n_th -1)]', modulus)
// Alternatively, just iterate over the whole thing and check the remainder.  Checked to make sure it's cheap
void Search::searchHeteroXLinkedIndexedDB(Run &run, Proteome &proteome, vector<pair< double,xLinkPosition> > &massModMap1, vector<pair< double,xLinkPosition> > &massModMap2, int cur_thread_num) {
    if (massModMap1.size() == 0 || massModMap2.size() == 0) {
        cerr << endl << "ERROR: most likely missing reverse sequences" << endl;
        exit(0);
    }
    zeroOutCnts(); ////
    
    linkerMass = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
    double max_query_mass = run.msnSpectra.rbegin()->first + MSMass::precursorMass2Diff(run.msnSpectra.rbegin()->first) + DigestOptions::precursorSteps.back();
    double min_query_mass = run.msnSpectra.begin()->first - MSMass::precursorMass2Diff(run.msnSpectra.begin()->first) + DigestOptions::precursorSteps.front();

    double min_mod_pep_mass = massModMap1.begin()->first;
    double max_mod_pep_mass2 = max_query_mass - min_mod_pep_mass + MSMass::protonMass - linkerMass;

    // remove modified peptides with masses too large to ever be considered
    vector<pair< double,xLinkPosition> >::iterator upperIt2 = upper_bound(massModMap2.begin(), massModMap2.end(), max_mod_pep_mass2, sortMHxLinkPositionPair());
    if (cur_thread_num == 0) cout << distance(massModMap2.begin(), upperIt2) << " mods considered in massModMap2" << endl;

    double min_mod_pep_mass2 = massModMap2.begin()->first;
    double max_mod_pep_mass = max_query_mass - min_mod_pep_mass2 + MSMass::protonMass - linkerMass;

    // remove modified peptides with masses too large to ever be considered
    vector<pair< double,xLinkPosition> >::iterator upperIt1 = upper_bound(massModMap1.begin(), massModMap1.end(), max_mod_pep_mass, sortMHxLinkPositionPair());
    if (cur_thread_num == 0) cout << distance(massModMap1.begin(), upperIt1) << " mods considered in massModMap1" << endl;

    int local_cnt = 0;
    for (vector<pair< double,xLinkPosition> >::iterator mmMIt = massModMap1.begin(); mmMIt != upperIt1; mmMIt++) {
        local_cnt++;
        if ((local_cnt - cur_thread_num - 1) % DigestOptions::maxNumThreads != 0)
            continue;
        vector<pair< double,xLinkPosition> >::iterator lo_vec = lower_bound(massModMap2.begin(), upperIt2, (min_query_mass - mmMIt->first + MSMass::protonMass - linkerMass), sortMHxLinkPositionPair());
        vector<pair< double,xLinkPosition> >::iterator hi_vec = upper_bound(lo_vec,              upperIt2, (max_query_mass - mmMIt->first + MSMass::protonMass - linkerMass), sortMHxLinkPositionPair());
        for(vector<pair< double,xLinkPosition> >::iterator mmMIt2 = lo_vec; mmMIt2 != hi_vec; mmMIt2++) {
            // make sure combined mods don't exceed the max
            // +1 for the xlink as a mod ???
            if ((mmMIt->second.mIt._modIndexVector.size() + mmMIt2->second.mIt._modIndexVector.size() + 1) <= DigestOptions::nMaxVariableMods) {
                matchXlinkedIndexedDBRange(run, proteome, mmMIt, mmMIt2, cur_thread_num);
            }
        }
    }
    return;
}

void Search::searchHomoXLinkedIndexedDB(Run &run, Proteome &proteome, vector<pair< double,xLinkPosition> > &massModMap, int cur_thread_num) {
    if (massModMap.size() == 0) {
        cerr << endl << "ERROR: most likely missing reverse sequences" << endl;
        exit(0);
    }
    zeroOutCnts(); ////
    
    linkerMass = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
    double max_query_mass = run.msnSpectra.rbegin()->first + MSMass::precursorMass2Diff(run.msnSpectra.rbegin()->first) + DigestOptions::precursorSteps.back();
    double min_query_mass = run.msnSpectra.begin()->first - MSMass::precursorMass2Diff(run.msnSpectra.begin()->first) + DigestOptions::precursorSteps.front();

    double min_mod_pep_mass = massModMap.begin()->first;
    double max_mod_pep_mass = max_query_mass - min_mod_pep_mass + MSMass::protonMass - linkerMass;

    // remove modified peptides with masses too large to ever be considered
    //cout << "... removing peptides too large to ever be considered" << endl;
    vector<pair< double,xLinkPosition> >::iterator upperIt = upper_bound(massModMap.begin(), massModMap.end(), max_mod_pep_mass, sortMHxLinkPositionPair());
    if (cur_thread_num == 0) cout << distance(massModMap.begin(), upperIt) << " mods stored in massModMap" << endl;

    int local_cnt = 0;
    for (vector<pair< double,xLinkPosition> >::iterator mmMIt = massModMap.begin(); mmMIt != upperIt; mmMIt++) {
        local_cnt++;
        if ((local_cnt - cur_thread_num - 1) % DigestOptions::maxNumThreads != 0)
            continue;
        vector<pair< double,xLinkPosition> >::iterator lo_vec = lower_bound(mmMIt, upperIt, (min_query_mass - mmMIt->first + MSMass::protonMass - linkerMass), sortMHxLinkPositionPair());
        vector<pair< double,xLinkPosition> >::iterator hi_vec = upper_bound(mmMIt, upperIt, (max_query_mass - mmMIt->first + MSMass::protonMass - linkerMass), sortMHxLinkPositionPair());
        for(vector<pair< double,xLinkPosition> >::iterator mmMIt2 = lo_vec; mmMIt2 != hi_vec; mmMIt2++) {
            // make sure combined mods don't exceed the max
            // +1 for the xlink as a mod ???
            if ((mmMIt->second.mIt._modIndexVector.size() + mmMIt2->second.mIt._modIndexVector.size() + 1) <= DigestOptions::nMaxVariableMods) {
                matchXlinkedIndexedDBRange(run, proteome, mmMIt, mmMIt2, cur_thread_num);
            }
        }
    }
    return;
}

// Helper function shared between Homo- and Hetero- xlinked searches.  Makes it easier/cleaner to multithread
void Search::matchXlinkedIndexedDBRange(Run &run, Proteome &proteome, vector<pair< double,xLinkPosition> >::iterator mmMIt, vector<pair< double,xLinkPosition> >::iterator mmMIt2, int cur_thread_num) {
    double query_mass = mmMIt->first + mmMIt2->first - MSMass::protonMass + linkerMass; // add linker mass
    double mass_diff = MSMass::precursorMass2Diff(query_mass);

   for (int step = 0; step < DigestOptions::precursorSteps.size(); step++) {
        multimap<double, Spectrum>::iterator specItHi = run.msnSpectra.lower_bound(query_mass + DigestOptions::precursorSteps[step] + mass_diff);
        multimap<double, Spectrum>::iterator specItLo = run.msnSpectra.upper_bound(query_mass + DigestOptions::precursorSteps[step] - mass_diff);
        if (specItLo != specItHi) {
            int cnt = mmMIt->second.mIt.numUnmodifiedCleavageSites + mmMIt2->second.mIt.numUnmodifiedCleavageSites -
                      DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].xlinkedCleavageSites;
            if (cnt > DigestOptions::xMaxMissedCleavages) {
//                cout << "JMDEBUG: omitting\t"
//                        << mmMIt->second.mIt.getFlankedSeq() << " - " << mmMIt2->second.mIt.getFlankedSeq() << "\t"
//                        << mmMIt->first << "-" << mmMIt2->first << "\t"
//                        << mmMIt->second.mIt.numUnmodifiedCleavageSites << "-" << mmMIt2->second.mIt.numUnmodifiedCleavageSites << "\t"
//                        << DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].xlinkedCleavageSites << "\t"
//                        << DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkType << endl;
                continue;
            }
            vector< pair<int, int> > localized_xlink_pairs = getXLinkPairs(mmMIt->second.pos, mmMIt2->second.pos);
            for (vector< pair<int, int> >::iterator vpIt = localized_xlink_pairs.begin(); vpIt != localized_xlink_pairs.end(); vpIt++) {
                /** here could filter out special cases with terminal protein xlinks that lead to underestimation of of missed cleavages
                 * but appears to not be worth the trouble
                 */
//                count_mc[(mmMIt->second.mIt.pep->_missedCleavPosArray.size()) + (mmMIt2->second.mIt.pep->_missedCleavPosArray.size())] ++;
                SpectralMatch xsm(specItLo, specItHi, proteome, mmMIt->second.mIt, mmMIt2->second.mIt, (*vpIt).first, (*vpIt).second);
                psmCntVec[cur_thread_num] += distance(specItLo, specItHi); ////
            }
        }
    }
    modCntVec[cur_thread_num] ++; ////
    if (fmod(getTotalModCnt(), 1000000.0) == 0.0) ////
        cout << "." << flush; ////
    return;
}

void Search::searchLoopLinkedPeptides(Run& run, Proteome& proteome) {
    cout << "Search::searchLooplinkedPeptides" << "\t"
         << DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].join2string() << endl;
    
    DigestOptions::filterVarModsForXlink();

    if (DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].enzymeSpecificity != DigestOptions::enzymeSpecificity ||
        DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].maxMissedCleavages != DigestOptions::maxMissedCleavages) {
        
        // re-init options
        cout << "Re-initializing [Digest] options: " << endl;
        DigestOptions::enzymeSpecificity = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].enzymeSpecificity;
        DigestOptions::maxMissedCleavages = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].maxMissedCleavages + 
                                            DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].xlinkedCleavageSites;

        proteome.reIndexPeptides();
    }
    
    double vm, rss;
    zeroOutCnts();
    if (DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].isHomoXLink) {
        if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
            searchHomoLoopLinkedMT(run, proteome);
        } else {
            boost::thread_group btg;
            for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
                if (DigestOptions::verbose) {
                    cout << "... starting thread # " << i << endl; ////
                }
                btg.add_thread(new boost::thread(boost::bind(&Search::searchHomoLoopLinkedMT, this, boost::ref(run), boost::ref(proteome), i)));
            }
            btg.join_all();
        }
    } else {
        if (DigestOptions::maxNumThreads == 1) { // Parallelization overhead is pretty small (~5-10%) but it's nice to avoid it
            searchHeteroLoopLinkedMT(run, proteome);
        } else {
            boost::thread_group btg;
            for (int i = 0; i < DigestOptions::maxNumThreads; i++) {
                if (DigestOptions::verbose) {
                    cout << "... starting thread # " << i << endl; ////
                }
                btg.add_thread(new boost::thread(boost::bind(&Search::searchHeteroLoopLinkedMT, this, boost::ref(run), boost::ref(proteome), i)));
            }
            btg.join_all();
        }
    }
    
    cout << endl << getTotalPsmCnt() << " matching operations for " << getTotalModCnt() << " mods of " << getTotalPepCnt() << " peptides" << endl << endl; ////
    JMUtil::procMemUsage(vm, rss); cout << "VM: " << vm << "; RSS: " << rss << endl; ////
    return;
}

void Search::searchHomoLoopLinkedMT(Run &run, Proteome &proteome, int cur_thread_num) {
    linkerMass = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
    int local_cnt = 0;
    for (tr1::unordered_map<Peptide,int>::iterator pIt = proteome.peptideSet.begin(); pIt != proteome.peptideSet.end(); pIt++) {
        local_cnt++;
        if ((local_cnt - cur_thread_num - 1) % DigestOptions::maxNumThreads != 0)
            continue;

        if (!(*pIt).first._liveIndex)
            continue;

        if (fmod(getTotalPepCnt(), 100000.0) == 0.0) ////
            cout << "." << flush; ////

        ModIterator mIt((*pIt).first), lastMod;
        while (mIt != lastMod) {
            modCntVec[cur_thread_num] ++; ////
            mIt.prepOutput();
            set<int> pos = Search::getXLinkPositions(mIt, DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkedRes1);
            if (pos.size() > 1) {
                double query_mass = (*pIt).first._mPlusH + mIt.modMass + linkerMass;
                double mass_diff = MSMass::precursorMass2Diff(query_mass);

                for (int step = 0; step < DigestOptions::precursorSteps.size(); step++) {
                    multimap<double, Spectrum>::iterator specItHi = run.msnSpectra.lower_bound(query_mass + DigestOptions::precursorSteps[step] + mass_diff);
                    multimap<double, Spectrum>::iterator specItLo = run.msnSpectra.upper_bound(query_mass + DigestOptions::precursorSteps[step] - mass_diff);
                    if (specItLo != specItHi) {
                        vector< pair<int, int> > localized_xlink_pairs = getXLinkPairs(pos);
                        for (vector< pair<int, int> >::iterator vpIt = localized_xlink_pairs.begin(); vpIt != localized_xlink_pairs.end(); vpIt++) {
                            //cout << query_mass << "\t" << mIt.modSequence << "\t" << (*pIt).first.getProteinReference() << "\t" << (*vpIt).first << " " << (*vpIt).second << endl; ////

                            /* Make sure there's a stretch of amino acids outside the cross-linked segment (on either end) that could 
                             * produce fragments available for matching since bonds within the segment don't break. Use the standard 
                             * length min limit instead of an arbitrary new one.  Need to worry about situation where there are only 
                             * 1 or 2 fragments possible and if they match by chance it skews the lda distribution because the fraction 
                             * ions matched is high
                             */
                            if (max(vpIt->first + 1, (int) (mIt.pep->_sequence.length() - vpIt->second - 1)) < DigestOptions::minLength) {
                                continue;
                            }                            
                           
                            SpectralMatch xsm(specItLo, specItHi, proteome, pIt, mIt, (*vpIt).first, (*vpIt).second);
                            psmCntVec[cur_thread_num] += distance(specItLo, specItHi); ////
                        }
                    }
                }

            }
            mIt++;
        }

        pepCntVec[cur_thread_num] ++; ////
    }
    return;

}

void Search::searchHeteroLoopLinkedMT(Run &run, Proteome &proteome, int cur_thread_num) {
    linkerMass = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
    int local_cnt = 0;
    for (tr1::unordered_map<Peptide,int>::iterator pIt = proteome.peptideSet.begin(); pIt != proteome.peptideSet.end(); pIt++) {
        local_cnt++;
        if ((local_cnt - cur_thread_num - 1) % DigestOptions::maxNumThreads != 0)
            continue;

        if (!(*pIt).first._liveIndex)
            continue;
        
        if (fmod(getTotalPepCnt(), 100000.0) == 0.0) ////
            cout << "." << flush; ////

        ModIterator mIt((*pIt).first), lastMod;
        while (mIt != lastMod) {
            modCntVec[cur_thread_num] ++; ////
            mIt.prepOutput();
            set<int> pos1 = Search::getXLinkPositions(mIt, DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkedRes1);
            set<int> pos2 = Search::getXLinkPositions(mIt, DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkedRes2);
            if (pos1.size() > 0 && pos2.size() > 0) {
                double query_mass = (*pIt).first._mPlusH + mIt.modMass + linkerMass;
                double mass_diff = MSMass::precursorMass2Diff(query_mass);

                for (int step = 0; step < DigestOptions::precursorSteps.size(); step++) {
                    multimap<double, Spectrum>::iterator specItHi = run.msnSpectra.lower_bound(query_mass + DigestOptions::precursorSteps[step] + mass_diff);
                    multimap<double, Spectrum>::iterator specItLo = run.msnSpectra.upper_bound(query_mass + DigestOptions::precursorSteps[step] - mass_diff);
                    if (specItLo != specItHi) {
                        vector< pair<int, int> > localized_xlink_pairs = getXLinkPairs(pos1, pos2);
                        for (vector< pair<int, int> >::iterator vpIt = localized_xlink_pairs.begin(); vpIt != localized_xlink_pairs.end(); vpIt++) {
                            //cout << query_mass << "\t" << mIt.modSequence << "\t" << (*pIt).first.getProteinReference() << "\t" << (*vpIt).first << " " << (*vpIt).second << endl; ////
                            
                            /* Make sure there's a stretch of amino acids outside the cross-linked segment (on either end) that could 
                             * produce fragments available for matching since bonds within the segment don't break. Use the standard 
                             * length min limit instead of an arbitrary new one.  Need to worry about situation where there are only 
                             * 1 or 2 fragments possible and if they match by chance it skews the lda distribution because the fraction 
                             * ions matched is high
                             */
                            int pos1 = min(vpIt->first, vpIt->second); // for hetero-looplinks can't assume that xlink positions are in increasing order
                            int pos2 = max(vpIt->first, vpIt->second);
                            if (max(pos1 + 1, (int) (mIt.pep->_sequence.length() - pos2 - 1)) < DigestOptions::minLength) {
                                continue;
                            }
                            
                            SpectralMatch xsm(specItLo, specItHi, proteome, pIt, mIt, (*vpIt).first, (*vpIt).second);
                            psmCntVec[cur_thread_num] += distance(specItLo, specItHi); ////
                        }
                    }
                }

            }
            mIt++;
        }

        pepCntVec[cur_thread_num] ++; ////
    }
    return;

}

// localize hetero-linked pair possibilities
vector< pair<int,int> > Search::getXLinkPairs(set<int> &pos1, set<int> &pos2) {
    vector< pair<int,int> > loc_pairs;
    for (set<int>::iterator setIt = pos1.begin(); setIt != pos1.end(); setIt++) {
        for (set<int>::iterator setIt2 = pos2.begin(); setIt2 != pos2.end(); setIt2++) {
            //cout << "\t" << *setIt << "-" << *setIt2 << endl; ////
            loc_pairs.push_back(pair<int,int>(*setIt, *setIt2));
        }
    }
    return loc_pairs;
}

// localize homo-linked pair possibilities
vector< pair<int,int> > Search::getXLinkPairs(set<int> &pos1) {
    vector< pair<int,int> > loc_pairs;
    if (pos1.size() <= 1)
        return loc_pairs;

    for (set<int>::iterator setIt = pos1.begin(); setIt != pos1.end(); setIt++) {
        set<int>::iterator setIt2 = setIt;
        setIt2++;
        for (; setIt2 != pos1.end(); setIt2++) {
            loc_pairs.push_back(pair<int,int>(*setIt, *setIt2));
        }
    }
    return loc_pairs;
}


// Find all possible pairs of xlink locations, handling tryptic ends for any enzyme, mods, termini,
// Handle the protein terminal mods separately to ensure they don't get removed from consideration
set<int> Search::getXLinkPositions(ModIterator& mIt, string linked_res) {
    set<int> pos;
    if (linked_res.compare(".") == 0) { // handle entangled peptides
        pos.insert(-1);
        return pos;
    }
    
    bool n_term = false;
    bool c_term = false;

    // match potential xlink positions
    const boost::regex linked_res_re("[" + linked_res + "]");
    boost::sregex_iterator reIt(mIt.pep->_sequence.begin(), mIt.pep->_sequence.end(), linked_res_re);
    boost::sregex_iterator rEnd;
    while (reIt != rEnd) {
        pos.insert((*reIt++).position());
    }
    // Add protein terminal positions if relevant
    if (mIt.pep->_flankResNterm == '+' && boost::regex_match("+", linked_res_re)) {
        n_term = true;
    }
    if (mIt.pep->_flankResCterm == '-' && boost::regex_match("-", linked_res_re)) {
        c_term = true;
    }
    
    // subtract modified positions
    for (vector<int>::iterator vIt = mIt._modIndexVector.begin(); vIt != mIt._modIndexVector.end(); vIt++) {
        if (mIt.pep->_modSites[*vIt - 1].second->_symbol == '[') {
            c_term = false;
        } else if (mIt.pep->_modSites[*vIt - 1].second->_symbol == ']') {
            n_term = false;
        } else {
            pos.erase(mIt.pep->_modSites[*vIt - 1].first);
        }
    }
    // subtract protease-specific terminal sites
    if (mIt.pep->_flankResNterm != '+') {
        const boost::regex before_sites_re("^" + DigestOptions::enzymeBeforeSites.first);
        boost::match_results<string::const_iterator> match_results_before_sites;
        if (boost::regex_search (mIt.pep->_sequence, match_results_before_sites, before_sites_re)) {
            //cout << "JMDEBUG 1: " << match_results_before_sites.position() << " " << match_results_before_sites.length() << endl; ////
            for (int p = 0; p < match_results_before_sites.length(); p++) {
                pos.erase(match_results_before_sites.position() + p);
            }
        }
    }
    if (mIt.pep->_flankResCterm != '-') {
        const boost::regex after_sites_re(DigestOptions::enzymeAfterSites.first + "$");
        boost::match_results<string::const_iterator> match_results_after_sites;
        if (boost::regex_search (mIt.pep->_sequence, match_results_after_sites, after_sites_re)) {
            //cout << "JMDEBUG 2: " << match_results_after_sites.position() << " " << match_results_after_sites.length() << endl; ////
            for (int p = 0; p < match_results_after_sites.length(); p++) {
                pos.erase(match_results_after_sites.position() + p);
            }
        }
    }

    if (n_term) {
        pos.insert(0);
    }
    if (c_term) {
        pos.insert(mIt.pep->_sequence.length() - 1);
    }

    return pos;
}

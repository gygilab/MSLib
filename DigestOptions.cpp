#include "MSMass.hpp"
#include "Modification.hpp"
#include "DigestOptions.hpp"

// Default General Options
string DigestOptions::param_file = "";
bool DigestOptions::verbose = false;
int DigestOptions::maxNumThreads = 1;

// Default Digest Options
string DigestOptions::fastaFile = "undef";
string DigestOptions::fastaHeaderRegEx = ">(.+?)\\s.*?";
string DigestOptions::reverseTag = "##";
double DigestOptions::minMass = 400.0;
double DigestOptions::maxMass = 6000.0;
int DigestOptions::minLength = 5;
int DigestOptions::maxLength = 100;
string DigestOptions::enzymeName = "Trypsin";
string DigestOptions::enzymeRegEx = "(?<=[KR])(?!P)"; //Trypsin
int DigestOptions::enzymeSpecificity = 2; // 0 = NoEnzyme; 1 = Partial; 2 = Fully Specific
pair<string,bool> DigestOptions::enzymeBeforeSites = make_pair("", false);
pair<string,bool> DigestOptions::enzymeAfterSites = make_pair("", true);
bool DigestOptions::swapReverseInitMet = false;
bool DigestOptions::swapReverseCutSites = false;
int DigestOptions::maxMissedCleavages = 2;
bool DigestOptions::cleaveInitMet = false;

// Default Mass Options
string DigestOptions::varModProteinSpecificDelimiter = ":";
map<char, double> DigestOptions::elementMassMod;
map<char, double> DigestOptions::staticAAMod;
multimap<char, Modification> DigestOptions::variableAAMod; // CURRENTLY ONLY SUPPORTS 1 MOD PER AA. NEED TO GET CLEVER TO DO MULTIPLES
string DigestOptions::variablyModAas = "";
string DigestOptions::orderedModSymbols = "*#@^~$"; // '[' and ']' are reserved for C- and N- termini
string DigestOptions::orderedNTerminalModSymbols = "]})";
string DigestOptions::orderedCTerminalModSymbols = "[{(";
bool DigestOptions::lowerCaseMods = 0;

// Default Spectrum Options
int DigestOptions::windowSize = 100;
int DigestOptions::peakDepth = 9;
bool DigestOptions::deIsotope = true;
bool DigestOptions::printDeisoSpectra = false;
double DigestOptions::filterLowMzMass = 0.0;
bool DigestOptions::filterCommonCorrel = false;
int DigestOptions::minRequiredPeaks = 5;
bool DigestOptions::removePrecursorNL = true;
bool DigestOptions::removeWaterNLoss = true;
bool DigestOptions::removePhosphoNLoss = false;
int DigestOptions::maxTheoFragCharge = 4;
double DigestOptions::maxMS2ItyCutoffForTIC = 2000;

// Default Search Options
double DigestOptions::precursorTol = 25;
bool DigestOptions::precursorTolRelative = 1;
double DigestOptions::fragmentTol = 0.4;
int DigestOptions::fragmentTolInt = MSMass::double2int(0.4);
bool DigestOptions::fragmentTolRelative = 0;
int DigestOptions::nMaxVariableMods = 3;
bool DigestOptions::noModSiteCleavage = true;
string DigestOptions::searchDir = "./";
string DigestOptions::outputFile = "xlink.out";
int DigestOptions::maxHitRank = 50;
int DigestOptions::maxHitOutputRank = 1;
string DigestOptions::scoringFunction = "Binomial";
vector<double> DigestOptions::precursorSteps {0.0};
bool DigestOptions::indexAllPeptides = 1;
bool DigestOptions::indexSeenSequences = true;

// Default X-link Search Options
string DigestOptions::searchType = "Peptide";
double DigestOptions::linkerMass = 0.0;
string DigestOptions::linkType = "";
string DigestOptions::linkedRes1 = "";
string DigestOptions::linkedRes2 = "";
vector<SearchOpt> DigestOptions::multiSearchOptions;
int DigestOptions::currentSearchOpt = 0;
int DigestOptions::xFilterTopProts = 50;
double DigestOptions::xFilterOMDiff = -1.0;
double DigestOptions::xFilterFragDevZ = 1000;
int DigestOptions::xFilterMinMatched = 0;
double DigestOptions::xFilterMinFracMatched = 0;
int DigestOptions::xMaxMissedCleavages = 3;
bool DigestOptions::skipREVREVsearch = true;
bool DigestOptions::skipREVFWDsearch = false;
bool DigestOptions::skipHalfREVFWDsearch = true;
int DigestOptions::loopFilterMinMatched = 6;

DigestOptions::DigestOptions(int ac, char *av[]) {
    parseProgramOptions(ac, av);
    printParams();
}

void DigestOptions::parseNumThreads (const int &input_num_th) {
    int numCPU = sysconf( _SC_NPROCESSORS_ONLN );
    if (input_num_th < 1)
        maxNumThreads = 1;
    else if (input_num_th > numCPU)
        maxNumThreads = numCPU;
    else
        maxNumThreads = input_num_th;
    return;
}

void DigestOptions::parsePrecursorTolerance(const string &tolerance_input) {
    string buf;
    stringstream ss(tolerance_input);
    vector<string> tokens; // vector to hold tokens
    while (ss >> buf) {
        tokens.push_back(buf);
    }
    if (tokens.size() != 2) {
        cerr << "ERROR: Precursor tolerance specified incorrectly" << endl;
        exit(0);
    } else {
        precursorTol = atof(tokens[0].c_str());

        std::transform(tokens[1].begin(), tokens[1].end(), tokens[1].begin(), ::toupper);
        if (tokens[1].compare("TH") == 0 || tokens[1].compare("DA") == 0) {
            precursorTolRelative = 0;
        } else if (tokens[1].compare("PPM") == 0) {
            precursorTolRelative = 1;
        } else {
            cerr << "ERROR: Precursor tolerance specified incorrectly" << endl;
            exit(0);
        }
    }

    if ((precursorSteps.front() < 0.0 || precursorSteps.back() > 0.0) &&
        precursorTolRelative == 0 && precursorTol > 1.0) {
        cerr << endl << "WARNING: precursorStepRange settings might not make sense when paired with current precursorTol settings" << endl;
        exit(0);
    }
}

void DigestOptions::parseFragmentTolerance(const string &tolerance_input) {
    string buf;
    stringstream ss(tolerance_input);
    vector<string> tokens; // vector to hold tokens
    while (ss >> buf) {
        tokens.push_back(buf);
    }
    if (tokens.size() != 2) {
        cerr << "ERROR: Fragment tolerance specified incorrectly" << endl;
        exit(0);
    } else {
        fragmentTol = atof(tokens[0].c_str());
        fragmentTolInt = MSMass::double2int(fragmentTol);

        std::transform(tokens[1].begin(), tokens[1].end(), tokens[1].begin(), ::toupper);
        if (tokens[1].compare("TH") == 0 || tokens[1].compare("DA") == 0) {
            fragmentTolRelative = 0;
            return;
        } else if (tokens[1].compare("PPM") == 0) {
            fragmentTolRelative = 1;
            fragmentTolInt /= 1000.0; // average
            return;
        } else {
            cerr << "ERROR: Fragment tolerance specified incorrectly" << endl;
            exit(0);
        }
    }
}

// Parse space-separated mass/residue:protein pairs like: 15.9949 M 79.9663 STY:ALBU_HUMAN:OVAL_CHICK
void DigestOptions::parseVariableMods(const vector<string> &var_mods) {
    //copy(var_mods.begin(), var_mods.end(), ostream_iterator<string>(cout, " ")); cout << endl; ////

    vector<string> parsed_var_mods;
    boost::regex re("\\s+");

    // first tokenize the whole string
    for (vector<string>::const_iterator vecIt = var_mods.begin(); vecIt != var_mods.end(); vecIt++) {
        boost::sregex_token_iterator i((*vecIt).begin(), (*vecIt).end(), re, -1);
        boost::sregex_token_iterator j;

        while (i != j) {
//            cout << *i << endl; ////
            parsed_var_mods.push_back(*i);
            i++;
        }
    }

    // there should be pairs
    if (parsed_var_mods.size() % 2 != 0) {
        cerr << "ERROR: Each variable mod has to contain residue(s) and mass information" << var_mods.size() << endl;
        exit(EXIT_FAILURE);
        return;
    }

    // sort out the residues from the masses
    vector<string>::const_iterator vecIt = parsed_var_mods.begin();
    int mod_cnt = 0;
    int n_term_mod_cnt = 0;
    int c_term_mod_cnt = 0;
    while (vecIt != parsed_var_mods.end()) {
        // 1
        string mass = *vecIt;
        // TODO: regex check
        double m = atof(mass.c_str());
        vecIt++;

        // 2
        string residues_prots = *vecIt;
        string residues;
        string prots = "";
        // TODO: regex check
        
        // 3
        size_t found = residues_prots.find_first_of(varModProteinSpecificDelimiter);
        if (found != string::npos) {
            vector<string> res_prot_vec;
            residues = residues_prots.substr(0, found);
            prots = residues_prots.substr(found + 1);
        } else {
            residues = residues_prots;
        }

        char s;
        if (mod_cnt <= (orderedModSymbols.size() - 1)) {
            s = orderedModSymbols[mod_cnt];
        } else {
            cerr << "ERROR: Hard maximum number of variable mods exceeded" << endl;
            exit(EXIT_FAILURE);
        }
        Modification mod(m, s, prots);

        for (string::const_iterator sIt = residues.begin(); sIt != residues.end(); sIt++) {
            /* 
             * Use separate set of mod symbols for termini. 
             * Protein terminal mods are keyed on '+' and '-' in the variableAAMod map but the actual modification
             * contains the '[' and ']' types of symbols ensuring that those are the ones that get displayed in output but
             * still allow addressing protein termini directly where needed
             * Number of allowed mods currently limited by reserved symbols
             */
            if (*sIt == '[' || *sIt == '-') {
                char c;

                if (c_term_mod_cnt <= (orderedCTerminalModSymbols.size() - 1)) {
                    c = orderedCTerminalModSymbols[c_term_mod_cnt];
                } else {
                    cerr << "ERROR: Hard maximum number of C-terminal variable mods exceeded" << endl;
                    exit(EXIT_FAILURE);
                }
                Modification c_terminal_mod(m, c, prots);
                variableAAMod.insert(pair<char,Modification>(*sIt, c_terminal_mod));
                c_term_mod_cnt++;
            } else if (*sIt == ']' || *sIt == '+') {
                char n;

                if (n_term_mod_cnt <= (orderedNTerminalModSymbols.size() - 1)) {
                    n = orderedNTerminalModSymbols[n_term_mod_cnt];
                } else {
                    cerr << "ERROR: Hard maximum number of N-terminal variable mods exceeded" << endl;
                    exit(EXIT_FAILURE);
                }
                Modification n_terminal_mod(m, n, prots);
                variableAAMod.insert(pair<char,Modification>(*sIt, n_terminal_mod));
                n_term_mod_cnt++;
            } else {
                variableAAMod.insert(pair<char,Modification>(*sIt, mod));
            }
        }
        mod_cnt++;
        
        vecIt++;
    }
   
    updateVariablyModAas();

    return;
}

// make a concatenated string of variablyModAas to check against for each peptide
void DigestOptions::updateVariablyModAas() {
    variablyModAas = "";
    for (multimap<char, Modification>::iterator mmapIt = variableAAMod.begin(); mmapIt != variableAAMod.end(); mmapIt++) {
        if (isalpha((*mmapIt).first)) // Non-alphabetic (terminal) residues casuse problems for regex downstream - don't need them here anyway
            variablyModAas += (*mmapIt).first;
    }
    return;
}

/* Remove mods from variableAAMod that are not xlink residues (or Met oxidation)
 * Allows adding mods for regular peptide search that will then be filtered out here
 */
void DigestOptions::filterVarModsForXlink() {
//    string linked_res = "M"; // don't remove Met oxidation mod from xlink consideration
//    for (int so_num = 0; so_num < DigestOptions::multiSearchOptions.size(); so_num++) {
//        linked_res += DigestOptions::multiSearchOptions[so_num].linkedRes1 + DigestOptions::multiSearchOptions[so_num].linkedRes2;
//    }
//    for (multimap<char, Modification>::iterator mmapIt = DigestOptions::variableAAMod.begin(); mmapIt != DigestOptions::variableAAMod.end();) {
//        size_t found = linked_res.find(mmapIt->first);
//        if (linked_res.find(mmapIt->first) == string::npos || mmapIt->second._mass < 15.9) { // hack to keep Met ox, but filter out hydroxamic acid
//            cout << "... removing variable mod from consideration: " << mmapIt->first << endl; ////
//            DigestOptions::variableAAMod.erase(mmapIt++);
//        } else {
//             ++mmapIt;
//        }
//    }
//    DigestOptions::updateVariablyModAas();
//    return;
}

// Parse space-separated mass/residue pairs like: 57.02146374 C
void DigestOptions::parseStaticMods(const vector<string> &static_mods) {
//    copy(static_mods.begin(), static_mods.end(), ostream_iterator<string>(cout, " ")); cout << endl; ////

    vector<string> parsed_static_mods;
    boost::regex re("\\s+");

    // first tokenize the whole string
    for (vector<string>::const_iterator vecIt = static_mods.begin(); vecIt != static_mods.end(); vecIt++) {
        boost::sregex_token_iterator i((*vecIt).begin(), (*vecIt).end(), re, -1);
        boost::sregex_token_iterator j;

        while (i != j) {
//            cout << *i << endl; ////
            parsed_static_mods.push_back(*i);
            i++;
        }
    }

    // there should be pairs
    if (parsed_static_mods.size() % 2 != 0) {
        cerr << "ERROR: Each variable mod has to contain residue(s) and mass information" << static_mods.size() << endl;
        exit(EXIT_FAILURE);
        return;
    }

    // sort out the residues from the masses
    vector<string>::const_iterator vecIt = parsed_static_mods.begin();
    while (vecIt != parsed_static_mods.end()) {
        // 1
        string mass = *vecIt;
        // TODO: regex check
        double m = atof(mass.c_str());
        vecIt++;

        // 2
        string residues = *vecIt;
        // TODO: regex check

        for (string::const_iterator sIt = residues.begin(); sIt != residues.end(); sIt++) {
            pair<map<char,double>::iterator, bool> ret_pair = staticAAMod.insert(pair<char,double>(*sIt, m));
            if (!ret_pair.second) {
                cerr << "WARNING: Multiple static modifications specified for '" << *sIt << "'" << endl;
                exit(EXIT_FAILURE);
            }
        }

        vecIt++;
    }
    return;
}

void DigestOptions::parseElementMods(const vector<string> &element_mods) {
//    copy(element_mods.begin(), element_mods.end(), ostream_iterator<string>(cout, " ")); cout << endl; ////

    vector<string> parsed_element_mods;
    boost::regex re("\\s+");

    // first tokenize the whole string
    for (vector<string>::const_iterator vecIt = element_mods.begin(); vecIt != element_mods.end(); vecIt++) {
        boost::sregex_token_iterator i((*vecIt).begin(), (*vecIt).end(), re, -1);
        boost::sregex_token_iterator j;

        while (i != j) {
//            cout << *i << endl; ////
            parsed_element_mods.push_back(*i);
            i++;
        }
    }

    // there should be pairs
    if (parsed_element_mods.size() % 2 != 0) {
        cerr << "ERROR: Each variable mod has to contain element(s) - mass pairs" << element_mods.size() << endl;
        exit(EXIT_FAILURE);
        return;
    }

    // sort out the elements from the masses
    vector<string>::const_iterator vecIt = parsed_element_mods.begin();
    while (vecIt != parsed_element_mods.end()) {
        // 1
        string mass = *vecIt;
        // TODO: regex check
        double m = atof(mass.c_str());
        vecIt++;

        // 2
        string elements = *vecIt;
        // TODO: regex check

        for (string::const_iterator sIt = elements.begin(); sIt != elements.end(); sIt++) {
            pair<map<char,double>::iterator, bool> ret_pair = elementMassMod.insert(pair<char,double>(*sIt, m));
            if (!ret_pair.second) {
                cerr << "WARNING: Multiple element modifications specified for '" << *sIt << "'" << endl;
                exit(EXIT_FAILURE);
            }
        }

        vecIt++;
    }
    return;
}

// Populate the vector<SearchOpt> multiSearchOptions
void DigestOptions::parseAddSearch(const vector<string> &add_search) {
    for (vector<string>::const_iterator vIt = add_search.begin(); vIt != add_search.end(); vIt++) {
        vector<string> tokens;
        istringstream iss(*vIt);
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tokens));

        if (tokens.size() == 1) {
            searchType = tokens[0].c_str();
            SearchOpt so = SearchOpt(maxMissedCleavages, enzymeSpecificity, searchType, linkType, linkerMass, linkedRes1, linkedRes2);
            DigestOptions::multiSearchOptions.push_back(so);
        } else if (tokens.size() == 7) {
            int max_mis_cleav = atoi(tokens[0].c_str());
            int enz_spec = atoi(tokens[1].c_str());
            string search_type = tokens[2].c_str();
            string link_type = tokens[3].c_str();
            double linker_mass = atof(tokens[4].c_str());
            string linked_res1 = tokens[5].c_str();
            string linked_res2 = tokens[6].c_str();

            // TODO: sanity check to make sure theat res1 and res2 don't actually share any residues
            // Also, sort them so that if this is a homo-link, a simple string comparison will determine that
            sort(linked_res1.begin(), linked_res1.end());
            sort(linked_res2.begin(), linked_res2.end());

            SearchOpt so = SearchOpt(max_mis_cleav, enz_spec, search_type, link_type, linker_mass, linked_res1, linked_res2);
            DigestOptions::multiSearchOptions.push_back(so);
        } else {
            cerr << "ERROR: searchType specified incorrectly" << endl;
            exit(0);
        }
    }
    // Check if no searchType was specified, add default
    if (add_search.size() == 0) {
        cerr << "WARNING: assuming default searchType ('Peptide')" << endl;
        SearchOpt so = SearchOpt(maxMissedCleavages, enzymeSpecificity, searchType, linkType, linkerMass, linkedRes1, linkedRes2);
        DigestOptions::multiSearchOptions.push_back(so);
    }
    return;
}

// Break down the regex into SEQUEST-style params to be used later in some cases
void DigestOptions::parseEnzymeRegEx(string enzyme_regex) {
    //cout << "DigestOptions::parseEnzymeRegEx" << endl; ////
    enzymeRegEx = enzyme_regex;
    string before_sites = "";
    string before_sites_regex = "";
    int before_sites_cut = 0;
    string after_sites = "";
    string after_sites_regex = "";
    int after_sites_cut = 0;
    boost::regex re_positive_look("\\?[<]{0,1}="); // positive lookahead/behind
    boost::regex re_negative_look("\\?[<]{0,1}!"); // negative lookahead/behind

    boost::regex re_after_sites("\\?[<][=!]([[]{0,1}[A-Z\\-]+[]]{0,1})");
    boost::match_results<string::const_iterator> regex_result_after_sites;
    if (boost::regex_search(enzyme_regex, regex_result_after_sites, re_after_sites)) {
        after_sites_regex = regex_result_after_sites[0];
        after_sites = regex_result_after_sites[1];
        
//        copy(regex_result_after_sites.begin(), regex_result_after_sites.end(), ostream_iterator<string>(cout, " ")); cout << endl; ////

        if (boost::regex_search(after_sites_regex, re_positive_look)) {
            after_sites_cut = 1;
        } else if (boost::regex_search(after_sites_regex, re_negative_look)) {
            after_sites_cut = 0;
        } else {
            cerr << "Error parsing P1 site cleavage rules" << endl;
            exit(0);
        }
        // TODO: need to use iterators here to make sure not more than one P1 site is defined
    }

    boost::regex re_before_sites("\\?[=!]([[]{0,1}[A-Z\\-]+[]]{0,1})");
    boost::match_results<string::const_iterator> regex_result_before_sites;
    if (boost::regex_search(enzyme_regex, regex_result_before_sites, re_before_sites)) {
        before_sites_regex = regex_result_before_sites[0];
        before_sites = regex_result_before_sites[1];

//        copy(regex_result_before_sites.begin(), regex_result_before_sites.end(), ostream_iterator<string>(cout, " ")); cout << endl; ////

        if (boost::regex_search(before_sites_regex, re_positive_look)) {
            before_sites_cut = 1;
        } else if (boost::regex_search(before_sites_regex, re_negative_look)) {
            before_sites_cut = 0;
        } else {
            cerr << "Error parsing P1' site cleavage rules" << endl;
            exit(0);
        }
        // TODO: need to use iterators here to make sure not more than one P1' site is defined
    }

//    cout << after_sites << " " << after_sites_cut << "\t\t" << before_sites << " " << before_sites_cut << endl; ////

    if ((before_sites.length() > 0 || after_sites.length() > 0) && (before_sites_cut > 0 || after_sites_cut > 0)) {
        enzymeBeforeSites = make_pair(before_sites, before_sites_cut);
        enzymeAfterSites = make_pair(after_sites, after_sites_cut);
    } else {
        cerr << "Error parsing proteolytic cleavage rules" << endl;
        exit(0);
    }
    return;
}

void DigestOptions::parseEnzymeSpecificity (const int &num_specific_ends) {
    if (num_specific_ends < 0 || num_specific_ends > 2) {
        enzymeSpecificity = 2;
    } else {
        enzymeSpecificity = num_specific_ends;
    }
    return;
}

void DigestOptions::parsePrecursorStepRange(const string &step_range) {
    double buf;
    stringstream ss(step_range);
    vector<double> tokens; // vector to hold tokens
    while (ss >> buf) {
        tokens.push_back(buf);
    }
    sort(tokens.begin(), tokens.end());
    if (tokens.front() > 0 || tokens.back() < 0) {
        cerr << endl << "ERROR: precursorStepRange specified incorrectly: range must include 0" << endl;
        exit(0);
    }
    if ((tokens.front() < 0 || tokens.back() > 0) &&
        precursorTolRelative == 0 && precursorTol > 1.0) {
        cerr << endl << "WARNING: these precursorStepRange settings might not make sense when paired with current precursorTol settings" << endl;
        exit(0);
    }
    precursorSteps = tokens;
}

void DigestOptions::printParams() {
    cout << "verbose                   = " << verbose << endl;
    cout << "maxNumThreads             = " << maxNumThreads << endl;
    cout << "param_file                = " << param_file << endl;

    cout << endl;

    cout << "Digest.FastaFile          = " << fastaFile << endl;
    cout << "Digest.FastaHeaderRegEx   = " << fastaHeaderRegEx << endl;
    cout << "Digest.MinMass            = " << minMass << endl;
    cout << "Digest.MaxMass            = " << maxMass << endl;
    cout << "Digest.MinLength          = " << minLength << endl;
    cout << "Digest.MaxLength          = " << maxLength << endl;
    cout << "Digest.EnzymeName         = " << enzymeName << endl;
    cout << "Digest.EnzymeRegEx        = " << enzymeRegEx << endl;
    cout << "                            " << enzymeName << " " << enzymeAfterSites.first << " " << enzymeAfterSites.second << "\t"
                                                                << enzymeBeforeSites.first << " " << enzymeBeforeSites.second << endl;
    
    cout << "Digest.swapReverseInitMet = " << swapReverseInitMet << endl;
    cout << "Digest.swapReverseCutSites= " << swapReverseCutSites << endl;
    cout << "Digest.MaxMissedCleavages = " << maxMissedCleavages << endl;
    cout << "Digest.enzymeSpecificity  = " << enzymeSpecificity << endl;
    cout << "Digest.cleaveInitMet      = " << cleaveInitMet << endl;

    cout << endl;

    for (map<char, double>::const_iterator eMMiT = elementMassMod.begin(); eMMiT != elementMassMod.end(); eMMiT++) {
        if ((*eMMiT).second != 0) {
            cout << "Mass.element " << (*eMMiT).first << "       = " << (*eMMiT).second << endl;
        }
    }

    for (map<char, double>::const_iterator sAAMiT = staticAAMod.begin(); sAAMiT != staticAAMod.end(); sAAMiT++) {
        if ((*sAAMiT).second != 0) {
            cout << "Mass.staticMod " << (*sAAMiT).first << " = " << (*sAAMiT).second << endl;
        }
    }

    for (multimap<char, Modification>::const_iterator vAAMiT = variableAAMod.begin(); vAAMiT != variableAAMod.end(); vAAMiT++) {
        if (((*vAAMiT).second)._mass != 0) {
            cout << "Mass.variableMod " << (*vAAMiT).first << " = " << ((*vAAMiT).second)._mass << "\t" << ((*vAAMiT).second)._symbol
                 << "\t" << boost::algorithm::join(vAAMiT->second._proteinIDs, ":") << endl;
        }
    }

    cout << endl;
    cout << "Spectrum.windowSize         = " << windowSize << endl;
    cout << "Spectrum.peakDepth          = " << peakDepth << endl;
    cout << "Spectrum.deIsotope          = " << deIsotope << endl;
    cout << "Spectrum.printDeisoSpectra  = " << printDeisoSpectra << endl;

    cout << "Spectrum.filterLowMzMass    = " << filterLowMzMass << endl;
    cout << "Spectrum.filterCommonCorrel = " << filterCommonCorrel << endl;
    cout << "Spectrum.removeWaterNLoss   = " << removeWaterNLoss << endl;
    cout << "Spectrum.removePhosphoNLoss = " << removePhosphoNLoss << endl;
    cout << "Spectrum.maxTheoFragCharge = "  << maxTheoFragCharge << endl;

    cout << endl;
    cout << "Search.precursorTol         = " << precursorTol << endl;
    cout << "Search.precursorTolRelative = " << precursorTolRelative << endl;
    cout << "Search.fragmentTol          = " << fragmentTol << endl;
    cout << "Search.fragmentTolInt       = " << fragmentTolInt << endl;
    cout << "Search.fragmentTolRelative  = " << fragmentTolRelative << endl;
    cout << "Search.nMaxVariableMods     = " << nMaxVariableMods << endl;
    cout << "Search.noModSiteCleavage    = " << noModSiteCleavage << endl;
    cout << "Search.searchDir            = " << searchDir << endl;
//    cout << "Search.searchType           = " << searchType << endl;
    cout << "Search.outputFile           = " << outputFile << endl;
    cout << "Search.maxHitOutputRank     = " << maxHitOutputRank << endl;
    
    cout << "Search.scoringFunction      = " << scoringFunction << endl;
    cout << fixed << setprecision(3) << 
            "Search.precursorStepRange   = " << precursorSteps.front() << " " << precursorSteps.back() << endl; // TODO: output the whole vector
    cout << "Search.indexSeenSequences   = " << indexSeenSequences << endl;
    cout << "Search.xFilterTopProts      = " << xFilterTopProts << endl;
//    cout << "Search.xFilterOMDiff        = " << xFilterOMDiff << endl;
//    cout << "Search.xFilterFragDevZ      = " << xFilterFragDevZ << endl;
//    cout << "Search.xFilterMinMatched    = " << xFilterMinMatched << endl;
//    cout << "Search.xFilterMinFracMatched= " << xFilterMinFracMatched << endl;
    cout << "Search.xMaxMissedCleavages  = " << xMaxMissedCleavages << endl;
    cout << "Search.skipREVREVsearch     = " << skipREVREVsearch << endl;
//    cout << "Search.skipREVFWDsearch     = " << skipREVFWDsearch << endl;
    cout << "Search.skipHalfREVFWDsearch = " << skipHalfREVFWDsearch << endl;

    cout << endl;
    for (vector<SearchOpt>::iterator vIt = DigestOptions::multiSearchOptions.begin(); vIt != DigestOptions::multiSearchOptions.end(); vIt++) {
        cout << "Search.addSearch        = " << (*vIt).join2string() << endl;
    }
}

void DigestOptions::parseProgramOptions(int ac, char *av[]) {
    try {
        // Declare a group of options that will be
        // allowed only on command line
        po::options_description general("General options");
        general.add_options()
                ("version",                                                                                            "print current version")
                ("verbose,v",                 po::value<bool>(&verbose)->default_value(verbose),                       "verbose output")
                ("help,h",                                                                                             "print this help message")
                ("param-file,p",              po::value<string> (&param_file),                                         "digest parameter file")
                ("maxNumThreads,t",           po::value<int>()->notifier(
                                                     boost::bind(&DigestOptions::parseNumThreads, this, _1) ),         "maximum threads to use [1]")
                ;

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description digest("Digest Options");
        digest.add_options()
                ("Digest.fastaFile,d",        po::value<string>(&fastaFile)       ->default_value(fastaFile),          "fasta file")
                ("Digest.fastaHeaderRegEx",   po::value<string>(&fastaHeaderRegEx)->default_value(fastaHeaderRegEx),   "fasta header id parse regex")
                ("Digest.minMass",            po::value<double>(&minMass)         ->default_value(minMass),            "minimum peptide mass")
                ("Digest.maxMass",            po::value<double>(&maxMass)         ->default_value(maxMass),            "maximum peptide mass")
                ("Digest.minLength",          po::value<int>(&minLength)          ->default_value(minLength),          "minimum peptide length")
                ("Digest.maxLength",          po::value<int>(&maxLength)          ->default_value(maxLength),          "maximum peptide length")
                ("Digest.enzymeName",         po::value<string>(&enzymeName)      ->default_value(enzymeName),         "enzyme name")
                ("Digest.enzymeRegEx,e",      po::value<string>()->notifier(
                                                     boost::bind(&DigestOptions::parseEnzymeRegEx, this, _1) ),        "enzyme regex")
                ("Digest.maxMissedCleavages", po::value<int>(&maxMissedCleavages) ->default_value(maxMissedCleavages), "maxiumum missed cleavages")
                ("Digest.enzymeSpecificity",  po::value<int>()->notifier(
                                                     boost::bind(&DigestOptions::parseEnzymeSpecificity, this, _1) ),  "enzyme specificity: 2 = full; 1 = partial; 0 = none")
                ("Digest.swapReverseInitMet", po::value<bool>(&swapReverseInitMet)->default_value(swapReverseInitMet), "swap N-terminal Met for reverse sequences to avoid FWD bias")
                ("Digest.swapReverseCutSites",po::value<bool>(&swapReverseCutSites)->default_value(swapReverseCutSites),"swap each enzyme cut site eg [KR] with previous residue to avoid exact mass matches for reverse sequences")
                ("Digest.cleaveInitMet",      po::value<bool>(&cleaveInitMet)     ->default_value(cleaveInitMet),      "cleave initiator Met when enzyme specificity == 2")
                ;

        po::options_description mass("Mass Options");
        mass.add_options()
                ("Mass.elementMods",     po::value< vector<string> >()->composing()->notifier(
                                            boost::bind(&DigestOptions::parseElementMods, this, _1) ),  "space-separated mass/element pairs, e.g. 1.0 N")
                ("Mass.staticMods",      po::value< vector<string> >()->composing()->notifier(
                                            boost::bind(&DigestOptions::parseStaticMods, this, _1) ),   "space-separated mass/residue pairs, e.g. 57.02146374 C")
                ("Mass.variableMods",    po::value< vector<string> >()->composing()->notifier(
                                            boost::bind(&DigestOptions::parseVariableMods, this, _1) ), "space-separated mass/residue pairs, e.g. 15.9949 M 79.9663 STY")
                ;

        po::options_description spectrum("Spectrum Options");
        spectrum.add_options()
                ("Spectrum.windowSize",         po::value<int>(&windowSize)          ->default_value(windowSize),         "window size for filtering MS/MS spectra")
                ("Spectrum.peakDepth",          po::value<int>(&peakDepth)           ->default_value(peakDepth),          "peak depth for filtering MS/MS spectra)")
                ("Spectrum.deIsotope",          po::value<bool>(&deIsotope)          ->default_value(deIsotope),          "de-isotope MS/MS spectra?")
                ("Spectrum.printDeisoSpectra",  po::value<bool>(&printDeisoSpectra)  ->default_value(printDeisoSpectra),  "print de-isotoped MS/MS spectra?")
                ("Spectrum.filterLowMzMass",    po::value<double>(&filterLowMzMass)  ->default_value(filterLowMzMass),    "filter low m/z below this mass")
                ("Spectrum.filterCommonCorrel", po::value<bool>(&filterCommonCorrel) ->default_value(filterCommonCorrel), "filter low m/z below this mass")
                ("Spectrum.removeWaterNLoss",   po::value<bool>(&removeWaterNLoss)   ->default_value(removeWaterNLoss),   "remove NL of waters from MS/MS spectra?")
                ("Spectrum.removePhosphoNLoss", po::value<bool>(&removePhosphoNLoss) ->default_value(removePhosphoNLoss), "remove phospho NL from MS/MS spectra?")
                ("Spectrum.maxTheoFragCharge",  po::value<int>(&maxTheoFragCharge)   ->default_value(maxTheoFragCharge),  "max theoretical MS/MS fragment charge")
                ;

        po::options_description search("Search Options");
        search.add_options()
                ("Search.precursorTol",     po::value<string>()->notifier(
                                                boost::bind(&DigestOptions::parsePrecursorTolerance, this, _1) ), "precursor ion tolerance ('X Th' or 'X ppm')")
                ("Search.fragmentTol",      po::value<string>()->notifier(
                                                boost::bind(&DigestOptions::parseFragmentTolerance, this, _1) ),  "fragment ion tolerance ('X Th' or 'X ppm')")
                ("Search.nMaxVariableMods", po::value<int>(&nMaxVariableMods)  ->default_value(nMaxVariableMods), "maximum number of variable mods")
                ("Search.noModSiteCleavage",po::value<bool>(&noModSiteCleavage)->default_value(noModSiteCleavage),"Terminal cleavage positions excluded from modifications. (default: True)")
                ("Search.searchDir,s",      po::value<string>(&searchDir)      ->default_value(searchDir),        "directory with dta files")
                ("Search.outputFile,o",     po::value<string>(&outputFile)     ->default_value(outputFile),       "file to output")
                ("Search.maxHitOutputRank", po::value<int>(&maxHitOutputRank)  ->default_value(maxHitOutputRank), "number of hits per spectrum to output")
                ("Search.precursorStepRange",po::value<string>()->notifier(
                                                boost::bind(&DigestOptions::parsePrecursorStepRange, this, _1) ), "precursor interval step range, e.g. '-1 1' or '-3 2' (default '0 0')")
                ("Search.scoringFunction",  po::value<string>(&scoringFunction)->default_value(scoringFunction),  "scoring function")
                // deal with multi-options
                ("Search.searchType",       po::value< vector<string> >()->composing()->notifier(
                                                boost::bind(&DigestOptions::parseAddSearch, this, _1) ),          "space-separated search options with the format: [maxMissedCleavages enzymeSpecificity SearchType linkType linkerMass linkedRes1 linkedRes2], e.g. 'XLink CC -2.0 C C', 'XLink 0_DE -18. DE- K+'")
                ("Search.indexSeenSequences",po::value<bool>(&indexSeenSequences)->default_value(true),           "index peptides (1-default/0), turn off to save memory (for very large databases), lose speed")
                ("Search.xFilterTopProts",   po::value<int>(&xFilterTopProts)->default_value(xFilterTopProts),    "number of first (forward or reverse) sequences to consider for xlink search")
//                ("Search.xFilterOMDiff",     po::value<double>(&xFilterOMDiff)->default_value(xFilterOMDiff),     "filter x-linked peptides by comparison to respective open mod scores")
//                ("Search.xFilterFragDevZ",   po::value<double>(&xFilterFragDevZ)->default_value(xFilterFragDevZ), "filter x-linked peptides by comparing matched peak deviations between the two")
//                ("Search.xFilterMinMatched", po::value<int>(&xFilterMinMatched)->default_value(xFilterMinMatched),"min matched peaks required for each x-linked peptide")
//                ("Search.xFilterMinFracMatched", po::value<double>(&xFilterMinFracMatched)->default_value(xFilterMinFracMatched),"min fraction matched peaks required for each x-linked peptide")
                ("Search.skipREVREVsearch",  po::value<bool>(&skipREVREVsearch)->default_value(skipREVREVsearch), "skip REV-REV search")
//                ("Search.skipREVFWDsearch",  po::value<bool>(&skipREVFWDsearch)->default_value(skipREVFWDsearch), "skip REV-FWD decoy search")
                ("Search.skipHalfREVFWDsearch",po::value<bool>(&skipHalfREVFWDsearch)->default_value(skipHalfREVFWDsearch), "cut number of REV-FWD decoys in half")
                ;

        po::options_description cmdline_options;
        cmdline_options.add(general).add(digest).add(mass).add(spectrum).add(search);

        po::options_description config_file_options;
        config_file_options.add(digest).add(mass).add(spectrum).add(search);

        po::options_description visible("Usage: XLsearch -p param_file [options] ...");
        visible.add(general).add(digest).add(mass).add(spectrum).add(search);

        po::positional_options_description p;
//        p.add("Digest.FastaFile", 1);

        po::variables_map vm;
       	po::parsed_options parsed = po::command_line_parser(ac, av).
                options(cmdline_options).run();
	po::store(parsed, vm);
        po::notify(vm);

        if (vm.count("help") | vm.count("param-file") == 0) { //// | vm.count("FastaFile") == 0) {
            cout << visible << "\n";
            exit(0);
            return;
        }

        if (vm.count("version")) {
            cout << "XL search, version 0.9" << endl;
            cout << "Copyright (C) 2018 Julian Mintseris & The President and Fellows of Harvard University" << endl;
            exit(0);
            return;
        }

        ifstream ifs(param_file.c_str());
        if (!ifs.is_open()) {
            cout << "Error opening parameter file " << param_file << endl;
            exit(0);
        }
        po::parsed_options parsed_cfg = po::parse_config_file(ifs, config_file_options, true);
        po::store(parsed_cfg, vm);
        po::notify(vm);
        
        // Check if no searchType was specified, add default
        if (DigestOptions::multiSearchOptions.size() == 0) {
            vector<string> empty;
            parseAddSearch(empty);
        }

    } catch (std::exception& e) {
        cout << e.what() << "\n";
        exit(0);
        return;
    }

    return;
}

DigestOptions::~DigestOptions() {
}


SearchOpt::SearchOpt (int max_mis_cleav, int enz_spec, string s_t, string l_t, double l_m, string l_r1, string l_r2)
    : maxMissedCleavages(max_mis_cleav), enzymeSpecificity(enz_spec), searchType(s_t), linkType(l_t), linkerMass(l_m), linkedRes1(l_r1), linkedRes2(l_r2) {
        if (enzymeSpecificity < 0 || enzymeSpecificity > 2) {
            enzymeSpecificity = 2;
        }

        if (linkedRes1.compare(linkedRes2) == 0) {
            isHomoXLink = true;
        } else {
            isHomoXLink = false;
        }
        
        string enzyme_aas = DigestOptions::enzymeAfterSites.first + DigestOptions::enzymeBeforeSites.first;
        // TODO: strip out non-alphanums
        effectiveMaxMissedCleavages1 = max_mis_cleav;
        effectiveMaxMissedCleavages2 = max_mis_cleav;
        xlinkedCleavageSites = 0;
        for (int i = 0; i < enzyme_aas.length(); i++) {
            if (linkedRes1.find(enzyme_aas[i]) != string::npos) {
                effectiveMaxMissedCleavages1 ++;
                xlinkedCleavageSites ++;
            }
            if (linkedRes2.find(enzyme_aas[i]) != string::npos) {
                xlinkedCleavageSites ++;
                effectiveMaxMissedCleavages2 ++;
            }
        }

    }
    
    string SearchOpt::join2string () {
        string effective_max_missed_cleavages = "[" + toStr(effectiveMaxMissedCleavages1) + "," + toStr(effectiveMaxMissedCleavages2) + "]";
        return (toStr(maxMissedCleavages) + effective_max_missed_cleavages + " " + toStr(enzymeSpecificity) + " " + searchType + " " + linkType + " " + toStr(linkerMass) + " " + linkedRes1 + " " + linkedRes2);
    }
#ifndef _DIGESTOPTIONS_HPP
#define	_DIGESTOPTIONS_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <exception>

#include "../JMUtil/JMUtil.hpp"
using namespace std;
using namespace JMUtil;

#include <tr1/unordered_map>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

class MSMass;
class Modification;

class SearchOpt {
public:
    int maxMissedCleavages;
    int enzymeSpecificity;
    string searchType;
    string linkType;
    double linkerMass;
    string linkedRes1;
    string linkedRes2;
    int effectiveMaxMissedCleavages1;
    int effectiveMaxMissedCleavages2;
    short int xlinkedCleavageSites;
    bool isHomoXLink;

    SearchOpt (int max_mis_cleav, int enz_spec, string s_t, string l_t, double l_m, string l_r1, string l_r2);
    
    string join2string ();
};

/**
 * Encapsulate the digest options into a class, reading them in from cmdline/config file
 * with boost::program_options.  The class helps keep the rest of the code independent
 * of the boost specifics and allows other option parsers to be subbed in later if desired
 */
class DigestOptions {
public:
    // General Options
    static string param_file;
    static bool verbose;
    static int maxNumThreads;

    // Digest Options
    static string fastaFile;
    static string fastaHeaderRegEx;
    static string reverseTag;
    static double minMass;
    static double maxMass;
    static int minLength;
    static int maxLength;
    static string enzymeName;
    static string enzymeRegEx;
    static int enzymeSpecificity;
    static pair<string,bool> enzymeBeforeSites;
    static pair<string,bool> enzymeAfterSites;
    static bool swapReverseInitMet;
    static bool swapReverseCutSites;
    static int maxMissedCleavages;
    static bool cleaveInitMet; // only when enzymeSpecificity == 2

    // Mass Options
    static string varModProteinSpecificDelimiter;
    static map<char, double> elementMassMod;
    static map<char, double> staticAAMod;
    static multimap<char, Modification> variableAAMod; // CURRENTLY ONLY SUPPORTS 1 MOD PER AA. NEED TO GET CLEVER TO DO MULTIPLES
    static string variablyModAas;
    static string orderedModSymbols;
    static string orderedNTerminalModSymbols;
    static string orderedCTerminalModSymbols;
    static bool lowerCaseMods;

    // Spectrum Options
    static int windowSize;
    static int peakDepth;
    static bool deIsotope;
    static bool printDeisoSpectra;
    static double filterLowMzMass;
    static bool filterCommonCorrel;
    static int minRequiredPeaks;
    static bool removePrecursorNL;
    static bool removeWaterNLoss;
    static bool removePhosphoNLoss;
    static int maxTheoFragCharge;
    static double maxMS2ItyCutoffForTIC;

    // Search Options
    static double precursorTol;
    static bool precursorTolRelative;
    static double fragmentTol;
    static int fragmentTolInt;
    static bool fragmentTolRelative;
    static int nMaxVariableMods;
    static bool noModSiteCleavage;
    static string searchDir;
    static string outputFile;
    static int maxHitRank;
    static int maxHitOutputRank;
    static string scoringFunction;
    static string fragmentType;
    static vector<double> precursorSteps;
    static bool indexAllPeptides; // internal flag; only 0 for spectrum-indexed 'Peptide' searches
    static bool indexSeenSequences; // for non-peptide-indexed searches

    // X-link Options
    static vector<SearchOpt> multiSearchOptions;
    static int currentSearchOpt;
    static string searchType;
    static double linkerMass;
    static string linkType;
    static string linkedRes1;
    static string linkedRes2;
    static int xFilterTopProts;
    static double xFilterOMDiff;
    static double xFilterFragDevZ;
    static int xFilterMinMatched;
    static double xFilterMinFracMatched;
    static int xMaxMissedCleavages;
    static bool skipREVREVsearch;
    static bool skipREVFWDsearch;
    static bool skipHalfREVFWDsearch;
    static int loopFilterMinMatched;

    DigestOptions(int ac, char *av[]);

    void printParams();
    void parseNumThreads (const int &input_num_th);
    void parsePrecursorTolerance(const string &tolerance_input);
    void parseFragmentTolerance(const string &tolerance_input);
    void parseVariableMods(const vector<string> &var_mods);
    static void updateVariablyModAas();
    static void filterVarModsForXlink();
    void parseStaticMods(const vector<string> &static_mods);
    void parseElementMods(const vector<string> &element_mods);
    void parseAddSearch(const vector<string> &add_search);
    void parseEnzymeRegEx(string enzyme_regex);
    void parsePrecursorStepRange(const string &step_range);
    void parseEnzymeSpecificity(const int &num_specific_ends);
    void parseMissedCleavages(const int &missed);
    virtual ~DigestOptions();

private:
    void parseProgramOptions(int ac, char *av[]);
};

#endif	/* _DIGESTOPTIONS_HPP */


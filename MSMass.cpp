#include "MSMass.hpp"

tr1::unordered_map<char,   double> MSMass::elementMass;
tr1::unordered_map<char,   string> MSMass::aaComp;
tr1::unordered_map<char,   double> MSMass::aa2mass;
tr1::unordered_map<string, double> MSMass::commonMass;
tr1::unordered_map<char,   double> MSMass::averagineAtomicComp;
double MSMass::numericMaxLimit = MSMass::int2double(INT_MAX); // switch to UINT_MAX here for m/z values > 2100
double MSMass::protonMass = 1.007276466771;
double MSMass::intDoubleConvertFactor = 1000000.0;

MSMass::MSMass(DigestOptions opt) {
    initElementMass();
    addElementMassMods(opt);

    initAAComp();

    initAAMass();
    addStaticAAMods(opt);
    initCommonMass();
    checkVariableMods(opt);
}

void MSMass::initElementMass() {
    MSMass::elementMass['H'] = 1.0078250321;
    MSMass::elementMass['C'] = 12.0000000000;
    MSMass::elementMass['N'] = 14.0030740052;
    MSMass::elementMass['O'] = 15.9949146221;
    MSMass::elementMass['P'] = 30.9737615120;
    MSMass::elementMass['S'] = 31.9720706912;
    return;
}

void MSMass::addElementMassMods(DigestOptions opt) {
    for (map<char, double>::const_iterator eMMiT = opt.elementMassMod.begin(); eMMiT != opt.elementMassMod.end(); eMMiT++) {
        tr1::unordered_map<char, double>::iterator umIt = elementMass.find((*eMMiT).first);
        if (umIt != elementMass.end()) {
            umIt->second += (*eMMiT).second;
        } else {
            cerr << "ERROR: Element '" << (*eMMiT).first << "' is not recognized" << endl;
            exit(EXIT_FAILURE);
        }
    }
    return;
}

void MSMass::initAAComp() {
    aaComp['A'] = "C3H5N1O1";
    aaComp['C'] = "C3H5N1O1S1";
    aaComp['D'] = "C4H5N1O3";
    aaComp['E'] = "C5H7N1O3";
    aaComp['F'] = "C9H9N1O1";
    aaComp['G'] = "C2H3N1O1";
    aaComp['H'] = "C6H7N3O1";
    aaComp['I'] = "C6H11N1O1";
    aaComp['K'] = "C6H12N2O1";
    aaComp['L'] = "C6H11N1O1";
    aaComp['M'] = "C5H9N1O1S1";
    aaComp['N'] = "C4H6N2O2";
    aaComp['P'] = "C5H7N1O1";
    aaComp['Q'] = "C5H8N2O2";
    aaComp['R'] = "C6H12N4O1";
    aaComp['S'] = "C3H5N1O2";
    aaComp['T'] = "C4H7N1O2";
    aaComp['U'] = "C5H5N1O2"; // Pyroglutamate - invitro = N-terminal Q -NH3 or N-terminal E -H2O
    aaComp['V'] = "C5H9N1O1";
    aaComp['W'] = "C11H10N2O1";
//    aaComp['X'] = "C6H11N1O1"; // Only SEQUEST does this.  Makes no sense because in FASTA format this means 'any aa'
    aaComp['Y'] = "C9H9N1O2";
    
    aaComp[']'] = ""; // peptide N-terminus - makes it easier to specify static/variable terminal peptide mods in params file
    aaComp['['] = ""; // peptide C-terminus - makes it easier to specify static/variable terminal peptide mods in params file
    aaComp['+'] = ""; // protein N-terminus - makes it easier to specify static/variable terminal protein mods in params file
    aaComp['-'] = ""; // protein C-terminus - makes it easier to specify static/variable terminal protein mods in params file
    return;
}

double MSMass::elComp2Mass (string el_comp) {
    char el;
    double numEl;
    double mass = 0.0;

    const boost::regex re("(\\w)(\\d+)");
    boost::sregex_iterator reIt(el_comp.begin(), el_comp.end(), re);
    boost::sregex_iterator end;
    boost::smatch m;
    while (reIt != end) {
        m = (*reIt++);
        el = ((string) m[1])[0];
        numEl = atof(((string) m[2]).c_str());
        mass += elementMass[el]*numEl;
    }
    return mass;
}

void MSMass::initAAMass() {
    for (tr1::unordered_map<char, string>::const_iterator umIt = aaComp.begin(); umIt != aaComp.end(); umIt++) {
        ////cout << umIt->first << "\t" << umIt->second << "\t" << elComp2Mass(umIt->second) << endl; ////
        aa2mass[umIt->first] = elComp2Mass(umIt->second);
    }
    // Ambiguous amino acids
    aa2mass['B'] = (aa2mass['D'] + aa2mass['N'])/2;
    aa2mass['Z'] = (aa2mass['E'] + aa2mass['Q'])/2;

    return;
}

void MSMass::initCommonMass() {
    commonMass["H2O"] = elComp2Mass("H2O1");
    commonMass["NH3"] = elComp2Mass("N1H3");
    commonMass["HPO3"] = elComp2Mass("H1P1O3");
    commonMass["Nterm"] = protonMass + commonMass["H2O"];
    commonMass["H2O_proton"] = commonMass["H2O"] + protonMass;
    commonMass["averagine"] = 111.1254;

    MSMass::averagineAtomicComp['C'] = 4.9384;
    MSMass::averagineAtomicComp['H'] = 7.7583;
    MSMass::averagineAtomicComp['N'] = 1.3577;
    MSMass::averagineAtomicComp['O'] = 1.4773;
    MSMass::averagineAtomicComp['S'] = 0.0417;
    return;
}

void MSMass::addStaticAAMods (DigestOptions opt) {
    for (map<char, double>::const_iterator sAAMiT = opt.staticAAMod.begin(); sAAMiT != opt.staticAAMod.end(); sAAMiT++) {
        tr1::unordered_map<char,double>::iterator umIt = aa2mass.find((*sAAMiT).first);
        if (umIt != aa2mass.end()) {
            umIt->second += (*sAAMiT).second;
        } else {
            cerr << "ERROR: Residue '" << (*sAAMiT).first << "' is not recognized" << endl;
            exit(EXIT_FAILURE);
        }
    }
    return;
}

void MSMass::checkVariableMods (DigestOptions opt) {
    for (multimap<char, Modification>::const_iterator vAAMiT = opt.variableAAMod.begin(); vAAMiT != opt.variableAAMod.end(); vAAMiT++) {
        tr1::unordered_map<char,double>::iterator umIt = aa2mass.find((*vAAMiT).first);
        if (umIt == aa2mass.end()) {
            cerr << "ERROR: Residue '" << (*vAAMiT).first << "' is not recognized" << endl;
            exit(EXIT_FAILURE);
        }
    }
    return;
}

/* Calculate the M+H mass of a bare peptide.  Bail and return 0.0 if encounter
 * an unexpected residue.  For now assume that protein mods cannot be static
 */
double MSMass::seq2MHMass (string seq) {
    return seq2MHMass('\0', seq, '\0');
}

/* Calculate the M+H mass of a bare peptide.  Bail and return 0.0 if encounter
 * an unexpected residue.  Assume protein terminal mods take precedence over peptide terminal mods
 */
double MSMass::seq2MHMass (char flank_res_nterm, string seq, char flank_res_cterm) {
    double mass = commonMass["Nterm"];
    if (flank_res_nterm == '+' && aa2mass[flank_res_nterm] > 0) {
        mass += aa2mass[flank_res_nterm]; // static protein terminal mods
    } else {
        mass += aa2mass[']']; // static peptide terminal mods
    }
    for (string::iterator i = seq.begin(); i != seq.end(); i++) {
        tr1::unordered_map<char,double>::iterator a2mIt = aa2mass.find(*i);
        if (a2mIt == aa2mass.end()) {
            mass -= numeric_limits<double>::max();
            //cerr << "WARNING: ... encountered unrecognized element: " << *i << " of indeterminate mass" << endl; ////
            return mass;
        } else {
            mass += (*a2mIt).second;
        }
    }
    if (flank_res_cterm == '-' && aa2mass[flank_res_cterm] > 0) {
        mass += aa2mass[flank_res_cterm]; // static protein terminal mods
    } else {
        mass += aa2mass['['];  // static peptide terminal mods
    }
    return mass;
}


void MSMass::seq2MassVec(char flank_res_nterm, string seq, char flank_res_cterm, vector<double> &massVec) {
    massVec.resize(seq.length());
    for (int i = 0; i < seq.length(); i++) {
        massVec[i] = aa2mass[ seq[i] ];
    }
    if (flank_res_nterm == '+' && aa2mass[flank_res_nterm] > 0) {
        massVec.front() += aa2mass[flank_res_nterm]; // static protein terminal mods
    } else {
        massVec.front() += aa2mass[']']; // static peptide terminal mods
    }
    if (flank_res_cterm == '-' && aa2mass[flank_res_cterm] > 0) {
        massVec.back() += aa2mass[flank_res_cterm]; // static protein terminal mods
    } else {
        massVec.back() += aa2mass['[']; // static peptide terminal mods
    }
    return;
}

// isobaric conversion functor
// maybe will add some other here depending on CID/HCD? M* == F?
// ideally, want to derive the conversion pairs directly from aa2mass and fragment tolernace info
char MSMass::i2l(char c) {
    if (c == 'I')
        return 'L';
    else
        return c;
}

string MSMass::makeIsobaric (const string &peptide){
    string isobaric = peptide;
    transform(isobaric.begin(), isobaric.end(), isobaric.begin(), i2l);
    return isobaric;
}

// Senko MW, Beu SC, McLafferty FW. J Am Soc Mass Spectrom. 1995;6:229–33.
string MSMass::averagineModel (double m) {
    double num_avrgn_res = m/commonMass["averagine"];
    tr1::unordered_map<char,int> avrgn_comp;
    string chnos = "CHNOS";
    double approx_mass = 0.0;
    for (int i = 0; i < chnos.length(); i++) {
        avrgn_comp[chnos[i]] = (JMUtil::roundTo( averagineAtomicComp[chnos[i]] * num_avrgn_res, 0 ));
        approx_mass += elementMass[chnos[i]] * (double) avrgn_comp[chnos[i]];
    }

    int H_diff = JMUtil::roundTo( (m - approx_mass)/elementMass['H'], 0);
    if ((avrgn_comp['H'] + H_diff) < 0) {
        avrgn_comp['C'] --;
        approx_mass -= elementMass['C'];
        H_diff = JMUtil::roundTo((m - approx_mass)/elementMass['H'], 0);
    }
    avrgn_comp['H'] += H_diff;

    string formula;
    for (int i = 0; i < chnos.length(); i++) {
        formula += (chnos[i] + toStr(avrgn_comp[ chnos[i] ]));
    }
    return formula;
}

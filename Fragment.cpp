#include "Fragment.hpp"

void PeptideFragmentTable::prettyPrint(string outfile) {
    std::ofstream outFile(outfile.c_str(), ios::out);
    outFile << "   = PeptideFragmentTable::prettyPrint =   " << endl; ////
    outFile << peptide << "\t" << "\t" << id << endl; ////
    /*     for (int i = 0; i < b[0].size(); i++) { */
    /*       outFile << "\t" << i << " " << b[0][i] << "\t" << y[0][i] << endl; */
    /*     } */
    outFile << endl;
    for (int c = 0; c < max_charge; c++) {
        outFile << "== +" << (c + 1) << " ==" << endl;
        int i = 0;
        // j here is the index into the pft.y vector
        int j = b[0].size() - i;
        outFile << (i + 1) << "  ";
        outFile << fixed << setprecision(3) << b[c][i] << "\t";
        outFile << fixed << setprecision(3) << "--";
        outFile << "  " << (j + 1) << endl;
        for (int i = 1; i < b[0].size(); i++) {
            int j = b[0].size() - i;
            outFile << (i + 1) << "  ";
            outFile << fixed << setprecision(3) << b[c][i] << "\t";
            outFile << fixed << setprecision(3) << y[c][j];
            outFile << "  " << (j + 1) << endl;
        }
        i = b[0].size();
        j = b[0].size() - i;
        outFile << (i + 1) << "  ";
        outFile << fixed << setprecision(3) << "--" << "\t";
        outFile << fixed << setprecision(3) << y[c][j];
        outFile << "  " << (j + 1) << endl;
    }
    outFile.close();
    return;
}

void XLinkedFragmentTable::prettyPrint(string outfile) {
    std::ofstream outFile(outfile.c_str(), ios::out);
    cout << "   = XLinkedFragmentTable::prettyPrint =   " << endl; ////
    outFile << peptide << "-" << peptide2 << "\t" << pos1 << "-" << pos2 << "; " << pos1_y << "-" << pos2_y << "\t" << id << endl; ////
    // TODO: the output gets messed up when symbols are used for mods.  Need to strip them out here and store new temp peptide and peptide2
    ////cout << "pos1: " << pos1 << "\t" << "pos1_y: " << pos1_y << endl; ////
    ////cout << "pos2: " << pos2 << "\t" << "pos2_y: " << pos2_y << endl; ////

    outFile << endl;
    for (int c = 0; c < max_charge; c++) {
        outFile << "== +" << (c + 1) << " ==" << endl;
        int i = 0;
        // j here is the index into the pft.y vector
        int j = peptide.length() - 1 - i;
        outFile << peptide[i] << "(1)" << (i + 1) << "  ";
        outFile << setw(9) << b[c][i] << "\t";
        outFile << setw(9) << "--";
        outFile << "  " << (j + 1) << "(1)" << endl;
        for (int i = 1; i < (peptide.length() - 1); i++) {
            int j = peptide.length() - 1 - i;
            outFile << peptide[i] << "(1)" << (i + 1) << "  ";
            outFile << setw(9) << b[c][i] << "\t";
            outFile << setw(9) << y[c][j];
            outFile << "  " << (j + 1) << "(1)" << endl;
        }
        i = peptide.length() - 1;
        j = peptide.length() - 1 - i;
        outFile << peptide[i] << "(1)" << (i + 1) << "  ";
        outFile << setw(9) << "--" << "\t";
        outFile << setw(9) << y[c][j];
        outFile << "  " << (j + 1) << "(1)" << endl;

        i = 0;
        // j here is the index into the pft.y vector
        j = peptide2.length() - 1 - i;
        outFile << peptide2[i] << "(2)" << (i + 1) << "  ";
        outFile << setw(9) << b2[c][i] << "\t";
        outFile << setw(9) << "-----";
        outFile << "  " << (j + 1) << "(2)" << endl;
        for (int i = 1; i < (peptide2.length() - 1); i++) {
            int j = peptide2.length() - 1 - i;
            outFile << peptide2[i] << "(2)" << (i + 1) << "  ";
            outFile << setw(9) << b2[c][i] << "\t";
            outFile << setw(9) << y2[c][j];
            outFile << "  " << (j + 1) << "(2)" << endl;
        }
        i = peptide2.length() - 1;
        j = peptide2.length() - 1 - i;
        outFile << peptide2[i] << "(2)" << (i + 1) << "  ";
        outFile << setw(9) << "-----" << "\t";
        outFile << setw(9) << y2[c][j];
        outFile << "  " << (j + 1) << "(2)" << endl;
        outFile << endl;
    }
    outFile << endl;
    outFile << "lf lf lf lf" << endl;
    if (lf.size() > 0) {
        for (int i = 0; i < lf.size(); i++) {
            for (int j = 0; j < lf[0].size(); j++) {
                outFile << lf[i][j] << " ";
            }
            outFile << endl;
        }
        ////outFile << "lf " << lf[0][0] << " " << lf[0][1] << endl; ////
    }
    outFile.close();
    return;
}

void LoopLinkedFragmentTable::prettyPrint(string outfile) {
    std::ofstream outFile(outfile.c_str(), ios::out);
    outFile << "   = PeptideFragmentTable::prettyPrint =   " << endl; ////
    outFile << peptide << "\t" << pos1 << "-" << pos2 << "\t" << id << endl; ////
    /*     for (int i = 0; i < b[0].size(); i++) { */
    /*       outFile << "\t" << i << " " << b[0][i] << "\t" << y[0][i] << endl; */
    /*     } */
    outFile << endl;
    for (int c = 0; c < max_charge; c++) {
        int i = 0;
        // j here is the index into the pft.y vector
        int j = b[0].size() - i;
        outFile << (i + 1) << "  ";
        outFile << fixed << setprecision(3) << b[c][i] << "\t";
        outFile << fixed << setprecision(3) << "--";
        outFile << "  " << (j + 1) << endl;
        for (int i = 1; i < b[0].size(); i++) {
            int j = b[0].size() - i;
            outFile << (i + 1) << "  ";
            outFile << fixed << setprecision(3) << b[c][i] << "\t";
            outFile << fixed << setprecision(3) << y[c][j];
            outFile << "  " << (j + 1) << endl;
        }
        i = b[0].size();
        j = b[0].size() - i;
        outFile << (i + 1) << "  ";
        outFile << fixed << setprecision(3) << "--" << "\t";
        outFile << fixed << setprecision(3) << y[c][j];
        outFile << "  " << (j + 1) << endl;
    }
    outFile.close();
    return;
}

void FragmentTable::printPeakList(string outfile) {
    std::ofstream outFile(outfile.c_str(), ios::out);
    for (vector<TheoreticalPeak>::iterator vecIt = fullPeakList.begin(); vecIt != fullPeakList.end(); vecIt++) { ////
        outFile << fixed << setprecision(5) << MSMass::int2double( (*vecIt).m_z )  << "\t";
        outFile << (*vecIt).charge << "\t";
        outFile << (*vecIt).ion_type << "\t";
        outFile << (*vecIt).pep_1_2 << endl;
    }
    outFile.close();
    return;
}


void PeptideFragmentTable::fragment() {
    double bprev = MSMass::protonMass;
    double yprev = MSMass::commonMass["H2O_proton"];

    for (int i = 0; i < massVec.size() - 1; i++) {
        // j here is the index into massVec
        int j = massVec.size() - 1 - i;
        b[0][i] = bprev + massVec[i];
        y[0][i] = yprev + massVec[j];
        //cout << i << " " << b[0][i] << "\t" << y[0][i] << " " << j << endl;
        bprev = b[0][i];
        yprev = y[0][i];
    }

    for (int c = 1; c < max_charge; c++) {
        for (int i = 0; i < massVec.size() - 1; i++) {
            b[c][i] = (b[0][i] + MSMass::protonMass * c) / (c + 1);
            y[c][i] = (y[0][i] + MSMass::protonMass * c) / (c + 1);
        }
    }

    return;
}


void XLinkedFragmentTable::fragment() {
    double bprev = MSMass::protonMass;
    double yprev = MSMass::commonMass["H2O_proton"];

    for (int i = 0; i < massVec.size() - 1; i++) {
        // j here is the index into massvec
        int j = massVec.size() - 1 - i;
        b[0][i] = bprev + massVec[i];
        if (i == pos1) {
            b[0][i] += m_h_mass2 - MSMass::protonMass + DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
        }
        y[0][i] = yprev + massVec[j];
        if (i == pos1_y) {
            y[0][i] += m_h_mass2 - MSMass::protonMass + DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
        }
        ////cout << i << " " << b[0][i] << "\t" << y[0][i] << " " << j << endl; ////
        bprev = b[0][i];
        yprev = y[0][i];
    }
    ////cout << endl; ////

    bprev = MSMass::protonMass;
    yprev = MSMass::commonMass["H2O_proton"];

    for (int i = 0; i < massVec2.size() - 1; i++) {
        // j here is the index into massvec
        int j = massVec2.size() - 1 - i;
        b2[0][i] = bprev + massVec2[i];
        if (i == pos2) {
            b2[0][i] += m_h_mass1 - MSMass::protonMass + DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
        }
        y2[0][i] = yprev + massVec2[j];
        if (i == pos2_y) {
            y2[0][i] += m_h_mass1 - MSMass::protonMass + DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
        }
        ////cout << i << " " << b2[0][i] << "\t" << y2[0][i] << " " << j << endl; ////
        bprev = b2[0][i];
        yprev = y2[0][i];
    }
    ////cout << endl; ////
//    prettyPrint("fragmentTable.out_" + peptide + "_" + peptide2 + "_" + toStr(pos1) + "-" + toStr(pos2)); ////
    for (int c = 1; c < max_charge; c++) {
        for (int i = 0; i < massVec.size() - 1; i++) {
            b[c][i] = (b[0][i] + MSMass::protonMass * c) / (c + 1);
            y[c][i] = (y[0][i] + MSMass::protonMass * c) / (c + 1);
            ////cout << "\t" << i << " " << b[c][i] << "\t" << y[c][i] << " " << c << endl;
        }
        for (int j = 0; j < massVec2.size() - 1; j++) {
            b2[c][j] = (b2[0][j] + MSMass::protonMass * c) / (c + 1);
            y2[c][j] = (y2[0][j] + MSMass::protonMass * c) / (c + 1);
            ////cout << "\t" << j << " " << b2[c][j] << "\t" << y2[c][j] << " " << c << endl;
        }
    }

//    /* Sometimes, when comparing a (bogus) end-to-end xlink that is really a missed cleavage
//     * the x-linked version gets a slightly higher score because it considers fewer theoretical
//     * peaks for equal number of matches.  Here, we try to consider theoretical x-link breaks
//     * (which sometimes do occur especially in 0-length x-links), primarily to increase that k
//     * and decrease the score relative to the (correct) regular peptide match.
//     * This does not seem to have the desired effect in all situations.
//     */
    
    if (pos1 < 0 && pos2 < 0) // don't bother with linker fragments if looking at entangled/mixed spectra
        return;
    
    int lf_max_charge = max_charge > 1 ? 2 : 1;
    if (DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkType[0] == '0') {

        int lf_0_num = 2;
        lf.resize(max_charge, vector<double> (lf_0_num));
        if        ((res1 == 'K' || (pos1 == 0 && pos1_term)) && (res2 == 'D' || res2 == 'E' || (pos2 == (massVec2.size() - 1) && pos2_term))) {  // peptide1 is the y-ion, peptide2 is the b-ion
            lf[0][0] = m_h_mass1;
            lf[0][1] = m_h_mass2 - MSMass::commonMass["H2O"];
        } else if ((res2 == 'K' || (pos2 == 0 && pos2_term)) && (res1 == 'D' || res1 == 'E' || (pos1 == (massVec.size() - 1) && pos1_term))) {  // peptide1 is the b-ion, peptide2 is the y-ion
            lf[0][0] = m_h_mass1 - MSMass::commonMass["H2O"];
            lf[0][1] = m_h_mass2;
        } else {
            ////cerr << "WARNING: Don't know how to assign " << peptide << " - " << peptide2 << "\t" << pos1 << "-" << pos2 << endl; ////
            lf.resize(max_charge, vector<double> (0));
            return;
        }
        /**/
        ////cout << "JMDEBUG lf orig " << lf[0][0] << " " << lf[0][1] << endl; ////
        for (int c = 1; c < lf_max_charge; c++) { // lf_max_charge ??? or max_charge ??? ???
            for (int i = 0; i < lf_0_num; i++) {
                lf[c][i] = (lf[0][i] + MSMass::protonMass * c) / (c + 1);
                ////cout << "JMDEBUG lf: " << lf[c][i] << " "; ////
            }
            ////cout << endl; ////
        }
        /**/
    } else if (DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkType.substr(0,2).compare("DE") == 0) {
        int lf_DE_num = 4;
        double linker_mass = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
        lf.resize(max_charge, vector<double> (lf_DE_num));
        lf[0][0] = m_h_mass1 - MSMass::commonMass["H2O"];
        lf[0][1] = m_h_mass2 - MSMass::commonMass["H2O"];
        /* want peptide + linker - H
         * but linkerMass already has 2 waters subtracted so it's really:
         * peptide + (- water + linker - water) + water
         */
        lf[0][2] = m_h_mass1 + linker_mass + MSMass::commonMass["H2O"];
        lf[0][3] = m_h_mass2 + linker_mass + MSMass::commonMass["H2O"];
        for (int c = 1; c < lf_max_charge; c++) {
            for (int i = 0; i < lf_DE_num; i++) {
                lf[c][i] = (lf[0][i] + MSMass::protonMass * c) / (c + 1);
                ////cout << "JMDEBUG lf: " << lf[c][i] << " "; ////
            }
            ////cout << endl; ////
        }
    } else if (DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkType.compare("BS3") == 0) {
        int lf_BS3_num = 4;
        double linker_mass = DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
        lf.resize(max_charge, vector<double> (lf_BS3_num));
        lf[0][0] = m_h_mass1;
        lf[0][1] = m_h_mass2;

        lf[0][2] = m_h_mass1 + linker_mass;
        lf[0][3] = m_h_mass2 + linker_mass;
        for (int c = 1; c < lf_max_charge; c++) {
            for (int i = 0; i < lf_BS3_num; i++) {
                lf[c][i] = (lf[0][i] + MSMass::protonMass * c) / (c + 1);
                ////cout << "JMDEBUG lf: " << lf[c][i] << " "; ////
            }
            ////cout << endl; ////
        }
    }
    return;
}

void LoopLinkedFragmentTable::fragment() {
    double bprev = MSMass::protonMass;
    double yprev = MSMass::commonMass["H2O_proton"];

    // For now just produce fragments that leave the cross-linker on the pos1 and pos2_y residue
    // Could have additional fragments based on the cross-linker staying with pos2 and pos1_y
    for (int i = 0; i < massVec.size() - 1; i++) {
        // j here is the index into massvec
        int j = massVec.size() - 1 - i;
        b[0][i] = bprev + massVec[i];
        if (i == pos1) {
            b[0][i] += DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
        }
        y[0][i] = yprev + massVec[j];
        if (i == pos2_y) {
            y[0][i] += DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
        }
        ////cout << i << " " << b[0][i] << "\t" << y[0][i] << " " << j << endl; ////
        bprev = b[0][i];
        yprev = y[0][i];
    }
    ////cout << endl; ////

    // make sure the linked positions are in the right sequence to simplify subsequent calculation
    int p1 = pos1;
    int p2 = pos2;
    int p1_y = pos1_y;
    int p2_y = pos2_y;
    if (pos1 > pos2) {
        swap(p1, p2);
        swap(p1_y, p2_y);
    }

    // remove from consideration fragments between loop-linked positions
    // They really don't appear to be there no matter how hard I look even with HCD. 
    // should just just erase them, but zero out for now for easier debugging
    for (int i = p1; i < p2; i++) {
        b[0][i] = 0;
    }
    for (int i = p2_y; i < p1_y; i++) {
        y[0][i] = 0;
    }

//    bprev = MSMass::protonMass;
//    yprev = MSMass::commonMass["H2O_proton"];
//
//    for (int i = 0; i < massVec2.size() - 1; i++) {
//        // j here is the index into massvec
//        int j = massVec2.size() - 1 - i;
//        b2[0][i] = bprev + massVec2[i];
//        if (i == pos2) {
//            b2[0][i] += m_h_mass1 - MSMass::protonMass + DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
//        }
//        y2[0][i] = yprev + massVec2[j];
//        if (i == pos2_y) {
//            y2[0][i] += m_h_mass1 - MSMass::protonMass + DigestOptions::multiSearchOptions[DigestOptions::currentSearchOpt].linkerMass;
//        }
//        ////cout << i << " " << b2[0][i] << "\t" << y2[0][i] << " " << j << endl; ////
//        bprev = b2[0][i];
//        yprev = y2[0][i];
//    }
//    ////cout << endl; ////

    for (int c = 1; c < max_charge; c++) {
        for (int i = 0; i < massVec.size() - 1; i++) {
            b[c][i] = (b[0][i] + MSMass::protonMass * c) / (c + 1);
            y[c][i] = (y[0][i] + MSMass::protonMass * c) / (c + 1);
            ////cout << "\t" << i << " " << b[c][i] << "\t" << y[c][i] << " " << c << endl;
        }
//        for (int j = 0; j < massVec2.size() - 1; j++) {
//            b2[c][j] = (b2[0][j] + MSMass::protonMass * c) / (c + 1);
//            y2[c][j] = (y2[0][j] + MSMass::protonMass * c) / (c + 1);
//            ////cout << "\t" << j << " " << b2[c][j] << "\t" << y2[c][j] << " " << c << endl;
//        }
    }
    return;
}

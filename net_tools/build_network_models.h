#ifndef BUILD_NETWORKS_MODELS_H
#define BUILD_NETWORKS_MODELS_H

#include "reaction_network.h"


void build_model_Yung_99(std::vector<species> &sp, std::vector<reaction> &re) {
	sp.push_back(species(sp.size(), "hv"));
	size_t ihv=0;
	sp.push_back(species(sp.size(), "O2"));
	size_t iO2=1;
	sp.push_back(species(sp.size(), "O"));
	size_t iO=2;
	sp.push_back(species(sp.size(), "O1D"));
	size_t iO1D=3;
	sp.push_back(species(sp.size(), "O3"));
	size_t iO3=4;
	sp.push_back(species(sp.size(), "HO2"));
	size_t iHO2=5;
	sp.push_back(species(sp.size(), "OH"));
	size_t iOH=6;
	sp.push_back(species(sp.size(), "H2O"));
	size_t iH2O=7;
	sp.push_back(species(sp.size(), "H2O2"));
	size_t iH2O2=8;
	sp.push_back(species(sp.size(), "NO"));
	size_t iNO=9;
	sp.push_back(species(sp.size(), "NO2"));
	size_t iNO2=10;	
	sp.push_back(species(sp.size(), "NO3"));
	size_t iNO3=11;
	sp.push_back(species(sp.size(), "N2O"));
	size_t iN2O=12;
	sp.push_back(species(sp.size(), "N2"));
	size_t iN2=13;
	sp.push_back(species(sp.size(), "N2O5"));
	size_t iN2O5=14;
	sp.push_back(species(sp.size(), "HNO2"));
	size_t iHNO2=15;
	sp.push_back(species(sp.size(), "HNO3"));
	size_t iHNO3=16;
	sp.push_back(species(sp.size(), "HO2NO2"));
	size_t iHO2NO2=17;
	sp.push_back(species(sp.size(), "CH3O2NO2"));
	size_t iCH3O2NO2=18;
	sp.push_back(species(sp.size(), "NH3"));
	size_t iNH3=19;
	sp.push_back(species(sp.size(), "Cl2"));
	size_t iCl2=20;
	sp.push_back(species(sp.size(), "ClO"));
	size_t iClO=21;
	sp.push_back(species(sp.size(), "ClOO"));
	size_t iClOO=22;
	sp.push_back(species(sp.size(), "OClO"));
	size_t iOClO=23;
	sp.push_back(species(sp.size(), "ClO3"));
	size_t iClO3=24;
	sp.push_back(species(sp.size(), "Cl2O"));
	size_t iCl2O=25;
	sp.push_back(species(sp.size(), "Cl2O2"));
	size_t iCl2O2=26;
	sp.push_back(species(sp.size(), "Cl2O3"));
	size_t iCl2O3=27;
	sp.push_back(species(sp.size(), "ClNO"));
	size_t iClNO=28;
	sp.push_back(species(sp.size(), "ClNO2"));
	size_t iClNO2=29;
	sp.push_back(species(sp.size(), "ClONO"));
	size_t iClONO=30;
	sp.push_back(species(sp.size(), "ClNO3"));
	size_t iClNO3=31;
	sp.push_back(species(sp.size(), "HCl"));
	size_t iHCl=32;
	sp.push_back(species(sp.size(), "HOCl"));
	size_t iHOCl=33;
	sp.push_back(species(sp.size(), "CH3Cl"));
	size_t iCH3Cl=34;
	sp.push_back(species(sp.size(), "CH2FCl"));
	size_t iCH2FCl=35;
	sp.push_back(species(sp.size(), "CHFCl2"));
	size_t iCHFCl2=36;
	sp.push_back(species(sp.size(), "CHF2Cl"));
	size_t iCHF2Cl=37;
	sp.push_back(species(sp.size(), "CCl4"));
	size_t iCCl4=38;
	sp.push_back(species(sp.size(), "CFCl3"));
	size_t iCFCl3=39;
	sp.push_back(species(sp.size(), "CF2Cl2"));
	size_t iCF2Cl2=40;
	sp.push_back(species(sp.size(), "CF2Cl"));
	size_t iCF2Cl=41;
	sp.push_back(species(sp.size(), "CFCl2"));
	size_t iCFCl2=42;
	sp.push_back(species(sp.size(), "CCl3"));
	size_t iCCl3=43;
	sp.push_back(species(sp.size(), "H"));
	size_t iH=44;
	sp.push_back(species(sp.size(), "Cl"));
	size_t iCl=45;
	sp.push_back(species(sp.size(), "CHF3"));
	size_t iCHF3=46;	
	sp.push_back(species(sp.size(), "H2"));
	size_t iH2=47;
	sp.push_back(species(sp.size(), "N"));
	size_t iN=48;
	sp.push_back(species(sp.size(), "CH3O2"));
	size_t iCH3O2=49;
	sp.push_back(species(sp.size(), "NH2"));
	size_t iNH2=50;
	sp.push_back(species(sp.size(), "Products"));
	size_t iProducts=51;
	sp.push_back(species(sp.size(), "CH3"));
	size_t iCH3=52;
	sp.push_back(species(sp.size(), "CHF2"));
	size_t iCHF2=53;
	sp.push_back(species(sp.size(), "CH3CFCl2"));
	size_t iCH3CFCl2=54;
	sp.push_back(species(sp.size(), "Br"));
	size_t iBr=55;
	sp.push_back(species(sp.size(), "Br2"));
	size_t iBr2=56;
	sp.push_back(species(sp.size(), "S"));
	size_t iS=57;
	sp.push_back(species(sp.size(), "OCS"));
	size_t iOCS=58;
	sp.push_back(species(sp.size(), "SH"));
	size_t iSH=59;
	sp.push_back(species(sp.size(), "H2S"));
	size_t iH2S=60;
	sp.push_back(species(sp.size(), "SO"));
	size_t iSO=61;
	sp.push_back(species(sp.size(), "SO2"));
	size_t iSO2=62;
	sp.push_back(species(sp.size(), "CH3CN"));
	size_t iCH3CN=63;
	sp.push_back(species(sp.size(), "CN"));
	size_t iCN=64;
	sp.push_back(species(sp.size(), "HCN"));
	size_t iHCN=65;
	sp.push_back(species(sp.size(), "CH3O"));
	size_t iCH3O=66;
	sp.push_back(species(sp.size(), "CH3CF2Cl"));
	size_t iCH3CF2Cl=67;
	sp.push_back(species(sp.size(), "CH2F2"));
	size_t iCH2F2=68;	
	sp.push_back(species(sp.size(), "C2Cl2F4"));
	size_t iC2Cl2F4=69;
	sp.push_back(species(sp.size(), "C2Cl3F3"));
	size_t iC2Cl3F3=70;
	sp.push_back(species(sp.size(), "CF3CHCl2"));
	size_t iCF3CHCl2=71;
	sp.push_back(species(sp.size(), "CF3CHFCl"));
	size_t iCF3CHFCl=72;
	sp.push_back(species(sp.size(), "CF3CCl3"));
	size_t iCF3CCl3=73;
	sp.push_back(species(sp.size(), "CF3CFCl2"));
	size_t iCF3CFCl2=74;
	sp.push_back(species(sp.size(), "CO2"));
	size_t iCO2=75;
	sp.push_back(species(sp.size(), "CH4"));
	size_t iCH4=76;
	sp.push_back(species(sp.size(), "H2CO"));
	size_t iH2CO=77;
	sp.push_back(species(sp.size(), "C2F4Br2"));
	size_t iC2F4Br2=78;
	sp.push_back(species(sp.size(), "CH3OOH"));
	size_t iCH3OOH=79;
	sp.push_back(species(sp.size(), "CF3Br"));
	size_t iCF3Br=80;
	sp.push_back(species(sp.size(), "BrCl"));
	size_t iBrCl=81;
	sp.push_back(species(sp.size(), "HOBr"));
	size_t iHOBr=82;
	sp.push_back(species(sp.size(), "CH3Br"));
	size_t iCH3Br=83;
	sp.push_back(species(sp.size(), "F"));
	size_t iF=84;
	sp.push_back(species(sp.size(), "HF"));
	size_t iHF=85;
	sp.push_back(species(sp.size(), "CO"));
	size_t iCO=86;
	sp.push_back(species(sp.size(), "HCO"));
	size_t iHCO=87;
	sp.push_back(species(sp.size(), "CF3"));
	size_t iCF3=88;
	sp.push_back(species(sp.size(), "CF4"));
	size_t iCF4=89;
	sp.push_back(species(sp.size(), "BrNO2"));
	size_t iBrNO2=90;
	sp.push_back(species(sp.size(), "SF6"));
	size_t iSF6=91;
	sp.push_back(species(sp.size(), "BrO"));
	size_t iBrO=92;
	sp.push_back(species(sp.size(), "CF2Br2"));
	size_t iCF2Br2=93;
	sp.push_back(species(sp.size(), "CF2ClBr"));
	size_t iCF2ClBr=94;
	sp.push_back(species(sp.size(), "HBr"));
	size_t iHBr=95;
	sp.push_back(species(sp.size(), "COF2"));
	size_t iCOF2=96;
	sp.push_back(species(sp.size(), "COFCl"));
	size_t iCOFCl=97;
	sp.push_back(species(sp.size(), "CHClO"));
	size_t iCHClO=98;
	sp.push_back(species(sp.size(), "C2F6"));
	size_t iC2F6=99;
	sp.push_back(species(sp.size(), "COCl2"));
	size_t iCOCl2=100;
	sp.push_back(species(sp.size(), "CHBr3"));
	size_t iCHBr3=101;
	sp.push_back(species(sp.size(), "BrNO3"));
	size_t iBrNO3=102;
	sp.push_back(species(sp.size(), "CF3CF2CHCl2"));
	size_t iCF3CF2CHCl2=103;
	sp.push_back(species(sp.size(), "CF2ClCF2CHFCl"));
	size_t iCF2ClCF2CHFCl=104;
	sp.push_back(species(sp.size(), "CH2O"));
	size_t iCH2O=105;
	sp.push_back(species(sp.size(), "O21d"));
	size_t iO21d=106;
	sp.push_back(species(sp.size(), "O21s"));
	size_t iO21s=107;
	sp.push_back(species(sp.size(), "HD"));
	size_t iHD=108;
	sp.push_back(species(sp.size(), "HONO"));
	size_t iHONO=109;
	sp.push_back(species(sp.size(), "NH"));
	size_t iNH=110;
	sp.push_back(species(sp.size(), "C2H2"));
	size_t iC2H2=111;
	sp.push_back(species(sp.size(), "CH3CHO"));
	size_t iCH3CHO=112;
	sp.push_back(species(sp.size(), "CH3CO"));
	size_t iCH3CO=113;
	sp.push_back(species(sp.size(), "C2H4"));
	size_t iC2H4=114;
	sp.push_back(species(sp.size(), "C3H6"));
	size_t iC3H6=115;
	sp.push_back(species(sp.size(), "C2H6"));
	size_t iC2H6=116;
	sp.push_back(species(sp.size(), "Adduct"));
	size_t iAdduct=117;
	sp.push_back(species(sp.size(), "CH3D"));
	size_t iCH3D=118;
	sp.push_back(species(sp.size(), "CH3OH"));
	size_t iCH3OH=119;
	sp.push_back(species(sp.size(), "C2H5"));
	size_t iC2H5=120;
	sp.push_back(species(sp.size(), "C3H7"));
	size_t iC3H7=121;
	sp.push_back(species(sp.size(), "CH3C_O_OH"));
	size_t iCH3C_O_OH=122;
	sp.push_back(species(sp.size(), "C2H5O2"));
	size_t iC2H5O2=123;
	sp.push_back(species(sp.size(), "HC_O_OH"));
	size_t iHC_O_OH=124;
	sp.push_back(species(sp.size(), "C3H8"));
	size_t iC3H8=125;
	sp.push_back(species(sp.size(), "C2H5OH"));
	size_t iC2H5OH=126;
	sp.push_back(species(sp.size(), "CH3C_O_CH3"));
	size_t iCH3C_O_CH3=127;
	sp.push_back(species(sp.size(), "CH3C_O_CH2"));
	size_t iCH3C_O_CH2=128;
	sp.push_back(species(sp.size(), "CH3CONO2"));
	size_t iCH3ONO2=129;
	sp.push_back(species(sp.size(), "CH3C_O_O2NO2_PAN_"));
	size_t iCH3C_O_O2NO2_PAN_=130;
	sp.push_back(species(sp.size(), "C2H5ONO2"));
	size_t iC2H5ONO2=131;
	sp.push_back(species(sp.size(), "C2H5OOH"));
	size_t iC2H5OOH=132;
	sp.push_back(species(sp.size(), "CH3C_O_O2"));
	size_t iCH3C_O_O2=133;
	sp.push_back(species(sp.size(), "CH2OH"));
	size_t iCH2OH=134;
	sp.push_back(species(sp.size(), "C2H5O"));
	size_t iC2H5O=135;
	sp.push_back(species(sp.size(), "FNO"));
	size_t iFNO=136;
	sp.push_back(species(sp.size(), "CH2F"));
	size_t iCH2F=137;
	sp.push_back(species(sp.size(), "CH3F"));
	size_t iCH3F=138;
	sp.push_back(species(sp.size(), "CF3OH"));
	size_t iCF3OH=139;
	sp.push_back(species(sp.size(), "CF3O"));
	size_t iCF3O=140;
	sp.push_back(species(sp.size(), "CH3CH2F"));
	size_t iCH3CH2F=141;
	sp.push_back(species(sp.size(), "CH3CF3"));
	size_t iCH3CF3=142;
	sp.push_back(species(sp.size(), "CHFCF3"));
	size_t iCHFCF3=143;
	sp.push_back(species(sp.size(), "CHF2CHF2"));
	size_t iCHF2CHF2=144;
	sp.push_back(species(sp.size(), "CH2FCH2F"));
	size_t iCH2FCH2F=145;
	sp.push_back(species(sp.size(), "CF2CF3"));
	size_t iCF2CF3=146;
	sp.push_back(species(sp.size(), "CHF2CF3"));
	size_t iCHF2CF3=147;
	sp.push_back(species(sp.size(), "CH3CHF2"));
	size_t iCH3CHF2=148;
	sp.push_back(species(sp.size(), "CH2FCF3"));
	size_t iCH2FCF3=149;
	sp.push_back(species(sp.size(), "CHFCH2F"));
	size_t iCHFCH2F=150;
	sp.push_back(species(sp.size(), "CF2CHF2"));
	size_t iCF2CHF2=151;
	sp.push_back(species(sp.size(), "CH2CF3"));
	size_t iCH2CF3=152;
	sp.push_back(species(sp.size(), "CH2FCHF2"));
	size_t iCH2FCHF2=153;
	sp.push_back(species(sp.size(), "CF2O"));
	size_t iCF2O=154;
	sp.push_back(species(sp.size(), "CF3O2"));
	size_t iCF3O2=155;
	sp.push_back(species(sp.size(), "CHCl2"));
	size_t iCHCl2=156;
	sp.push_back(species(sp.size(), "CHCl3"));
	size_t iCHCl3=157;
	sp.push_back(species(sp.size(), "HOClO"));
	size_t iHOClO=158;
	sp.push_back(species(sp.size(), "CH2Cl2"));
	size_t iCH2Cl2=159;
	sp.push_back(species(sp.size(), "ClONO2"));
	size_t iClONO2=160;
	sp.push_back(species(sp.size(), "CH2Cl"));
	size_t iCH2Cl=161;
	sp.push_back(species(sp.size(), "CH2ClF"));
	size_t iCH2ClF=162;
	sp.push_back(species(sp.size(), "CHClF"));
	size_t iCHClF=163;
	sp.push_back(species(sp.size(), "C2Cl4"));
	size_t iC2Cl4=164;
	sp.push_back(species(sp.size(), "CH3OCl"));
	size_t iCH3OCl=165;
	sp.push_back(species(sp.size(), "CH2CFCl2"));
	size_t iCH2CFCl2=166;
	sp.push_back(species(sp.size(), "CH3CCl3"));
	size_t iCH3CCl3=167;
	sp.push_back(species(sp.size(), "CH2CCl3"));
	size_t iCH2CCl3=168;
	sp.push_back(species(sp.size(), "C2HCl3"));
	size_t iC2HCl3=169;
	sp.push_back(species(sp.size(), "CCl3CHO"));
	size_t iCCl3CHO=170;
	sp.push_back(species(sp.size(), "CHCl2CF2Cl"));
	size_t iCHCl2CF2Cl=171;
	sp.push_back(species(sp.size(), "CCl2CF2Cl"));
	size_t iCCl2CF2Cl=172;
	sp.push_back(species(sp.size(), "CCl3CO"));
	size_t iCCl3CO=173;
	sp.push_back(species(sp.size(), "CH2CF2Cl"));
	size_t iCH2CF2Cl=174;
	sp.push_back(species(sp.size(), "CH2ClCF2Cl"));
	size_t iCH2ClCF2Cl=175;
	sp.push_back(species(sp.size(), "CHClCF2Cl"));
	size_t iCHClCF2Cl=176;
	sp.push_back(species(sp.size(), "CHFClCFCl2"));
	size_t iCHFClCFCl2=177;
	sp.push_back(species(sp.size(), "CHCl2CF3"));
	size_t iCHCl2CF3=178;
	sp.push_back(species(sp.size(), "C2H4O2"));
	size_t iC2H4O2=179;
	sp.push_back(species(sp.size(), "CFClCF3"));
	size_t iCFClCF3=180;
	sp.push_back(species(sp.size(), "CCl2CF3"));
	size_t iCCl2CF3=181;
	sp.push_back(species(sp.size(), "CFClCFCl2"));
	size_t iCFClCFCl2=182;
	sp.push_back(species(sp.size(), "CHFClCF3"));
	size_t iCHFClCF3=183;
	sp.push_back(species(sp.size(), "CH2ClCF3"));
	size_t iCH2ClCF3=184;
	sp.push_back(species(sp.size(), "CFClCF2Cl"));
	size_t iCFClCF2Cl=185;
	sp.push_back(species(sp.size(), "CHClCF3"));
	size_t iCHClCF3=186;
	sp.push_back(species(sp.size(), "CHFClCF2Cl"));
	size_t iCHFClCF2Cl=187;
	sp.push_back(species(sp.size(), "CF3H"));
	size_t iCF3H=188;
	sp.push_back(species(sp.size(), "CHFCl"));
	size_t iCHFCl=189;
	sp.push_back(species(sp.size(), "CBr3"));
	size_t iCBr3=190;
	sp.push_back(species(sp.size(), "CH3CHF"));
	size_t iCH3CHF=191;
	sp.push_back(species(sp.size(), "CHBr2"));
	size_t iCHBr2=192;
	sp.push_back(species(sp.size(), "CHF2Br"));
	size_t iCHF2Br=193;
	sp.push_back(species(sp.size(), "CH2Br"));
	size_t iCH2Br=194;
	sp.push_back(species(sp.size(), "CF2Br"));
	size_t iCF2Br=195;
	sp.push_back(species(sp.size(), "CHClBr"));
	size_t iCHClBr=196;
	sp.push_back(species(sp.size(), "CH3CO3NO2"));
	size_t iCH3CO3NO2=197;
	sp.push_back(species(sp.size(), "CH2CH2F"));
	size_t iCH2CH2F=198;
	sp.push_back(species(sp.size(), "CH3CF2"));
	size_t iCH3CF2=199;
	sp.push_back(species(sp.size(), "CH2CHF2"));
	size_t iCH2CHF2=200;
	sp.push_back(species(sp.size(), "CH2ClBr"));
	size_t iCH2ClBr=201;
	sp.push_back(species(sp.size(), "CHFCHF2"));
	size_t iCHFCHF2=202;
	sp.push_back(species(sp.size(), "CHFBrCF3"));
	size_t iCHFBrCF3=203;
	sp.push_back(species(sp.size(), "CH2Br2"));
	size_t iCH2Br2=204;
	sp.push_back(species(sp.size(), "CClBrCF3"));
	size_t iCClBrCF3=205;
	sp.push_back(species(sp.size(), "CH2FCF2"));
	size_t iCH2FCF2=206;
	sp.push_back(species(sp.size(), "CH2BrCF3"));
	size_t iCH2BrCF3=207;
	sp.push_back(species(sp.size(), "CFBrCF3"));
	size_t iCFBrCF3=208;
	sp.push_back(species(sp.size(), "CF2BrCFCl"));
	size_t iCF2BrCFCl=209;
	sp.push_back(species(sp.size(), "CHBrCF3"));
	size_t iCHBrCF3=210;
	sp.push_back(species(sp.size(), "CHClBrCF3"));
	size_t iCHClBrCF3=211;
	sp.push_back(species(sp.size(), "CF2BrCHFCl"));
	size_t iCF2BrCHFCl=212;
	sp.push_back(species(sp.size(), "CF2BrCF2Br"));
	size_t iCF2BrCF2Br=213;
	sp.push_back(species(sp.size(), "CH2BrO2"));
	size_t iCH2BrO2=214;
	sp.push_back(species(sp.size(), "I2"));
	size_t iI2=215;
	sp.push_back(species(sp.size(), "IO"));
	size_t iIO=216;
	sp.push_back(species(sp.size(), "I"));
	size_t iI=217;
	sp.push_back(species(sp.size(), "HOI"));
	size_t iHOI=218;
	sp.push_back(species(sp.size(), "INO"));
	size_t iINO=219;
	sp.push_back(species(sp.size(), "CS"));
	size_t iCS=220;
	sp.push_back(species(sp.size(), "HI"));
	size_t iHI=221;
	sp.push_back(species(sp.size(), "CH3I"));
	size_t iCH3I=222;
	sp.push_back(species(sp.size(), "CH2I"));
	size_t iCH2I=223;
	sp.push_back(species(sp.size(), "CF3I"));
	size_t iCF3I=224;
	sp.push_back(species(sp.size(), "CS2"));
	size_t iCS2=225;
	sp.push_back(species(sp.size(), "CH3SCH3"));
	size_t iCH3SCH3=226;
	sp.push_back(species(sp.size(), "CH3SO"));
	size_t iCH3SO=227;
	sp.push_back(species(sp.size(), "FSH"));
	size_t iFSH=228;
	sp.push_back(species(sp.size(), "F2"));
	size_t iF2=229;
	sp.push_back(species(sp.size(), "BrSH"));
	size_t iBrSH=230;
	sp.push_back(species(sp.size(), "ClSH"));
	size_t iClSH=231;
	sp.push_back(species(sp.size(), "CH2SCH3"));
	size_t iCH2SCH3=232;
	sp.push_back(species(sp.size(), "HSO"));
	size_t iHSO=233;
	sp.push_back(species(sp.size(), "SO3"));
	size_t iSO3=234;
	sp.push_back(species(sp.size(), "CH3S"));
	size_t iCH3S=235;
	sp.push_back(species(sp.size(), "CH3SH"));
	size_t iCH3SH=236;
	sp.push_back(species(sp.size(), "CH3SCH2"));
	size_t iCH3SCH2=237;
	sp.push_back(species(sp.size(), "CH3SSCH3"));
	size_t iCH3SSCH3=238;
	sp.push_back(species(sp.size(), "Na"));
	size_t iNa=239;
	sp.push_back(species(sp.size(), "NaO"));
	size_t iNaO=240;
	sp.push_back(species(sp.size(), "NaO2"));
	size_t iNaO2=241;
	sp.push_back(species(sp.size(), "NaCl"));
	size_t iNaCl=242;
	sp.push_back(species(sp.size(), "HSO2"));
	size_t iHSO2=243;
	sp.push_back(species(sp.size(), "HOSO2"));
	size_t iHOSO2=244;
	sp.push_back(species(sp.size(), "CH2SH"));
	size_t iCH2SH=245;
	sp.push_back(species(sp.size(), "HOCH2CH2"));
	size_t iHOCH2CH2=246;
	sp.push_back(species(sp.size(), "NaOH"));
	size_t iNaOH=247;
	sp.push_back(species(sp.size(), "CH3SS"));
	size_t iCH3SS=248;
	sp.push_back(species(sp.size(), "CH3SSO"));
	size_t iCH3SSO=249;
	sp.push_back(species(sp.size(), "CH3SO2"));
	size_t iCH3SO2=250;
	sp.push_back(species(sp.size(), "CH3SOO"));
	size_t iCH3SOO=251;
	sp.push_back(species(sp.size(), "CH3SCH2O2"));
	size_t iCH3SCH2O2=252;
	sp.push_back(species(sp.size(), "HOCHCH"));
	size_t iHOCHCH=253;
	sp.push_back(species(sp.size(), "CH3SCH2O"));
	size_t iCH3SCH2O=254;
	sp.push_back(species(sp.size(), "CH3ONO"));
	size_t iCH3ONO=255;
	sp.push_back(species(sp.size(), "C2H5ONO"));
	size_t iC2H5ONO=256;
	sp.push_back(species(sp.size(), "FNO2"));
	size_t iFNO2=257;
	sp.push_back(species(sp.size(), "FONO2"));
	size_t iFONO2=258;
	sp.push_back(species(sp.size(), "ClCO"));
	size_t iClCO=259;
	sp.push_back(species(sp.size(), "ClC2H2"));
	size_t iClC2H2=260;
	sp.push_back(species(sp.size(), "INO2"));
	size_t iINO2=261;
	sp.push_back(species(sp.size(), "IONO2"));
	size_t iIONO2=262;
	sp.push_back(species(sp.size(), "NaO3"));
	size_t iNaO3=263;
	sp.push_back(species(sp.size(), "NaCO3"));
	size_t iNaCO3=264;
	sp.push_back(species(sp.size(), "BrONO2"));
	size_t iBrONO2=265;
	sp.push_back(species(sp.size(), "NaHCO3"));
	size_t iNaHCO3=266;
	sp.push_back(species(sp.size(), "CF2ClO2"));
	size_t iCF2ClO2=267;
	sp.push_back(species(sp.size(), "C2Cl5"));
	size_t iC2Cl5=268;
	sp.push_back(species(sp.size(), "ClC2H4"));
	size_t iClC2H4=269;
	sp.push_back(species(sp.size(), "CH2ClO2"));
	size_t iCH2ClO2=270;
	sp.push_back(species(sp.size(), "CFCl2O2"));
	size_t iCFCl2O2=271;
	sp.push_back(species(sp.size(), "CHCl2O2"));
	size_t iCHCl2O2=272;
	sp.push_back(species(sp.size(), "CCl3O2"));
	size_t iCCl3O2=273;
	sp.push_back(species(sp.size(), "N2O3"));
	size_t iN2O3=274;
	sp.push_back(species(sp.size(), "N2O4"));
	size_t iN2O4=275;
	sp.push_back(species(sp.size(), "FOO"));
	size_t iFOO=276;
	sp.push_back(species(sp.size(), "ClO_O2"));
	size_t iClO_O2=277;
	sp.push_back(species(sp.size(), "FO"));
	size_t iFO=278;
	sp.push_back(species(sp.size(), "FO2"));
	size_t iFO2=279;
	
	
	// R1
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 2);
	re.back().photodiss(0.0, 1.2e-12, 3.0e-9);

	// R2
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iO1D, 1);
	re.back().photodiss(0.0, 0.0, 4.7e-10);

	// R3
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iO, 1); 
	re.back().photodiss(2.1e-4, 2.1e-4, 6.2e-4);

	// R4
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iO1D, 1); 
	re.back().photodiss(1.2e-5, 1.5e-5, 3.7e-3);

	// R5
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO, 1); 
	re.back().photodiss(0.0, 5.5e-7, 2.8e-4);

	// R6
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iOH, 1); 
	re.back().photodiss(0.0, 2.4e-12, 1.9e-6);

	// R7
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iO1D, 1); 
	re.back().photodiss(0.0, 0.0, 2.1e-7);

	// R8
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 2);
	re.back().add_product(iO, 1); 
	re.back().photodiss(0.0, 0.0, 2.5e-7);

	// R9
	re.push_back(reaction()); 
	re.back().add_educt(iH2O2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 2); 
	re.back().photodiss(3.4e-6, 3.3e-6, 4.8e-5);

	// R10
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iN, 1);
	re.back().add_product(iO, 1); 
	re.back().photodiss(0.0, 7.4e-11, 2.4e-6);

	// R11
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO, 1); 
	re.back().photodiss(4.6e-3, 4.7e-3, 4.7e-3);

	// R12
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO, 1); 
	re.back().photodiss(7.7e-2, 7.8e-2, 7.8e-2);

	// R13
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1); 
	re.back().photodiss(6.2e-3, 6.3e-3, 6.3e-3);

	// R14
	re.push_back(reaction()); 
	re.back().add_educt(iN2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO1D, 1); 
	re.back().photodiss(0.0, 1.8e-9, 5.0e-7);

	// R15
	re.push_back(reaction()); 
	re.back().add_educt(iN2O5, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iNO3, 1); 
	re.back().photodiss(2.1e-5, 1.0e-6, 9.4e-5);

	// R16
	re.push_back(reaction()); 
	re.back().add_educt(iN2O5, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iNO3, 1);
	re.back().add_product(iO, 1); 
	re.back().photodiss(1.2e-6, 1.4e-6, 2.2e-4);

	// R17
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iNO2, 1); 
	re.back().photodiss(1.2e-6, 2.0e-6, 1.5e-4);

	// R18
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iNO3, 1); 
	re.back().photodiss(6.5e-7, 1.1e-6, 8.0e-5);

	// R19
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2NO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCH3O2, 1);
	re.back().add_product(iNO2, 1); 
	re.back().photodiss(5.3e-6, 8.5e-6, 2.3e-4);

	// R20
	re.push_back(reaction()); 
	re.back().add_educt(iNH3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNH2, 1);
	re.back().add_product(iH, 1); 
	re.back().photodiss(0.0, 4.1e-7, 5.4e-5);

	// R21
	re.push_back(reaction()); 
	re.back().add_educt(iCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 2); 
	re.back().photodiss(1.2e-3, 1.3e-3, 1.4e-3);

	// R22
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iO, 1); 
	re.back().photodiss(6.2e-5, 8.9e-5, 3.0e-3);

	// R23
	re.push_back(reaction()); 
	re.back().add_educt(iClOO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iO, 1); 
	re.back().photodiss(0.0, 4.3e-10, 4.6e-3);

	// R24
	re.push_back(reaction()); 
	re.back().add_educt(iOClO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iO, 1); 
	re.back().photodiss(4.2e-2, 4.5e-2, 4.3e-2);

	// R25
	re.push_back(reaction()); 
	re.back().add_educt(iClO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iO2, 1); 
	re.back().photodiss(2.8e-4, 2.8e-4, 2.8e-4);

	// R26
	re.push_back(reaction()); 
	re.back().add_educt(iCl2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iClO, 1); 
	re.back().photodiss(5.1e-4, 5.6e-4, 1.9e-3);

	// R27
	re.push_back(reaction()); 
	re.back().add_educt(iCl2O2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iClOO, 1); 
	re.back().photodiss(7.1e-4, 7.8e-4, 3.1e-3);

	// R28
	re.push_back(reaction()); 
	re.back().add_educt(iCl2O3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iClOO, 1); 
	re.back().photodiss(2.7e-4, 4.0e-4, 1.0e-2);

	// R29
	re.push_back(reaction()); 
	re.back().add_educt(iClNO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iNO, 1); 
	re.back().photodiss(1.6e-3, 1.6e-3, 2.6e-3);

	// R30
	re.push_back(reaction()); 
	re.back().add_educt(iClNO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iNO2, 1); 
	re.back().photodiss(1.8e-4, 2.0e-4, 7.8e-4);

	// R31
	re.push_back(reaction()); 
	re.back().add_educt(iClONO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iNO2, 1); 
	re.back().photodiss(2.2e-3, 2.4e-3, 4.0e-3);

	// R32
	re.push_back(reaction()); 
	re.back().add_educt(iClNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iNO3, 1); 
	re.back().photodiss(4.0e-6, 4.0e-6, 2.1e-4);

	// R33
	re.push_back(reaction()); 
	re.back().add_educt(iClNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iNO2, 1); 
	re.back().photodiss(2.2e-5, 2.2e-5, 2.1e-4);

	// R34
	re.push_back(reaction()); 
	re.back().add_educt(iClNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iClONO, 1); 
	re.back().photodiss(2.5e-6, 2.1e-6, 3.6e-5);

	// R35
	re.push_back(reaction()); 
	re.back().add_educt(iHCl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iCl, 1); 
	re.back().photodiss(0.0, 1.7e-9, 6.1e-7);

	// R36
	re.push_back(reaction()); 
	re.back().add_educt(iHOCl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iCl, 1); 
	re.back().photodiss(1.3e-4, 1.4e-4, 2.3e-4);

	// R37
	re.push_back(reaction()); 
	re.back().add_educt(iCH3Cl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iCH3, 1); 
	re.back().photodiss(0.0, 1.1e-9, 5.1e-7);

	// R38
	re.push_back(reaction()); 
	re.back().add_educt(iCH2FCl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 1.2e-10, 7.8e-8);

	// R39
	re.push_back(reaction()); 
	re.back().add_educt(iCHFCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCFCl2, 1);
	re.back().add_product(iH, 1); 
	re.back().photodiss(0.0, 3.5e-9, 1.4e-6);

	// R40
	re.push_back(reaction()); 
	re.back().add_educt(iCHF2Cl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCF2Cl, 1);
	re.back().add_product(iH, 1); 
	re.back().photodiss(0.0, 1.1e-11, 1.1e-8);

	// R41
	re.push_back(reaction()); 
	re.back().add_educt(iCCl4, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCCl3, 1);
	re.back().add_product(iCl, 1); 
	re.back().photodiss(0.0, 6.7e-8, 1.8e-5);

	// R42
	re.push_back(reaction()); 
	re.back().add_educt(iCFCl3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCFCl2, 1);
	re.back().add_product(iCl, 1); 
	re.back().photodiss(0.0, 3.8e-8, 8.0e-6);

	// R43
	re.push_back(reaction()); 
	re.back().add_educt(iCF2Cl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCF2Cl, 1);
	re.back().add_product(iCl, 1); 
	re.back().photodiss(0.0, 3.1e-9, 1.9e-6);

	// R44
	re.push_back(reaction()); 
	re.back().add_educt(iCH3CCl3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 5.7e-8, 9.2e-6);

	// R45
	re.push_back(reaction()); 
	re.back().add_educt(iCH3CFCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 9.8e-9, 1.4e-6);

	// R46
	re.push_back(reaction()); 
	re.back().add_educt(iCH3CF2Cl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 8.6e-11, 1.3e-8);

	// R47
	re.push_back(reaction()); 
	re.back().add_educt(iC2Cl3F3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 6.3e-9, 1.6e-6);

	// R48
	re.push_back(reaction()); 
	re.back().add_educt(iC2Cl2F4, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 3.5e-10, 2.1e-7);

	// R49
	re.push_back(reaction()); 
	re.back().add_educt(iCF3CHCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 7.5e-9, 1.1e-6);

	// R50
	re.push_back(reaction()); 
	re.back().add_educt(iCF3CHFCl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 7.9e-11, 1.2e-8);

	// R51
	re.push_back(reaction()); 
	re.back().add_educt(iCF3CCl3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 7.2e-9, 2.6e-6);

	// R52
	re.push_back(reaction()); 
	re.back().add_educt(iCF3CFCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 4.9e-10, 2.1e-7);

	// R53
	re.push_back(reaction()); 
	re.back().add_educt(iCF3CF2CHCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 0.0, 0.0);

	// R54
	re.push_back(reaction()); 
	re.back().add_educt(iCF2ClCF2CHFCl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 9.8e-9, 1.4e-6);

	// R55
	re.push_back(reaction()); 
	re.back().add_educt(iCHClO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iHCO, 1);
	re.back().add_product(iCl, 1);
	re.back().photodiss(2.8e-4, 2.8e-4, 2.8e-4);

	// R56
	re.push_back(reaction()); 
	re.back().add_educt(iCOCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 2);
	re.back().add_product(iCO, 1);
	re.back().photodiss(8.4e-19, 2.5e-8, 2.5e-5);

	// R57
	re.push_back(reaction()); 
	re.back().add_educt(iCHClO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iHCO, 1);
	re.back().add_product(iCl, 1);
	re.back().photodiss(2.8e-4, 2.8e-4, 2.8e-4);

	// R58
	re.push_back(reaction()); 
	re.back().add_educt(iCOFCl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 1.3e-8, 3.2e-6);

	// R59
	re.push_back(reaction()); 
	re.back().add_educt(iHF, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iF, 1);
	re.back().photodiss(0.0, 0.0, 1.1e-7);

	// R60
	re.push_back(reaction()); 
	re.back().add_educt(iCF4, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCF3, 1);
	re.back().add_product(iF, 1);
	re.back().photodiss(0.0, 0.0, 2.0e-10);

	// R61
	re.push_back(reaction()); 
	re.back().add_educt(iC2F6, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCF3, 2);
	re.back().photodiss(0.0, 0.0, 4.6e-20);

	// R62
	re.push_back(reaction()); 
	re.back().add_educt(iCOF2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 9.9e-10, 1.6e-7);

	// R63
	re.push_back(reaction()); 
	re.back().add_educt(iSF6, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 0.0, 0.0);

	// R64
	re.push_back(reaction()); 
	re.back().add_educt(iBr2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iBr, 2);
	re.back().photodiss(1.7e-2, 1.6e-2, 1.6e-2);

	// R65
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iBr, 1);
	re.back().add_product(iO, 1);
	re.back().photodiss(2.1e-2, 2.3e-2, 2.3e-2);

	// R66
	re.push_back(reaction()); 
	re.back().add_educt(iBrNO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iBr, 1);
	re.back().add_product(iNO2, 1);
	re.back().photodiss(0.0, 0.0, 0.0);

	// R67
	re.push_back(reaction()); 
	re.back().add_educt(iBrNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iBr, 1);
	re.back().add_product(iNO3, 1);
	re.back().photodiss(5.2e-4, 5.5e-4, 1.1e-3);

	// R68
	re.push_back(reaction()); 
	re.back().add_educt(iHBr, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iBr, 1);
	re.back().photodiss(0.0, 1.1e-7, 2.2e-5);

	// R69
	re.push_back(reaction()); 
	re.back().add_educt(iHOBr, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iBr, 1);
	re.back().photodiss(3.7e-4, 3.8e-4, 6.3e-4);

	// R70
	re.push_back(reaction()); 
	re.back().add_educt(iCH3Br, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCH3, 1);
	re.back().add_product(iBr, 1);
	re.back().photodiss(0.0, 8.9e-8, 2.6e-5);

	// R71
	re.push_back(reaction()); 
	re.back().add_educt(iCHBr3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(1.6e-7, 8.8e-7, 6.2e-4);

	// R72
	re.push_back(reaction()); 
	re.back().add_educt(iBrCl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iBr, 1);
	re.back().add_product(iCl, 1);
	re.back().photodiss(4.8e-3, 4.8e-3, 4.7e-3);

	// R73
	re.push_back(reaction()); 
	re.back().add_educt(iCF2ClBr, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCF2Cl, 1);
	re.back().add_product(iBr, 1);
	re.back().photodiss(1.5e-15, 1.4e-7, 5.9e-5);

	// R74
	re.push_back(reaction()); 
	re.back().add_educt(iCF2Br2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(6.4e-10, 1.2e-7, 2.3e-4);

	// R75
	re.push_back(reaction()); 
	re.back().add_educt(iCF3Br, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCF3, 1);
	re.back().add_product(iBr, 1);
	re.back().photodiss(0.0, 1.4e-8, 4.8e-6);

	// R76
	re.push_back(reaction()); 
	re.back().add_educt(iC2F4Br2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 1.4e-7, 5.0e-5);

	// R77
	re.push_back(reaction()); 
	re.back().add_educt(iCH4, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iProducts, 1);
	re.back().photodiss(0.0, 0.0, 6.6e-7);

	// R78
	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iO, 1);
	re.back().photodiss(0.0, 5.3e-14, 3.2e-9);

	// R79
	re.push_back(reaction()); 
	re.back().add_educt(iH2CO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iHCO, 1);
	re.back().add_product(iH, 1);
	re.back().photodiss(1.5e-5, 1.8e-5, 3.8e-5);

	// R80
	re.push_back(reaction()); 
	re.back().add_educt(iH2CO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iCO, 1);
	re.back().photodiss(2.6e-5, 3.8e-5, 5.6e-5);

	// R81
	re.push_back(reaction()); 
	re.back().add_educt(iCH3OOH, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCH3O, 1);
	re.back().add_product(iOH, 1);
	re.back().photodiss(2.6e-6, 3.0e-6, 2.7e-5);

	// R82
	re.push_back(reaction()); 
	re.back().add_educt(iHCN, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iCN, 1);
	re.back().photodiss(0.0, 8.3e-13, 1.8e-6);

	// R83
	re.push_back(reaction()); 
	re.back().add_educt(iCH3CN, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCH3, 1);
	re.back().add_product(iCN, 1);
	re.back().photodiss(0.0, 0.0, 7.0e-7);

	// R84
	re.push_back(reaction()); 
	re.back().add_educt(iSO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iO, 1);
	re.back().photodiss(0.0, 6.6e-7, 9.4e-5);

	// R85
	re.push_back(reaction()); 
	re.back().add_educt(iH2S, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iSH, 1);
	re.back().photodiss(0.0, 5.9e-7, 1.6e-4);

	// R86
	re.push_back(reaction()); 
	re.back().add_educt(iOCS, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iS, 1);
	re.back().photodiss(5.4e-13, 9.0e-9, 1.3e-5);





	// R87
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO2, 2);
	re.back().sorder(8.0e-12, 2060,250, 8.0e-15, 1.15);

	// R88
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(3.2e-11, -70,100, 4.0e-11, 1.2);

	// R89
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO2, 2);
	re.back().sorder(1.2e-10, 0,100, 1.2e-10, 1.3);

	// R90
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iO, 2);
	re.back().sorder(1.2e-10, 0,100, 1.2e-10, 1.3);

	// R91
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iH, 1);
	re.back().sorder(1.1e-10, 0,100, 1.1e-10, 1.1);

	// R92
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iOH, 2);
	re.back().sorder(2.2e-10, 0,100, 2.2e-10, 1.2);

	// R93
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iN2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iN2, 1);
	re.back().sorder(1.8e-11, -110,100, 2.6e-11, 1.2);

	// R94
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iN2O, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(4.9e-11, 0,100, 4.9e-11, 1.3);

	// R95
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iN2O, 1);
	re.back().add_product(iNO, 2);
	re.back().sorder(6.7e-11, 0,100, 6.7e-11, 1.3);

	// R96
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iNH3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iNH2, 1);
	re.back().sorder(2.5e-10, 0,100, 2.5e-10, 1.3);

	// R97
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iCO2, 1);
	re.back().sorder(7.4e-11, -120,100, 1.1e-10, 1.2);

	// R98
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iCH3, 1);
	re.back().sorder(1.2e-10, 0,100, 1.2e-10, 1.2);

	// R99
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iCH2O, 1);
	re.back().sorder(3.0e-11, 0,100, 3.0e-11, 1.5);

	// R100
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.5e-10, 0,100, 1.5e-10, 1.2);

	// R101
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iHF, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iF, 1);
	re.back().sorder(1.4e-10, 0,100, 1.4e-10, 2.0);

	// R102
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iHBr, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.5e-10, 0,100, 1.5e-10, 2.0);
	
	// R103
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iSF6, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0,100, 1.8e-14, 1.5);
	
	// R104
	re.push_back(reaction()); 
	re.back().add_educt(iO21d, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-16, 0);
	
	// R105
	re.push_back(reaction()); 
	re.back().add_educt(iO21d, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(3.6e-18, 220,100, 1.7e-18, 1.2);
	
	// R106
	re.push_back(reaction()); 
	re.back().add_educt(iO21d, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iO2, 2);
	re.back().sorder(5.2e-11, 2840,500, 3.8e-15, 1.2);

	// R107
	re.push_back(reaction()); 
	re.back().add_educt(iO21d, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 4.8e-18, 1.5);

	// R108
	re.push_back(reaction()); 
	re.back().add_educt(iO21d, 1);
	re.back().add_educt(iN2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1e-20, 0);

	// R109
	re.push_back(reaction()); 
	re.back().add_educt(iO21d, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2e-20, 0);

	// R110
	re.push_back(reaction()); 
	re.back().add_educt(iO21s, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 8e-14, 5.0);

	// R111
	re.push_back(reaction()); 
	re.back().add_educt(iO21s, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.9e-17, 1.5);

	// R112
	re.push_back(reaction()); 
	re.back().add_educt(iO21s, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.2e-11, 0, 200, 2.2e-11, 1.2);

	// R113
	re.push_back(reaction()); 
	re.back().add_educt(iO21s, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 5.4e-12, 1.3);

	// R114
	re.push_back(reaction()); 
	re.back().add_educt(iO21s, 1);
	re.back().add_educt(iN, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1e-13, 0);

	// R115
	re.push_back(reaction()); 
	re.back().add_educt(iO21s, 1);
	re.back().add_educt(iN2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.1e-15, 0, 200, 2.1e-15, 1.2);

	// R116
	re.push_back(reaction()); 
	re.back().add_educt(iO21s, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.2e-13, 0, 200, 4.2e-13, 1.2);

	// R117
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iH, 1);
	re.back().sorder(2.2e-11, -120, 100, 3.3e-11, 1.2);

	// R118
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(3.0e-11, -200, 100, 5.9e-11, 1.2);
	
	// R119
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iH2O2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iHO2, 1);
	re.back().sorder(1.4e-12, 2000, 1000, 1.7e-15, 2.0);

	// R120
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.4e-10, 470, 200, 2.9e-11, 1.25);

	// R121
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(8.1e-11, 0, 100, 8.1e-11, 1.3);

	// R122
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.6e-11, 940, 300, 6.8e-14, 1.3);

	// R123
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iH, 1);
	re.back().sorder(5.5e-12, 2000, 100, 6.7e-15, 1.1);

	// R124
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHD, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(5.0e-12, 2130, 200, 4.0e-15, 1.2);

	// R125
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 2);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO, 1);
	re.back().sorder(4.2e-12, 240, 240, 1.9e-12, 1.4);
	
	// R126
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(4.8e-11, -250, 200, 1.1e-10, 1.3);

	// R127
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iH2O2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iHO2, 1);
	re.back().sorder(2.9e-12, 160, 100, 1.7e-12, 1.2);

	// R128
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 2);
	re.back().sorder(1.1e-14, 500, 500, 2.0e-15, 1.3);

	// R129
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 2);
	re.back().add_product(iH2O2, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.3e-13, -600, 200, 1.7e-12, 1.3);

	// R130
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(6.5e-12, -120, 120, 9.7e-12, 1.1);

	// R131
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(1.0e-11, 0, 150, 1.0e-11, 1.5);

	// R132
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iN2O5, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.0e-16, 0);

	// R133
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHNO3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iNO3, 1);
	re.back().sorder(0, 0, 0, 3.0e-17, 0);

	// R134
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(7.8e-11, 3400, 750, 8.6e-16, 3.0);

	// R135
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iNO, 1);
	re.back().sorder(4.0e-10, 340, 300, 1.3e-10, 1.3);

	// R136
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.2e-11, 1.5);

	// R137
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHONO, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(1.8e-11, 390, 200, 4.5e-12, 1.5);

	// R138
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHNO3, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iNO3, 1);
	re.back().sorder(0, 0, 0, 0, 1.3);

	// R139
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.3e-12, -380, 270, 4.6e-12, 1.5);

	// R140
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iNH3, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iNH2, 1);
	re.back().sorder(1.7e-12, 710, 200, 1.6e-13, 1.2);

	// R141
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iOH, 1);
	re.back().sorder(3.5e-12, -250, 50, 8.1e-12, 1.15);

	// R142
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iHNO3, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(0, 0, 0, 3.5e-12, 1.5);

	// R143
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO, 1);
	re.back().sorder(1.5e-11, 3600, 400, 8.5e-17, 1.25);

	// R144
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(0, 0, 0, 2.0e-16, 0);

	// R145
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO, 1);
	re.back().sorder(2.1e-11, -100, 100, 3.0e-11, 1.3);

	// R146
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iN2O, 1);
	re.back().add_product(iO, 1);
	re.back().sorder(5.8e-12, -220, 100, 1.2e-11, 1.5);

	// R147
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.0e-12, 1400, 200, 1.8e-14, 1.1);

	// R148
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iNO2, 2);
	re.back().sorder(1.5e-11, -170, 100, 2.6e-11, 1.3);

	// R149
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO3, 2);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.2e-13, 2450, 150, 3.2e-17, 1.15);

	// R150
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 2);
	re.back().add_product(iNO2, 2);
	re.back().add_product(iO2, 1);
	re.back().sorder(8.5e-13, 2450, 500, 2.3e-16, 1.5);

	// R151
	re.push_back(reaction()); 
	re.back().add_educt(iNH2, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 6.0e-21, 0);

	// R152
	re.push_back(reaction()); 
	re.back().add_educt(iNH2, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.3e-12, 930, 500, 1.9e-13, 3.0);

	// R153
	re.push_back(reaction()); 
	re.back().add_educt(iNH2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.0e-12, -450, 150, 1.8e-11, 1.3);

	// R154
	re.push_back(reaction()); 
	re.back().add_educt(iNH2, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.1e-12, -650, 250, 1.9e-11, 3.0);

	// R155
	re.push_back(reaction()); 
	re.back().add_educt(iNH, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.9e-11, 0, 300, 4.9e-11, 1.5);

	// R156
	re.push_back(reaction()); 
	re.back().add_educt(iNH, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(3.5e-13, -1140, 500, 1.6e-11, 2.0);

	// R157
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iHNO2, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iHNO3, 1);
	re.back().sorder(0, 0, 0, 5.0e-19, 0);

	// R158
	re.push_back(reaction()); 
	re.back().add_educt(iN2O5, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iHNO3, 2);
	re.back().sorder(0, 0, 0, 2.0e-21, 0);

	// R159
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.1e-10, 0, 250, 1.1e-10, 1.3);
	
	// R160
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHCN, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.1e-11, 4000, 1000, 1.5e-17, 10);
	
	// R161
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iC2H2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(3.0e-11, 1600, 250, 1.4e-13, 1.3);
	
	// R162
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iH2CO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(3.4e-11, 1600, 250, 1.6e-13, 1.25);
	
	// R163
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCH3CHO, 1);
	re.back().add_product(iCH3CO, 1);
	re.back().add_product(iOH, 1);
	re.back().sorder(1.8e-11, 1100, 200, 4.5e-13, 1.25);
	
	// R164
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iC2H2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-14, 4100, 500, 1e-20, 3.0);

	// R165
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iC2H4, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.2e-14, 2630, 100, 1.7e-18, 1.25);

	// R166
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iC3H6, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(6.5e-15, 1900, 200, 1.1e-17, 1.2);
	
	// R167
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.5e-13, 0, 300, 1.5e-13, 1.3);   // WARNING
	
	// R168
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iCH3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(4.6e-12, 1965, 200, 6.3e-15, 1.1);
	
	// R169
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3D, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.1e-12, 2000, 200, 5.0e-15, 1.15);
	
	// R170
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iH2CO, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iHCO, 1);
	re.back().sorder(2.0e-12, -480, 200, 1.0e-11, 1.25);
	
	// R171
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3OH, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(6.7e-12, 600, 300, 8.9e-13, 1.2);
	
	// R172
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3OOH, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(3.8e-12, -200, 200, 7.4e-12, 1.5);
	
	// R173
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHC_O_OH, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.0e-13, 0, 200, 4.0e-13, 1.3);
	
	// R174
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHCN, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.2e-13, 400, 150, 3.1e-14, 3.0);

	// R175
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iC2H6, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iC2H5, 1);
	re.back().sorder(8.7e-12, 1070, 100, 2.4e-13, 1.1);
	
	// R176
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iC3H8, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iC3H7, 1);
	re.back().sorder(1.0e-11, 660, 100, 1.1e-12, 1.2);
	
	// R177
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3CHO, 1);
	re.back().add_product(iCH3CO, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(5.6e-12, -270, 200, 1.4e-11, 1.2);
	
	// R178
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iC2H5OH, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(7.0e-12, 235, 100, 3.2e-12, 1.3);
	
	// R179
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3C_O_OH, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.0e-13, -200, 400, 8.0e-13, 1.3);
	
	// R180
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3C_O_CH3, 1);
	re.back().add_product(iCH3C_O_CH2, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(2.2e-12, 685, 300, 2.2e-13, 1.15);
	
	// R181
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3CN, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.4e-12, 1664, 200, 1.7e-14, 1.5);
	
	// R182
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3ONO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(5.0e-13, 890, 500, 2.4e-14, 3.0);
	
	// R183
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3C_O_O2NO2_PAN_, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 4.0e-14, 0);
	
	// R184
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iC2H5ONO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(8.2e-13, 450, 300, 1.8e-13, 3.0);
	
	// R185
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iCH2O, 1);
	re.back().add_product(iAdduct, 1);    // WARNING
	re.back().sorder(6.7e-15, -600, 600, 5.0e-14, 5.0);
	
	// R186
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iCH3O2, 1);
	re.back().add_product(iCH3OOH, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(3.8e-13, -800, 400, 5.6e-12, 2.0);
	
	// R187
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iC2H5O2, 1);
	re.back().add_product(iC2H5OOH, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(7.5e-13, -700, 250, 8.0e-12, 1.5);
	
	// R188
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iCH3C_O_O2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.5e-13, -1000, 600, 1.3e-11, 2.0);
	
	// R189
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 4.0e-19, 0);
	
	// R190
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iCH2O, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 5.8e-16, 0);
	
	// R191
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iCH3CHO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.4e-12, 1900, 300, 2.4e-15, 1.3);
	
	// R192
	re.push_back(reaction()); 
	re.back().add_educt(iCH3, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.0e-16, 0);
	
	// R193
	re.push_back(reaction()); 
	re.back().add_educt(iCH3, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(5.4e-12, 220, 150, 2.6e-12, 2.0);
	
	// R194
	re.push_back(reaction()); 
	re.back().add_educt(iHCO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iHO2, 1);
	re.back().sorder(3.5e-12, -140, 140, 5.5e-12, 1.3);
	
	// R195
	re.push_back(reaction()); 
	re.back().add_educt(iCH2OH, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCH2O, 1);
	re.back().add_product(iHO2, 1);
	re.back().sorder(9.1e-12, 0, 200, 9.1e-12, 1.3);
	
	// R196
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCH2O, 1);
	re.back().add_product(iHO2, 1);
	re.back().sorder(3.9e-14, 900, 300, 1.9e-15, 1.5);
	
	// R197
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iCH2O, 1);
	re.back().add_product(iHONO, 1);
	re.back().sorder(1.1e-11, 1200, 600, 2.0e-13, 5.0);
	
	// R198
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.0e-17, 0);
	
	// R199
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 2);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.5e-13, -190, 190, 4.7e-13, 1.5);
	
	// R200
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iCH3O, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(3.0e-12, -280, 60, 7.7e-12, 1.15);
	
	// R201
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 1);
	re.back().add_educt(iCH3C_O_O2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.3e-12, -640, 200, 1.1e-11, 1.5);
	
	// R202
	re.push_back(reaction()); 
	re.back().add_educt(iC2H5, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iC2H4, 1);
	re.back().add_product(iHO2, 1);
	re.back().sorder(0, 0, 0, 2.0e-14, 0);
	
	// R203
	re.push_back(reaction()); 
	re.back().add_educt(iC2H5O, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCH3CHO, 1);
	re.back().add_product(iHO2, 1);
	re.back().sorder(6.3e-14, 550, 200, 1.0e-14, 1.5);
	
	// R204
	re.push_back(reaction()); 
	re.back().add_educt(iC2H5O2, 2);
	re.back().add_product(iProducts, 1);
	re.back().sorder(6.8e-14, 0, 300, 6.8e-14, 2.0);
	
	// R205
	re.push_back(reaction()); 
	re.back().add_educt(iC2H5O2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.6e-12, -365, 150, 8.7e-12, 1.2);
	
	// R206
	re.push_back(reaction()); 
	re.back().add_educt(iCH3C_O_O2, 2);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.9e-12, -500, 150, 1.5e-11, 1.5);
	
	// R207
	re.push_back(reaction()); 
	re.back().add_educt(iCH3C_O_O2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(5.3e-12, -360, 150, 1.8e-11, 1.4);
	
	// R208
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iFO, 1);
	re.back().add_product(iF, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.7e-11, 0, 250, 2.7e-11, 3.0);
	
	// R209
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iFO2, 1);
	re.back().add_product(iFO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(5.0e-11, 0, 250, 5.0e-11, 5.0);
	
	// R210
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3F, 1);
	re.back().add_product(iCH2F, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(3.0e-12, 1500, 300, 2.0e-14, 1.1);
	
	// R211
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2F2, 1);
	re.back().add_product(iCHF2, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.9e-12, 1550, 200, 1.0e-14, 1.2);
	
	// R212
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHF3, 1);
	re.back().add_product(iCHF2, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.0e-12, 2440, 200, 2.8e-16, 1.3);
	
	// R213
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCF3OH, 1);
	re.back().add_product(iCF3O, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(0, 0, 0, 2.0e-17, 0);
	
	// R214
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3CH2F, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(7.0e-12, 1100, 300, 1.7e-13, 1.4);
	
	// R215
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3CHF2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.4e-12, 1260, 200, 3.5e-14, 1.2);
	
	// R216
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2FCH2F, 1);
	re.back().add_product(iCHFCH2F, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(4.0e-12, 1070, 500, 1.1e-13, 2.0);
	
	// R217
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3CF3, 1);
	re.back().add_product(iCH2CF3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.8e-12, 2170, 150, 1.2e-15, 1.1);
		
	// R218
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2FCHF2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.0e-12, 1650, 300, 1.6e-14, 1.5);
	
	// R219
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2FCF3, 1);
	re.back().add_product(iCHFCF3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.5e-12, 1750, 200, 4.2e-15, 1.1);
	
	// R220
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHF2CHF2, 1);
	re.back().add_product(iCF2CHF2, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.6e-12, 1680, 300, 5.7e-15, 2.0);
	
	// R221
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHF2CF3, 1);
	re.back().add_product(iCF2CF3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(5.6e-13, 1700, 300, 1.9e-15, 1.3);
	
	// R222
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iFO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.2e-11, 230, 200, 1.0e-11, 1.5);
	
	// R223
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iHF, 1);
	re.back().add_product(iH, 1);
	re.back().sorder(1.4e-10, 500, 200, 2.6e-11, 1.2);
	
	// R224
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iHF, 1);
	re.back().add_product(iOH, 1);
	re.back().sorder(1.4e-11, 0, 200, 1.4e-11, 1.3);
	
	// R225
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iHNO3, 1);
	re.back().add_product(iHF, 1);
	re.back().add_product(iNO3, 1);
	re.back().sorder(6.0e-12, -400, 200, 2.3e-11, 1.3);
	
	// R226
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iHF, 1);
	re.back().add_product(iCH3, 1);
	re.back().sorder(1.6e-10, 260, 200, 6.7e-11, 1.4);
	
	// R227
	re.push_back(reaction()); 
	re.back().add_educt(iFO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-14, 0);
	
	// R228
	re.push_back(reaction()); 
	re.back().add_educt(iFO, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iF, 1);
	re.back().sorder(8.2e-12, -300, 200, 2.2e-11, 1.5);
	
	// R229
	re.push_back(reaction()); 
	re.back().add_educt(iFO, 2);
	re.back().add_product(iF, 2);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.0e-11, 0, 250, 1.0e-11, 1.5);
	
	// R230
	re.push_back(reaction()); 
	re.back().add_educt(iFO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.4e-16, 0);
	
	// R231
	re.push_back(reaction()); 
	re.back().add_educt(iFO2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iFNO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(7.5e-12, 690, 400, 7.5e-13, 2.0);
	
	// R232
	re.push_back(reaction()); 
	re.back().add_educt(iFO2, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(3.8e-11, 2040, 500, 4.0e-14, 2.0);
	
	// R233
	re.push_back(reaction()); 
	re.back().add_educt(iFO2, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 5.1e-16, 0);
	
	// R234
	re.push_back(reaction()); 
	re.back().add_educt(iFO2, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-16, 0);
	
	// R235
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iFO2, 1);
	re.back().add_product(iCF2O, 1);
	re.back().sorder(3.0e-11, 5000, 0, 1.5e-18, 0);
	
	// R236
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iCF3O2, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.0e-12, 1400, 600, 1.8e-14, 1.3);
	
	// R237
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iCF3OH, 1);
	re.back().sorder(3.0e-12, 3600, 0, 2.0e-17, 0);
	
	// R238
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iCF2O, 1);
	re.back().add_product(iFNO, 1);
	re.back().sorder(3.7e-11, -110, 70, 5.4e-11, 1.2);
	
	// R239
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-15, 0);
	
	// R240
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iCH3, 1);
	re.back().add_product(iCF3OH, 1);
	re.back().sorder(2.6e-12, 1420, 200, 2.2e-14, 1.1);
	
	// R241
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O, 1);
	re.back().add_educt(iC2H6, 1);
	re.back().add_product(iC2H5, 1);
	re.back().add_product(iCF3OH, 1);
	re.back().sorder(4.9e-12, 400, 100, 1.3e-12, 1.2);
	
	// R242
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O2, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iCF3O, 1);
	re.back().add_product(iO2, 2);
	re.back().sorder(0, 0, 0, 3.0e-15, 0);
	
	// R243
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O2, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iCF3O, 1);
	re.back().add_product(iCO2, 1);
	re.back().sorder(0, 0, 0, 5.0e-16, 0);
	
	// R244
	re.push_back(reaction()); 
	re.back().add_educt(iCF3O2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iCF3O, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(5.4e-12, -320, 150, 1.6e-11, 1.1);

	// R245
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(3.0e-11, -70, 70, 3.8e-11, 1.2);

	// R246
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iOClO, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.4e-12, 960, 300, 1.0e-13, 2.0);
	
	// R247
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCl2O, 1);
	re.back().add_product(iClO, 2);
	re.back().sorder(2.7e-11, 530, 150, 4.5e-12, 1.3);

	// R248
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iCl, 1);
	re.back().sorder(1.0e-11, 3300, 350, 1.5e-16, 2.0);

	// R249
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHOCl, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iClO, 1);
	re.back().sorder(1.7e-13, 0, 300, 1.7e-13, 3.0);

	// R250
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iClONO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.9e-12, 800, 200, 2.0e-13, 1.5);
	
	// R251
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iOClO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.1e-12, 4700, 1000, 3.0e-19, 2.5);

	// R252
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iCl2O2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-19, 0);

	// R253
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iHOCl, 1);
	re.back().add_product(iCl, 1);
	re.back().sorder(1.4e-12, 900, 400, 6.7e-14, 1.2);

	// R254
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.1e-11, -120, 150, 1.7e-11, 1.5);
	
	// R255
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHOClO, 1);
	re.back().add_product(iHOCl, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(4.5e-13, -800, 200, 6.8e-12, 2.0);
	
	// R256
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iCl, 1);
	re.back().sorder(2.6e-12, 350, 100, 8.0e-13, 1.2);
	
	// R257
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHOCl, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iClO, 1);
	re.back().sorder(3.0e-12, 500, 500, 5.0e-13, 3.0);

	// R258
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iClNO2, 1);
	re.back().add_product(iHOCl, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(2.4e-12, 1250, 300, 3.6e-14, 2.0);

	// R259
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iClONO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.2e-12, 330, 200, 3.9e-13, 1.5);
	
	// R260
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3Cl, 1);
	re.back().add_product(iCH2Cl, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(4.0e-12, 1400, 250, 3.6e-14, 1.2);
	
	// R261
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2Cl2, 1);
	re.back().add_product(iCHCl2, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(2.1e-12, 900, 100, 1.0e-13, 1.1);
	
	// R262
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHCl3, 1);
	re.back().add_product(iCCl3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.1e-12, 715, 150, 1.0e-13, 1.2);

	// R263
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCCl4, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-12, 2300, 0, 5.0e-16, 0);
	
	// R264
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCFCl3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-12, 3700, 0, 5.0e-18, 0);
	
	// R265
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCF2Cl2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-12, 3600, 0, 6.0e-18, 0);
	
	// R266
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2ClF, 1);
	re.back().add_product(iCHClF, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(2.8e-12, 1270, 200, 3.9e-14, 1.2);
	
	// R267
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHFCl2, 1);
	re.back().add_product(iCFCl2, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.7e-12, 1250, 150, 2.6e-14, 1.2);
	
	// R268
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHF2Cl, 1);
	re.back().add_product(iCF2Cl, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.0e-12, 1600, 150, 4.7e-15, 1.1);
	
	// R269
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3OCl, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.4e-12, 360, 200, 7.2e-13, 3.0);
	
	// R270
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3CCl3, 1);
	re.back().add_product(iCH2CCl3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.8e-12, 1550, 150, 1.0e-14, 1.1);
	
	// R271
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iC2HCl3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.9e-13, -450, 200, 2.2e-12, 1.25);
	
	// R272
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iC2Cl4, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(9.4e-12, 1200, 200, 1.7e-13, 1.25);
	
	// R273
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCCl3CHO, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iCCl3CO, 1);
	re.back().sorder(8.2e-12, 600, 300, 1.1e-12, 1.5);
		
	// R274
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3CFCl2, 1);
	re.back().add_product(iCH2CFCl2, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.7e-12, 1700, 150, 5.7e-15, 1.2);

	// R275
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3CF2Cl, 1);
	re.back().add_product(iCH2CF2Cl, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.3e-12, 1800, 150, 3.1e-15, 1.2);
	
	// R276
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2ClCF2Cl, 1);
	re.back().add_product(iCHClCF2Cl, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(3.6e-12, 1600, 400, 1.7e-14, 2.0);
		
	// R277
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHCl2CF2Cl, 1);
	re.back().add_product(iCCl2CF2Cl, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.0e-12, 900, 150, 4.9e-14, 1.2);
		
	// R278
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHFClCFCl2, 1);
	re.back().add_product(iCFClCFCl2, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.0e-12, 1250, 150, 1.5e-14, 1.1);
		
	// R279
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2ClCF3, 1);
	re.back().add_product(iCHClCF3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(2.3e-12, 1640, 200, 9.4e-15, 1.2);
		
	// R280
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHCl2CF3, 1);
	re.back().add_product(iCCl2CF3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(7.0e-13, 900, 150, 3.4e-14, 1.2);
		
	// R281
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHFClCF2Cl, 1);
	re.back().add_product(iCFClCF2Cl, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(9.2e-13, 1280, 150, 1.3e-14, 1.2);
		
	// R282
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHFClCF3, 1);
	re.back().add_product(iCFClCF3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(8.0e-13, 1350, 150, 8.6e-15, 1.2);
		
	// R283
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.8e-11, -170, 200, 3.2e-11, 1.5);
		
	// R284
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iClO, 1);
	re.back().sorder(4.1e-11, 450, 200, 9.1e-12, 2.0);
		
	// R285
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iHOCl, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(4.8e-13, -700, 250, 5.0e-12, 1.4);
			
	// R286
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(iClONO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-21, 0);
		
	// R287
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iOClO, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iClO, 1);
	re.back().sorder(2.5e-12, 600, 300, 3.4e-13, 2.0);
			
	// R288
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iCl2O2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-14, 0);
			
	// R289
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iHNO3, 1);
	re.back().add_product(iCl, 1);
	re.back().sorder(0, 0, 0, 5.0e-17, 0);
		
	// R290
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-21, 0);
	
	// R291
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.9e-11, 260, 100, 1.2e-11, 1.15);
	
	// R292
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iH, 1);
	re.back().sorder(3.7e-11, 2300, 200, 1.6e-14, 1.25);
		
	// R293
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iH2O2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iHO2, 1);
	re.back().sorder(1.1e-11, 980, 500, 4.1e-13, 1.5);
		
	// R294
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(2.4e-11, 0, 400, 2.4e-11, 1.5);
		
	// R295
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iHNO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-16, 0);
		
	// R296
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH3, 1);
	re.back().sorder(1.1e-11, 980, 500, 4.1e-13, 1.1);
	
	// R297
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3D, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 7.4e-14, 2.0);
		
	// R298
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iH2CO, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iHCO, 1);
	re.back().sorder(8.1e-11, 30, 100, 7.3e-11, 1.15);
	
	// R299
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3O2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.6e-10, 1.5);
		
	// R300
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3OH, 1);
	re.back().add_product(iCH2OH, 1);
	re.back().add_product(iHCl, 1);
	re.back().sorder(5.4e-11, 0, 250, 5.4e-11, 1.5);
	
	// R301
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iC2H6, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iC2H5, 1);
	re.back().sorder(7.7e-11, 90, 90, 5.7e-11, 1.1);
	
	// R302
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iC2H5O2, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iC2H5O, 1);
	re.back().sorder(0, 0, 0, 7.4e-11, 2.0);
	
	// R303
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iC2H5O2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iC2H4O2, 1);
	re.back().sorder(0, 0, 0, 7.7e-11, 2.0);
	
	// R304
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CN, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.6e-11, 2140, 300, 1.2e-14, 2.0);

	// R305
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CO3NO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1e-14, 0);

	// R306
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iC3H8, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iC3H7, 1);
	re.back().sorder(1.2e-11, -40, 250, 1.4e-10, 1.3);
	
	// R307
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iOClO, 1);
	re.back().add_product(iClO, 2);
	re.back().sorder(03.4e-11, -160, 200, 5.8e-11, 1.25);
		
	// R308
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iClOO, 1);
	re.back().add_product(iCl2, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.3e-10, 0, 250, 2.3e-10, 3.0);
	
	// R309
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iClOO, 1);
	re.back().add_product(iClO, 2);
	re.back().sorder(01.2e-11, 0, 250, 1.2e-11, 3.0);

	// R310
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCl2O, 1);
	re.back().add_product(iCl2, 1);
	re.back().add_product(iClO, 1);
	re.back().sorder(6.2e-11, -130, 130, 9.6e-11, 1.2);
		
	// R311
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCl2O2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-10, 2.0);

	// R312
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iHOCl, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.5e-12, 130, 250, 1.6e-12, 1.5);
	
	// R313
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iClNO, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iCl2, 1);
	re.back().sorder(5.8e-12, -100, 200, 8.1e-11, 1.5);
	
	// R314
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iClONO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(6.5e-12, -135, 50, 1.0e-11, 1.2);
	
	// R315
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3Cl, 1);
	re.back().add_product(iCH2Cl, 1);
	re.back().add_product(iHCl, 1);
	re.back().sorder(3.2e-11, 1250, 200, 4.8e-13, 1.2);
	
	// R316
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2Cl2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCHCl2, 1);
	re.back().sorder(3.1e-11, 1350, 500, 3.3e-13, 1.5);
	
	// R317
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCHCl3, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCCl3, 1);
	re.back().sorder(8.2e-12, 1325, 300, 9.6e-14, 1.3);
		
	// R318
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3F, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH2F, 1);
	re.back().sorder(2.0e-11, 1200, 500, 3.5e-13, 1.3);

	// R319
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2F2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCHF2, 1);
	re.back().sorder(1.2e-11, 1630, 500, 5.0e-14, 1.5);

	// R320
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCF3H, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCF3, 1);
	re.back().sorder(0, 0, 0, 3.0e-18, 5.0);
	
	// R321
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2FCl, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCHFCl, 1);
	re.back().sorder(1.2e-11, 1390, 500, 1.1e-13, 2.0);

	// R322
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCHFCl2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCFCl2, 1);
	re.back().sorder(5.5e-12, 1675, 200, 2.0e-14, 1.3);

	// R323
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCHF2Cl, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCF2Cl, 1);
	re.back().sorder(5.9e-12, 2430, 200, 1.7e-15, 1.3);
	
	// R324
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CCl3, 1);
	re.back().add_product(iCH2CCl3, 1);
	re.back().add_product(iHCl, 1);
	re.back().sorder(2.8e-12, 1790, 400, 7.0e-15, 2.0);
	
	// R325
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CH2F, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH3CHF, 1);
	re.back().sorder(1.8e-11, 290, 500, 6.8e-12, 3.0);
	
	// R326
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CH2F, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH2CH2F, 1);
	re.back().sorder(1.4e-11, 880, 500, 7.3e-13, 3.0);

	// R327
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CH2F, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH2CH2F, 1);
	re.back().sorder(1.4e-11, 880, 500, 7.3e-13, 3.0);

	// R328
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CHF2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH3CF2, 1);
	re.back().sorder(6.4e-12, 950, 500, 2.6e-13, 1.3);
	
	// R329
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CHF2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH2CHF2, 1);
	re.back().sorder(7.2e-12, 2390, 500, 2.4e-15, 3.0);
	
	// R330
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2FCH2F, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCHFCH2F, 1);
	re.back().sorder(2.6e-11, 1060, 500, 7.5e-13, 3.0);
	
	// R331
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CFCl2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH2CFCl2, 1);
	re.back().sorder(1.8e-12, 2000, 300, 2.2e-15, 1.2);

	// R332
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CF2Cl, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH2CF2Cl, 1);
	re.back().sorder(1.4e-12, 2420, 500, 4.2e-16, 1.2);

	// R333
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3CF3, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH2CF3, 1);
	re.back().sorder(1.2e-11, 3880, 500, 2.6e-17, 5.0);
		
	// R334
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2FCHF2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH2FCF2, 1);
	re.back().sorder(5.5e-12, 1610, 500, 2.5e-14, 3.0);
	
	// R335
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2FCHF2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCHFCHF2, 1);
	re.back().sorder(7.7e-12, 1720, 500, 2.4e-14, 3.0);

	// R336
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2ClCF3, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCHClCF3, 1);
	re.back().sorder(1.8e-12, 1710, 500, 5.9e-15, 3.0);
	
	// R337
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2FCF3, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCHFCF3, 1);
	re.back().sorder(0, 0, 0, 1.5e-15, 1.2);
	
	// R338
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCHF2CHF2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCF2CHF2, 1);
	re.back().sorder(78.5e-12, 2430, 500, 2.2e-15, 1.5);
	
	// R339
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCHCl2CF3, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCCl2CF3, 1);
	re.back().sorder(4.4e-12, 1750, 500, 1.2e-14, 1.3);

	// R340
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCHFClCF3, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCFClCF3, 1);
	re.back().sorder(1.1e-12, 1800, 500, 2.7e-15, 1.3);

	// R341
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCHF2CF3, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCF2CF3, 1);
	re.back().sorder(0, 0, 0, 2.4e-16, 1.3);

	// R342
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iClOO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(0, 0, 0, 1.4e-17, 0);

	// R343
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iOClO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.0e-12, 4000, 0, 1.0e-18, 0);
	
	// R344
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-12, 4800, 0, 1.0e-19, 0);

	// R345
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iCl, 1);
	re.back().sorder(6.4e-12, -290, 100, 1.7e-11, 1.15);
	
	// R346
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iClOO, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(4.7e-13, 0, 400, 4.7e-13, 1.5);
	
	// R347
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iN2O, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-12, 4300, 0, 6.0e-19, 0);

	// R348
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-12, 3700, 0, 4.0e-18, 0);
	
	// R349
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-12, 3700, 0, 4.0e-18, 0);
	
	// R350
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iH2CO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-12, 2100, 0, 1.0e-15, 0);
	
	// R351
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iCH3O2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(3.3e-12, 115, 115, 2.2e-12, 1.5);

	// R352
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 2);
	re.back().add_product(iCl2, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.0e-12, 1590, 300, 4.8e-15, 1.5);
	
	// R353
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 2);
	re.back().add_product(iClOO, 1);
	re.back().add_product(iCl, 1);
	re.back().sorder(3.0e-11, 2450, 500, 8.0e-15, 1.5);
	
	// R354
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 2);
	re.back().add_product(iOClO, 1);
	re.back().add_product(iCl, 1);
	re.back().sorder(3.5e-13, 1370, 300, 3.5e-15, 1.5);
	
	// R355
	re.push_back(reaction()); 
	re.back().add_educt(iHCl, 1);
	re.back().add_educt(iClONO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-20, 0);

	// R356
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iBrO, 1);
	re.back().add_product(iBr, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.9e-11, -230, 150, 4.1e-11, 1.5);
	
	// R357
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHBr, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iBr, 1);
	re.back().sorder(5.8e-12, 1500, 200, 3.8e-14, 1.3);
	
	// R358
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHOBr, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iBrO, 1);
	re.back().sorder(1.2e-10, 430, 300, 2.8e-11, 3.0);

	// R359
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iBr2, 1);
	re.back().add_product(iHOBr, 1);
	re.back().add_product(iBr, 1);
	re.back().sorder(4.2e-11, 0, 600, 4.2e-11, 1.3);
	
	// R360
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iBrO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 7.5e-11,  3.0);
	
	// R361
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHBr, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iBr, 1);
	re.back().sorder(11e-11, 0, 250, 1.1e-11, 1.2);
	
	// R362
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3Br, 1);
	re.back().add_product(iCH2Br, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(4.0e-12, 1470, 150, 2.9e-14, 1.1);
	
	// R363
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2Br2, 1);
	re.back().add_product(iCHBr2, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(2.4e-12, 900, 300, 1.2e-13, 1.1);
	
	// R364
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHBr3, 1);
	re.back().add_product(iCBr3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.6e-12, 710, 200, 1.5e-13, 2.0);

	// R365
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHF2Br, 1);
	re.back().add_product(iCF2Br, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.1e-12, 1400, 200, 1.0e-14, 1.1);
	
	// R366
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2ClBr, 1);
	re.back().add_product(iCHClBr, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(2.3e-12, 930, 150, 1.0e-13, 1.2);
	
	// R367
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCF2ClBr, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.5e-16, 0);

	// R368
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCF2Br2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 5.0e-16, 0);
	
	// R369
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCF3Br, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.2e-16, 0);
	
	// R370
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH2BrCF3, 1);
	re.back().add_product(iCHBrCF3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.4e-12, 1340, 200, 1.6e-14, 1.3);
	
	// R371
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHFBrCF3, 1);
	re.back().add_product(iCFBrCF3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(7.2e-13, 1110, 150, 1.8e-14, 1.5);
	
	// R372
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCHClBrCF3, 1);
	re.back().add_product(iCClBrCF3, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(1.3e-12, 995, 150, 4.5e-14, 1.5);

	// R373
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCF2BrCHFCl, 1);
	re.back().add_product(iCF2BrCFCl, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(9.3e-13, 1250, 150, 1.4e-14, 1.5);

	// R374
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCF2BrCF2Br, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.5e-16, 0);
	
	// R375
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iBr, 1);
	re.back().add_product(iHBr, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.5e-11, 600, 600, 2.0e-12, 2.0);

	// R376
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iBrO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(3.4e-12, -540, 200, 2.1e-11, 1.5);

	// R377
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iHBr, 1);
	re.back().add_product(iHNO3, 1);
	re.back().add_product(iBr, 1);
	re.back().sorder(0, 0, 0, 1.0e-16, 0);
	
	// R378
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2ClBr, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCHClBr, 1);
	re.back().sorder(4.3e-11, 1370, 500, 4.3e-13, 3.0);
	
	// R379
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3Br, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCH2Br, 1);
	re.back().sorder(1.5e-11, 1060, 100, 4.3e-13, 1.2);
	
	// R380
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH2Br2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCHBr2, 1);
	re.back().sorder(6.4e-12, 810, 100, 4.2e-13, 1.2);
	
	// R381
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iBrO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.7e-11, 800, 200, 1.2e-12, 1.2);

	// R382
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iH2O2, 1);
	re.back().add_product(iHBr, 1);
	re.back().add_product(iHO2, 1);
	re.back().sorder(1.0e-11, 3000, 0, 5.0e-16, 0);
	
	// R383
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iBrO, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(0, 0, 0, 1.6e-11, 2.0);
	
	// R384
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iH2CO, 1);
	re.back().add_product(iHBr, 1);
	re.back().add_product(iHCO, 1);
	re.back().sorder(1.7e-11, 800, 200, 1.1e-12, 1.3);
	
	// R385
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iOClO, 1);
	re.back().add_product(iBrO, 1);
	re.back().add_product(iClO, 1);
	re.back().sorder(2.6e-11, 1300, 300, 3.4e-13, 2.0);
	
	// R386
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iCl2O, 1);
	re.back().add_product(iBrCl, 1);
	re.back().add_product(iClO, 1);
	re.back().sorder(2.1e-11, 470, 150, 4.3e-12, 1.3);
	
	// R387
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iCl2O2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.0e-12, 2.0);

	// R388
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.0e-12, 3200, 0, 2.0e-17, 0);
	
	// R389
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iBr, 1);
	re.back().sorder(8.8e-12, -260, 130, 2.1e-11, 1.15);
	
	// R390
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-12, 3.0);
	
	// R391
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iBr, 1);
	re.back().add_product(iOClO, 1);
	re.back().sorder(1.6e-12, -430, 200, 6.8e-12, 1.25);
	
	// R393
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iBr, 1);
	re.back().add_product(iClOO, 1);
	re.back().sorder(2.9e-12, -220, 200, 6.1e-12, 1.25);
	
	// R394
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iBrCl, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(5.8e-13, -170, 200, 1.0e-12, 1.25);
	
	// R395
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 2);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.5e-12, -230, 150, 3.2e-12, 1.15);
	
	// R396
	re.push_back(reaction()); 
	re.back().add_educt(iCH2BrO2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iCH2O, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iBr, 1);
	re.back().sorder(4.0e-12, -300, 200, 1.1e-11, 1.5);
	
	// R397
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iI2, 1);
	re.back().add_product(iIO, 1);
	re.back().add_product(iI, 1);
	re.back().sorder(1.4e-10, 0, 250, 1.4e-10, 1.4);
	
	// R398
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iIO, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iI, 1);
	re.back().sorder(0, 0, 0, 1.2e-10, 2.0);
	
	// R399
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iI2, 1);
	re.back().add_product(iHOI, 1);
	re.back().add_product(iI, 1);
	re.back().sorder(0, 0, 0, 1.8e-10, 2.0);
	
	// R400
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHI, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iI, 1);
	re.back().sorder(0, 0, 0, 3.0e-11, 2.0);
	
	// R401
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3I, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iCH2I, 1);
	re.back().sorder(3.1e-12, 1120, 500, 7.2e-14, 3.0);
	
	// R402
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCF3I, 1);
	re.back().add_product(iHOI, 1);
	re.back().add_product(iCF3, 1);
	re.back().sorder(0, 0, 0, 3.1e-14, 5.0);
	
	// R403
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iI, 1);
	re.back().add_product(iHI, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.5e-11, 1090, 500, 3.8e-13, 2.0);
	
	// R404
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iIO, 1);
	re.back().add_product(iHOI, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(0, 0, 0, 8.4e-11, 1.5);
	
	// R405
	re.push_back(reaction()); 
	re.back().add_educt(iI, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iIO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.3e-11, 870, 200, 1.2e-12, 1.2);
	
	// R406
	re.push_back(reaction()); 
	re.back().add_educt(iI, 1);
	re.back().add_educt(iBrO, 1);
	re.back().add_product(iIO, 1);
	re.back().add_product(iBr, 1);
	re.back().sorder(0, 0, 0, 1.2e-11, 2.0);
	
	// R407
	re.push_back(reaction()); 
	re.back().add_educt(iIO, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iI, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(9.1e-12, -240, 150, 2.0e-11, 1.2);
	
	// R408
	re.push_back(reaction()); 
	re.back().add_educt(iIO, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(5.1e-12, -280, 200, 1.3e-11, 2.0);
	
	// R409
	re.push_back(reaction()); 
	re.back().add_educt(iIO, 1);
	re.back().add_educt(iBrO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 6.9e-11, 1.5);
	
	// R409
	re.push_back(reaction()); 
	re.back().add_educt(iIO, 2);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.5e-11, -500, 500, 8.0e-11, 1.5);
	
	// R410
	re.push_back(reaction()); 
	re.back().add_educt(iINO, 2);
	re.back().add_product(iI2, 1);
	re.back().add_product(iNO, 2);
	re.back().sorder(8.4e-11, 2620, 600, 1.3e-14, 2.5);
	
	// R412
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iSH, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iH, 1);
	re.back().sorder(0, 0, 0, 1.6e-10, 5.0);
	
	// R413
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCS, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iS, 1);
	re.back().sorder(2.7e-10, 760, 250, 2.1e-11, 1.1);

	// R414
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iSH, 1);
	re.back().sorder(9.2e-12, 1800, 550, 2.2e-14, 1.7);

	// R415
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iSO, 1);
	re.back().sorder(2.1e-11, 2200, 150, 1.3e-14, 1.2);

	// R416
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCS2, 1);
	re.back().add_product(iCS, 1);
	re.back().add_product(iSO, 1);
	re.back().sorder(3.2e-11, 650, 150, 3.6e-12, 1.2);

	// R417
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iCH3SO, 1);
	re.back().add_product(iCH3, 1);
	re.back().sorder(1.3e-11, -410, 100, 5.0e-11, 1.1);
	
	// R418
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCH3SSCH3, 1);
	re.back().add_product(iCH3SO, 1);
	re.back().add_product(iCH3S, 1);
	re.back().sorder(5.5e-11, -250, 100, 1.3e-10, 1.3);

	// R419
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-20, 0);
	
	// R420
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-18, 0);
	
	// R421
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iSO3, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(3.0e-12, 7000, 0, 2.0e-22, 0);
	
	// R422
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iSH, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(6.0e-12, 75, 75, 4.7e-12, 1.2);
	
	// R423
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.1e-13, 1200, 500, 1.9e-15, 2.0);
	
	// R424
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3SH, 1);
	re.back().add_product(iCH3S, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(9.9e-12, -360, 100, 3.3e-11, 1.2);
	
	// R425
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iCH2SCH3, 1);
	re.back().sorder(1.2e-11, 260, 100, 5.0e-12, 1.15);
	
	// R426
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCH3SSCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(6.0e-11, -400, 200, 2.3e-10, 1.2);
  
	// R427
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iS, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iSO, 1);
	re.back().sorder(0, 0, 0, 6.6e-11, 3.0);
	
	// R428
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iSO, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iSO2, 1);
	re.back().sorder(0, 0, 0, 8.6e-11, 2.0);
	
	// R429
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.0e-15, 0);
	
	// R430
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iCH3SH, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 4.0e-15, 0);
		
	// R431
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 5.0e-15, 0);
		
	// R432
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-18, 0);
		
	// R433
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-26, 0);
		
	// R434
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 8.0e-16, 0);
		
	// R435
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-16, 0);
	
	// R436
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iCS2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 4.0e-16, 0);
	
	// R437
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iCH3SH, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(4.4e-13, -210, 210, 8.9e-13, 1.25);

	// R438
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iCH3SCH2, 1);
	re.back().add_product(iHNO3, 1);
	re.back().sorder(1.9e-13, -500, 200, 1.0e-12, 1.2);
	
	// R439
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iCH3SSCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.3e-12, 270, 270, 5.3e-15, 1.4);
	
	// R440
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 7.0e-21, 0);
		
	// R441
	re.push_back(reaction()); 
	re.back().add_educt(iN2O5, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-17, 0);
	
	// R442
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-21, 0);
		
	// R443
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.4e-10, 2.0);
	
	// R444
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iSH, 1);
	re.back().sorder(3.7e-11, -210, 100, 7.4e-11, 1.25);
	
	// R445
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-16, 0);
	
	// R446
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCS2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 4.0e-15, 0);
	
	// R447
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCH3SH, 1);
	re.back().add_product(iCH3S, 1);
	re.back().add_product(iHCl, 1);
	re.back().sorder(1.2e-10, -150, 50, 2.0e-10, 1.25);
	
	// R448
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-16, 0);
	
	// R449
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 9.5e-15, 2.0);
	
	// R450
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iSO, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iSO2, 1);
	re.back().sorder(2.8e-11, 0, 50, 2.8e-11, 1.3);
	
	// R451
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iSO3, 1);
	re.back().sorder(0, 0, 0, 4.0e-18, 0);
	
	// R452
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iHBr, 1);
	re.back().add_product(iSH, 1);
	re.back().sorder(1.4e-11, 2750, 300, 1.4e-15, 2.0);
	
	// R453
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iCH3SH, 1);
	re.back().add_product(iCH3S, 1);
	re.back().add_product(iHBr, 1);
	re.back().sorder(9.2e-12, 390, 100, 2.5e-12, 2.0);
	
	// R454
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.5e-14, -850, 200, 2.6e-13, 1.3);
	
	// R455
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(iSO, 1);
	re.back().add_product(iBr, 1);
	re.back().add_product(iSO2, 1);
	re.back().sorder(0, 0, 0, 5.7e-11, 1.4);
	
	// R456
	re.push_back(reaction()); 
	re.back().add_educt(iIO, 1);
	re.back().add_educt(iCH3SH, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 6.6e-16, 2.0);
		
	// R457
	re.push_back(reaction()); 
	re.back().add_educt(iIO, 1);
	re.back().add_educt(iCH3SCH3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.2e-14, 1.5);
	
	// R458
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iO, 1);
	re.back().sorder(2.3e-12, 0, 200, 2.3e-12, 1.2);
	
	// R459
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(0, 0, 0, 1.2e-11, 2.0);
	
	// R460
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iO, 1);
	re.back().sorder(2.6e-13, 2400, 500, 8.4e-17, 2.0);

	// R461
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(3.6e-12, 1100, 200, 9.0e-14, 1.2);
	
	// R462
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iNO, 1);
	re.back().sorder(1.4e-11, 0, 50, 1.4e-11, 1.2);
	
	// R463
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iOClO, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iClO, 1);
	re.back().sorder(0, 0, 0, 1.9e-12, 3.0);
	
	// R464
	re.push_back(reaction()); 
	re.back().add_educt(iSO3, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-19, 10.0);
	
	// R465
	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iSO, 1);
	re.back().sorder(0, 0, 0, 4.0e-19, 0);
	
	// R466
	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iHSO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(9.0e-12, 280, 200, 3.5e-12, 1.3);
	
	// R467
	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iH2O2, 1);
	re.back().add_product(iProducts, 1);;
	re.back().sorder(0, 0, 0, 5.0e-15, 0);
	
	// R468
	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iHSO, 1);
	re.back().add_product(iNO, 1);
	re.back().sorder(2.9e-11, -240, 50, 6.5e-11, 1.2);
	
	// R469
	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iClSH, 1);
	re.back().add_product(iCl, 1);
	re.back().sorder(1.7e-11, 690, 200, 1.7e-12, 2.0);
	
	// R470
	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iBrCl, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.3e-11, -350, 200, 7.4e-11, 2.0);
	
	// R471
	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iBr2, 1);
	re.back().add_product(iBrSH, 1);
	re.back().add_product(iBr, 1);
	re.back().sorder(6.0e-11, -160, 160, 1.0e-10, 2.0);
	
	// R472
	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iF2, 1);
	re.back().add_product(iFSH, 1);
	re.back().add_product(iF, 1);
	re.back().sorder(4.3e-11, 1390, 200, 4.0e-13, 2.0);
	
	// R473
	re.push_back(reaction()); 
	re.back().add_educt(iHSO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 2.0e-17, 0);
	
	// R474
	re.push_back(reaction()); 
	re.back().add_educt(iHSO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-13, 1.3);
	
	// R475
	re.push_back(reaction()); 
	re.back().add_educt(iHSO, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-15, 0);
	
	// R476
	re.push_back(reaction()); 
	re.back().add_educt(iHSO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iHSO2, 1);
	re.back().add_product(iNO, 1);
	re.back().sorder(0, 0, 0, 9.6e-12, 2.0);
	
	// R477
	re.push_back(reaction()); 
	re.back().add_educt(iHSO2, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iSO2, 1);
	re.back().sorder(0, 0, 0, 3.0e-13, 3.0);
	
	// R478
	re.push_back(reaction()); 
	re.back().add_educt(iHOSO2, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iSO3, 1);
	re.back().sorder(1.3e-12, 330, 200, 4.4e-13, 1.2);
	
	// R479
	re.push_back(reaction()); 
	re.back().add_educt(iCS, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iOCS, 1);
	re.back().add_product(iO, 1);
	re.back().sorder(0, 0, 0, 2.9e-19, 2.0);
	
	// R480
	re.push_back(reaction()); 
	re.back().add_educt(iCS, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iOCS, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(0, 0, 0, 3.0e-16, 3.0);
	
	// R481
	re.push_back(reaction()); 
	re.back().add_educt(iCS, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iOCS, 1);
	re.back().add_product(iNO, 1);
	re.back().sorder(0, 0, 0, 7.6e-17, 3.0);
	
	// R482
	re.push_back(reaction()); 
	re.back().add_educt(iCH3S, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.0e-18, 0);
	
	// R483
	re.push_back(reaction()); 
	re.back().add_educt(iCH3S, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.0e-12, -290, 100, 5.3e-12, 1.15);
	
	// R484
	re.push_back(reaction()); 
	re.back().add_educt(iCH3S, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.0e-13, 0);
	
	// R485
	re.push_back(reaction()); 
	re.back().add_educt(iCH3S, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iCH3SO, 1);
	re.back().add_product(iNO, 1);
	re.back().sorder(2.1e-11, 320, 100, 6.1e-11, 1.15);
	
	// R486
	re.push_back(reaction()); 
	re.back().add_educt(iCH2SH, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 6.5e-12, 2.0);
	
	// R487
	re.push_back(reaction()); 
	re.back().add_educt(iCH2SH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.5e-11, 2.0);
	
	// R488
	re.push_back(reaction()); 
	re.back().add_educt(iCH2SH, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.9e-11, 2.0);
	
	// R489
	re.push_back(reaction()); 
	re.back().add_educt(iCH2SH, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 5.2e-11, 2.0);
	
	// R490
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 6.0e-13, 1.5);
	
	// R491
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iCH3SO2, 1);
	re.back().add_product(iNO, 1);
	re.back().sorder(0, 0, 0, 1.2e-11, 1.4);
	
	// R492
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SOO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 8.0e-13, 0);
	
	// R493
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SOO, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(1.1e-11, 0, 100, 1.1e-11, 2.0);
	
	// R494
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SO2, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.2e-11, 0, 100, 2.2e-11, 2.0);
	
	// R495
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SCH2, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 3.0e-10, 2.0);
	
	// R496
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SCH2O2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iCH3SCH2O, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(0, 0, 0, 1.9e-11, 2.0);
	
	// R497
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SS, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 4.6e-13, 2.0);
	
	// R498
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SS, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 1.8e-11, 2.0);
	
	// R499
	re.push_back(reaction()); 
	re.back().add_educt(iCH3SSO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(0, 0, 0, 4.5e-12, 2.0);
	
	// R500
	re.push_back(reaction()); 
	re.back().add_educt(iNa, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNaO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.0e-9, 95, 50, 7.3e-10, 1.2);
	
	// R501
	re.push_back(reaction()); 
	re.back().add_educt(iNa, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNaO2, 1);
	re.back().add_product(iO, 1);
	re.back().sorder(0, 0, 0, 4.0e-11, 0);
	
	// R502
	re.push_back(reaction()); 
	re.back().add_educt(iNa, 1);
	re.back().add_educt(iN2O, 1);
	re.back().add_product(iNaO, 1);
	re.back().add_product(iN2, 1);
	re.back().sorder(02.8e-10, 1600, 400, 1.2e-12, 1.2);
	
	// R503
	re.push_back(reaction()); 
	re.back().add_educt(iNa, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iNaCl, 1);
	re.back().add_product(iCl, 1);
	re.back().sorder(7.3e-10, 0, 200, 7.3e-10, 1.3);
	
	// R504
	re.push_back(reaction()); 
	re.back().add_educt(iNaO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iNa, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(3.7e-10, 0, 400, 3.7e-10, 3.0);
	
	// R505
	re.push_back(reaction()); 
	re.back().add_educt(iNaO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNaO2, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(1.1e-9, 570, 300, 1.6e-10, 1.5);
	
	// R506
	re.push_back(reaction()); 
	re.back().add_educt(iNaO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNa, 1);
	re.back().add_product(iO2, 2);
	re.back().sorder(6.0e-11, 0, 800, 6.0e-11, 3.0);
	
	// R507
	re.push_back(reaction()); 
	re.back().add_educt(iNaO, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iNaOH, 1);
	re.back().add_product(iH, 1);
	re.back().sorder(2.6e-11, 0, 600, 2.6e-11, 2.0);
	
	// R508
	re.push_back(reaction()); 
	re.back().add_educt(iNaO, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iNaOH, 1);
	re.back().add_product(iOH, 1);
	re.back().sorder(2.2e-10, 0, 400, 2.2e-10, 2.0);
	
	// R509
	re.push_back(reaction()); 
	re.back().add_educt(iNaO, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iNa, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(1.5e-10, 0, 400, 1.5e-10, 4.0);
	
	// R510
	re.push_back(reaction()); 
	re.back().add_educt(iNaO, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.8e-10, 0, 400, 2.8e-10, 3.0);
	
	// R511
	re.push_back(reaction()); 
	re.back().add_educt(iNaO2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iNaO, 1);
	re.back().add_product(iO2, 1);
	re.back().sorder(2.2e-11, 0, 600, 2.2e-11, 5.0);
	
	// R512
	re.push_back(reaction()); 
	re.back().add_educt(iNaO2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iNaO, 1);
	re.back().add_product(iNO2, 1);
	re.back().sorder(0, 0, 0, 1.0e-14, 0);
		
	// R513
	re.push_back(reaction()); 
	re.back().add_educt(iNaO2, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iProducts, 1);
	re.back().sorder(2.3e-10, 0, 400, 2.3e-10, 3.0);
	
	// R514
	re.push_back(reaction()); 
	re.back().add_educt(iNaOH, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iNaCl, 1);
	re.back().add_product(iH2O, 1);
	re.back().sorder(2.8e-10, 0, 400, 2.8e-10, 3.0);
	
	// R515
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iO3, 1);
	re.back().assrea(6.0e-34, 2.3, 0, 0);
	
	// R516
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iHO2, 1);
	re.back().assrea(5.7e-32, 1.6, 7.5e-11, 0.0);	
	
	// R517
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 2);
	re.back().add_product(iH2O2, 1);
	re.back().assrea(6.2e-31, 1.0, 2.6e-11, 0.0);
	
	// R518
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iNO2, 1);
	re.back().assrea(9.0e-32, 1.5, 3.0e-11, 0.0);
	
	// R519
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iNO3, 1);
	re.back().assrea(9.0e-32, 2.0, 2.2e-11, 0.0);

	// R520
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iHONO, 1);
	re.back().assrea(7.0e-31, 2.6, 3.6e-11, 0.1);
	
	// R521
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iHNO3, 1);
	re.back().assrea(2.5e-30, 4.4, 1.6e-11, 1.7);
	
	// R522
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iHO2NO2, 1);
	re.back().assrea(1.8e-31, 3.2, 4.7e-12, 1.4);
  
	// R523
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iN2O5, 1);
	re.back().assrea(2.2e-30, 3.9, 1.5e-12, 0.7);
	
	// R524
	re.push_back(reaction()); 
	re.back().add_educt(iCH3, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCH3O2, 1);
	re.back().assrea(4.5e-31, 3.0, 1.8e-12, 1.7);

	// R525
	re.push_back(reaction()); 
	re.back().add_educt(iC2H5, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iC2H5O2, 1);
	re.back().assrea(1.5e-28, 3.0, 8.0e-12, 0.0);
	
	// R526
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iC2H2, 1);
	re.back().add_product(iHOCHCH, 1);
	re.back().assrea(5.5e-30, 0.0, 8.3e-13, -2.0);
	
	// R527
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iC2H4, 1);
	re.back().add_product(iHOCH2CH2, 1);
	re.back().assrea(1.0e-28, 0.8, 8.8e-12, 0.0);
	
	// R528
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iCH3ONO, 1);
	re.back().assrea(1.4e-29, 3.8, 3.6e-11, 0.6);
		
	// R529
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iCH3ONO2, 1);
	re.back().assrea(1.1e-28, 4.0, 1.6e-11, 1.0);
	
	// R530
	re.push_back(reaction()); 
	re.back().add_educt(iC2H5O, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iC2H5ONO, 1);
	re.back().assrea(2.8e-27, 4.0, 5.0e-11, 1.0);
	
	// R531
	re.push_back(reaction()); 
	re.back().add_educt(iC2H5O, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iC2H5ONO2, 1);
	re.back().assrea(2.0e-27, 4.0, 6.5e-12, 2.0);
		
	// R532
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iCH3O2NO2, 1);
	re.back().assrea(1.5e-30, 4.0, 6.5e-12, 2.0);
		
	// R533
	re.push_back(reaction()); 
	re.back().add_educt(iCH3C_O_O2, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iCH3O2NO2, 1);
	re.back().assrea(9.7e-29, 5.6, 9.3e-12, 1.5);
	
       // R534
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iFO2, 1);
	re.back().assrea(4.4e-33, 1.2, 0, 0);
	
        // R535
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iFNO, 1);
	re.back().assrea(1.8e-31, 1.0, 2.8e-10, 0.0);
	
        // R536
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iFNO2, 1);
	re.back().assrea(6.3e-32, 2.0, 2.6e-10, 0.0);
	
        // R537
	re.push_back(reaction()); 
	re.back().add_educt(iFO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iFONO2, 1);
	re.back().assrea(2.6e-31, 1.3, 2.0e-11, 1.5);
	
	// R538
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iClOO, 1);
	re.back().assrea(2.7e-33, 1.5, 0, 0);
	
        // R539
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iClNO, 1);
	re.back().assrea(9.0e-32, 1.6, 0, 0);
	
        // R540
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iClONO, 1);
	re.back().assrea(1.3e-30, 2.0, 1.0e-10, 1.0);
	
        // R541
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iClNO2, 1);
	re.back().assrea(1.8e-31, 2.0, 1.0e-10, 1.0);
	
        // R542
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iClCO, 1);
	re.back().assrea(1.3e-33, 3.8, 0, 0);

	// R543
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iC2H2, 1);
	re.back().add_product(iClC2H2, 1);
	re.back().assrea(5.9e-30, 2.1, 2.1e-10, 1.0);
	
	// R544
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iC2H4, 1);
	re.back().add_product(iClC2H4, 1);
	re.back().assrea(1.6e-29, 3.3, 3.1e-10, 1.0);

	// R545
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iC2Cl4, 1);
	re.back().add_product(iC2Cl5, 1);
	re.back().assrea(1.4e-28, 8.5, 4.0e-11, 1.2);
	
	// R546
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iClONO2, 1);
	re.back().assrea(1.8e-31, 3.4, 1.5e-11, 1.9);
	
	// R547
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 2);
	re.back().add_product(iCl2O2, 1);
	re.back().assrea(2.2e-32, 3.1, 3.5e-12, 1.0);
	
        // R548
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iOClO, 1);
	re.back().add_product(iCl2O3, 1);
	re.back().assrea(6.2e-32, 4.7, 2.4e-11, 0.0);
	
	// R549
	re.push_back(reaction()); 
	re.back().add_educt(iOClO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iClO3, 1);
	re.back().assrea(1.9e-31, 1.1, 3.1e-11, 0.0);
	
        // R550
	re.push_back(reaction()); 
	re.back().add_educt(iCH2Cl, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCH2ClO2, 1);
	re.back().assrea(1.9e-30, 3.2, 2.9e-12, 1.2);
	
	// R551
	re.push_back(reaction()); 
	re.back().add_educt(iCHCl2, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCHCl2O2, 1);
	re.back().assrea(1.3e-30, 4.0, 2.8e-12, 1.4);
	
        // R552
	re.push_back(reaction()); 
	re.back().add_educt(iCCl3, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCCl3O2, 1);
	re.back().assrea(6.9e-31, 6.4, 2.4e-12, 2.1);
	
	// R553
	re.push_back(reaction()); 
	re.back().add_educt(iCFCl2, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCFCl2O2, 1);
	re.back().assrea(5.0e-30, 4.0, 6.0e-12, 1.0);
	
        // R554
	re.push_back(reaction()); 
	re.back().add_educt(iCF2Cl, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iCF2ClO2, 1);
	re.back().assrea(3.0e-30, 4.0, 3.0e-12, 1.0);
	
	// R555
	re.push_back(reaction()); 
	re.back().add_educt(iBr, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iBrNO2, 1);
	re.back().assrea(4.2e-31, 2.4, 2.7e-11, 0.0);
	
	// R556
	re.push_back(reaction()); 
	re.back().add_educt(iBrO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iBrONO2, 1);
	re.back().assrea(5.2e-31, 3.2, 6.9e-12, 2.9);
	
	// R557
	re.push_back(reaction()); 
	re.back().add_educt(iI, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iINO, 1);
	re.back().assrea(1.8e-32, 1.0, 1.7e-11, 0.0);
	
        // R558
	re.push_back(reaction()); 
	re.back().add_educt(iI, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iINO2, 1);
	re.back().assrea(3.0e-31, 1.0, 6.6e-11, 0.0);
	
        // R559
	re.push_back(reaction()); 
	re.back().add_educt(iIO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iIONO2, 1);
	re.back().assrea(5.9e-31, 3.5, 9.0e-12, 1.5);
	
        // R560
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iSO3, 1);
	re.back().assrea(1.3e-33, -3.6, 0, 0);
	
	// R561
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iHOSO2, 1);
	re.back().assrea(3.0e-31, 3.3, 1.5e-12, 0.0);
	
        // R562
	re.push_back(reaction()); 
	re.back().add_educt(iNa, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iNaO2, 1);
	re.back().assrea(3.2e-30, 1.4, 6.0e-10, 0.0);
	
	// R563
	re.push_back(reaction()); 
	re.back().add_educt(iNaO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iNaO3, 1);
	re.back().assrea(3.5e-30, 2.0, 5.7e-10, 0.0);
	
	// R564
	re.push_back(reaction()); 
	re.back().add_educt(iNaO, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iNaCO3, 1);
	re.back().assrea(8.7e-28, 2.0, 6.5e-10, 0.0);
	
	// R565
	re.push_back(reaction()); 
	re.back().add_educt(iNaOH, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iNaHCO3, 1);
	re.back().assrea(1.3e-28, 2.0, 6.8e-10, 0.0);
	
	
	// R566
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iHO2NO2, 1);
	re.back().eqrea(2.1e-27, 10900, 1.6e-11, 5);
	
	// R567
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iN2O3, 1);
	re.back().eqrea(3.3e-27, 4667, 2.1e-20, 2);
	
	// R568
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 2);
	re.back().add_product(iN2O4, 1);
	re.back().eqrea(5.2e-29, 6643, 2.5e-19, 2);
	
	// R569
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iN2O5, 1);
	re.back().eqrea(2.7e-27, 11000, 2.9e-11, 1.3);
		
	// R570
	re.push_back(reaction()); 
	re.back().add_educt(iCH3C_O_O2, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iCH3C_O_O2NO2_PAN_, 1);
	re.back().eqrea(9.0e-29, 14000, 2.3e-8, 2);     // WARNING
	
	// R571
	re.push_back(reaction()); 
	re.back().add_educt(iF, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iFOO, 1);
	re.back().eqrea(3.2e-25, 6100, 2.5e-16, 1);      // WARNING
	
	// R572
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iClOO, 1);
	re.back().eqrea(5.7e-25, 2500, 2.5e-21, 2);
	
	// R573
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iClCO, 1);
	re.back().eqrea(1.6e-26, 4000, 1.1e-19, 5);
	
	
	// R574
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iClO_O2, 1);
	re.back().eqrea(2.9e-26, 3700, 7.2e-21, 0);
	
	// R575
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 2);
	re.back().add_product(iCl2O2, 1);
	re.back().eqrea(1.3e-27, 8744, 7.2e-15, 1.5);
	
	// R576
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iOClO, 1);
	re.back().add_product(iCl2O3, 1);
	re.back().eqrea(1.1e-24, 5455, 9.8e-17, 3);
	

}



void build_model_Krasnopolsky_2007(std::vector<species> &sp, std::vector<reaction> &re) {
	sp.push_back(species(sp.size(), "SO3"));
        sp.push_back(species(sp.size(), "H2O"));
        sp.push_back(species(sp.size(), "H2SO4"));
        sp.push_back(species(sp.size(), "CO"));
        sp.push_back(species(sp.size(), "CO2"));
        sp.push_back(species(sp.size(), "SO2"));
        sp.push_back(species(sp.size(), "OCS"));
        sp.push_back(species(sp.size(), "S2"));
        sp.push_back(species(sp.size(), "M"));
        sp.push_back(species(sp.size(), "S"));
        sp.push_back(species(sp.size(), "S3"));
        sp.push_back(species(sp.size(), "NO"));
        sp.push_back(species(sp.size(), "SNO"));
        sp.push_back(species(sp.size(), "SH"));
        sp.push_back(species(sp.size(), "H2S"));
        sp.push_back(species(sp.size(), "H2"));
        sp.push_back(species(sp.size(), "H"));
        sp.push_back(species(sp.size(), "HCl"));
        sp.push_back(species(sp.size(), "Cl"));
        sp.push_back(species(sp.size(), "ClSO2"));
        sp.push_back(species(sp.size(), "SCl"));
        sp.push_back(species(sp.size(), "SO"));
        sp.push_back(species(sp.size(), "OH"));
        sp.push_back(species(sp.size(), "rH"));
        sp.push_back(species(sp.size(), "HSCl"));
        sp.push_back(species(sp.size(), "SO2Cl2"));
        sp.push_back(species(sp.size(), "(SO)2"));
        sp.push_back(species(sp.size(), "hv"));
        sp.push_back(species(sp.size(), "Cl2"));
	
	// IDs
	size_t iSO3=0;
	size_t iH2O=1;
	size_t iH2SO4=2;
	size_t iCO=3;
	size_t iCO2=4;
	size_t iSO2=5;
	size_t iOCS=6;
	size_t iS2=7;
	size_t iM=8;
	size_t iS=9;
	size_t iS3=10;
	size_t iNO=11;
	size_t iSNO=12;
	size_t iSH=13;
	size_t iH2S=14;
	size_t iH2=15;
	size_t iH=16;
	size_t iHCl=17;
	size_t iCl=18;
	size_t iClSO2=19;
	size_t iSCl=20;
	size_t iSO=21;
	size_t iOH=22;
	size_t irH=23;
	size_t iHSCl=24;
	size_t iSO2Cl2=25;
	size_t iS2O2=26;
	size_t ihv=27;
	size_t iCl2=28;
	
	
	// R1
	re.push_back(reaction()); 
	re.back().add_educt(iSO3, 1);
	re.back().add_educt(iH2O, 2);
	re.back().add_product(iH2SO4, 1);
	re.back().add_product(iH2O, 1);  
	
        // R2
	re.push_back(reaction()); 
	re.back().add_educt(iH2SO4, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iSO3, 1);
	re.back().add_product(iH2O, 2);  

	// R3
	re.push_back(reaction()); 
	re.back().add_educt(iSO3, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iSO2, 1);
	
	// R4
 	re.push_back(reaction()); 
	re.back().add_educt(iSO3, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iS2O2, 1);  
        
	// R5
 	re.push_back(reaction()); 
	re.back().add_educt(iS2O2, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iS2, 1);  

	// R6
 	re.push_back(reaction()); 
	re.back().add_educt(iSO, 2);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iS, 1);
	
	// R7
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iSO, 2);
	
	// R8
 	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iSO, 1);
	
	// R9
 	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(iSO, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iSO2, 1);

	// R10
 	re.push_back(reaction()); 
	re.back().add_educt(iS3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iS, 1);
	
	// R11
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iM, 1);
	
	// R12
 	re.push_back(reaction()); 
	re.back().add_educt(iS2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iS, 2);
	re.back().add_product(iM, 1);
		
	// R13
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iS2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iS3, 1);
	re.back().add_product(iM, 1);
		
	// R14
 	re.push_back(reaction()); 
	re.back().add_educt(iS3, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iS, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iM, 1);		

	// R15
 	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iS, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iOCS, 1);
	re.back().add_product(iM, 1);
	
	// R16
 	re.push_back(reaction()); 
	re.back().add_educt(iOCS, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iS, 1);
	re.back().add_product(iM, 1);	
	
	// R17
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iS2, 1);	

	// R18
 	re.push_back(reaction()); 
	re.back().add_educt(iS2, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iOCS, 1);
	re.back().add_product(iS, 1);	

	// R19
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iS3, 1);
	re.back().add_product(iS2, 2);
	
	// R20
 	re.push_back(reaction()); 
	re.back().add_educt(iS2, 2);
	re.back().add_product(iS3, 1);
	re.back().add_product(iS, 1);
	
	// R21
 	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iS3, 1);
	re.back().add_product(iOCS, 1);
	re.back().add_product(iS2, 1);
	
	// R22
 	re.push_back(reaction()); 
	re.back().add_educt(iOCS, 1);
	re.back().add_educt(iS2, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iS3, 1);
	
	// R23
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iSNO, 1);
	re.back().add_product(iM, 1);
		
	// R24
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iSNO, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iNO, 1);
		
	// R25
 	re.push_back(reaction()); 
	re.back().add_educt(iS2, 1);
	re.back().add_educt(iSNO, 1);
	re.back().add_product(iS3, 1);
	re.back().add_product(iNO, 1);		
	
	// R26
 	re.push_back(reaction()); 
	re.back().add_educt(iSH, 2);
	re.back().add_product(iH2S, 1);
	re.back().add_product(iS, 1);
	
	// R27
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iSH, 2);
	
	// R28
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iSH, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iH, 1);

    	// R29
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iS2, 1);
	re.back().add_product(iS, 1);
	re.back().add_product(iSH, 1);
	
	// R30
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iSH, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iS, 1);
	
	// R31
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iSH, 1);
	re.back().add_product(iH, 1);
	
	// R32
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iSH, 1);
	
	// R33
 	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iSH, 1);
	re.back().add_product(iOCS, 1);
	re.back().add_product(iH, 1);
	
	// R34
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iCl, 1);
	
	// R35
 	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iH, 1);
	
	// R36
 	re.push_back(reaction()); 
	re.back().add_educt(iH2S, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iSH, 1);
	
	// R37
 	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iH2S, 1);
	re.back().add_product(iCl, 1);
	
	// R38
 	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iSH, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iS, 1);
	
	// R39
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iSH, 1);
	re.back().add_product(iCl, 1);
	
	// R40
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iSNO, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iSH, 1);
	
	// R41
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iH2S, 1);
	re.back().add_product(iM, 1);
	
	// R42
 	re.push_back(reaction()); 
	re.back().add_educt(iH2S, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iSH, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iM, 1);    

   	// R43
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iS3, 1);
	re.back().add_product(iSH, 1);
	re.back().add_product(iS2, 1);

	// R44
 	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iClSO2, 1);
	re.back().add_product(iM, 1);
	
	// R45
 	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iM, 1);
	
	// R46
 	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 2);
	re.back().add_product(iSO2Cl2, 1);
	re.back().add_product(iSO2, 1);
	
	// R47
 	re.push_back(reaction()); 
	re.back().add_educt(iSO2Cl2, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iClSO2, 2);
	
	// R48
 	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iCl2, 1);
	
	// R49
 	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 1);
	re.back().add_educt(iS, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iSCl, 1);
	
	// R50
 	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iHCl, 1);

	// R51
 	re.push_back(reaction()); 
	re.back().add_educt(iSO2Cl2, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iClSO2, 1);
	re.back().add_product(iCl2, 1);
	
	// R52
 	re.push_back(reaction()); 
	re.back().add_educt(iSO2Cl2, 1);
	re.back().add_educt(iS, 1);
	re.back().add_product(iClSO2, 1);
	re.back().add_product(iSCl, 1);
	
	// R53
 	re.push_back(reaction()); 
	re.back().add_educt(iSO2Cl2, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iClSO2, 1);
	re.back().add_product(iHCl, 1);

	// R54
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCl, 1);
	
	// R55
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iSCl, 1);
	re.back().add_product(iCl, 1);
	
	// R56
 	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_product(iS, 1);
	re.back().add_product(iCl2, 1);
	
	// R57
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iCl, 1);
	
	// R58
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iS, 1);
	
	// R59
 	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iHSCl, 1);
	re.back().add_product(iCl, 1);
	
	// R60
 	re.push_back(reaction()); 
	re.back().add_educt(iHSCl, 1);
	re.back().add_educt(iSH, 1);
	re.back().add_product(iH2S, 1);
	re.back().add_product(iSCl, 1);
	
	// R61
 	re.push_back(reaction()); 
	re.back().add_educt(iHSCl, 1);
	re.back().add_educt(iS, 1);
	re.back().add_product(iSH, 1);
	re.back().add_product(iSCl, 1);
	
	// R62
 	re.push_back(reaction()); 
	re.back().add_educt(iHSCl, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iSCl, 1);
	
	// R63
 	re.push_back(reaction()); 
	re.back().add_educt(iHSCl, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iSCl, 1);
	re.back().add_product(iHCl, 1);
	
	// R64
 	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iHCl, 1);
	
	// R65
 	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iS, 1);
	
	// R66
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iSH, 1);
	
	// R67
 	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iH, 1);
	
	// R68
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iH2, 1);
	
	// R69
 	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iCl, 1);
	
	// R70
 	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iHCl, 1);
	
	// R71
 	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iSH, 1);
	
	// R72
 	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iH2S, 1);
	re.back().add_product(iH, 1);
	
	// R73
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iH2S, 1);
	re.back().add_product(iSH, 1);
	re.back().add_product(iH2, 1);
	
	// R74
 	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iH2S, 1);
	re.back().add_product(iH, 1);
	
	// R75
 	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iH, 1);
	
	// R76
 	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iOH, 1);
	
	// R77
 	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iSH, 1);
	
	// R78
 	re.push_back(reaction()); 
	re.back().add_educt(iSH, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iOCS, 1);
	
	// R79
 	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iH, 1);
	
	// R80
 	re.push_back(reaction()); 
	re.back().add_educt(iSO2, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iOH, 1);
	
	// R81
 	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(irH, 1);
	
	// R82
 	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iS, 1);
	re.back().add_product(iOH, 1);  
}


void build_model_Krasnopolsky_2012(std::vector<species> &sp, std::vector<reaction> &re) {
  
	sp.push_back(species(sp.size(), "CO2"));
        sp.push_back(species(sp.size(), "CO"));
        sp.push_back(species(sp.size(), "O"));
        sp.push_back(species(sp.size(), "Os"));
        sp.push_back(species(sp.size(), "O2"));
        sp.push_back(species(sp.size(), "O2s"));
        sp.push_back(species(sp.size(), "O3"));
        sp.push_back(species(sp.size(), "H2O"));
        sp.push_back(species(sp.size(), "H"));
        sp.push_back(species(sp.size(), "OH"));
        sp.push_back(species(sp.size(), "HO2"));
        sp.push_back(species(sp.size(), "H2"));
        sp.push_back(species(sp.size(), "SO2"));
        sp.push_back(species(sp.size(), "OCS"));
        sp.push_back(species(sp.size(), "SO"));
        sp.push_back(species(sp.size(), "SO3"));
        sp.push_back(species(sp.size(), "HSO3"));
        sp.push_back(species(sp.size(), "H2SO4"));
        sp.push_back(species(sp.size(), "S"));
        sp.push_back(species(sp.size(), "S2"));
        sp.push_back(species(sp.size(), "S3"));
        sp.push_back(species(sp.size(), "Sa"));
        sp.push_back(species(sp.size(), "S2O"));
        sp.push_back(species(sp.size(), "S2O2"));
        sp.push_back(species(sp.size(), "HCl"));
        sp.push_back(species(sp.size(), "Cl"));
        sp.push_back(species(sp.size(), "ClO"));
        sp.push_back(species(sp.size(), "ClO2"));
        sp.push_back(species(sp.size(), "Cl2"));
        sp.push_back(species(sp.size(), "ClCO"));
        sp.push_back(species(sp.size(), "ClCO3"));
        sp.push_back(species(sp.size(), "COCl2"));
        sp.push_back(species(sp.size(), "OSCl"));
        sp.push_back(species(sp.size(), "SCl"));
        sp.push_back(species(sp.size(), "SCl2"));
        sp.push_back(species(sp.size(), "ClSO2"));
        sp.push_back(species(sp.size(), "SO2Cl2"));
        sp.push_back(species(sp.size(), "N"));
        sp.push_back(species(sp.size(), "NO"));
        sp.push_back(species(sp.size(), "NO2"));
        sp.push_back(species(sp.size(), "NO3"));
        sp.push_back(species(sp.size(), "ClNO"));
        sp.push_back(species(sp.size(), "HO2NO2"));
        sp.push_back(species(sp.size(), "SNO"));
        sp.push_back(species(sp.size(), "M"));
        sp.push_back(species(sp.size(), "hv"));
        sp.push_back(species(sp.size(), "N2"));
	
	// IDs
	size_t iCO2=0;
        size_t iCO=1;
        size_t iO=2;
        size_t iOs=3;
        size_t iO2=4;
        size_t iO2s=5;
        size_t iO3=6;
        size_t iH2O=7;
        size_t iH=8;
        size_t iOH=9;
        size_t iHO2=10;
        size_t iH2=11;
        size_t iSO2=12;
        size_t iOCS=13;
        size_t iSO=14;
        size_t iSO3=15;
        size_t iHSO3=16;
        size_t iH2SO4=17;
        size_t iS=18;
        size_t iS2=19;
        size_t iS3=20;
        size_t iSa=21;
        size_t iS2O=22;
        size_t iS2O2=23;
        size_t iHCl=24;
        size_t iCl=25;
        size_t iClO=26;
        size_t iClO2=27;
        size_t iCl2=28;
        size_t iClCO=29;
        size_t iClCO3=30;
        size_t iCOCl2=31;
        size_t iOSCl=32;
        size_t iSCl=33;
        size_t iSCl2=34;
        size_t iClSO2=35;
        size_t iSO2Cl2=36;
        size_t iN=37;
        size_t iNO=38;
        size_t iNO2=39;
        size_t iNO3=40;
        size_t iClNO=41;
        size_t iHO2NO2=42;
        size_t iSNO=43;
        size_t iM=44;
        size_t ihv=45;
        size_t iN2=46;
	
	
	// R1
	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iO, 1);  
	
	// R2
	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iOs, 1);
	
	// R3
	re.push_back(reaction()); 
	re.back().add_educt(iSO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iO, 1);
	
	// R4
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 2);
	
	// R5
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iOs, 1);
	
	// R6
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO2s, 1);
	re.back().add_product(iOs, 1);
	
	// R7
	re.push_back(reaction()); 
	re.back().add_educt(iSO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iO, 1);
	
	// R8
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iS, 1);
	re.back().add_product(iO, 1);
	
	// R9
	re.push_back(reaction()); 
	re.back().add_educt(iOCS, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iS, 1);
	
	// R10
	re.push_back(reaction()); 
	re.back().add_educt(iS2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iS, 1);
	
	// R11
	re.push_back(reaction()); 
	re.back().add_educt(iS2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iS, 2);
	
	// R12
	re.push_back(reaction()); 
	re.back().add_educt(iS3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iS, 1);
	
	// R13
	re.push_back(reaction()); 
	re.back().add_educt(iS2O2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iS, 1);
	
	// R14
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iH, 1);
	
	// R15
	re.push_back(reaction()); 
	re.back().add_educt(iHCl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iCl, 1);
	
	// R16
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iO, 1);
	
	// R17
	re.push_back(reaction()); 
	re.back().add_educt(iClO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iO, 1);
	
	// R18
	re.push_back(reaction()); 
	re.back().add_educt(iClCO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iClO, 1);
	
	// R19
	re.push_back(reaction()); 
	re.back().add_educt(iCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 2);
	
	// R20
	re.push_back(reaction()); 
	re.back().add_educt(iSCl, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iS, 1);
	re.back().add_product(iCl, 1);
	
	// R21
	re.push_back(reaction()); 
	re.back().add_educt(iCOCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iClCO, 1);
	re.back().add_product(iCl, 1);
	
	// R22
	re.push_back(reaction()); 
	re.back().add_educt(iCOCl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iClCO, 1);
	re.back().add_product(iCl, 1);
	
	// R23
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iN, 1);
	re.back().add_product(iO, 1);
	
	// R24
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO, 1);
	
	// R25
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO, 1);
	
	// R26
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);
	
	// R27
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iNO2, 1);
	
	// R28
	re.push_back(reaction()); 
	re.back().add_educt(iClNO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iNO, 1);
	
	// R29
	re.push_back(reaction()); 
	re.back().add_educt(iSNO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iS, 1);
	re.back().add_product(iNO, 1);
	
	// R30
	re.push_back(reaction()); 
	re.back().add_educt(iSO2Cl2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iClSO2, 1);
	re.back().add_product(iCl, 1);
	
	// R31
	re.push_back(reaction()); 
	re.back().add_educt(iOs, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iCO2, 1);
	
	// R32
	re.push_back(reaction()); 
	re.back().add_educt(iOs, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iOH, 2);
	
	// R33
	re.push_back(reaction()); 
	re.back().add_educt(iOs, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iCl, 1);
	
	// R34
	re.push_back(reaction()); 
	re.back().add_educt(iO2s, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(ihv, 1);
	
	// R35
	re.push_back(reaction()); 
	re.back().add_educt(iO2s, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iCO2, 1);
	
	// R36
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iM, 1);
	
	// R37
	re.push_back(reaction()); 
	re.back().add_educt(iO, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iO2s, 1);
	re.back().add_product(iM, 1);
	
	// R38
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iO3, 1);
	re.back().add_product(iM, 1);

	// R39
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO2, 2);
	
	// R40
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iM, 1);
	
	// R41
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 1);
	
	// R42
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iH, 1);
	
	// R43
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iH, 1);
	
	// R44
	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iH, 1);
	
	// R45
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iOH, 2);
	
	// R46
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iO2, 1);
	
	// R47
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO, 1);
	
	// R48
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO2, 1);
	
	// R49
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iO2, 1);
	
	// R50
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iCl, 1);
	
	// R51
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iCl, 1);

	// R52
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHCl, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iCl, 1);
	
	// R53
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iH, 1);
	
	// R54
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iO2, 1);
	
	// R55
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iClO, 1);
	
	// R56
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iO2, 1);
	
	// R57
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iCl, 1);
	
	// R58
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iHCl, 1);
	re.back().add_product(iCl, 1);
	
	// R59
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iCl2, 1);
	re.back().add_product(iM, 1);
	
	// R60
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iO2, 1);
	
	// R61
	re.push_back(reaction()); 
	re.back().add_educt(iClO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iO2, 1);
	
	// R62
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iClCO, 1);
	re.back().add_product(iM, 1);
	
	// R63
	re.push_back(reaction()); 
	re.back().add_educt(iClCO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iM, 1);

	// R64
	re.push_back(reaction()); 
	re.back().add_educt(iClCO, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iCl2, 1);
	
	// R65
	re.push_back(reaction()); 
	re.back().add_educt(iClCO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iCl, 1);

	// R66
	re.push_back(reaction()); 
	re.back().add_educt(iClCO, 1);
	re.back().add_educt(iClNO, 1);
	re.back().add_product(iCOCl2, 1);
	re.back().add_product(iNO, 1);
	
	// R67
	re.push_back(reaction()); 
	re.back().add_educt(iClCO, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iCOCl2, 1);
	re.back().add_product(iCl, 1);
	
	// R68
	re.push_back(reaction()); 
	re.back().add_educt(iClCO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iCl, 1);

	// R69
	re.push_back(reaction()); 
	re.back().add_educt(iClCO, 1);
	re.back().add_educt(iS, 1);
	re.back().add_product(iOCS, 1);
	re.back().add_product(iCl, 1);	

	// R70
	re.push_back(reaction()); 
	re.back().add_educt(iClCO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iClCO3, 1);
	re.back().add_product(iM, 1);

	// R71
	re.push_back(reaction()); 
	re.back().add_educt(iClCO3, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iClO, 1);

	// R72
	re.push_back(reaction()); 
	re.back().add_educt(iClCO3, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iCl, 1);	
	re.back().add_product(iO2, 1);
	
	// R73
	re.push_back(reaction()); 
	re.back().add_educt(iClCO3, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iCl, 1);	
	re.back().add_product(iOH, 1);	
	
	// R74
	re.push_back(reaction()); 
	re.back().add_educt(iCOCl2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iClCO, 1);
	re.back().add_product(iClO, 1);
	
	// R75
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(ihv, 1);
	
	// R76
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iM, 1);	
	
	// R77
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO, 1);
	
	// R78
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iM, 1);
	
	// R79
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);
	
	// R80
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iOH, 1);
	
	// R81
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);
	
	// R82
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iNO3, 1);
	re.back().add_product(iM, 1);
	
	// R83
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iNO2, 2);
	
	// R84
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);
	
	// R85
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO, 1);
	
	// R86
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iClNO, 1);
	re.back().add_product(iM, 1);
	
	// R87
	re.push_back(reaction()); 
	re.back().add_educt(iClNO, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iCl2, 1);
	re.back().add_product(iNO, 1);
	
	// R88
	re.push_back(reaction()); 
	re.back().add_educt(iClNO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iClO, 1);
	
	// R89
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iCl, 1);
	
	// R90
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iClO2, 1);
	re.back().add_product(iClNO, 1);
	re.back().add_product(iO2, 1);
	
	// R91
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iOSCl, 1);
	
	// R92
	re.push_back(reaction()); 
	re.back().add_educt(iSO2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iSO3, 1);
	re.back().add_product(iM, 1);
	
	// R93
	re.push_back(reaction()); 
	re.back().add_educt(iSO2, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHSO3, 1);
	re.back().add_product(iM, 1);

	// R94
	re.push_back(reaction()); 
	re.back().add_educt(iSO2, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iClSO2, 1);
	re.back().add_product(iM, 1);

	// R95
	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iHCl, 1);
	
	// R96
	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iCl2, 1);
	
	// R97
	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iClO, 1);
	
 	// R98
	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iClO, 1);
	
	// R99
	re.push_back(reaction()); 
	re.back().add_educt(iHSO3, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iSO3, 1);
	re.back().add_product(iHO2, 1);
	
	// R100
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iO, 1);
	
	// R101
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iO2, 1);
	
	// R102
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iSNO, 1);
	re.back().add_product(iM, 1);
	
	// R103
	re.push_back(reaction()); 
	re.back().add_educt(iSNO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iSO, 1);

	// R104
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iM, 1);

	// R105
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iO, 1);
	
	// R106
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iO2, 1);
	
	// R107
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iH, 1);
	
	// R108
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iClO, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iCl, 1);
	
	// R109
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iClO2, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iClO, 1);
	
	// R110
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iOH, 1);

	// R111
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 2);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iS, 1);
	
	// R112
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iNO, 1);
	
	// R113
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iSO3, 1);
	re.back().add_product(iSO2, 2);
	
	// R114
	re.push_back(reaction()); 
	re.back().add_educt(iSO3, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iO2, 1);

 	// R115
	re.push_back(reaction()); 
	re.back().add_educt(iSO3, 1);
	re.back().add_educt(iH2O, 2);
	re.back().add_product(iH2SO4, 1);
	re.back().add_product(iH2O, 1);

	// R116
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iClO2, 1);
	re.back().add_product(iM, 1);
	
	// R117
	re.push_back(reaction()); 
	re.back().add_educt(iClO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iM, 1);
	
	// R118
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iCl, 1);
	
	// R119
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iM, 1);
	
	// R120
	re.push_back(reaction()); 
	re.back().add_educt(iS, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iM, 1);
	
	// R121
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iS2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iS3, 1);
	re.back().add_product(iM, 1);
	
	// R122
	re.push_back(reaction()); 
	re.back().add_educt(iS2, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iSa, 1);
	re.back().add_product(iM, 1);
	
	// R123
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iS3, 1);
	re.back().add_product(iS2, 2);
	
	// R124
	re.push_back(reaction()); 
	re.back().add_educt(iS2, 1);
	re.back().add_educt(iS3, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iSa, 1);
	re.back().add_product(iM, 1);
	
	// R125
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iS2, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iS, 1);
	
	// R126
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iS3, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iS2, 1);
	
	// R127
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iSa, 1);
	re.back().add_product(iSO, 1);
	
	// R128
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iS2O2, 1);
	re.back().add_product(iM, 1);
	
	// R129
	re.push_back(reaction()); 
	re.back().add_educt(iS2O2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iS2O, 1);
	re.back().add_product(iO2, 1);
	
	// R130
	re.push_back(reaction()); 
	re.back().add_educt(iS2O2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iSO, 1);
	
	// R131
	re.push_back(reaction()); 
	re.back().add_educt(iS2O2, 1);
	re.back().add_educt(iSO, 1);
	re.back().add_product(iS2O, 1);
	re.back().add_product(iSO2, 1);
	
	// R132
	re.push_back(reaction()); 
	re.back().add_educt(iS2O, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iSO, 2);
	
	// R133
	re.push_back(reaction()); 
	re.back().add_educt(iS2O, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iSCl, 1);
	re.back().add_product(iSO, 1);
	
	// R134
	re.push_back(reaction()); 
	re.back().add_educt(iClCO3, 1);
	re.back().add_educt(iSO, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iCO2, 1);
	
	// R135
	re.push_back(reaction()); 
	re.back().add_educt(iClCO3, 1);
	re.back().add_educt(iSO2, 1);
	re.back().add_product(iCl, 1);
	re.back().add_product(iSO3, 1);
	re.back().add_product(iCO2, 1);

	// R136
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iCO, 1);
	
	// R137
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iOCS, 1);
	re.back().add_product(iS2, 1);
	re.back().add_product(iCO, 1);
	
	// R138
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iOCS, 1);
	re.back().add_product(iM, 1);
	
	// R139
	re.push_back(reaction()); 
	re.back().add_educt(iS, 1);
	re.back().add_educt(iCl2, 1);
	re.back().add_product(iSCl, 1);
	re.back().add_product(iCl, 1);
	
	// R140
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iCl, 1);
	
	// R141
	re.push_back(reaction()); 
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iSCl2, 1);
	re.back().add_product(iM, 1);
	
	// R142
	re.push_back(reaction()); 
	re.back().add_educt(iCl2, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_product(iSCl2, 1);
	re.back().add_product(iCl, 1);
	
	// R143
	re.push_back(reaction()); 
	re.back().add_educt(iClO2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iClO, 1);
	re.back().add_product(iO2, 1);
	
	// R144
	re.push_back(reaction()); 
	re.back().add_educt(iClO2, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iCl2, 1);
	re.back().add_product(iO2, 1);
	
	// R145
	re.push_back(reaction()); 
	re.back().add_educt(iClO2, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iClO, 2);
	
	// R146
	re.push_back(reaction()); 
	re.back().add_educt(iSO, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iOSCl, 1);
	re.back().add_product(iM, 1);
	
	// R147
	re.push_back(reaction()); 
	re.back().add_educt(iOSCl, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iSO2, 1);
	re.back().add_product(iCl, 1);
	
	// R148
	re.push_back(reaction()); 
	re.back().add_educt(iOSCl, 1);
	re.back().add_educt(iCl, 1);
	re.back().add_product(iSO, 1);
	re.back().add_product(iCl2, 1);
	
	// R149
	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 2);
	re.back().add_product(iSO2Cl2, 1);
	re.back().add_product(iSO2, 1);
	
	// R150
	re.push_back(reaction()); 
	re.back().add_educt(iClSO2, 1);
	re.back().add_educt(iSCl, 1);
	re.back().add_product(iSO2Cl2, 1);
	re.back().add_product(iS, 1);
	
	// R151
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2NO2, 1);
	re.back().add_product(iM, 1);
	
	// R152
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iM, 1);
	
	// R153
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);
    re.back().add_product(iH2O, 1);
}


void build_model_Kasting_1985(std::vector<species> &sp, std::vector<reaction> &re) {
    	sp.push_back(species(sp.size(), "H2O"));
	sp.push_back(species(sp.size(), "hv"));
	sp.push_back(species(sp.size(), "H"));
	sp.push_back(species(sp.size(), "OH"));
	sp.push_back(species(sp.size(), "O2"));
	sp.push_back(species(sp.size(), "M"));
	sp.push_back(species(sp.size(), "HO2"));
	sp.push_back(species(sp.size(), "O3"));
	sp.push_back(species(sp.size(), "O1D"));
	sp.push_back(species(sp.size(), "CO"));
	sp.push_back(species(sp.size(), "CO2"));
	sp.push_back(species(sp.size(), "NO"));
	sp.push_back(species(sp.size(), "NO2"));
	sp.push_back(species(sp.size(), "H2O2"));
	sp.push_back(species(sp.size(), "HCO"));
	sp.push_back(species(sp.size(), "H2CO"));
	sp.push_back(species(sp.size(), "CH4"));
	sp.push_back(species(sp.size(), "CH3O2"));
	sp.push_back(species(sp.size(), "CH3OOH"));
	sp.push_back(species(sp.size(), "N2O"));
	sp.push_back(species(sp.size(), "HNO2"));
	sp.push_back(species(sp.size(), "HNO3"));
	sp.push_back(species(sp.size(), "H2"));
	sp.push_back(species(sp.size(), "O"));
	sp.push_back(species(sp.size(), "N2"));
	sp.push_back(species(sp.size(), "CH3"));
	sp.push_back(species(sp.size(), "H3CO"));
	sp.push_back(species(sp.size(), "HO2NO2"));
	sp.push_back(species(sp.size(), "N"));
	sp.push_back(species(sp.size(), "CH2_1"));
	sp.push_back(species(sp.size(), "CH2_3"));
	
	// IDs
	size_t iH2O=0;
	size_t ihv=1;
	size_t iH=2;
	size_t iOH=3;
	size_t iO2=4;
	size_t iM=5;
	size_t iHO2=6;
	size_t iO3=7;
	size_t iO1D=8;
	size_t iCO=9;
	size_t iCO2=10;
	size_t iNO=11;
	size_t iNO2=12;
	size_t iH2O2=13;
	size_t iHCO=14;
	size_t iH2CO=15;
	size_t iCH4=16;
	size_t iCH3O2=17;
	size_t iCH3OOH=18;
	size_t iN2O=19;
	size_t iHNO2=20;
	size_t iHNO3=21;
	size_t iH2=22;
	size_t iO=23;
	size_t iN2=24;
	size_t iCH3=25;
	size_t iH3CO=26;
	size_t iHO2NO2=27;
	size_t iN=28;
	size_t iCH2_1=29;
	size_t iCH2_3=30;
	
	
	
	// R1
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(iO1D, 1);
	re.back().add_product(iOH, 2);  
	
	// R2
	re.push_back(reaction()); 
	re.back().add_educt(iH2, 1);
	re.back().add_educt(iO1D, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iH, 1);	
	
	// R3
	re.push_back(reaction()); 
	re.back().add_educt(iH2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iH, 1);	
	
	// R4
	re.push_back(reaction()); 
	re.back().add_educt(iH2, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iH, 1);
	
	// R5
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 1);
	
	// R6
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iM, 1);
		
	// R7
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iO2, 1);
		
	// R8
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO, 1);
		
	// R9
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iOH, 2);
	
	// R10
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iO2, 1);
	
	// R11
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO2, 1);
	
	// R12
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iO2, 1);
	
	// R13
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 1);
	
	// R14
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 2);
	
	// R15
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 2);
	re.back().add_product(iH2O2, 1);
	re.back().add_product(iO2, 1);
	
	// R16
	re.push_back(reaction()); 
	re.back().add_educt(iH2O2, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iH2O, 1);
	
	// R17
	re.push_back(reaction()); 
	re.back().add_educt(iH2O2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iOH, 1);
	
	// R18
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iH2O2, 1);
	re.back().add_product(iM, 1);
	
	// R19
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO2, 2);
	
	// R20
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 2);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO, 1);
	
	// R21
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iN2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iN2, 1);
	
	// R22
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iO2, 1);
	
	// R23
	re.push_back(reaction()); 
	re.back().add_educt(iO, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iM, 1);
	
	// R24
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iO3, 1);
	re.back().add_product(iM, 1);
	
	
	// R25
	re.push_back(reaction()); 
	re.back().add_educt(iH, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iM, 1);
	
	// R26
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iM, 1);

	// R27
	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iH, 1);
	
	// R28
	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iM, 1);

	// R29
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHCO, 1);
	re.back().add_product(iM, 1);
	
	// R30
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHCO, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iCO, 1);
	
	// R31
	re.push_back(reaction()); 
	re.back().add_educt(iHCO, 2);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iCO, 1);
	
	// R32
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHCO, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iCO, 1);
	
	// R33
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHCO, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iCO2, 1);
	
	// R34
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHCO, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iCO, 1);
	
	// R35
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iHCO, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iCO, 1);
	
	// R36
	re.push_back(reaction()); 
	re.back().add_educt(iH2CO, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iHCO, 1);
	
	// R37
	re.push_back(reaction()); 
	re.back().add_educt(iH2CO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iHCO, 1);
	
	// R38
	re.push_back(reaction()); 
	re.back().add_educt(iH2CO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iHCO, 1);
	re.back().add_product(iOH, 1);

	// R39
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iO1D, 1);
	
	// R40
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 2);
	
	// R41
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iOH, 1);

	// R42
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iO1D, 1);
	
	// R43
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iO, 1);
	
	// R44
	re.push_back(reaction()); 
	re.back().add_educt(iH2O2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 2);
	
	// R45
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO, 1);
	
	// R46
	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iO1D, 1);
	
	// R47
	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iO, 1);
	
	// R48
	re.push_back(reaction()); 
	re.back().add_educt(iH2CO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iCO, 1);
	
	// R49
	re.push_back(reaction()); 
	re.back().add_educt(iH2CO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iHCO, 1);
	re.back().add_product(iH, 1);

	// R50
	re.push_back(reaction()); 
	re.back().add_educt(iHCO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iCO, 1);
	
	// R51
	re.push_back(reaction()); 
	re.back().add_educt(iCH4, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCH2_1, 1);
	re.back().add_product(iH2, 1);
	
	// R52
	re.push_back(reaction()); 
	re.back().add_educt(iCH3OOH, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH3CO, 1);
	re.back().add_product(iOH, 1);

	// R53
	re.push_back(reaction()); 
	re.back().add_educt(iN2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO, 1);
	
	// R54
	re.push_back(reaction()); 
	re.back().add_educt(iHNO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iOH, 1);
	
	// R55
	re.push_back(reaction()); 
	re.back().add_educt(iHNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iOH, 1);
	
	// R56
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iNO2, 1);
	
	// R57
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iN, 1);
	re.back().add_product(iO, 1);
	
	// R58
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO, 1);
	
	// R59
	re.push_back(reaction()); 
	re.back().add_educt(iCH4, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iCH3, 1);
	re.back().add_product(iH2O, 1);
	
	// R60
	re.push_back(reaction()); 
	re.back().add_educt(iCH4, 1);
	re.back().add_educt(iO1D, 1);
	re.back().add_product(iCH3, 1);
	re.back().add_product(iOH, 1);
	
	// R61
	re.push_back(reaction()); 
	re.back().add_educt(iCH4, 1);
	re.back().add_educt(iO1D, 1);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iH2, 1);

	// R62
	re.push_back(reaction()); 
	re.back().add_educt(iCH2_1, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iCH3, 2);
	
	// R63
	re.push_back(reaction()); 
	re.back().add_educt(iCH2_1, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iO, 1);
	
	// R64
	re.push_back(reaction()); 
	re.back().add_educt(iCH2_1, 1);
	re.back().add_educt(iN2, 1);
	re.back().add_product(iCH2_3, 1);
	re.back().add_product(iN2, 1);
	
	// R65
	re.push_back(reaction()); 
	re.back().add_educt(iCH2_3, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iCH3, 1);
	re.back().add_product(iH, 1);
	
	// R66
	re.push_back(reaction()); 
	re.back().add_educt(iCH2_3, 1);
	re.back().add_educt(iCH4, 1);
	re.back().add_product(iCH3, 2);
	
	// R67
	re.push_back(reaction()); 
	re.back().add_educt(iCH2_3, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iO, 1);
	
	// R68
	re.push_back(reaction()); 
	re.back().add_educt(iCH3, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iCH3O2, 1);
	re.back().add_product(iM, 1);
	
	// R69
	re.push_back(reaction()); 
	re.back().add_educt(iCH3, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iH, 1);
	
	// R70
	re.push_back(reaction()); 
	re.back().add_educt(iCH3, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iH2, 1);
		
	// R71
	re.push_back(reaction()); 
	re.back().add_educt(iCH3, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iHO2, 1);
		
	// R72
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iCH3OOH, 1);
	re.back().add_product(iO2, 1);
	
	// R73
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 2);
	re.back().add_product(iH3CO, 2);
	re.back().add_product(iO2, 1);
		
	// R74
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iH3CO, 1);
	re.back().add_product(iNO2, 1);
		
	// R75
	re.push_back(reaction()); 
	re.back().add_educt(iCH3O2, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iH3CO, 1);
	re.back().add_product(iHO2, 1);
	
	// R76
	re.push_back(reaction()); 
	re.back().add_educt(iCH3OOH, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iCH3O2, 1);
	re.back().add_product(iH2O, 1);
	
	// R77
	re.push_back(reaction()); 
	re.back().add_educt(iH3CO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iHO2, 1);
	
	// R78
	re.push_back(reaction()); 
	re.back().add_educt(iH3CO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iOH, 1);
	
	// R79
	re.push_back(reaction()); 
	re.back().add_educt(iH3CO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iH2CO, 1);
	re.back().add_product(iH2O, 1);
	
	// R80
	re.push_back(reaction()); 
	re.back().add_educt(iN2O, 1);
	re.back().add_educt(iO1D, 1);
	re.back().add_product(iNO, 2);
	
	// R81
	re.push_back(reaction()); 
	re.back().add_educt(iN2O, 1);
	re.back().add_educt(iO1D, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO2, 1);
	
	// R82
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO, 1);
	
	// R83
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);
	
	// R84
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iH, 1);
		
	// R84
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO, 1);
		
	// R85
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);
		
	// R86
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);
		
	// R87
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iM, 1);
		
	// R88
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iOH, 1);
		
	// R89
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHNO2, 1);
	re.back().add_product(iM, 1);
		
	// R90
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);
		
	// R91
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHNO3, 1);
	re.back().add_product(iM, 1);
		
	// R92
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2NO2, 1);
	re.back().add_product(iM, 1);
		
	// R93
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iOH, 1);
	
	// R94
	re.push_back(reaction()); 
	re.back().add_educt(iHNO3, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO, 1);
	
	// R95
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);

	// R96
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 1);
		
	// R97
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iM, 1);
}



void build_model_Nair_1994(std::vector<species> &sp, std::vector<reaction> &re) {
    	sp.push_back(species(sp.size(), "O2"));
	sp.push_back(species(sp.size(), "hv"));
	sp.push_back(species(sp.size(), "O"));
	sp.push_back(species(sp.size(), "O1D"));
	sp.push_back(species(sp.size(), "O3"));
	sp.push_back(species(sp.size(), "H2"));
	sp.push_back(species(sp.size(), "H"));
	sp.push_back(species(sp.size(), "OH"));
	sp.push_back(species(sp.size(), "HO2"));
	sp.push_back(species(sp.size(), "H2O"));
	sp.push_back(species(sp.size(), "H2O2"));
	sp.push_back(species(sp.size(), "CO2"));
	sp.push_back(species(sp.size(), "CO"));
	sp.push_back(species(sp.size(), "M"));
	sp.push_back(species(sp.size(), "N2"));
	sp.push_back(species(sp.size(), "N"));
	sp.push_back(species(sp.size(), "N2D"));
	sp.push_back(species(sp.size(), "NO"));
	sp.push_back(species(sp.size(), "NO2"));
	sp.push_back(species(sp.size(), "NO3"));
	sp.push_back(species(sp.size(), "N2O"));
	sp.push_back(species(sp.size(), "N2O5"));
	sp.push_back(species(sp.size(), "HNO2"));
	sp.push_back(species(sp.size(), "HNO3"));
	sp.push_back(species(sp.size(), "HO2NO2"));
	sp.push_back(species(sp.size(), "Op"));
	sp.push_back(species(sp.size(), "O2p"));
	sp.push_back(species(sp.size(), "em"));
	sp.push_back(species(sp.size(), "CO2p"));
	sp.push_back(species(sp.size(), "CO2Hp"));
	sp.push_back(species(sp.size(), "igrain"));
	sp.push_back(species(sp.size(), "H2Ograin"));
	
	// IDs
	size_t iO2=0;
    	size_t ihv=1;
    	size_t iO=2;
    	size_t iO1D=3;
    	size_t iO3=4;
    	size_t iH2=5;
    	size_t iH=6;
    	size_t iOH=7;
    	size_t iHO2=8;
    	size_t iH2O=9;
    	size_t iH2O2=10;
    	size_t iCO2=11;
    	size_t iCO=12;
    	size_t iM=13;
    	size_t iN2=14;
    	size_t iN=15;
    	size_t iN2D=16;
    	size_t iNO=17;
    	size_t iNO2=18;
    	size_t iNO3=19;
    	size_t iN2O=20;
    	size_t iN2O5=21;
    	size_t iHNO2=22;
    	size_t iHNO3=23;
    	size_t iHO2NO2=24;
    	size_t iOp=25;
    	size_t iO2p=26;
    	size_t iem=27;
    	size_t iCO2p=28;
	size_t iCO2Hp=29;
	size_t igrain=30;
	size_t iHO2grain=31;
    	
	
	// R1
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 2);  

  	// R2
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iO1D, 1);
	
  	// R3
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iO, 1);
	
  	// R4
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iO1D, 1);

	// R5
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iO1D, 1);
	
	// R6
	re.push_back(reaction()); 
	re.back().add_educt(iH2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 2);
	
  	// R7
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iH, 1);

	// R8
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO, 1);
	
  	// R9
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 1);
	re.back().add_product(iOH, 1);
	
  	// R10
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iO1D, 1);
	
  	// R11
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iH, 2);
	re.back().add_product(iO, 1);
	
  	// R12
	re.push_back(reaction()); 
	re.back().add_educt(iH2O2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 2);
	
  	// R13
	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iO, 1);

	// R14
	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iO1D, 1);

	// R15
	re.push_back(reaction()); 
	re.back().add_educt(iO, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iM, 1);

	// R16
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iN2, 1);
	re.back().add_product(iO3, 1);
	re.back().add_product(iN2, 1);
	
	// R17
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iO3, 1);
	re.back().add_product(iCO2, 1);

	// R18
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO2, 2);

	// R19
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iM, 1);
	
	// R20
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iO2, 1);

	// R21
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO2, 2);
	
	// R22
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iO, 2);
	
	// R23
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iCO2, 1);

	// R24
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iCO2, 1);
	
	// R25
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iH2O, 1);
	re.back().add_product(iOH, 2);
	
	// R26
	re.push_back(reaction()); 
	re.back().add_educt(iH, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iM, 1);
	
	// R27
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iM, 1);

	// R28
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 1);
	
	// R29
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iOH, 2);
	
	// R30
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iH2, 1);
	re.back().add_product(iO2, 1);
	
	// R31
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO, 1);
	
	// R32
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iH, 1);
	
	// R33
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iH, 1);
	
	// R34
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 1);
	
	// R35
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iH2O2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iHO2, 1);
	
	// R36
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 2);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO, 1);
	
	// R37
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iH2O2, 1);
	re.back().add_product(iM, 1);
	
	// R38
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iO2, 1);
	
	// R39
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iH, 1);
	
	// R40
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO2, 1);
	
	// R41
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iH2O2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iHO2, 1);
	
	// R42
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iCO, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iH, 1);
	
	// R43
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iO2, 2);
	
	// R44
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 2);
	re.back().add_product(iH2O2, 1);
	re.back().add_product(iO2, 1);

	// R45
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iH2O2, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iM, 1);
	
	// R46
	re.push_back(reaction()); 
	re.back().add_educt(iN2, 1);
	re.back().add_product(iN, 2);

	// R47
	re.push_back(reaction()); 
	re.back().add_educt(iN2, 1);
	re.back().add_product(iN2D, 2);

	// R48
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iN, 1);
	re.back().add_product(iO, 1);
	
	// R49
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO, 2);
	
	// R50
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO, 1);
	
	// R51
	re.push_back(reaction()); 
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);
	
	// R52
	re.push_back(reaction()); 
	re.back().add_educt(iN2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO1D, 1);

	// R53
	re.push_back(reaction()); 
	re.back().add_educt(iN2O5, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iNO3, 1);
	
	// R54
	re.push_back(reaction()); 
	re.back().add_educt(iHNO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iNO, 1);

	// R55
	re.push_back(reaction()); 
	re.back().add_educt(iHNO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iOH, 1);

	// R56
	re.push_back(reaction()); 
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iNO2, 1);
	
	// R57
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO, 1);
	
	// R58
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);

	// R59
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iH, 1);

	// R60
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iOH, 1);

	// R61
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO, 1);
	
	// R62
	re.push_back(reaction()); 
	re.back().add_educt(iN, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iN2O, 1);
	re.back().add_product(iO, 1);
	
	// R63
	re.push_back(reaction()); 
	re.back().add_educt(iN2D, 1);
	re.back().add_educt(iO, 1);
	re.back().add_product(iN, 1);
	re.back().add_product(iO, 1);

	// R64
	re.push_back(reaction()); 
	re.back().add_educt(iN2D, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iCO, 1);
	
	// R65
	re.push_back(reaction()); 
	re.back().add_educt(iN2D, 1);
	re.back().add_educt(iN2, 1);
	re.back().add_product(iN, 1);
	re.back().add_product(iN2, 1);
	
	// R66
	re.push_back(reaction()); 
	re.back().add_educt(iN2D, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO, 1);
	
	// R67
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iM, 1);

	// R68
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iO2, 1);
	
	// R69
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iNO3, 1);
	re.back().add_product(iM, 1);

	// R70
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);
	
	// R71
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);

	// R72
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iN2, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iN2, 1);

	// R73
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iN2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iN2O, 1);
	re.back().add_product(iM, 1);
	
	// R74
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iN2O, 1);
	re.back().add_product(iNO, 2);
	
	// R75
	re.push_back(reaction()); 
	re.back().add_educt(iO1D, 1);
	re.back().add_educt(iN2O, 1);
	re.back().add_product(iN2, 1);
	re.back().add_product(iO2, 1);

	// R76
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);
	
	// R77
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iOH, 1);

	// R78
	re.push_back(reaction()); 
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iNO2, 2);
	
	// R79
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iNO, 1);
	
	// R80
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iOH, 1);
	re.back().add_product(iNO2, 1);

	// R81
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iNO, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHNO2, 1);
	re.back().add_product(iM, 1);
	
	// R82
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHNO3, 1);
	re.back().add_product(iM, 1);

	// R83
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iNO2, 1);

	// R84
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHNO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iNO2, 1);

	// R85
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHNO3, 1);
	re.back().add_product(iNO3, 1);
	re.back().add_product(iH2O, 1);

 	// R86
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHO2NO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);

	// R87
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2NO2, 1);
	re.back().add_product(iM, 1);

	// R88
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iHNO3, 1);
	
	// R89
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iNO3, 1);
	re.back().add_product(iO2, 1);

	// R90
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iN2O5, 1);
	re.back().add_product(iM, 1);
	
	// R91
	re.push_back(reaction()); 
	re.back().add_educt(iNO2, 1);
	re.back().add_educt(iNO3, 1);
	re.back().add_product(iNO, 1);
	re.back().add_product(iNO2, 1);
	re.back().add_product(iO2, 1);	

	// R92
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOp, 1);
	re.back().add_product(iem, 1);

	// R92
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iOp, 1);
	re.back().add_product(iem, 1);
	
	// R93
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO2p, 1);
	re.back().add_product(iem, 1);
	
	// R94
	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO2p, 1);
	re.back().add_product(iem, 1);
	
	// R95
	re.push_back(reaction()); 
	re.back().add_educt(iCO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iOp, 1);
	re.back().add_product(iem, 1);
	
	// R96
	re.push_back(reaction()); 
	re.back().add_educt(iO2p, 1);
	re.back().add_educt(iem, 1);
	re.back().add_product(iO, 2);
	
	// R97
	re.push_back(reaction()); 
	re.back().add_educt(iCO2p, 1);
	re.back().add_educt(iem, 1);
	re.back().add_product(iCO, 1);
	re.back().add_product(iO, 1);
	
	// R98
	re.push_back(reaction()); 
	re.back().add_educt(iOp, 1);
	re.back().add_educt(iCO2, 1);
	re.back().add_product(iO2p, 1);
	re.back().add_product(iCO, 1);
	
	// R99
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCO2p, 1);
	re.back().add_product(iO2p, 1);
	re.back().add_product(iCO, 1);

	// R100
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iCO2p, 1);
	re.back().add_product(iOp, 1);
	re.back().add_product(iCO2, 1);
	
	// R101
	re.push_back(reaction()); 
	re.back().add_educt(iCO2p, 1);
	re.back().add_educt(iH2, 1);
	re.back().add_product(iCO2Hp, 1);
	re.back().add_product(iH, 1);
	
	// R102
	re.push_back(reaction()); 
	re.back().add_educt(iCO2Hp, 1);
	re.back().add_educt(iem, 1);
	re.back().add_product(iCO2, 1);
	re.back().add_product(iH, 1);
	
	// R103
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 1);
	re.back().add_educt(igrain, 1);
	re.back().add_product(iHO2grain, 1);
	
	// R104
	re.push_back(reaction()); 
	re.back().add_educt(iHO2grain, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iO2, 1);	
}



// TODO Vervollständigen und überprüfen
void build_model_Pavlov_2001(std::vector<species> &sp, std::vector<reaction> &re) {
    sp.push_back(species(sp.size(), "hv"));
    sp.push_back(species(sp.size(), "O"));
    sp.push_back(species(sp.size(), "O1D"));
    sp.push_back(species(sp.size(), "O2"));
    sp.push_back(species(sp.size(), "O3"));
    sp.push_back(species(sp.size(), "C"));
    sp.push_back(species(sp.size(), "H"));
    sp.push_back(species(sp.size(), "CO"));
    sp.push_back(species(sp.size(), "CO2"));
    sp.push_back(species(sp.size(), "HCO"));
    sp.push_back(species(sp.size(), "H2CO"));
    sp.push_back(species(sp.size(), "CH"));
    sp.push_back(species(sp.size(), "1CH2"));
    sp.push_back(species(sp.size(), "3CH2"));
    sp.push_back(species(sp.size(), "CH3"));
    sp.push_back(species(sp.size(), "CH4"));
    sp.push_back(species(sp.size(), "C2H2"));
    sp.push_back(species(sp.size(), "C2H4"));
    sp.push_back(species(sp.size(), "C2H6"));
    sp.push_back(species(sp.size(), "C3H3"));
    sp.push_back(species(sp.size(), "C3H6"));
    sp.push_back(species(sp.size(), "C3H8"));
    sp.push_back(species(sp.size(), "CH2CCH2"));
    sp.push_back(species(sp.size(), "CH3C2H"));
    sp.push_back(species(sp.size(), "H2"));
    sp.push_back(species(sp.size(), "C2H"));
    sp.push_back(species(sp.size(), "C2"));
    sp.push_back(species(sp.size(), "S"));
    sp.push_back(species(sp.size(), "S2"));
    sp.push_back(species(sp.size(), "S3"));
    sp.push_back(species(sp.size(), "S4"));
    sp.push_back(species(sp.size(), "HSO"));
    sp.push_back(species(sp.size(), "H2SO4"));
    sp.push_back(species(sp.size(), "HS"));
    sp.push_back(species(sp.size(), "M"));
    sp.push_back(species(sp.size(), "H2S"));
    sp.push_back(species(sp.size(), "1SO2"));
    sp.push_back(species(sp.size(), "3SO2"));
    sp.push_back(species(sp.size(), "SO"));
    sp.push_back(species(sp.size(), "SO2"));
    sp.push_back(species(sp.size(), "HNO"));
    sp.push_back(species(sp.size(), "HNO2"));
    sp.push_back(species(sp.size(), "HNO3"));
    sp.push_back(species(sp.size(), "C3H2"));
    sp.push_back(species(sp.size(), "C2H5"));
    sp.push_back(species(sp.size(), "CH2CO"));
    sp.push_back(species(sp.size(), "CH3CHO"));
    sp.push_back(species(sp.size(), "C2H5CHO"));
    sp.push_back(species(sp.size(), "H2O"));
    sp.push_back(species(sp.size(), "OH"));
    sp.push_back(species(sp.size(), "H2O2"));
    sp.push_back(species(sp.size(), "HO2"));
    sp.push_back(species(sp.size(), "NO"));
    sp.push_back(species(sp.size(), "N"));
    sp.push_back(species(sp.size(), "NH"));
    sp.push_back(species(sp.size(), "NH2"));
    sp.push_back(species(sp.size(), "NH2s"));
    sp.push_back(species(sp.size(), "NH3"));
    sp.push_back(species(sp.size(), "NO2"));
    sp.push_back(species(sp.size(), "N2H3"));
    sp.push_back(species(sp.size(), "N2H4"));

    // IDs
    size_t ihv=0;
    size_t iO=1;
    size_t iO1D=2;
    size_t iO2=3;
    size_t iO3=4;
    size_t iC=5;
    size_t iH=6;
    size_t iCO=7;
    size_t iCO2=8;
    size_t iHCO=9;
    size_t iH2CO=10;
    size_t iCH=11;
    size_t i1CH2=12;
    size_t i3CH2=13;
    size_t iCH3=14;
    size_t iCH4=15;
    size_t iC2H2=16;
    size_t iC2H4=17;
    size_t iC2H6=18;
    size_t iC3H3=19;
    size_t iC3H6=20;
    size_t iC3H8=21;
    size_t iCH2CCH2=22;
    size_t iCH3C2H=23;
    size_t iH2=24;
    size_t iC2H=25;
    size_t iC2=26;
    size_t iS=27;
    size_t iS2=28;
    size_t iS3=29;
    size_t iS4=30;
    size_t iHSO=31;
    size_t iH2SO4=32;
    size_t iHS=33;
    size_t iM=34;
    size_t iH2S=35;
    size_t i1SO2=36;
    size_t i3SO2=37;
    size_t iSO=38;
    size_t iSO2=39;
    size_t iHNO=40;
    size_t iHNO2=41;
    size_t iHNO3=42;
    size_t iC3H2=43;
    size_t iC2H5=44;
    size_t iCH2CO=45;
    size_t iCH3CHO=46;
    size_t iC2H5CHO=47;
    size_t iH2O=48;
    size_t iOH=49;
    size_t iH2O2=50;
    size_t iHO2=51;
    size_t iNO=52;
    size_t iN=53;
    size_t iNH=54;
    size_t iNH2=55;
    size_t iNH2s=56;
    size_t iNH3=57;
    size_t iNO2=58;
    size_t iN2H3=59;
    size_t iN2H4=60;


    // R1
    re.push_back(reaction()); 
    re.back().add_educt(iCH, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC, 1);
    re.back().add_product(iH, 1);  

    // R2
    re.push_back(reaction()); 
    re.back().add_educt(iCH3, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(i1CH2, 1);
    re.back().add_product(iH, 1); 
 
    // R3
    re.push_back(reaction()); 
    re.back().add_educt(iCH4, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(i1CH2, 1);
    re.back().add_product(iH2, 1); 

    // R4
    re.push_back(reaction()); 
    re.back().add_educt(iCH4, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCH3, 1);
    re.back().add_product(iH, 1); 

    // R5
    re.push_back(reaction()); 
    re.back().add_educt(iCH4, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCH3, 1);
    re.back().add_product(iH, 1); 

    // R6
    re.push_back(reaction()); 
    re.back().add_educt(iC2H2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H, 1);
    re.back().add_product(iH, 1); 

    // R7
    re.push_back(reaction()); 
    re.back().add_educt(iC2H2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2, 1);
    re.back().add_product(iH2, 1); 

    // R8
    re.push_back(reaction()); 
    re.back().add_educt(iC2H4, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H2, 1);
    re.back().add_product(iH2, 1); 

    // R9
    re.push_back(reaction()); 
    re.back().add_educt(iC2H4, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H2, 1);
    re.back().add_product(iH, 2); 

    // R10
    re.push_back(reaction()); 
    re.back().add_educt(iC2H6, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H2, 1);
    re.back().add_product(iH2, 2); 

    // R11
    re.push_back(reaction()); 
    re.back().add_educt(iC2H6, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H4, 1);
    re.back().add_product(iH, 2); 

    // R12
    re.push_back(reaction()); 
    re.back().add_educt(iC2H6, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H4, 1);
    re.back().add_product(iH2, 1); 

    // R13
    re.push_back(reaction()); 
    re.back().add_educt(iC2H6, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCH3, 2);

    // R14
    re.push_back(reaction()); 
    re.back().add_educt(iC2H6, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCH4, 1);
    re.back().add_product(i1CH2, 1); 

    // R15
    re.push_back(reaction()); 
    re.back().add_educt(iC3H3, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC3H2, 1);
    re.back().add_product(iH, 1); 

    // R16
    re.push_back(reaction()); 
    re.back().add_educt(iCH2CCH2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC3H3, 1);
    re.back().add_product(iH, 1); 

    // R17
    re.push_back(reaction()); 
    re.back().add_educt(iCH2CCH2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC3H2, 1);
    re.back().add_product(iH2, 1); 

    // R18
    re.push_back(reaction()); 
    re.back().add_educt(iCH2CCH2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H2, 1);
    re.back().add_product(i3CH2, 1); 

    // R19
    re.push_back(reaction()); 
    re.back().add_educt(iCH3C2H, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC3H3, 1);
    re.back().add_product(iH, 1); 

    // R20
    re.push_back(reaction()); 
    re.back().add_educt(iCH3C2H, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC3H2, 1);
    re.back().add_product(iH2, 1); 

    // R21
    re.push_back(reaction()); 
    re.back().add_educt(iCH3C2H, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCH3, 1);
    re.back().add_product(iC2H, 1); 

    // R22
    re.push_back(reaction()); 
    re.back().add_educt(iC3H6, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H2, 1);
    re.back().add_product(iCH3, 1);
    re.back().add_product(iH, 1); 

    // R23
    re.push_back(reaction()); 
    re.back().add_educt(iC3H6, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCH2CCH2, 1);
    re.back().add_product(iH2, 1); 

    // R24
    re.push_back(reaction()); 
    re.back().add_educt(iC3H6, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H4, 1);
    re.back().add_product(i3CH2, 1); 

    // R25
    re.push_back(reaction()); 
    re.back().add_educt(iC3H6, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H, 1);
    re.back().add_product(iCH4, 1); 
    re.back().add_product(iH, 1);

    // R26
    re.push_back(reaction()); 
    re.back().add_educt(iC3H8, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC3H6, 1);
    re.back().add_product(iH2, 1); 

    // R27
    re.push_back(reaction()); 
    re.back().add_educt(iC3H8, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H6, 1);
    re.back().add_product(i1CH2, 1); 

    // R28
    re.push_back(reaction()); 
    re.back().add_educt(iC3H8, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H4, 1);
    re.back().add_product(iCH4, 1); 

    // R29
    re.push_back(reaction()); 
    re.back().add_educt(iC3H8, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCH3, 1);
    re.back().add_product(iC2H5, 1); 

    // R30
    re.push_back(reaction()); 
    re.back().add_educt(iCO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCO, 1);
    re.back().add_product(iO, 1); 

    // R31
    re.push_back(reaction()); 
    re.back().add_educt(iCO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCO, 1);
    re.back().add_product(iO1D, 1); 

    // R32
    re.push_back(reaction()); 
    re.back().add_educt(iH2CO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iH2, 1);
    re.back().add_product(iCO, 1); 

    // R33
    re.push_back(reaction()); 
    re.back().add_educt(iH2CO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iHCO, 1);
    re.back().add_product(iH, 1); 

    // R34
    re.push_back(reaction()); 
    re.back().add_educt(iHCO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iH, 1);
    re.back().add_product(iCO, 1); 

    // R35
    re.push_back(reaction()); 
    re.back().add_educt(iCH2CO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(i3CH2, 1);
    re.back().add_product(iCO, 1); 

    // R36
    re.push_back(reaction()); 
    re.back().add_educt(iCH3CHO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCH4, 1);
    re.back().add_product(iCO, 1); 

    // R37
    re.push_back(reaction()); 
    re.back().add_educt(iCH3CHO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iCH4, 1);
    re.back().add_product(iCO, 1); 

    // R38
    re.push_back(reaction()); 
    re.back().add_educt(iC2H5CHO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iC2H5, 1);
    re.back().add_product(iHCO, 1); 

    // R39
    re.push_back(reaction()); 
    re.back().add_educt(iH2O, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iH, 1);
    re.back().add_product(iOH, 1); 

    // R40
    re.push_back(reaction()); 
    re.back().add_educt(iH2O2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iOH, 2); 

    // R41
    re.push_back(reaction()); 
    re.back().add_educt(iHO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iOH, 1);
    re.back().add_product(iO, 1); 

    // R42
    re.push_back(reaction()); 
    re.back().add_educt(iO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iO, 1);
    re.back().add_product(iO1D, 1); 

    // R43
    re.push_back(reaction()); 
    re.back().add_educt(iO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iO, 2); 

    // R44
    re.push_back(reaction()); 
    re.back().add_educt(iO3, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iO2, 1);
    re.back().add_product(iO1D, 1); 

    // R45
    re.push_back(reaction()); 
    re.back().add_educt(iO3, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iO2, 1);
    re.back().add_product(iO, 1); 

    // R46
    re.push_back(reaction()); 
    re.back().add_educt(iNO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iN, 1);
    re.back().add_product(iO, 1); 

    // R47
    re.push_back(reaction()); 
    re.back().add_educt(iNO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iNO, 1);
    re.back().add_product(iO, 1); 

    // R48
    re.push_back(reaction()); 
    re.back().add_educt(iHNO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iNO, 1);
    re.back().add_product(iOH, 1); 

    // R49
    re.push_back(reaction()); 
    re.back().add_educt(iHNO3, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iNO2, 1);
    re.back().add_product(iOH, 1); 

    // R50
    re.push_back(reaction()); 
    re.back().add_educt(iHNO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iNO, 1);
    re.back().add_product(iH, 1); 

    // R51
    re.push_back(reaction()); 
    re.back().add_educt(iNH3, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iNH2, 1);
    re.back().add_product(iH, 1); 
 

    // R52
    re.push_back(reaction()); 
    re.back().add_educt(iN2H4, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iN2H3, 1);
    re.back().add_product(iH, 1); 
 
    // R53
    re.push_back(reaction()); 
    re.back().add_educt(iNH, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iN, 1);
    re.back().add_product(iH, 1); 
 
    // R54
    re.push_back(reaction()); 
    re.back().add_educt(iNH2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iNH, 1);
    re.back().add_product(iH, 1); 
 
    // R55
    re.push_back(reaction()); 
    re.back().add_educt(iNH2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iNH2s, 1); 
 
    // R56
    re.push_back(reaction()); 
    re.back().add_educt(iSO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iSO, 1);
    re.back().add_product(iO, 1); 
 
    // R57
    re.push_back(reaction()); 
    re.back().add_educt(iSO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(i1SO2, 1); 
 
    // R58
    re.push_back(reaction()); 
    re.back().add_educt(iSO2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(i3SO2, 1); 
 
    // R59
    re.push_back(reaction()); 
    re.back().add_educt(iH2S, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iHS, 1);
    re.back().add_product(iH, 1); 
 
    // R60
    re.push_back(reaction()); 
    re.back().add_educt(iH2SO4, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iSO2, 1);
    re.back().add_product(iOH, 2); 
 
    // R61
    re.push_back(reaction()); 
    re.back().add_educt(iHSO, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iHS, 1);
    re.back().add_product(iO, 1); 
 
    // R62
    re.push_back(reaction()); 
    re.back().add_educt(iS2, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iS, 2); 
 
    // R63
    re.push_back(reaction()); 
    re.back().add_educt(iS3, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iS2, 1);
    re.back().add_product(iS, 1); 
 
    // R64
    re.push_back(reaction()); 
    re.back().add_educt(iS4, 1);
    re.back().add_educt(ihv, 1);
    re.back().add_product(iS2, 2); 
 



    // R65
    re.push_back(reaction()); 
    re.back().add_educt(iC, 1);
    re.back().add_educt(iOH, 1);
    re.back().add_product(iCO, 1);
    re.back().add_product(iH, 1); 

    // R66
    re.push_back(reaction()); 
    re.back().add_educt(iC, 1);
    re.back().add_educt(iH2, 1);
    re.back().add_educt(iM, 1);
    re.back().add_product(i3CH2, 1);
    re.back().add_product(iM, 1); 
 
    // R67
    re.push_back(reaction()); 
    re.back().add_educt(iC, 1);
    re.back().add_educt(iO2, 1);
    re.back().add_product(iCO, 1);
    re.back().add_product(iO, 1); 
 
 

}



// TODO Die erstellung des Netzwerks verfollständigen...
void build_model_Wang_1998(std::vector<species> &sp, std::vector<reaction> &re) {
	sp.push_back(species(sp.size(), "hv"));
	sp.push_back(species(sp.size(), "O3"));
	sp.push_back(species(sp.size(), "Os"));
	sp.push_back(species(sp.size(), "O2"));
	sp.push_back(species(sp.size(), "H2O"));
	sp.push_back(species(sp.size(), "OH"));
	sp.push_back(species(sp.size(), "N2"));
	sp.push_back(species(sp.size(), "O"));
	sp.push_back(species(sp.size(), "CO"));
	sp.push_back(species(sp.size(), "CO2"));
	sp.push_back(species(sp.size(), "H"));
	sp.push_back(species(sp.size(), "HO2"));
	sp.push_back(species(sp.size(), "NO"));
	sp.push_back(species(sp.size(), "NO2"));
	sp.push_back(species(sp.size(), "HNO3"));
	sp.push_back(species(sp.size(), "NO3"));
	sp.push_back(species(sp.size(), "N2O5"));
	sp.push_back(species(sp.size(), "H2O2"));
	sp.push_back(species(sp.size(), "CH4"));
	sp.push_back(species(sp.size(), "CH3"));
	sp.push_back(species(sp.size(), "CH3O2"));
	sp.push_back(species(sp.size(), "CH3O"));
	sp.push_back(species(sp.size(), "CH2O"));
	sp.push_back(species(sp.size(), "CH3O2H"));
	sp.push_back(species(sp.size(), "CHO"));
	sp.push_back(species(sp.size(), "N2O"));
	sp.push_back(species(sp.size(), "CF2Cl2"));
	sp.push_back(species(sp.size(), "CFCl3"));
	sp.push_back(species(sp.size(), "SO3"));
	sp.push_back(species(sp.size(), "H2SO4"));
	sp.push_back(species(sp.size(), "SO2"));
	sp.push_back(species(sp.size(), "HOSO2"));
	sp.push_back(species(sp.size(), "H2SO4_aq"));
	sp.push_back(species(sp.size(), "HSO4_m"));
	sp.push_back(species(sp.size(), "H_p"));
	sp.push_back(species(sp.size(), "HNO3_aq"));
	sp.push_back(species(sp.size(), "NO3_m"));
	sp.push_back(species(sp.size(), "CH2O_aq"));
	sp.push_back(species(sp.size(), "SO2_aq"));
	sp.push_back(species(sp.size(), "HSO3_m"));
	sp.push_back(species(sp.size(), "SO3_2m"));
	sp.push_back(species(sp.size(), "H2O2_aq"));
	sp.push_back(species(sp.size(), "OH_aq"));
	sp.push_back(species(sp.size(), "HO2_aq"));
	sp.push_back(species(sp.size(), "O2_m"));
	sp.push_back(species(sp.size(), "products"));
	
	// IDs
	size_t ihv=0;
	size_t iO3=1;
	size_t iOs=2;
	size_t iO2=3;
	size_t iH2O=4;
	size_t iOH=5;
	size_t iN2=6;
	size_t iO=7;
	size_t iCO=8;
	size_t iCO2=9;
	size_t iH=10;
	size_t iHO2=11;
	size_t iNO=12;
	size_t iNO2=13;
	size_t iHNO3=14;
	size_t iNO3=15;
	size_t iN2O5=16;
	size_t iH2O2=17;
	size_t iCH4=18;
	size_t iCH3=19;
	size_t iCH3O2=20;
	size_t iCH3O=21;
	size_t iCH2O=22;
	size_t iCH3O2H=23;
	size_t iCHO=24;
	size_t iN2O=25;
	size_t iCF2Cl2=26;
	size_t iCFCl3=27;
	size_t iSO3=28;
	size_t iH2SO4=29;
	size_t iSO2=30;
	size_t iHOSO2=31;
	size_t iH2SO4_aq=32;
	size_t iHSO4_m=33;
	size_t iH_p=34;
	size_t iHNO3_aq=35;
	size_t iNO3_m=36;
	size_t iCH2O_aq=37;
	size_t iSO2_aq=38;
	size_t iHSO3_m=39;
	size_t iSO3_2m=40;
	size_t iH2O2_aq=41;
	size_t iOH_aq=42;
	size_t iHO2_aq=43;
	size_t iO2_m=44;
	size_t iproducts=45;
	
	
	// R1
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 2);
	re.back().add_product(iOs, 1);
	re.back().add_product(iO2, 1);  
	
	// R2
	re.push_back(reaction()); 
	re.back().add_educt(iOs, 1);
	re.back().add_educt(iH2O, 2);
	re.back().add_product(iOH, 2);

	// R3
	re.push_back(reaction()); 
	re.back().add_educt(iOs, 1);
	re.back().add_educt(iN2, 2);
	re.back().add_product(iO, 1);
	re.back().add_product(iN2, 1);
	
	// R4
	re.push_back(reaction()); 
	re.back().add_educt(iOs, 1);
	re.back().add_educt(iO2, 2);
	re.back().add_product(iO, 1);
	re.back().add_product(iO2, 1); 
	
	// R5
	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iOH, 2);
	re.back().add_product(iH, 1);
	re.back().add_product(iCO2, 1);
	
	// R6
	re.push_back(reaction()); 
	re.back().add_educt(iCO, 1);
	re.back().add_educt(iOH, 2);
	re.back().add_product(iH, 1);
	re.back().add_product(iCO2, 1);  
}




void build_model_Hadac_2011(std::vector<species> &sp, std::vector<reaction> &re) {
	sp.push_back(species(sp.size(), "O"));
	size_t iO=0;
	sp.push_back(species(sp.size(), "O2"));
	size_t iO2=1;
	sp.push_back(species(sp.size(), "O3"));
	size_t iO3=2;
	sp.push_back(species(sp.size(), "H"));
	size_t iH=3;
	sp.push_back(species(sp.size(), "OH"));
	size_t iOH=4;
	sp.push_back(species(sp.size(), "HO2"));
	size_t iHO2=5;
	sp.push_back(species(sp.size(), "H2"));
	size_t iH2=6;
	sp.push_back(species(sp.size(), "H2O2"));
	size_t iH2O2=7;
	sp.push_back(species(sp.size(), "M"));
	size_t iM=8;	
	sp.push_back(species(sp.size(), "H2O"));
	size_t iH2O=9;
	sp.push_back(species(sp.size(), "hv"));
	size_t ihv=10;

	
	// R1
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO3, 1);
	re.back().add_product(iO2, 2);
	re.back().set_k(9.08e-16);
	
	// R2
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iH, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iOH, 1);
	re.back().set_k(1.78e-11);
	
	// R3
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iH, 1);
	re.back().set_k(4.4e-11);
	
	// R4
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iHO2, 1);
	re.back().set_k(2.37e-14);
	
	// R5
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iOH, 1);
	re.back().set_k(3.5e-11);
	
	// R6
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iH2, 1);
	re.back().set_k(5.4e-12);
	
	// R7
	re.push_back(reaction()); 
	re.back().add_educt(iH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iOH, 2);
	re.back().set_k(2.26e-12);
	
	// R8
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 2);
	re.back().add_educt(iO, 1);
	re.back().add_product(iH2O, 1);
	re.back().set_k(1.12e-12);
	
	// R9
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iHO2, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iH2O, 1);
	re.back().set_k(4.0e-11);
	
	// R10
	re.push_back(reaction()); 
	re.back().add_educt(iHO2, 2);
	re.back().add_product(iO2, 1);
	re.back().add_product(iH2O2, 1);
	re.back().set_k(2.5e-12);
	
	// R11
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iO3, 1);
	re.back().add_product(iM, 1);
	re.back().set_k(1.05e-34);
	
	// R12
	re.push_back(reaction()); 
	re.back().add_educt(iO, 2);
	re.back().add_educt(iM, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iM, 1);
	re.back().set_k(8.23e-33);
	
	// R13
	re.push_back(reaction()); 
	re.back().add_educt(iO, 1);
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iM, 1);
	re.back().set_k(1.4e-31);

	// R14
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(iH, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iHO2, 1);
	re.back().add_product(iM, 1);
	re.back().set_k(8.08e-32);

	// R15
	re.push_back(reaction()); 
	re.back().add_educt(iOH, 1);
	re.back().add_educt(iH, 1);
	re.back().add_educt(iM, 1);
	re.back().add_product(iH2O, 1);
	re.back().add_product(iM, 1);
	re.back().set_k(3.31e-27);
	
	// R16
	re.push_back(reaction()); 
	re.back().add_educt(iO2, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 2);
	re.back().set_k(1.26e-8);
	
	// R17
	re.push_back(reaction()); 
	re.back().add_educt(iH2O, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO, 1);
	re.back().add_product(iOH, 1);
	re.back().set_k(3.4e-6);

	// R18
	re.push_back(reaction()); 
	re.back().add_educt(iO3, 1);
	re.back().add_educt(ihv, 1);
	re.back().add_product(iO2, 1);
	re.back().add_product(iO, 1);
	re.back().set_k(7.1e-8);	
	
	
}



#endif // BUILD_NETWORKS_MODELS_H

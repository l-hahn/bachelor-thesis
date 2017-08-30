#include "alphabet.hpp"

namespace alphabet{
	std::map<unsigned char, unsigned char> ProteinBit = {{'A',6},{'C',7}
		,{'D',8},{'E',9},{'F',10},{'G',11},{'H',12},{'I',13},{'K',14},{'L',15}
		,{'M',16},{'N',17},{'P',18},{'Q',19},{'R',20},{'S',21},{'T',22},{'V',23}
		,{'W',24},{'X',25},{'Y',26},{'*',27}};
	std::map<unsigned char, unsigned char> BitProtein = {{6,'A'},{7,'C'},{8,'D'}
		,{9,'E'},{10,'F'},{11,'G'},{12,'H'},{13,'I'},{14,'K'},{15,'L'},{16,'M'}
		,{17,'N'},{18,'P'},{19,'Q'},{20,'R'},{21,'S'},{22,'T'},{23,'V'},{24,'W'}
		,{25,'X'},{26,'Y'},{27,'*'}};

	std::map<unsigned char, unsigned char> AlphaReduct = {{'A','S'},{'C','C'}
		,{'D','K'},{'E','K'},{'F','F'},{'G','G'},{'H','H'},{'I','I'},{'K','K'}
		,{'L','I'},{'M','M'},{'N','K'},{'P','P'},{'Q','K'},{'R','K'},{'S','S'}
		,{'T','S'},{'V','I'},{'W','W'},{'X','X'},{'Y','Y'},{'*','*'}};
	std::map<unsigned char, unsigned char> AlphaBitReduct = {{'A',21},{'C',7}
		,{'D',14},{'E',14},{'F',10},{'G',11},{'H',12},{'I',13},{'K',14},{'L',13}
		,{'M',16},{'N',14},{'P',18},{'Q',14},{'R',14},{'S',21},{'T',21},{'V',13}
		,{'W',24},{'X',25},{'Y',26},{'*',27}};
	std::map<unsigned char, unsigned char> BitReduct = {{6,21},{7,7},{8,14},{9,14}
		,{10,10},{11,11},{12,12},{13,13},{14,14},{15,13},{16,16},{17,14},{18,18}
		,{19,14},{20,14},{21,21},{22,21},{23,13},{24,24},{25,25},{26,26},{27,27}};

	std::map<std::string, unsigned char> CodonProtein = {{"AAA",'K'},{"AAC",'N'}
		,{"AAG",'K'},{"AAT",'N'},{"ACA",'T'},{"ACC",'T'},{"ACG",'T'},{"ACT",'T'}
		,{"AGA",'R'},{"AGC",'S'},{"AGG",'R'},{"AGT",'S'},{"ATA",'I'},{"ATC",'I'}
		,{"ATG",'M'},{"ATT",'I'},{"CAA",'Q'},{"CAC",'H'},{"CAG",'Q'},{"CAT",'H'}
		,{"CCA",'P'},{"CCC",'P'},{"CCG",'P'},{"CCT",'P'},{"CGA",'R'},{"CGC",'R'}
		,{"CGG",'R'},{"CGT",'R'},{"CTA",'L'},{"CTC",'L'},{"CTG",'L'},{"CTT",'L'}
		,{"GAA",'E'},{"GAC",'D'},{"GAG",'E'},{"GAT",'D'},{"GCA",'A'},{"GCC",'A'}
		,{"GCG",'A'},{"GCT",'A'},{"GGA",'G'},{"GGC",'G'},{"GGG",'G'},{"GGT",'G'}
		,{"GTA",'V'},{"GTC",'V'},{"GTG",'V'},{"GTT",'V'},{"TAA",'*'},{"TAC",'Y'}
		,{"TAG",'*'},{"TAT",'Y'},{"TCA",'S'},{"TCC",'S'},{"TCG",'S'},{"TCT",'S'}
		,{"TGA",'*'},{"TGC",'C'},{"TGG",'W'},{"TGT",'C'},{"TTA",'L'},{"TTC",'F'}
		,{"TTG",'L'},{"TTT",'F'}};
	std::map<unsigned short,unsigned char> CodonProteinBit = {{73,14},{74,17}
		,{75,14},{76,17},{81,22},{82,22},{83,22},{84,22},{89,20},{90,21},{91,20}
		,{92,21},{97,13},{98,13},{99,16},{100,13},{137,19},{138,12},{139,19}
		,{140,12},{145,18},{146,18},{147,18},{148,18},{153,20},{154,20},{155,20}
		,{156,20},{161,15},{162,15},{163,15},{164,15},{201,9},{202,8},{203,9}
		,{204,8},{209,6},{210,6},{211,6},{212,6},{217,11},{218,11},{219,11}
		,{220,11},{225,23},{226,23},{227,23},{228,23},{265,27},{266,26},{267,27}
		,{268,26},{273,21},{274,21},{275,21},{276,21},{281,27},{282,7},{283,24}
		,{284,7},{289,15},{290,10},{291,15},{292,10}};

	std::vector<unsigned char> ProteinAlphabet = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','*'};
	std::vector<unsigned char> DnaAlphabet = {'A','C','G','T','N'};
	std::vector<unsigned char> ProteinUniqAlphabet = {'D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','X','Y','*'};
	std::vector<unsigned char> ProteinBitAlphabet = {6,7,8,9,11,12,13,14,16,17,18,19,20,22,23,24,25,27};
	std::vector<unsigned char> DnaBitAlphabet = {1,2,3,4,5};

	unsigned DnaBitLength = 3;
	unsigned ProteinBitLength = 5;

	std::map<unsigned char, unsigned char> BitDna ={{1,'A'},{2,'C'},{3,'G'},{4,'T'},{5,'N'}};
	std::map<unsigned char, unsigned char> DnaBit ={{'A',1},{'C',2},{'G',3},{'T',4},{'N',5}};
	unsigned char DnaInvert(unsigned char Letter){
		if(Letter < 5){
			switch(Letter){
				case 1:
					Letter = 4;
					break;
				case 2:
					Letter = 3;
					break;
				case 3:
					Letter = 2;
					break;
				case 4:
					Letter = 1;
					break;
			}
		}
		else{
			switch(Letter){
				case 'A':
					Letter = 'T';
					break;
				case 'C':
					Letter = 'G';
					break;
				case 'G':
					Letter = 'C';
					break;
				case 'T':
					Letter = 'A';
					break;
			}
		}
		return Letter;
	}

	std::map<unsigned char, unsigned char> & AlphaMap(int BitLeng){
		if(BitLeng == 5){
			return ProteinBit;
		}
		else{
			return DnaBit;
		}
	}
	std::map<unsigned char, unsigned char> & BitMap(int BitLeng){
		if(BitLeng == 5){
			return BitProtein;
		}
		else{
			return BitDna;
		}
	}
};

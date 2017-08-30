#include "dboptions.hpp"

namespace dboptions{
	std::string InSeqFam = "";
	std::string OutSwDb = "Database";
	std::string PatternFile = "";
	double SeqCovThreshold = 1;
	unsigned WordPerSeq = 1;
	unsigned PatternNumber = 1;
	unsigned PatternWeight = 6;
	unsigned PatternMaxDC = 6;
	unsigned PatternMinDC = 6;
	unsigned FamilySize = 0;
	unsigned FamilyOffset = 10;
	unsigned ThreadNumber = 1;
	unsigned BitLength = 5;
	bool AlphaReduction = false;
	bool Greedy = true;
	bool Overlap = false;
	bool SeqCov = false;
	bool SplitFile = false;

	void parse_length(const char* Str){
		std::string Length = std::string(Str);
		unsigned Pos = (unsigned) Length.size();
		unsigned Leng = Pos;

		for(unsigned i = 0; i < Leng; i++){
			if(std::isdigit(Length[i]) == false){
				Pos = i;	
			}
		}

		if(Pos < Leng){
			char *Tmp1, *Tmp2;
			Tmp1 = new char[Pos];
			Tmp2 = new char[Leng-Pos-1];
			for(uint32_t i = 0; i < Leng; i++){
				if(i < Pos){
					Tmp1[i] = Length[i];
				}
				if(i > Pos){
					Tmp2[i-Pos-1] = Length[i];
				}
			}
			PatternMinDC = std::atoi(Tmp1);
			PatternMaxDC = std::atoi(Tmp2);
			delete Tmp1;
			delete Tmp2;
			
		}
		else{
			PatternMinDC = std::atoi(Str);
			PatternMaxDC = std::atoi(Str);
		}
	}
};
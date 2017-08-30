#ifndef DBOPTIONS_HPP_
#define DBOPTIONS_HPP_

#include <cstdlib>
#include <string>

namespace dboptions{
	extern std::string InSeqFam;
	extern std::string OutSwDb;
	extern std::string PatternFile;
	extern double SeqCovThreshold;
	extern unsigned WordPerSeq;
	extern unsigned PatternNumber;
	extern unsigned PatternWeight;
	extern unsigned PatternMaxDC;
	extern unsigned PatternMinDC;
	extern unsigned FamilySize;
	extern unsigned FamilyOffset;
	extern unsigned ThreadNumber;
	extern unsigned BitLength;
	extern bool AlphaReduction;
	extern bool Greedy;
	extern bool Overlap;
	extern bool SeqCov;
	extern bool SplitFile;

	void parse_length(const char* Length);
};

#endif
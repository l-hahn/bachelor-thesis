#ifndef DSOPTIONS_HPP_
#define DSOPTIONS_HPP_

#include <string>

namespace dsoptions{
	extern std::string InSeq;
	extern std::string InSwDb;
	extern std::string OutDetect;
	extern double HitThreshold;
	extern signed BlockSize;
	extern signed WordPerSeq;
	extern unsigned FamilySize;
	extern unsigned ThreadNumber;
	extern unsigned BitLength;
	extern bool AlphaReduction;
	extern bool Translated;
};

#endif
#ifndef SPACEDFAMILY_HPP_
#define SPACEDFAMILY_HPP_

#include <cstdint>
#include <cstdlib>
#include <vector>
#include "familyscore.hpp"

class spacedword_family{
	public:
		spacedword_family(int64_t BitSpacedWord);
		spacedword_family(int64_t BitSpacedWord, unsigned Family, unsigned Position, double FamilyScore);
		spacedword_family(int64_t BitSpacedWord, family_score & FamilyScore);

		int64_t bits() const;
		unsigned size() const;
		void push_back(unsigned Family, unsigned Position, double Score);
		void push_back(family_score & FamilyScore);

		family_score operator[](size_t Idx);
		bool operator==(const spacedword_family & Spw) const;
		bool operator<(const spacedword_family & Spw) const;
		bool operator>(const spacedword_family & Spw) const;

		typedef std::vector<family_score>::iterator iterator;
		iterator begin();
		iterator end();

	private:
		int64_t _BitSpacedWord;
		std::vector<family_score> _Family;
};
#endif
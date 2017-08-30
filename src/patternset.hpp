#ifndef PATTERNSET_HPP_
#define PATTERNSET_HPP_

#include <algorithm>
#include <limits>
#include <iostream>
#include <vector>
#include "pattern.hpp"

class patternset{
	public:
		patternset();
		patternset(unsigned Numb, unsigned MaxW, unsigned MinW, unsigned MaxDc, unsigned MinDc);

		void push_back(pattern & Pat);
		void push_back(pattern && Pat);
		void push_back(std::string & Pat);
		void push_back(std::string && Pat);
		void push_back(std::vector<char> & Pat);
		void push_back(std::vector<unsigned char> & Pat);
		void random(unsigned Size, unsigned Length, unsigned Weight);
		void random(unsigned Size, unsigned MinLength, unsigned MaxLength, unsigned Weight);
		void random(unsigned Size, unsigned MinLength, unsigned MaxLength, unsigned MinWeight, unsigned MaxWeight);
		void sort();
		bool is_uniq(const pattern & Pat) const;


		pattern & operator[](size_t Idx);
		bool operator<(const patternset & P) const;
		bool operator>(const patternset & P) const;

		unsigned max_weight() const;
		unsigned min_weight() const;
		unsigned weight() const;
		unsigned max_dontcare() const;
		unsigned min_dontcare() const;
		unsigned dontcare() const;
		unsigned max_length() const;
		unsigned min_length() const;
		unsigned length() const;
		unsigned size() const;

		typedef std::vector<pattern>::iterator iterator;
		iterator begin();
		iterator end();

	private:
		std::vector<pattern> _PatternSet;
		double _Score;
		unsigned _MaxWeight;
		unsigned _MinWeight;
		unsigned _MaxDontCare;
		unsigned _MinDontCare;
};
#endif
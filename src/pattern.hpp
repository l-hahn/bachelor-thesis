#ifndef PATTERN_HPP_
#define PATTERN_HPP_

#include <ctime>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

class pattern{
	public:
		pattern();
		pattern(unsigned DontCare, unsigned Weight);
		pattern(unsigned DontCare, unsigned Weight, uint64_t Seed);
		pattern(std::string & StrPat);
		pattern(std::string && StrPat);
		pattern(std::vector<char> & StrPat);
		pattern(std::vector<unsigned char> & StrPat);

		std::string to_string();
		bool is_match(unsigned Pos) const;
		void bit_swap(unsigned PosA, unsigned PosB);
		void random_swap();
		void random_swap(uint64_t Seed);

		void random();
		void random(unsigned DontCare, unsigned Weight);
		void random(unsigned DontCare, unsigned Weight, uint64_t Seed);

		unsigned char operator[](size_t Idx) const;
		bool operator<(const pattern & P) const;
		bool operator>(const pattern & P) const;
		bool operator==(const pattern & P) const;

		void set_score(double Scr);
		double score() const;
		unsigned weight() const;
		unsigned length() const;
		unsigned dontcare() const;
		unsigned get_overlap(const pattern & P, int Shift) const;

		typedef std::vector<unsigned char>::iterator iterator;
		iterator begin();
		iterator end();


	private:
		void _check_pattern();

		std::vector<unsigned char> _VectorPattern;
		std::vector<unsigned char> _MatchPos;
		std::vector<unsigned char> _DCPos;
		double _Score;
		uint64_t _BitPattern;
		bool _IsBit;
};
#endif
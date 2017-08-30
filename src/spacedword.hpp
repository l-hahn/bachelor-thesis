#ifndef SPACEDWORD_HPP_
#define SPACEDWORD_HPP_

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include "alphabet.hpp"
#include "pattern.hpp"

template<int _BitLength>
class spacedword{
	public:
		spacedword();
		spacedword(std::vector<char> & SeqVec, pattern & Pat, int64_t SeqPos, signed SeqNo = -1);
		spacedword(std::vector<unsigned char> & SeqVec, pattern & Pat, int64_t SeqPos, signed SeqNo = -1);
		spacedword(int64_t Spw, int64_t SeqPos = -1, signed SeqNo = -1);
		
		std::string to_string() const;
		void push_back(unsigned char Letter);

		int64_t position() const;
		int64_t bits() const;
		int64_t bits_comp() const;
		unsigned size() const;
		signed sequence() const;
		signed counter() const;
		void set_position(int64_t Pos);
		void set_sequence(signed No);
		void set_counter(signed Ctr);

		unsigned char operator[](unsigned Idx) const;
		
		bool operator==(const spacedword & Spw) const;
		bool operator!=(const spacedword & Spw) const;
		bool operator<(const spacedword & Spw) const;
		bool operator>(const spacedword & Spw) const;
		
	private:
		template<typename T>
		void create_word(T & Vec, pattern & Pat);


		int64_t _BitSpacedWord;
		int64_t _SequencePos;
		signed _SequenceNo;
		signed _SequenceCtr;
};
typedef spacedword<5> protein_spacedword;
typedef spacedword<3> dna_spacedword;
#include "spacedword.cpp"
#endif

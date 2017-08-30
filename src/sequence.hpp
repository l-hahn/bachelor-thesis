#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include <string>
#include <vector>
#include "alphabet.hpp"
#include "pattern.hpp"
#include "spacedword.hpp"

class sequence{
	public:
		sequence();
		sequence(std::string & Seq, std::string SeqName = "");
		sequence(std::string && Seq, std::string SeqName = "");
		sequence(std::vector<char> & Seq, std::string SeqName = "");
		sequence(std::vector<char> && Seq, std::string SeqName = "");
		sequence(std::vector<unsigned char> & Seq, std::string SeqName = "");
		sequence(std::vector<unsigned char> && Seq, std::string SeqName = "");

		void push_back(unsigned char Letter);
		void to_bit(bool IsProtein = true);
		std::string to_string() const;
		bool is_bit() const;
		bool is_protein() const;
		bool is_dna() const;

		template<typename T>
		void spaced_words(std::vector<T> && Spw, pattern & Pat, signed SeqNo = -1);
		template<typename T>
		void spaced_words(std::vector<T> & Spw, pattern & Pat, signed SeqNo = -1);

		std::string name() const;
		size_t size() const;
		void set_name(std::string & Name);
		void set_name(std::string && Name);

		// template<typename T>
		// void set_sequence(T && Seq);
		template<typename T>
		void set_sequence(T & Seq);

		unsigned char operator[](size_t Idx) const;
		bool operator<(const sequence & Seq) const;
		bool operator>(const sequence & Seq) const;
		typedef std::vector<unsigned char>::iterator iterator;
		iterator begin();
		iterator end();

	private:

		std::vector<unsigned char> _Sequence;
		std::string _SequenceName;
};

// template<typename T>
// void sequence::set_sequence(T && Seq){
// 	set_sequence(Seq);	
// }
template<typename T>
void sequence::set_sequence(T & Seq){
	for(auto & Letter : Seq){
		if(std::isalpha(Letter)){
			_Sequence.push_back(std::toupper(Letter));
		}
		else if(Letter <= alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1]){
			_Sequence.push_back(Letter);
		}
		else if(Letter == '*'){
			_Sequence.push_back(Letter);
		}
	}	
}


template<typename T>
void sequence::spaced_words(std::vector<T> && Spw, pattern & Pat, signed SeqNo){
	spaced_words(Spw, Pat, SeqNo);
}
template<typename T>
void sequence::spaced_words(std::vector<T> & Spw, pattern & Pat, signed SeqNo){
	for(signed i = 0; i < (signed)_Sequence.size() - (signed)Pat.length() + 1; i++){
			Spw.push_back(T(_Sequence, Pat, i, SeqNo));
	}
}

#endif
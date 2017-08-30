#ifndef SPACEDWORDDB_HPP_
#define SPACEDWORDDB_HPP_

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <string>
#include <sstream>
#include <sys/stat.h> 
#include <vector>
#include "alphabet.hpp"
#include "spacedfamily.hpp"
#include "pattern.hpp"

class spacedword_db{
	public:
		spacedword_db(int64_t BitLength, unsigned WordPerSeq = 1);
		spacedword_db(std::string & FileName);
		spacedword_db(std::vector<spacedword_family> & SpwFam, pattern & Pat, int64_t BitLength, unsigned WordPerSec = 1);

		void push_back(std::vector<spacedword_family> & SpwFam, pattern & Pat);
		void to_file(std::string & FileName, bool SplitData = false);
		void from_file(std::string & FileName);
		void merge(const spacedword_db & DBAddition);

		void add_family_id(std::string && FamID, unsigned ID);
		void add_family_id(std::string & FamID, unsigned ID);
		std::string id_fam_name(unsigned ID) const;
		signed fam_name_id(std::string && Str) const;
		signed fam_name_id(std::string & Str) const;

		std::vector<spacedword_family> spaced_families(size_t Idx) const;
		pattern pattern_families(size_t Idx) const;
		std::string to_string(int64_t BitWord) const;
		int64_t to_bit(std::string && StrWord) const;
		int64_t to_bit(std::string & StrWord) const;
		size_t size() const;
		size_t families_size() const;
		unsigned bit_length() const;
		bool is_protein() const;
		unsigned words_per_seq() const;

	private:
		std::string file_part_name(std::string & FileName, int Part);
		std::vector<std::string> get_name_id() const;
		unsigned new_id();

		std::vector< std::vector<spacedword_family> > _SpacedWordsList;
		std::vector<pattern> _Patterns;
		std::vector<std::string> _FamilyName;
		int64_t _BitLength;
		unsigned _WordPerSeq;
};
#endif
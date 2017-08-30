#ifndef SEQUENCEFAMILY_HPP_
#define SEQUENCEFAMILY_HPP_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <vector>
#include <utility>
#include "alphabet.hpp"
#include "pattern.hpp"
#include "sequence.hpp"
#include "spacedword.hpp"

enum FileType {FALSETYPE, STOCKHOLM, FASTA};
void parse_fasta_header(std::string & Line, std::string & SeqName, std::string & FamID, std::string & FamName);
FileType file_format(std::string & SeqFile);
unsigned info_start(std::string & Str);
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

class sequence_family{
	public:
		sequence_family();
		sequence_family(std::string & SeqFile);
		sequence_family(std::string && SeqFile);
		sequence_family(std::vector<sequence> & SeqVec);
		sequence_family(std::string & SeqFile, std::string & ID, std::string & Name, std::string & Description);
		sequence_family(std::string && SeqFile, std::string & ID, std::string & Name, std::string & Description);
		sequence_family(std::vector<sequence> & SeqVec, std::string & ID, std::string & Name, std::string & Description);

		void push_back(sequence & Seq);
		void translate();
		void to_bit();
		void sort();

		std::string id() const;
		std::string name() const;
		std::string description() const;
		bool is_translated() const;
		bool is_protein() const;
		bool is_bit() const;
		size_t size(size_t MinSize = 0)const;
		void set_id(std::string && ID);
		void set_id(std::string & ID);
		void set_name(std::string && Name);
		void set_name(std::string & Name);
		void set_description(std::string && Description);
		void set_description(std::string & Description);
		void print() const;

		sequence & operator[](size_t Idx);

		typedef std::vector<sequence>::iterator iterator;
		iterator begin();
		iterator end();
		void erase(iterator Idx);

		static std::vector<sequence_family> families_file(std::string & SeqFile);

	private:
		std::vector<sequence> translate_sequence(sequence & Seq, signed SeqPos = -1);
		void read_family(std::string & SeqFile);
		static bool read_file_part(std::ifstream & Input,std::vector<sequence> & Sequences,
			std::string & ID, std::string & Name, std::string & Description, FileType SeqType);
		static bool read_fasta(std::ifstream & Input, std::vector<sequence> & Sequences,
			std::string & FamID, std::string & FamName, std::string & Description, FileType SeqType);
		static bool read_stockholm(std::ifstream & Input, std::vector<sequence> & Sequences,
			std::string & FamID, std::string & FamName, std::string & Description, FileType SeqType);

		std::vector<sequence> _Sequences;
		std::string _ID;
		std::string _Name;
		std::string _Description;
		bool _IsTranslated;
		bool _IsProtein;
		bool _IsBit;
};
#endif
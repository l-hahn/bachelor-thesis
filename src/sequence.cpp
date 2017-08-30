#include "sequence.hpp"

sequence::sequence(){
}

sequence::sequence(std::string & Seq, std::string SeqName): _SequenceName(SeqName){
	set_sequence(Seq);
}

sequence::sequence(std::string && Seq, std::string SeqName): _SequenceName(SeqName){
	set_sequence(Seq);
}

sequence::sequence(std::vector<char> & Seq, std::string SeqName): _SequenceName(SeqName){
	set_sequence(Seq);
}

sequence::sequence(std::vector<char> && Seq, std::string SeqName): _SequenceName(SeqName){
	set_sequence(Seq);
}

sequence::sequence(std::vector<unsigned char> & Seq, std::string SeqName): _SequenceName(SeqName){
	set_sequence(Seq);
}

sequence::sequence(std::vector<unsigned char> && Seq, std::string SeqName): _SequenceName(SeqName){
	set_sequence(Seq);
}

void sequence::push_back(unsigned char Letter){
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

std::string sequence::to_string() const{
	if(is_bit()){
		std::string StringSeq;
		if(is_protein()){
			for(auto & Letter : _Sequence){
				StringSeq.push_back(alphabet::BitProtein[Letter]);
			}
		}
		else{
			for(auto & Letter : _Sequence){
				StringSeq.push_back(alphabet::BitDna[Letter]);
			}
		}
		return StringSeq;
	}
	return std::string(_Sequence.begin(), _Sequence.end());
}

void sequence::to_bit(bool IsProtein){
	if(!IsProtein){
		for(auto & Letter : _Sequence){
			Letter = alphabet::DnaBit[Letter];
		}
		return;
	}
	for(auto & Letter : _Sequence){
		Letter = alphabet::ProteinBit[Letter];
	}
}

bool sequence::is_bit() const{
	if(_Sequence[0] <= alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1]){
		return true;
	}
	return false;
}

bool sequence::is_protein() const{
	if(is_bit()){
		if(_Sequence[0] > alphabet::DnaBitAlphabet[alphabet::DnaBitAlphabet.size()-1]){
			return true;
		}
		return false;
	}
	else{
		auto FoundProtein  = std::find_first_of(alphabet::ProteinUniqAlphabet.begin(),
			alphabet::ProteinUniqAlphabet.end(),_Sequence.begin(),_Sequence.end());
		if(FoundProtein == alphabet::ProteinUniqAlphabet.end()){
			return false;
		}
		return true;
	}
}

bool sequence::is_dna() const{
	return is_protein() == false;
}

std::string sequence::name() const{
	return _SequenceName;
}


unsigned char sequence::operator[](size_t Idx) const{
	return _Sequence[Idx];
}

bool sequence::operator<(const sequence & Seq) const{
	return _Sequence.size() < Seq._Sequence.size();
}

bool sequence::operator>(const sequence & Seq) const{
	return _Sequence.size() < Seq._Sequence.size();
}

sequence::iterator sequence::begin(){
	return _Sequence.begin();
}

sequence::iterator sequence::end(){
	return _Sequence.end();
}


size_t sequence::size() const{
	return _Sequence.size();
}

void sequence::set_name(std::string & Name){
	_SequenceName = Name;
}

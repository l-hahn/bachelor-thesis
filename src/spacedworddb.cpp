#include "spacedworddb.hpp"
spacedword_db::spacedword_db(int64_t BitLength, unsigned WordPerSeq):_BitLength(BitLength), _WordPerSeq(WordPerSeq){
}

spacedword_db::spacedword_db(std::string & FileName){
	from_file(FileName);
}

spacedword_db::spacedword_db(std::vector<spacedword_family> & SpwFam, pattern & Pat, int64_t BitLength, unsigned WordPerSeq):_BitLength(BitLength), _WordPerSeq(WordPerSeq){
	push_back(SpwFam, Pat);
}

void spacedword_db::push_back(std::vector<spacedword_family> & SpwFam, pattern & Pat){
	_SpacedWordsList.push_back(SpwFam);
	_Patterns.push_back(Pat);
	std::sort(_SpacedWordsList[_SpacedWordsList.size()-1].begin(), _SpacedWordsList[_SpacedWordsList.size()-1].end());
}

void spacedword_db::to_file(std::string & FileName, bool SplitData){
	std::string Type = "Protein";
	if(_BitLength != 5){
		Type = "DNA";
	}
	if(SplitData){
		std::string FolderName = FileName+"-DataBase";
		mkdir(FolderName.c_str(),0744);
		FolderName += "/" + FileName;
		for(size_t i = 0; i < size(); i++){
			std::string FilePart = file_part_name(FolderName, i+1);
			std::ofstream Output(FilePart);
			Output << "#<" << Type << " Spaced Word Data Base WPS " << _WordPerSeq << ">#" << std::endl;
			for(auto & SpwFam : _SpacedWordsList[i]){
				Output << "<" << to_string(SpwFam.bits()) << ">" << std::endl;
				for(auto & Fam : SpwFam){
					Output << id_fam_name(Fam.id()) << "\t" << Fam.score() << "\t" << Fam.position() << std::endl;
				}
			}
			Output.close();
		}
	}
	else{
		std::ofstream Output(FileName+".swdb");
		Output << "#<" << Type << " Spaced Word Data Base WPS " << _WordPerSeq << ">#" << std::endl;
		for(size_t i = 0; i < size(); i++){
			Output << "#" << _Patterns[i].to_string() << "#" << std::endl;
			for(auto & SpwFam : _SpacedWordsList[i]){
				Output << "<" << to_string(SpwFam.bits()) << ">" << std::endl;
				for(auto & Fam : SpwFam){
					Output << id_fam_name(Fam.id()) << "\t" << Fam.score() << "\t" << Fam.position() << std::endl;
				}
			}
		}
		Output.close();
	}
}

void spacedword_db::from_file(std::string & FileName){
	std::ifstream Input(FileName);
	std::string Parser;
	std::getline(Input, Parser);
	if(Parser.substr(0,2) == "#<"){
		if(Parser[2] == 'P'){
			_BitLength = alphabet::ProteinBitLength;
		}
		else if(Parser[2] == 'D'){
			_BitLength = alphabet::DnaBitLength;
		}
		else{
			std::cerr << "FATAL ERROR!\nCannot parse Sequence Type!" <<std::endl;
			std::cerr <<"Header-Format: #<[DNA|PROTEIN] ... [WPS:<INT>]>#" << std::endl;
			std::exit(-1);
		}
		std::string WPS;
		for(unsigned Pos = Parser.size()-3; Pos >= 0; Pos--){
			if(std::isdigit(Parser[Pos])){
				WPS.push_back(Parser[Pos]);
			}
			else{
				break;
			}
		}
		_WordPerSeq = std::stoi(std::string(WPS.rbegin(), WPS.rend()));
		if(_WordPerSeq == 0){
			std::cerr << "FATAL ERROR!\nCannot parse Words per Sequence(WPS) parameter!" <<std::endl;
			std::cerr <<"Header-Format: #<[DNA|PROTEIN] ... [WPS:<INT>]>#" << std::endl;
			std::exit(-1);
		}
		std::vector<spacedword_family> TmpSpwFam;
		while(!Input.eof()){
			Input >> Parser;
			if(Parser.size() != 0){
				switch(Parser[0]){
					case '<':
						Parser = Parser.substr(1,Parser.size()-2);
						TmpSpwFam.push_back(spacedword_family(to_bit(Parser)));
						break;
					case '#':
						if(TmpSpwFam.size() != 0){
							_SpacedWordsList.push_back(TmpSpwFam);
							TmpSpwFam.clear();
						}
						Parser = Parser.substr(1,Parser.size()-2);
						_Patterns.push_back(pattern(Parser));
						break;
					default:
						signed ID = fam_name_id(Parser);
						if(ID < 0){
							ID = new_id();
							add_family_id(Parser, ID);
						}
						double Score;
						unsigned Position;
						Input >> Score >> Position;
						TmpSpwFam[TmpSpwFam.size()-1].push_back(ID, Position, Score);
				}
				Parser.clear();
			}
		}
		if(TmpSpwFam.size() != 0){
			_SpacedWordsList.push_back(TmpSpwFam);
			TmpSpwFam.clear();
		}
	}
	else{
		std::cerr << "FATAL ERROR!\nCannot parse File!" <<std::endl;
		std::cerr << "Header-Format: #<[DNA|PROTEIN] ... [WPS:<INT>]>#" << std::endl;
		std::cerr << "Data-Format:\n#[Pattern]#\n<[SpacedWord]>\n[FamID]\t[Score]\t[Order Number]\n..." << std::endl;
		std::exit(-1);
	}
	Input.close();
}

void spacedword_db::merge(const spacedword_db & DBAddition){
	std::vector<std::string> NameVec = DBAddition.get_name_id();
	if(_FamilyName.size() < NameVec.size()){
		_FamilyName.resize(NameVec.size(),"");
	}
	for(unsigned i = 0; i < NameVec.size(); i++){
		if(NameVec[i].size() != 0){
			_FamilyName[i] = NameVec[i];
		}
	}

	for(size_t i = 0; i < DBAddition.size(); i++){
		unsigned PatPos;
		auto PatIt = std::find(_Patterns.begin(), _Patterns.end(), DBAddition.pattern_families(i));
		if(PatIt != _Patterns.end()){
			PatPos = std::distance(_Patterns.begin(), PatIt);
		}
		else{
			_Patterns.push_back(DBAddition.pattern_families(i));
			_SpacedWordsList.resize(_SpacedWordsList.size()+1);
			PatPos = _Patterns.size()-1;
		}
		std::vector<spacedword_family> DBEntry = DBAddition.spaced_families(i);
		for(auto & SwFam : DBEntry){
			auto SwIt = std::find(_SpacedWordsList[PatPos].begin(), _SpacedWordsList[PatPos].end(), SwFam);
			if(SwIt != _SpacedWordsList[PatPos].end()){
				for(auto & FamScore : SwFam){
					SwIt->push_back(FamScore);
				}
			}
			else{
				_SpacedWordsList[PatPos].push_back(SwFam);
			}
		}
		std::sort(_SpacedWordsList[PatPos].begin(),_SpacedWordsList[PatPos].end());
	}
}


void spacedword_db::add_family_id(std::string && FamID, unsigned ID){
	add_family_id(FamID, ID);
}
void spacedword_db::add_family_id(std::string & FamID, unsigned ID){
	if(ID >= _FamilyName.size()){
		_FamilyName.resize(ID+1,"");
	}
	_FamilyName[ID] = FamID;
}

std::vector<std::string> spacedword_db::get_name_id() const{
	return _FamilyName;
}

unsigned spacedword_db::new_id(){
	_FamilyName.push_back("");
	return _FamilyName.size()-1;
}

std::string spacedword_db::id_fam_name(unsigned ID) const{
	return _FamilyName[ID];
}


signed spacedword_db::fam_name_id(std::string && Str) const{
	return fam_name_id(Str);
}
signed spacedword_db::fam_name_id(std::string & Str) const{
	auto Pos = std::find(_FamilyName.begin(), _FamilyName.end(), Str);
	if(Pos != _FamilyName.end()){
		return std::distance(_FamilyName.begin(), Pos);
	}
	return -1;
}

std::string spacedword_db::to_string(int64_t BitWord) const{
	std::map<unsigned char, unsigned char> & BitMap = alphabet::BitMap(_BitLength);
	std::string AlphaWord = "";
	char BitMask = (1 << _BitLength)-1; 
	while(BitWord > 0 && _BitLength > 0){
		AlphaWord.push_back(BitMap[BitWord & BitMask]);
		BitWord >>= _BitLength;
	}
	std::reverse(AlphaWord.begin(),AlphaWord.end());
	return AlphaWord;
}

int64_t spacedword_db::to_bit(std::string && StrWord) const{
	return to_bit(StrWord);
}
int64_t spacedword_db::to_bit(std::string & StrWord) const{
	std::map<unsigned char, unsigned char> & AlphaMap = alphabet::AlphaMap(_BitLength);
	int64_t BitWord = 0;
	for(auto & Letter : StrWord){
		BitWord <<= _BitLength;
		if(Letter <= alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1]){
			BitWord |= Letter;
		}
		else{
			BitWord |= AlphaMap[Letter];
		}
	}
	return BitWord;
}


std::vector<spacedword_family> spacedword_db::spaced_families(size_t Idx) const{
	return _SpacedWordsList[Idx];
}

pattern spacedword_db::pattern_families(size_t Idx) const{
	return _Patterns[Idx];
}

size_t spacedword_db::size() const{
	return _Patterns.size();
}

size_t spacedword_db::families_size() const{
	return _FamilyName.size();
}

unsigned spacedword_db::bit_length() const{
	return _BitLength;
}

unsigned spacedword_db::words_per_seq() const{
	return _WordPerSeq;
}

bool spacedword_db::is_protein() const{
	return _BitLength == alphabet::ProteinBitLength;
}

std::string spacedword_db::file_part_name(std::string & FileName, int Part){
	std::ostringstream FilePart;
	int Width = std::log10(size())+1;
	FilePart << std::setfill('0') << std::setw(Width) << Part;
	return (FileName + "." + FilePart.str() + ".swdb");
}
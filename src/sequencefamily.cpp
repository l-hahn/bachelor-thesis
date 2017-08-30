#include "sequencefamily.hpp"

sequence_family::sequence_family(){
}

sequence_family::sequence_family(std::string & SeqFile){
	_IsTranslated = false;
	_IsBit = false;
	read_family(SeqFile);
}

sequence_family::sequence_family(std::string && SeqFile){
	_IsTranslated = false;
	_IsBit = false;
	read_family(SeqFile);
}

sequence_family::sequence_family(std::vector<sequence> & SeqVec): _Sequences(SeqVec){
	_IsProtein = _Sequences[0].is_protein();
	_IsTranslated = false;
	_IsBit = false;
}

sequence_family::sequence_family(std::string & SeqFile, std::string & ID, std::string & Name, std::string & Description){
	_IsTranslated = false;
	_IsBit = false;
	read_family(SeqFile);
	_ID = ID;
	_Name = Name;
	_Description = Description;
}

sequence_family::sequence_family(std::string && SeqFile, std::string & ID, std::string & Name, std::string & Description){
	_IsTranslated = false;
	_IsBit = false;
	read_family(SeqFile);
	_ID = ID;
	_Name = Name;
	_Description = Description;
}

sequence_family::sequence_family(std::vector<sequence> & SeqVec, std::string & ID, std::string & Name, std::string & Description): _Sequences(SeqVec), _ID(ID), _Name(Name), _Description(Description){
	_IsProtein = _Sequences[0].is_protein();
	_IsTranslated = false;
	_IsBit = false;
}


void sequence_family::push_back(sequence & Seq){
	if(_Sequences.size() == 0){
		_IsProtein = Seq.is_protein();
		_IsTranslated = _IsProtein;
		_IsBit = Seq.is_bit();
	}
	if(_IsBit && !Seq.is_bit()){
		Seq.to_bit();
	}
	else if(!_IsBit && Seq.is_bit()){
		Seq = sequence(Seq.to_string(), Seq.name());
	}
	if(_IsProtein && Seq.is_dna()){
		translate_sequence(Seq);
	}
	_Sequences.push_back(Seq);
	sort();
}

void sequence_family::translate(){
	if(!_IsTranslated && !_IsProtein){
		_IsTranslated = true;
		_IsProtein = true;
		if(!_IsBit){
			to_bit();
		}
		std::vector<sequence> NewFamily;
		for(auto & Seq : _Sequences){
			std::vector<sequence> NewSequences = translate_sequence(Seq);
			NewFamily.insert(NewFamily.end(), NewSequences.begin(), NewSequences.end());
		}
		_Sequences = std::move(NewFamily);
	}
}

std::vector<sequence> sequence_family::translate_sequence(sequence & Seq, signed SeqPos){
	pattern Codon("111");
	std::vector<sequence> OrfSequences(6); 
	std::vector<dna_spacedword> AminoAcid;
	Seq.spaced_words(AminoAcid, Codon, SeqPos);
	size_t Ctr = 0;
	for(auto & Spw : AminoAcid){
		OrfSequences[2*(Ctr%3)].push_back(alphabet::CodonProteinBit[Spw.bits()]);
		Ctr++;
	}
	std::reverse(Seq.begin(), Seq.end());
	AminoAcid.clear();
	Seq.spaced_words(AminoAcid, Codon, SeqPos);
	Ctr = 0;
	for(auto & Spw : AminoAcid){
		OrfSequences[2*(Ctr%3)+1].push_back(alphabet::CodonProteinBit[Spw.bits_comp()]);
		Ctr++;
	}
	Ctr = 0;
	std::string Name = Seq.name();
	for(auto & NewSeq : OrfSequences){
		NewSeq.set_name(Name);
		Ctr++;
	}
	return OrfSequences;
}

void sequence_family::to_bit(){
	if(!_IsBit){
		_IsBit = true;
		for(auto & Seq : _Sequences){
			Seq.to_bit(_IsProtein);
		}
	}
}

void sequence_family::sort(){
	std::sort(_Sequences.begin(), _Sequences.end());
}


std::string sequence_family::id() const{
	return _ID;
}

std::string sequence_family::name() const{
	return _Name;
}

std::string sequence_family::description() const{
	return _Description;
}

size_t sequence_family::size(size_t MinSize) const{
	size_t SeqSize = _Sequences.size();
	for(auto & Seq : _Sequences){
		if(Seq.size() < MinSize){
			SeqSize--;
		}
		else{
			break;
		}
	}
	return SeqSize;
}

bool sequence_family::is_translated() const{
	return _IsTranslated;
}
bool sequence_family::is_protein() const{
	return _IsProtein;
}
bool sequence_family::is_bit() const{
	return _IsBit;
}

void sequence_family::set_id(std::string && ID){
	set_id(ID);
}
void sequence_family::set_id(std::string & ID){
	_ID = ID;
}

void sequence_family::set_name(std::string && Name){
	set_name(Name);
}
void sequence_family::set_name(std::string & Name){
	_Name = Name;
}

void sequence_family::set_description(std::string && Description){
	set_description(Description);
}
void sequence_family::set_description(std::string & Description){
	_Description = Description;
}


void sequence_family::read_family(std::string & SeqFile){
	FileType SeqType = file_format(SeqFile);
	std::ifstream Input(SeqFile);
 	_IsProtein = read_file_part(Input, _Sequences, _ID, _Name, _Description, SeqType);
 	sort();
 	to_bit();
	Input.close();
}

void sequence_family::print() const{
	std::cout << std::endl;
	std::cout << "ID\t\t" << _ID << std::endl;
	std::cout << "Name\t\t" << _Name << std::endl;
	std::cout << "Description\t" << _Description << std::endl;
	std::cout << "IsProtein\t" << _IsProtein << std::endl;
	std::cout << "IsTranslated\t" << _IsTranslated << std::endl;
	std::cout << "IsBit\t\t" << _IsBit << std::endl;

	unsigned ctr = 0;
	for(auto & Seq : _Sequences){
		std::cout << ++ctr << "\t" << Seq.name() << "\t" << Seq.to_string() << std::endl;
	}
	std::cout << std::endl;
}


sequence & sequence_family::operator[](size_t Idx){
	return _Sequences[Idx];
}

sequence_family::iterator sequence_family::begin(){
	return _Sequences.begin();
}

sequence_family::iterator sequence_family::end(){
	return _Sequences.end();
}

void sequence_family::erase(sequence_family::iterator Idx){
	_Sequences.erase(Idx);
}


std::vector<sequence_family> sequence_family::families_file(std::string & SeqFile){
	FileType SeqType = file_format(SeqFile);
	std::ifstream Input(SeqFile);
	std::vector<sequence_family> FamVec;
	size_t Ctr = 0;
	while(!Input.eof()){
		std::vector<sequence> SeqVec;
		std::string ID, Name, Description;
		bool IsProtein = read_file_part(Input, SeqVec, ID, Name, Description, SeqType);
		if(SeqVec.size() > 0){
			std::cout << "\r" << std::string(80,' ') << "\r";
			std::cout << "\rProcessed family " << ++Ctr << " " << ID << ", Size " << SeqVec.size() << std::flush;
			FamVec.push_back(sequence_family(SeqVec, ID, Name, Description));
			FamVec[FamVec.size()-1]._IsProtein = IsProtein;
			FamVec[FamVec.size()-1].to_bit();
			FamVec[FamVec.size()-1].sort();
		}
	}
	std::cout << "\r" << std::string(80,' ') << "\r";
	return FamVec;
}

bool sequence_family::read_file_part(std::ifstream & Input, std::vector<sequence> & Sequences,
	std::string & FamID, std::string & FamName, std::string & Description, FileType SeqType){
	if(SeqType == FASTA){
		return read_fasta(Input, Sequences, FamID, FamName, Description, SeqType);
	}
	else if(SeqType == STOCKHOLM){
		return read_stockholm(Input, Sequences, FamID, FamName, Description, SeqType);
	}
	else{
		Input.close();
		std::cerr << "Fatal error! Wrong file type!\nMust be FASTA or STOCKHOLM file!" << std::endl;
		std::exit(-1);
	}
}

bool sequence_family::read_fasta(std::ifstream & Input, std::vector<sequence> & Sequences,
	std::string & FamID, std::string & FamName, std::string & Description, FileType SeqType){
	std::string Line, SeqName, TmpSeq = "";
	size_t ProtCount = 0, ReadPos = Input.tellg();
	Description = "UNKNOWN";
	std::getline(Input,Line);
	parse_fasta_header(Line, SeqName, FamID, FamName);
	while(std::getline(Input,Line)){
		if(Line[0] != '>' && Line.size() != 0){
			TmpSeq += Line;
		}
		else if(Line.size() != 0){
			sequence NewSequence(TmpSeq,SeqName);
			if(NewSequence.is_protein()){
				ProtCount++;
			}
			Sequences.push_back(NewSequence);
			TmpSeq = "";
			std::string TmpName, TmpFamID;
			parse_fasta_header(Line, SeqName, TmpFamID, TmpName);
			if(TmpFamID != FamID){
				Input.seekg(ReadPos);
				break;
			}
		}
		ReadPos = Input.tellg();
	}
	if(TmpSeq.size() != 0 && SeqName.size() != 0){
		sequence NewSequence(TmpSeq,SeqName);
		if(NewSequence.is_protein()){
				ProtCount++;
		}
		Sequences.push_back(NewSequence);
	}
	if(ProtCount == 0){
		return false;
	}
	return true;
}

bool sequence_family::read_stockholm(std::ifstream & Input, std::vector<sequence> & Sequences,
	std::string & FamID, std::string & FamName, std::string & Description, FileType SeqType){
	std::string Line;
	size_t ReadPos = Input.tellg(), ProtCount = 0;
	bool ReadOnce = false;
	while(std::getline(Input,Line)){
		if(Line.size() !=0 && Line[0] != '#'){
			std::vector<std::string> SeqVec = split(Line, ' ');
			if(SeqVec.size() == 2){
				sequence NewSequence(SeqVec[1],SeqVec[0]);
				if(NewSequence.is_protein()){
					ProtCount++;
				}
				Sequences.push_back(NewSequence);
			}
			ReadPos = Input.tellg();
		}
		else if(Line.substr(0,8) == "#=GF ID " || Line.substr(0,8) == "#=GF ID\t"){
			if(ReadOnce && Sequences.size() != 0){
				Input.seekg(ReadPos);
				break;
			}
			ReadOnce = true;
			FamName = Line.substr(info_start(Line));
		}
		else if(Line.substr(0,8) == "#=GF DE " || Line.substr(0,8) == "#=GF DE\t"){
			if(ReadOnce && Sequences.size() != 0){
				Input.seekg(ReadPos);
				break;
			}
			ReadOnce = true;
			Description = Line.substr(info_start(Line));
		}
		else if(Line.substr(0,8) == "#=GF AC " || Line.substr(0,8) == "#=GF AC\t"){
			if(ReadOnce && Sequences.size() != 0){
				Input.seekg(ReadPos);
				break;
			}
			ReadOnce = true;
			FamID = Line.substr(info_start(Line));
		}
	}
	if(Sequences.size() != 0){
		return ProtCount/Sequences.size();
	}
	return true;
}


void parse_fasta_header(std::string & Line, std::string & SeqName, std::string & FamID, std::string & FamName){
	std::replace(Line.begin(), Line.end(), '>', ' ');
	std::replace(Line.begin(), Line.end(), ';', ' ');
	std::vector<std::string> InfoVec = split(Line, ' ');
	if(InfoVec.size()  == 4){
		SeqName = std::move(InfoVec[0]);
		FamID = std::move(InfoVec[2]);
		FamName = std::move(InfoVec[3]);
	}
	else{
		SeqName = std::move(InfoVec[0]);
		FamID = "UNKNOWN";
		FamName = "UNKNOWN";
	}
}


FileType file_format(std::string & SeqFile){
	std::ifstream Input(SeqFile);
	assert(Input.is_open());

	std::string Head;
	std::getline(Input, Head);
	Input.close();
	
	if(Head.size() == 0){
		std::cerr << "Fatal Error! EMPTY FILE!" << std::endl;
		std::exit(-1);
	}

	if(Head[0] == '>')
		return FASTA;
	else if(Head[0] == '#' && Head.find("STOCKHOLM") != std::string::npos){
		return STOCKHOLM;
	}
	else{
		return FALSETYPE;
	}
}


unsigned info_start(std::string & Str){
	unsigned pos = 7;
	while(Str[pos] == ' '){
		pos++;
	}
	return pos;
}

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
    	if(item.empty() == false){

        	elems.push_back(item);
    	}
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

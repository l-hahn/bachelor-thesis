#include "pattern.hpp"

pattern::pattern():  _Score(1), _BitPattern(0), _IsBit(false){
	random();
}

pattern::pattern(unsigned DontCare, unsigned Weight): _Score(1), _BitPattern(0),_IsBit(false){
	random(DontCare, Weight);
}

pattern::pattern(unsigned DontCare, unsigned Weight, uint64_t Seed): _Score(1), _BitPattern(0),_IsBit(false){
	random(DontCare, Weight, Seed);
}

pattern::pattern(std::string & StrPat) : _Score(1), _BitPattern(0),_IsBit(false){
	_VectorPattern = std::vector<unsigned char>(StrPat.begin(), StrPat.end());
	_check_pattern();
}


pattern::pattern(std::string && StrPat) : _Score(1), _BitPattern(0),_IsBit(false){
	_VectorPattern = std::vector<unsigned char>(StrPat.begin(), StrPat.end());
	_check_pattern();
}


pattern::pattern(std::vector<char> & StrPat) : _Score(1), _BitPattern(0),_IsBit(false){
	_VectorPattern = std::vector<unsigned char>(StrPat.begin(), StrPat.end());
	_check_pattern();
}


pattern::pattern(std::vector<unsigned char> & StrPat) : _Score(1), _BitPattern(0),_IsBit(false){
	_VectorPattern = std::vector<unsigned char>(StrPat.begin(), StrPat.end());
	_check_pattern();
}


void pattern::_check_pattern(){
	unsigned Ctr = 0;

	if(_VectorPattern.size() < 64){
		_IsBit = true;
	}

	for(auto & C : _VectorPattern){
		switch(C){
			case '1':
			case 'X':
			case 'x':
			case '#':
			case '*':
				C = 1;
				break;
			case '0':
			case 'O':
			case 'o':
			case '-':
				C = 0;
				break;
		}
		if(C != 1 && C != 0){
			std::cerr << "Illegal character in pattern!\nFormat: match = {1,X,x,#,*,} | don't care = {0,O,o,-,}" << std::endl;
			std::exit(-1);
		}
		if(C == 1){
			_MatchPos.push_back(Ctr);
		}
		else{
			_DCPos.push_back(Ctr);
		}
		_BitPattern <<= 1;
		_BitPattern |= C;
		Ctr++;
	}
}


void pattern::random(){
	std::uniform_int_distribution<unsigned> RandomLength(2, 63);
	std::random_device RandBit;
	std::default_random_engine BitRandom(RandBit());
	std::mt19937 BitGenerator(BitRandom());
	unsigned RandLength = RandomLength(BitGenerator);
	std::uniform_int_distribution<unsigned> RandomWeight(2, RandLength);
	unsigned RandWeight = RandomWeight(BitGenerator);
	uint64_t Seed = RandBit();
	random(RandLength - RandWeight, RandWeight, Seed);
	// unsigned RandLength = rand()%(63-2+1)+2;
	// unsigned RandWeight = rand()%(RandLength-2 + 1) + 2;
	// random(RandLength - RandWeight, RandWeight, rand());
}

void pattern::random(unsigned DontCare, unsigned Weight){
	std::random_device RandomBit;
	random(DontCare, Weight, RandomBit());
	// random(DontCare, Weight, rand());
}

void pattern::random(unsigned DontCare, unsigned Weight, uint64_t Seed){
	unsigned Length = DontCare + Weight;
	if(Length < 64){
		_IsBit = true;
	}
	if(Weight < 2){
		std::cerr << "Illegal value for weight!\nMinimum allowed weight is 2!" << std::endl;
		std::exit(-1);
	}
	std::default_random_engine BitRandom(Seed);
	std::mt19937 BitGenerator(BitRandom());
	std::uniform_int_distribution<unsigned> BitPosition(1, Length-2);

	_MatchPos.clear();
	_DCPos.clear();
	_BitPattern = 0;
	_VectorPattern = std::vector<unsigned char>(Length,0);
	
	_VectorPattern[0] = 1;
	_VectorPattern[Length-1] = 1;
	_BitPattern = ((uint64_t) 1 << (Length-1)) | (uint64_t)1;
	Weight -= 2;

	while(Weight > 0){
		unsigned Pos = BitPosition(BitGenerator);
		// unsigned Pos = rand()%(Length-2 -1 + 1) + 1;
		if(_VectorPattern[Pos] == 0){
			Weight--;
			_VectorPattern[Pos] = 1;
			_BitPattern |= (uint64_t)1 << Pos;
		}
	}
	for(unsigned i = 0; i < _VectorPattern.size(); i++){
		if(_VectorPattern[i] == 1){
			_MatchPos.push_back(i);
		}
		else{
			_DCPos.push_back(i);
		}
	}
}

std::string pattern::to_string(){
	std::string Pattern(_VectorPattern.begin(), _VectorPattern.end());
	std::transform(Pattern.begin(), Pattern.end(), Pattern.begin(),[](unsigned char C) { return (C+'0'); });
	return Pattern;
}


bool pattern::is_match(unsigned Pos) const{
	return _VectorPattern[Pos] == 1;
}


void pattern::bit_swap(unsigned PosA, unsigned PosB){
	std::swap(_VectorPattern[PosA],_VectorPattern[PosB]);
	if(_IsBit && _VectorPattern[PosA] !=_VectorPattern[PosB]){
		_BitPattern ^= ((uint64_t)1 << (_VectorPattern.size()-PosA-1)) | ((uint64_t)1 << (_VectorPattern.size()-PosB-1));
			if(_VectorPattern[PosA] != 0){
				std::swap(PosA,PosB);
			}
			_MatchPos.erase(std::find(_MatchPos.begin(), _MatchPos.end(), PosA));
			_DCPos.erase(std::find(_DCPos.begin(), _DCPos.end(), PosB));
			_MatchPos.push_back(PosB);
			_DCPos.push_back(PosA);
			std::sort(_MatchPos.begin(), _MatchPos.end());
			std::sort(_DCPos.begin(), _DCPos.end());
	}
}

void pattern::random_swap(){
	std::random_device RandomBit;
	random_swap(RandomBit());
	// random_swap(rand());
}

void pattern::random_swap(uint64_t Seed){
	if(_DCPos.size() != 0){
		std::default_random_engine BitRandom(Seed);
		std::mt19937 BitGenerator(BitRandom());
		std::uniform_int_distribution<unsigned> SwapPosition(1, _MatchPos.size()-2);
		unsigned MatchPosition = _MatchPos[SwapPosition(BitGenerator)];
		SwapPosition = std::uniform_int_distribution<unsigned>(0, _DCPos.size()-1);
		unsigned DCPosition = _DCPos[SwapPosition(BitGenerator)];
		// unsigned MatchPosition = rand()%(_MatchPos.size() - 2 - 1 + 1) + 1;
		// unsigned DCPosition = rand()%(_DCPos.size() - 1 - 0 + 1) + 0;
		bit_swap(MatchPosition, DCPosition);
	}
}


unsigned char pattern::operator[](size_t Idx) const{
	return _VectorPattern[Idx];
}

bool pattern::operator<(const pattern & P) const{
	return _Score < P._Score;
}

bool pattern::operator>(const pattern & P) const{
	return _Score > P._Score;
}

bool pattern::operator==(const pattern & P) const{
	if(_IsBit){
		return _BitPattern == P._BitPattern;
	}
	return _VectorPattern == P._VectorPattern;
}




void pattern::set_score(double Scr){
	_Score = Scr;
}

double pattern::score() const{
	return _Score;
}

unsigned pattern::weight() const{
	return _MatchPos.size();
}

unsigned pattern::length() const{
	return _VectorPattern.size();
}
unsigned pattern::dontcare() const{
	return _VectorPattern.size() - _MatchPos.size();
}

unsigned pattern::get_overlap(const pattern & P, int Shift) const{
	unsigned Overlap = 0;
	if(P._IsBit && _IsBit){
		uint64_t ShiftPattern = _BitPattern;
		uint64_t BitPat = P._BitPattern;

		if(Shift < 0){
			std::swap(ShiftPattern, BitPat);
			Shift = std::abs(Shift);
		}
		
		ShiftPattern >>= Shift;
		ShiftPattern &= BitPat;
		unsigned ShiftSize = std::log2(ShiftPattern)+1;
		
		for(unsigned i = 0; i < ShiftSize; i++){
			Overlap += (ShiftPattern & (uint64_t)1);
			ShiftPattern >>= 1;
		}
	}
	else{
		auto PatAIt = _VectorPattern.begin();
		auto PatBIt = P._VectorPattern.begin();
		
		Shift += (int)P._VectorPattern.size()-(int)_VectorPattern.size();
		
		if(Shift > 0){
			PatBIt += Shift;
		}
		else{
			PatAIt -= Shift;
		}
		
		for(;PatAIt != _VectorPattern.end() && PatBIt != P._VectorPattern.end(); PatBIt++, PatAIt++){
			Overlap += *PatAIt * *PatBIt;
		}
	}
	return Overlap;
}

pattern::iterator pattern::begin(){
	return _MatchPos.begin();
}

pattern::iterator pattern::end(){
	return _MatchPos.end();
}
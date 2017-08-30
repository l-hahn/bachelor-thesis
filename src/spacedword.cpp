template<int _BitLength>
spacedword<_BitLength>::spacedword(): _BitSpacedWord(0), _SequencePos(-1), _SequenceNo(-1), _SequenceCtr(1){}

template<int _BitLength>
spacedword<_BitLength>::spacedword(std::vector<char> & SeqVec, pattern & Pat, int64_t SeqPos, signed SeqNo): _SequencePos(SeqPos), _SequenceNo(SeqNo), _SequenceCtr(1){
	create_word(SeqVec, Pat);
}

template<int _BitLength>
spacedword<_BitLength>::spacedword(std::vector<unsigned char> & SeqVec, pattern & Pat, int64_t SeqPos, signed SeqNo): _SequencePos(SeqPos), _SequenceNo(SeqNo), _SequenceCtr(1){
	create_word(SeqVec, Pat);
}

template<int _BitLength>
spacedword<_BitLength>::spacedword(int64_t Spw, int64_t SeqPos, signed SeqNo): _BitSpacedWord(Spw), _SequencePos(SeqPos), _SequenceNo(SeqNo), _SequenceCtr(1){

}


template<int _BitLength>
std::string spacedword<_BitLength>::to_string() const{
	std::map<unsigned char, unsigned char> & BitMap = alphabet::BitMap(_BitLength);
	std::string AlphaWord = "";
	int64_t BitWord = _BitSpacedWord;
	char BitMask = (1 << _BitLength)-1; 
	while(BitWord > 0 && _BitLength > 0){
		AlphaWord.push_back(BitMap[BitWord & BitMask]);
		BitWord >>= _BitLength;
	}
	std::reverse(AlphaWord.begin(),AlphaWord.end());
	return AlphaWord;
}

template<int _BitLength>
template<typename T>
void spacedword<_BitLength>::create_word(T & Vec, pattern & Pat){
	auto SeqIt = Vec.begin() + _SequencePos;
	_BitSpacedWord = 0;
	for(auto & PatIt : Pat){
		push_back(*(SeqIt+PatIt));
	}
}

template<int _BitLength>
void spacedword<_BitLength>::push_back(unsigned char Letter){
	std::map<unsigned char, unsigned char> & AlphaMap = alphabet::AlphaMap(_BitLength);
	_BitSpacedWord <<= _BitLength;
	if(Letter > alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1]){
		_BitSpacedWord |= AlphaMap[Letter];
	}
	else{
		_BitSpacedWord |= Letter;
	}
}

template<int _BitLength>
unsigned spacedword<_BitLength>::size() const{
	return (std::log2(_BitSpacedWord)+1)/_BitLength;
}

template<int _BitLength>
int64_t spacedword<_BitLength>::position() const{
	return _SequencePos;
}

template<int _BitLength>
int64_t spacedword<_BitLength>::bits() const{
	return _BitSpacedWord;
}

template<int _BitLength>
int64_t spacedword<_BitLength>::bits_comp() const{
	int64_t RevCompBits=0, TmpBits = _BitSpacedWord;
	unsigned char BitMask = (1 << _BitLength)-1;
	while(TmpBits > 0){
		RevCompBits <<= _BitLength;
		RevCompBits |= alphabet::DnaInvert(TmpBits & BitMask);
		TmpBits >>= _BitLength;
	}
	while(RevCompBits > 0){
		TmpBits <<= _BitLength;
		TmpBits |= (RevCompBits & BitMask);
		RevCompBits >>=_BitLength;
	}
	RevCompBits = TmpBits;
	return RevCompBits;
}

template<int _BitLength>
signed spacedword<_BitLength>::sequence() const{
	return _SequenceNo;
}

template<int _BitLength>
signed spacedword<_BitLength>::counter() const{
	return _SequenceCtr;
}

template<int _BitLength>
void spacedword<_BitLength>::set_position(int64_t Pos){
	_SequencePos = Pos;
}

template<int _BitLength>
void spacedword<_BitLength>::set_sequence(signed No){
	_SequenceNo = No;
}

template<int _BitLength>
void spacedword<_BitLength>::set_counter(signed Ctr){
	_SequenceCtr = Ctr;
}

template<int _BitLength>
unsigned char spacedword<_BitLength>::operator[](unsigned Idx) const{
	std::map<unsigned char, unsigned char> & BitMap = alphabet::BitMap(_BitLength);
	char BitMask = (1 << _BitLength)-1;
	int64_t BitLetter = _BitSpacedWord;
	for(unsigned i = 0; i < size()-Idx; i++){
		BitLetter >>= _BitLength;
	}
	return BitMap[BitLetter & BitMask];
}


template<int _BitLength>
bool spacedword<_BitLength>::operator==(const spacedword<_BitLength> & Spw) const{
	return _BitSpacedWord == Spw._BitSpacedWord;
}

template<int _BitLength>
bool spacedword<_BitLength>::operator!=(const spacedword<_BitLength> & Spw) const{
	return _BitSpacedWord != Spw._BitSpacedWord;
}

template<int _BitLength>
bool spacedword<_BitLength>::operator<(const spacedword<_BitLength> & Spw) const{
	if(_BitSpacedWord == Spw._BitSpacedWord){
		return _SequenceNo < Spw._SequenceNo;
	}
	return _BitSpacedWord < Spw._BitSpacedWord;
}

template<int _BitLength>
bool spacedword<_BitLength>::operator>(const spacedword<_BitLength> & Spw) const{
	if(_BitSpacedWord == Spw._BitSpacedWord){
		return _SequenceNo > Spw._SequenceNo;
	}
	return _BitSpacedWord > Spw._BitSpacedWord;
}

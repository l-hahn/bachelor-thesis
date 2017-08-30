#include "spacedfamily.hpp"

spacedword_family::spacedword_family(int64_t BitSpacedWord):_BitSpacedWord(BitSpacedWord){	
}

spacedword_family::spacedword_family(int64_t BitSpacedWord, unsigned Family, unsigned Position, double FamilyScore):_BitSpacedWord(BitSpacedWord){
	push_back(Family, Position, FamilyScore);	
}

spacedword_family::spacedword_family(int64_t BitSpacedWord, family_score & FamilyScore):_BitSpacedWord(BitSpacedWord){
	push_back(FamilyScore);	
}


int64_t spacedword_family::bits() const{
	return _BitSpacedWord;
}

unsigned spacedword_family::size() const{
	return _Family.size();
}

void spacedword_family::push_back(unsigned Family, unsigned Position, double FamilyScore){
	_Family.push_back(family_score(Family, Position, FamilyScore));
}

void spacedword_family::push_back(family_score & FamilyScore){
	_Family.push_back(FamilyScore);
}


family_score spacedword_family::operator[](size_t Idx){
	return _Family[Idx];
}

bool spacedword_family::operator==(const spacedword_family & Spw) const{
	return _BitSpacedWord == Spw._BitSpacedWord;
}

bool spacedword_family::operator<(const spacedword_family & Spw) const{
	return _BitSpacedWord < Spw._BitSpacedWord;
}

bool spacedword_family::operator>(const spacedword_family & Spw) const{
	return _BitSpacedWord > Spw._BitSpacedWord;
}


spacedword_family::iterator spacedword_family::begin(){
	return _Family.begin();
}

spacedword_family::iterator spacedword_family::end(){
	return _Family.end();
}
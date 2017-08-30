#include "patternset.hpp"

patternset::patternset():_Score(1),_MaxWeight(0),_MinWeight(std::numeric_limits<unsigned>::max()),
	_MaxDontCare(0),_MinDontCare(std::numeric_limits<unsigned>::max()){	
}

patternset::patternset(unsigned Size, unsigned MaxW, unsigned MinW, 
	unsigned MaxD, unsigned MinD):_Score(1),_MaxWeight(MaxW),
	_MinWeight(MinW),_MaxDontCare(MaxD),_MinDontCare(MinD){
	random(Size, MaxW, MinW, MaxD, MinD);
}


void patternset::sort(){
	std::sort(_PatternSet.begin(), _PatternSet.end());
}

void patternset::push_back(pattern & Pat){
	_MaxDontCare = std::max(_MaxDontCare, Pat.dontcare());
	_MinDontCare = std::min(_MinDontCare, Pat.dontcare());
	_MaxWeight = std::max(_MaxWeight, Pat.weight());
	_MinWeight = std::min(_MinWeight, Pat.weight());
	_PatternSet.push_back(Pat);
}

void patternset::push_back(pattern && Pat){
	push_back(Pat);
}

void patternset::push_back(std::string & Pat){
	push_back(pattern(Pat));
}

void patternset::push_back(std::string && Pat){
	push_back(pattern(Pat));
}

void patternset::push_back(std::vector<char> & Pat){
	push_back(pattern(Pat));
}

void patternset::push_back(std::vector<unsigned char> & Pat){
	push_back(pattern(Pat));
}

void patternset::random(unsigned Size, unsigned Length, unsigned Weight){
	random(Size, Length, Length, Weight, Weight);
}

void patternset::random(unsigned Size, unsigned MinLength, unsigned MaxLength,
	unsigned Weight){
	random(Size, MinLength, MaxLength, Weight, Weight);
}

void patternset::random(unsigned Size, unsigned MinLength, unsigned MaxLength,
	unsigned MinWeight, unsigned MaxWeight){
	_PatternSet.clear();
	if(Size >= 2){
		double CurLeng = MinLength;
		double CurWeight = MinWeight;
		double StepLeng = (MaxLength - MinLength+1)/(double)(Size);
		double StepWeight = (MaxWeight - MinWeight+1)/(double)(Size);
		for(unsigned i = 0; i < Size-1; i++){
			_PatternSet.push_back(pattern((unsigned)CurLeng,(unsigned)CurWeight));
			CurLeng += StepLeng;
			CurWeight += StepWeight;
		}
		_PatternSet.push_back(pattern(MaxLength,MaxWeight));
	}
	else{
		_PatternSet.push_back(pattern((MaxLength + MinLength)/2,(MaxWeight + MinWeight)/2));
	}
}

bool patternset::is_uniq(const pattern & Pat) const{
	return (std::find_if(_PatternSet.begin(),_PatternSet.end(), [& Pat](const pattern & P){
		if( & Pat != & P){
			return (Pat == P);
		}
		return false;
	}) == _PatternSet.end());
}


pattern & patternset::operator[](size_t Idx){
	return _PatternSet[Idx];
}

bool patternset::operator<(const patternset & P) const{
	return _Score < P._Score;
}

bool patternset::operator>(const patternset & P) const{
	return _Score > P._Score;
}

unsigned patternset::max_weight() const{
	return _MaxWeight;
}
unsigned patternset::min_weight() const{
	return _MinWeight;
}
unsigned patternset::weight() const{
	unsigned MeanWeight = 0;
	for(auto & Pat : _PatternSet){
		MeanWeight += Pat.weight();
	}
	return MeanWeight/_PatternSet.size();
}

unsigned patternset::max_dontcare() const{
	return _MaxDontCare;
}
unsigned patternset::min_dontcare() const{
	return _MinDontCare;
}
unsigned patternset::dontcare() const{
	unsigned MeanDC = 0;
	for(auto & Pat : _PatternSet){
		MeanDC += Pat.dontcare();
	}
	return MeanDC/_PatternSet.size();
}

unsigned patternset::max_length() const{
	return _MaxWeight + _MaxDontCare;
}

unsigned patternset::min_length() const{
	return _MinWeight + _MinDontCare;
}
unsigned patternset::length() const{
	unsigned MeanLength = 0;
	for(auto & Pat : _PatternSet){
		MeanLength += Pat.length();
	}
	return MeanLength/_PatternSet.size();
}

unsigned patternset::size() const{
	return _PatternSet.size();
}


patternset::iterator patternset::begin(){
	return _PatternSet.begin();
}

patternset::iterator patternset::end(){
	return _PatternSet.end();
}
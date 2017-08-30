#include "spacedhit.hpp"

spacedhit::spacedhit(size_t FamilyID, size_t SeqNumber, unsigned short Orf, size_t DBFamID, size_t Position, double Score): _Score(Score), _FamilyID(FamilyID), _SeqNumber(SeqNumber), _DBFamID(DBFamID), _Position(Position), _Orf(Orf){
}


bool spacedhit::operator<(const spacedhit & SwHit) const{
	if(_FamilyID == SwHit._FamilyID){
		if(_SeqNumber == SwHit._SeqNumber){
			if(_DBFamID == SwHit._DBFamID){
				if(_Position == SwHit._Position){
					return _Orf < SwHit._Orf;
				}
				return _Position < SwHit._Position;
			}
			return _DBFamID < SwHit._DBFamID;
		}
		return _SeqNumber < SwHit._SeqNumber;
	}
	return _FamilyID < SwHit._FamilyID;
}

bool spacedhit::operator>(const spacedhit & SwHit) const{
	if(_FamilyID == SwHit._FamilyID){
		if(_SeqNumber == SwHit._SeqNumber){
			if(_DBFamID == SwHit._DBFamID){
				if(_Position == SwHit._Position){
					return _Orf > SwHit._Orf;
				}
				return _Position > SwHit._Position;
			}
			return _DBFamID > SwHit._DBFamID;
		}
		return _SeqNumber > SwHit._SeqNumber;
	}
	return _FamilyID > SwHit._FamilyID;
}


size_t spacedhit::fam_id() const{
	return _FamilyID;
}

size_t spacedhit::seq_num() const{
	return _SeqNumber;
}

size_t spacedhit::db_fam_id() const{
	return _DBFamID;
}

size_t spacedhit::position() const{
	return _Position;
}

unsigned short spacedhit::orf() const{
	return _Orf;
}

double spacedhit::score() const{
	return _Score;
}

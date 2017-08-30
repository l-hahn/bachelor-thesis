#ifndef SPACEDHIT_HPP_
#define SPACEDHIT_HPP_

#include <cstdlib>

class spacedhit{
	public:
		spacedhit(size_t FamilyID, size_t SeqNumber, unsigned short Orf, size_t DBFamID, size_t Position, double Score);
		bool operator<(const spacedhit & SwHit) const;
		bool operator>(const spacedhit & SwHit) const;

		size_t fam_id() const;
		size_t seq_num() const;
		size_t db_fam_id() const;
		size_t position() const;
		unsigned short orf() const;
		double score() const;

	private:
		double _Score;
		size_t _FamilyID;
		size_t _SeqNumber;
		size_t _DBFamID;
		size_t _Position;
		unsigned short _Orf;
};
#endif
#ifndef FAMILYSCORE_HPP_
#define FAMILYSCORE_HPP_

class family_score{
	public:
		family_score(unsigned ID, unsigned Position, double Score);
		
		unsigned id() const;
		unsigned position() const;
		double score() const;

		bool operator<(const family_score & F) const;
		bool operator>(const family_score & F) const;
		bool operator==(const family_score & F) const;


	private:
		unsigned _ID;
		unsigned _Position;
		double _Score;
};
#endif
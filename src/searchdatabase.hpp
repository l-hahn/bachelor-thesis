#ifndef SEARCHDATABASE_HPP_
#define SEARCHDATABASE_HPP_

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <omp.h>
#include <vector>
#include "alphabet.hpp"
#include "dsoptions.hpp"
#include "familyscore.hpp"
#include "pattern.hpp"
#include "sequence.hpp"
#include "sequencefamily.hpp"
#include "spacedfamily.hpp"
#include "spacedhit.hpp"
#include "spacedword.hpp"
#include "spacedworddb.hpp"

void search_database();
std::vector<spacedhit> family_hit(sequence_family & Family, unsigned FamID, spacedword_db & SwDataBase);
std::vector<spacedhit> sequence_hit(sequence_family & Family, unsigned FamID, pattern & Pat, std::vector<spacedword_family> & DbSw);

template<typename T>
std::vector<spacedhit> sw_sequence_hit(T & NOUSE, sequence_family & Family, unsigned FamID, pattern & Pat, std::vector<spacedword_family> & DbSw, size_t & CtrSeq);

void init_omp();
void eval_classify(spacedword_db & SwDataBase, std::vector<sequence_family> & Families, std::vector<spacedhit> & DataBaseHit);

template<typename T>
std::vector<spacedhit> sw_sequence_hit(T & NOUSE, sequence_family & Family, unsigned FamID, pattern & Pat, std::vector<spacedword_family> & DbSw, size_t & CtrSeq){
	std::vector<spacedhit> SwMatches;
	std::vector<T> SpacedWords;
	size_t Ctr = 0;
	for(auto & Seq : Family){
		std::vector<T> TmpSw;
		Seq.spaced_words(TmpSw, Pat, Ctr);
		if(TmpSw.size() > 0){
			SpacedWords.insert(SpacedWords.end(), TmpSw.begin(), TmpSw.end());
			Ctr++;
		}
	}

//	std::cout << "Family " << Family.id() << " " << Family.size() << " -> " << Ctr << std::endl;

	CtrSeq = Ctr;

	std::sort(SpacedWords.begin(), SpacedWords.end(), [](const T & SwA, const T & SwB){
		return SwA.bits() < SwB.bits();
	});

	auto DbSwIt = DbSw.begin();
	auto SwIt = SpacedWords.begin();
	while(DbSwIt != DbSw.end() && SwIt != SpacedWords.end()){
		if(SwIt->bits() < DbSwIt->bits()){
			while(SwIt->bits() < DbSwIt->bits()){
				SwIt++;
				if(SwIt == SpacedWords.end()){
					break;
				}
			}
		}
		else if(SwIt->bits() == DbSwIt->bits()){
			while(SwIt->bits() == DbSwIt->bits()){
				for(auto & FamScore : *DbSwIt){
					int64_t Pos = SwIt->position();
					if(dsoptions::Translated){
						Pos *= 3; //Protein -> DNA retranslation for positions, if origin was DNA!
						Pos += 1;
					}
					SwMatches.push_back(spacedhit(FamScore.position(), SwIt->sequence(), 0, FamScore.id(), Pos, FamScore.score()));
				}
				SwIt++;
			}
		}
		else{
			while(SwIt->bits() > DbSwIt->bits()){
				DbSwIt++;
				if(DbSwIt == DbSw.end()){
					break;
				}
			}
		}
	}

	std::sort(SwMatches.begin(), SwMatches.end(), [](const spacedhit & SwHA, const spacedhit & SwHB){
		if(SwHA.db_fam_id() == SwHB.db_fam_id()){
			if(SwHA.seq_num() == SwHB.seq_num()){
				return SwHA.position() < SwHB.position();
			}
			return SwHA.seq_num() < SwHB.seq_num();
		}
		return SwHA.db_fam_id() < SwHB.db_fam_id();
	});

	std::vector<spacedhit> TrueSwMatches;
	auto BucketStart = SwMatches.begin(), BucketIter = SwMatches.begin();
	bool NotEnd = true, BlockDist = false;
	if(SwMatches.size() != 0){
		while(NotEnd){
			if(BucketIter == SwMatches.end()){
				NotEnd = false;
				BlockDist = true;
			}
			else{
				if(BucketStart->db_fam_id() != BucketIter->db_fam_id()){
					BlockDist = true;
				}
				else if(BucketStart->seq_num() != BucketIter->seq_num()){
					BlockDist = true;
				}
				else if((BucketIter->position() - BucketStart->position()) > (unsigned)dsoptions::BlockSize){
					BlockDist = true;
				}
			}
			if(BlockDist){
				double BlockScore = BucketStart->score();
				unsigned MatchCount = 1;
				for(auto BuckIt = BucketStart + 1; BuckIt < BucketIter; BuckIt++){
					if(BuckIt->fam_id() > (BuckIt-1)->fam_id()){
						BlockScore += BuckIt->score();
						MatchCount++;
					}
				}
				if(MatchCount >= (unsigned)dsoptions::WordPerSeq && BlockScore >= dsoptions::HitThreshold){
					unsigned short SeqOrf = 0;
					size_t SeqID = BucketStart->seq_num();
					if(dsoptions::Translated){
						SeqOrf = SeqID % 6;
						SeqID /= 6;
					}
					size_t Pos = ((BucketIter - 1)->position() - BucketStart->position())/2 + BucketStart->position();
					TrueSwMatches.push_back(spacedhit(FamID, SeqID, SeqOrf, BucketStart->db_fam_id(), Pos , BlockScore));
				}
				BucketStart = BucketIter;
				BlockDist = false;
			}
			if(NotEnd){
				BucketIter++;
			}
		}
		SwMatches = std::move(TrueSwMatches);
	}
	return SwMatches;
}
#endif
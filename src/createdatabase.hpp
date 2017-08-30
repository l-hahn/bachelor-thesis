#ifndef CREATEDATABASE_HPP_
#define CREATEDATABASE_HPP_

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <omp.h>
#include <string>
#include <vector>
#include "alphabet.hpp"
#include "dboptions.hpp"
#include "sequencefamily.hpp"
#include "spacedword.hpp"
#include "spacedworddb.hpp"
#include "pattern.hpp"
#include "patternset.hpp"

void create_database();
void create_family_spacedwords(spacedword_db & SpacedWordFamily, sequence_family & Family, unsigned FamID, patternset & PatSet,  double & SeqCover);
std::vector<spacedword_family> create_spacedword_buckets(patternset & UniqSet, pattern & Pat, sequence_family & Family, unsigned FamID,  double & SeqCover);
void get_seqcov(std::vector<sequence_family> & Families, std::vector<double> & SequenceCoverage);
void get_overlap(spacedword_db & SwDataBase);

template<typename T>
std::vector<spacedword_family> make_buckets(T & NOUSE, patternset & UniqSet, pattern & Pat, sequence_family & Family, unsigned FamID,  double & SeqCover);

void clean_line();
void init_omp();

template<typename T>
std::vector<spacedword_family> make_buckets(T & NOUSE, patternset & UniqSet, pattern & Pat, sequence_family & Family, unsigned FamID, double & SeqCover){
	std::vector<spacedword_family> SpacedBuckets;
	size_t Ctr = 0, FamSize = Family.size(Pat.length());
	while(!UniqSet.is_uniq(Pat) && Ctr++ < 50){
		Pat.random_swap();
	}
	if(Ctr == 50){
		return SpacedBuckets;
	}
	std::vector<T> SpacedWords;
	Ctr = 0; 
	for(auto & Seq : Family){
		std::vector<T> TmpSpw;
		Seq.spaced_words(TmpSpw, Pat, Ctr);
		if(TmpSpw.size() > 0){
			Ctr++;
			std::sort(TmpSpw.begin(), TmpSpw.end());
			TmpSpw.erase(std::unique(TmpSpw.begin(), TmpSpw.end(), [](const T & SpwA, const T & SpwB){ if(SpwA == SpwB){return SpwA.sequence() == SpwB.sequence();}return false;}), TmpSpw.end());
			SpacedWords.insert(SpacedWords.end(), TmpSpw.begin(), TmpSpw.end());
		}
	}
	std::sort(SpacedWords.begin(), SpacedWords.end());

	auto BucketStart = SpacedWords.begin();
	for(auto BucketIter = SpacedWords.begin(); BucketIter != SpacedWords.end(); BucketIter++){
		if(*BucketIter != *BucketStart){
			signed Length = std::distance(BucketStart, BucketIter);
			for(auto BucketIt = BucketStart; BucketIt != BucketIter; BucketIt++){
				BucketIt->set_counter(Length);
			}
			BucketStart = BucketIter;
		}
	}
	if(BucketStart != SpacedWords.end()){
		signed Length = std::distance(BucketStart, SpacedWords.end());
		for(auto BucketIt = BucketStart; BucketIt != SpacedWords.end(); BucketIt++){
			BucketIt->set_counter(Length);
		}
	}

	std::sort(SpacedWords.begin(), SpacedWords.end(),[](const T & SpwA, const T & SpwB){
		if(SpwA.counter() == SpwB.counter()){
			return SpwA < SpwB;
		}
		return SpwA.counter() > SpwB.counter();
	});

	size_t ReSize = 0;
	for(auto & Seq : SpacedWords){
		if(Seq.counter() == 1){
			break;
		}
		ReSize++;
	}
	SpacedWords.resize(ReSize);

	if(dboptions::Greedy && SpacedWords.size() != 0){
		std::vector< T > SpacedPositions;
		std::vector<int> UsedSeq(FamSize,0);
		size_t MeanPos = 0, SeqCounter = 0, MaxBucketLeng = SpacedWords[0].counter();
		bool Uncovered = false;
		BucketStart = SpacedWords.begin();
		for(auto BucketIter = SpacedWords.begin(); BucketIter != SpacedWords.end(); BucketIter++){
			if(*BucketIter != *BucketStart){
				if(Uncovered){
					Uncovered = false;
					MeanPos /= std::distance(BucketStart, BucketIter)+1;
					SpacedPositions.push_back(T(BucketStart->bits(), MeanPos, BucketStart->counter()));
					if((double)SeqCounter >= (double)(dboptions::SeqCovThreshold*FamSize)){
						break;
					}
				}
				BucketStart = BucketIter;
				MeanPos = 0;
			}
			if(UsedSeq[BucketIter->sequence()] == 0){
				Uncovered = true;
			}
			UsedSeq[BucketIter->sequence()]++;
			if((size_t)UsedSeq[BucketIter->sequence()] == dboptions::WordPerSeq){
				SeqCounter++;
			}
			MeanPos += BucketIter->position();
		}
		if(BucketStart != SpacedWords.end() && Uncovered){
			MeanPos /= std::distance(BucketStart, SpacedWords.end())+1;
			SpacedPositions.push_back(T(BucketStart->bits(), MeanPos, BucketStart->counter()));
		}
		SpacedWords = std::move(SpacedPositions);

		std::sort(SpacedWords.begin(), SpacedWords.end(), [](const T & SpwA, const T & SpwB){
			if(SpwA.position() == SpwB.position()){
				return SpwA < SpwB;
			}
			return SpwA.position() < SpwB.position();
		});

		Ctr = 0;
		for(auto & Spw : SpacedWords){
			SpacedBuckets.push_back(spacedword_family(Spw.bits(), FamID, Ctr++, Spw.sequence()/(double)MaxBucketLeng));
		}
		SeqCover = 0;
		for(auto & Val : UsedSeq){
			if((size_t)Val >= dboptions::WordPerSeq){
				SeqCover++;
			}
		}
		SeqCover /= FamSize;
	}
	else if(SpacedWords.size() != 0){
		// TODO: NOT GREEDY, look for best serving SW in remaining List.
	}
	return SpacedBuckets;
}

#endif
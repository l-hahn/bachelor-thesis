#include "createdatabase.hpp"


void create_database(){
	patternset PatSet;
	if(dboptions::PatternFile.size() != 0){
		//TODO READ PATTERN FILE
	}
	std::vector<sequence_family> Families = sequence_family::families_file(dboptions::InSeqFam);
	clean_line();
	if(Families.size() == 1){
		std::cout << Families.size() << " Family was read. Inputtype: ";
	}
	else{
		std::cout << Families.size() << " Families were read. Inputtype: ";
	}
	dboptions::FamilySize = Families.size();
	if(Families[0].is_protein()){
		dboptions::BitLength = alphabet::ProteinBitLength;
		std::cout << "Protein " << std::endl;
	}
	else{
		std::cout << "DNA" << std::endl;
		dboptions::BitLength = alphabet::DnaBitLength;
	}
	// std::ofstream FamLengOut("pfamleng.dat");
	// FamLengOut << "#ID\tSize\tMinLeng\tMeanLeng\tMaxLeng" << std::endl;
	// for(auto & Fam : Families){
	// 	size_t MaxLeng = 0, MinLeng = Fam[0].size();
	// 	double MeanLeng = 0;
	// 	for(auto & Seq : Fam){
	// 		MinLeng = std::min(MinLeng, Seq.size());
	// 		MaxLeng = std::max(MaxLeng, Seq.size());
	// 		MeanLeng += Seq.size();
	// 	}
	// 	MeanLeng /= Fam.size();
	// 	FamLengOut << Fam.id() << "\t" << Fam.size() << "\t" << MinLeng << "\t" << size_t(MeanLeng) << "\t" << MaxLeng << std::endl;
	// }
	// FamLengOut.close();
	// std::exit(0);

	std::vector<double> SequenceCoverage(Families.size(), -1);

	std::cout << "\nCreating SpacedWordDataBase ..." << std::flush;
	init_omp();
	spacedword_db FamiliesSpacedWords(dboptions::BitLength, dboptions::WordPerSeq);
	#pragma omp parallel
	{
		spacedword_db local_db(dboptions::BitLength, dboptions::WordPerSeq);
		#pragma omp for schedule(dynamic) nowait
		for(unsigned i = 0; i < dboptions::FamilySize; i++){
			if(Families[i].size() >= dboptions::FamilyOffset){
				create_family_spacedwords(local_db, Families[i], i, PatSet, SequenceCoverage[i]);
			}
		}
		if(omp_get_num_threads() == 1){
			std::swap(local_db, FamiliesSpacedWords);
		}
		else{
			#pragma omp critical
			FamiliesSpacedWords.merge(local_db);
		}
	}
	std::cout << " Done!" << std::endl;
	if(dboptions::SplitFile == false){
		std::cout << "Writing database to file '" << dboptions::OutSwDb << ".swdb' ..." << std::flush;
	}
	else{
		std::cout << "Writing database to multiple files into '" << dboptions::OutSwDb << "-DataBase/' ..." << std::flush;
	}
	FamiliesSpacedWords.to_file(dboptions::OutSwDb, dboptions::SplitFile);
	std::cout << " Done!\n" << std::endl;
	if(dboptions::SeqCov){
		get_seqcov(Families, SequenceCoverage);
	}
	if(dboptions::Overlap && FamiliesSpacedWords.size() != 0){
		get_overlap(FamiliesSpacedWords);
	}
}

void create_family_spacedwords(spacedword_db & SpacedWordDatabase, sequence_family & Family, unsigned FamID, patternset & PatSet, double & SeqCover){
	SpacedWordDatabase.add_family_id(Family.id(), FamID);
	if(PatSet.size() != 0){
		//TODO READ PATTERN FILE
	}
	else{
		patternset UniqSet;//(dboptions::PatternNumber, dboptions::PatternWeight, dboptions::PatternWeight, dboptions::PatternMaxDC, dboptions::PatternMinDC);
		for(unsigned i = 0; i < dboptions::PatternNumber; i++){
			pattern Pat(dboptions::PatternMaxDC, dboptions::PatternWeight);
			std::vector<spacedword_family> SpacedWordBuckets;
			SpacedWordBuckets = create_spacedword_buckets(UniqSet, Pat, Family, FamID, SeqCover);
			if(SpacedWordBuckets.size() != 0){
				UniqSet.push_back(Pat);
				SpacedWordDatabase.push_back(SpacedWordBuckets, Pat);
			}
		}
	}
}

std::vector<spacedword_family> create_spacedword_buckets(patternset & UniqSet, pattern & Pat, sequence_family & Family, unsigned FamID, double & SeqCover){
	if(dboptions::BitLength == 5){
		protein_spacedword NOUSE;
		return make_buckets(NOUSE, UniqSet, Pat, Family, FamID, SeqCover);
	}
	dna_spacedword NOUSE;
	return make_buckets(NOUSE, UniqSet, Pat, Family, FamID, SeqCover);
}

void get_seqcov(std::vector<sequence_family> & Families, std::vector<double> & SequenceCoverage){
	std::cout << "Writing results into file '" << dboptions::OutSwDb << ".sqcv' ... " << std::flush; 
	std::ofstream Output(dboptions::OutSwDb + ".sqcv");
	for(size_t i = 0; i < SequenceCoverage.size(); i++){
		if(SequenceCoverage[i] >= 0){
			Output << Families[i].id() << "\t" << SequenceCoverage[i] << std::endl;
		}
	}
	Output.close();
	std::cout << "Done!\n" << std::endl;
}

void get_overlap(spacedword_db & SwDataBase){
	std::cout << "\nCalculating Overlap among Database ... " << std::flush;
	std::vector<size_t> FamSwCtr(SwDataBase.families_size(),0);
	std::vector< std::vector<double> > Overlap(SwDataBase.families_size(), std::vector<double>(SwDataBase.families_size(), 0));
	for(size_t p = 0; p < SwDataBase.size(); p++){
		std::vector<spacedword_family> FamSw = SwDataBase.spaced_families(p);
		for(auto & Sw : FamSw){
			size_t PfSize = Sw.size() - 1;
			for(size_t i = 0; i < Sw.size()-1; i++){
				FamSwCtr[Sw[i].id()]++;
				for(size_t j = i+1; j < Sw.size(); j++){
					Overlap[Sw[i].id()][Sw[j].id()] += PfSize;
					Overlap[Sw[j].id()][Sw[i].id()] += PfSize;
				}
			}
			FamSwCtr[Sw[Sw.size()-1].id()]++;
		}
	}

	for(size_t i = 0; i < Overlap.size(); i++){
		for(size_t j = 0; j < Overlap[i].size(); j++){
			if(FamSwCtr[i] != 0){
				Overlap[i][j] /= FamSwCtr[i];
			}
		}
	}

	std::vector<double> Intervals = {0, 0.001, 0.01, 0.1, 1};
	std::vector< std::vector<double> > OvlIntervals(SwDataBase.families_size(), std::vector<double>(Intervals.size(),0));
	for(size_t i = 0; i < Overlap.size(); i++){
		for(size_t j = 0; j < Overlap[i].size(); j++){
			for(unsigned k = 0; k < Intervals.size(); k++){
				if(Overlap[i][j] <= Intervals[k] && i != j){
					OvlIntervals[i][k]++;
				}
			}
		}
		if(SwDataBase.families_size() > 1){
			for(auto & Val : OvlIntervals[i]){
				Val /= SwDataBase.families_size()-1;
			}
		}
	}
	std::cout << "Done!" << std::endl;

	std::cout << "Writing results into file '" << dboptions::OutSwDb << ".ovlp' ... " << std::flush; 
	std::ofstream Output(dboptions::OutSwDb + ".ovlp");
	Output << "#'PFamID'";
	for(auto & Val : Intervals){
		Output << "\t'" << Val << "'";
	}
	Output << std::endl;
	size_t Pfam = 0;
	for(auto Ovl : OvlIntervals){
		Output << SwDataBase.id_fam_name(Pfam++);
		for(size_t i = 0; i < Intervals.size(); i++){
			Output << "\t" << Ovl[i];
		}
		Output << std::endl;
	}
	Output.close();
	std::cout << "Done!\n" << std::endl;
}


void clean_line(){
	std::vector<char> EmptyVector(80,' ');
	std::cout << "\r" << std::string(EmptyVector.begin(), EmptyVector.end()) << "\r";
}


void init_omp(){
	#ifdef _OPENMP
	omp_set_dynamic(0);
	omp_set_num_threads(dboptions::ThreadNumber);
	#endif	
}
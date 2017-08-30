#include "searchdatabase.hpp"

std::vector<size_t> FamSizeUsed;

void search_database(){
	if(dsoptions::InSwDb.size() == 0 || dsoptions::InSeq.size() == 0){
		std::cerr << "FATAL ERROR!\nINPUT FILES missing!" << std::endl;
		std::exit(-1);
	}
	std::vector<sequence_family> Families = sequence_family::families_file(dsoptions::InSeq);
	dsoptions::FamilySize = Families.size();
	if(Families.size() < 2){
		if(Families[0].size() < 2){
			std::cout << Families[0].size() << " Sequence was read. Inputtype: ";
		}
		else{
			std::cout << Families[0].size() << " Sequences were read. Inputtype: ";
		}
	}
	else if(Families.size() == 0){
		std::cerr << "FATAL ERROR!\nNo sequence input!" << std::endl;
		std::exit(-1);
	}
	else{
		std::cout << Families.size() << " Families were read. Inputtype: ";
	}
	if(Families[0].is_protein()){
		std::cout << "Protein " << std::endl;
	}
	else{
		std::cout << "DNA" << std::endl;
	}
	std::cout << std::endl;
	std::cout << "Reading database '" << dsoptions::InSwDb << "' ... " << std::flush;
	spacedword_db SwDataBase(dsoptions::InSwDb);
	std::cout << "\n--> Done!" << std::endl;
	dsoptions::BitLength = SwDataBase.bit_length();
	if(dsoptions::WordPerSeq < 0){
		dsoptions::WordPerSeq = SwDataBase.words_per_seq();
	}
	unsigned SeqType = 0;
	if(Families[0].is_protein()){
		SeqType = alphabet::ProteinBitLength;
		if(dsoptions::BlockSize == -1){
			dsoptions::BlockSize = 1167; //DNA default BlockSize divided by 3, codon size.
		}
	}
	else{
		SeqType = alphabet::DnaBitLength;
		if(dsoptions::BlockSize == -1){
			dsoptions::BlockSize = 3500; //Typical DNA gene size.
		}
	}
	std::cout << "\n\tBlockSize = " << dsoptions::BlockSize << std::endl;
	std::cout << "\tThreshold = " << dsoptions::HitThreshold << std::endl;
	std::cout << "\tWords per Seq = " << dsoptions::WordPerSeq << std::endl << std::endl;

	if(dsoptions::BitLength < SeqType){
		std::cerr << "TYPE ERROR!\nDataBase Type is DNA, Sequences Type is Protein! Cannot convert SpacedWord to Protein!" << std::endl;
		std::exit(-1);
	}
	if(dsoptions::BitLength == alphabet::ProteinBitLength && SeqType == alphabet::DnaBitLength){
		std::cout << "DataBase Type is Protein, Sequences Type is DNA, translating sequences to protein ... " << std::flush;
			for(auto & Fam : Families){
				Fam.translate();
			}
		std::cout << "\n--> Done!" << std::endl;
		dsoptions::Translated = true;
	}

	FamSizeUsed = std::vector<size_t>(Families.size(),0);

	init_omp();
	std::vector<spacedhit> DataBaseHit;
	std::cout << "Matching spaced words database to input sequences  ... " << std::endl;
	for(unsigned i = 0; i < Families.size(); i++){
		if(Families[i].size() > 1){
			std::cout << "\rProcessing Family " << i+1 << "/" << Families.size() <<  " " << Families[i].id() << std::flush;
			std::vector<spacedhit> FamilyHit =  family_hit(Families[i], i, SwDataBase);
			DataBaseHit.insert(DataBaseHit.end(), FamilyHit.begin(), FamilyHit.end());
		}
	}
	std::cout << "\r" << std::string(80,' ') << "\r--> Done!" << std::endl;
	//TODO: SPLIT OUTPUT, MULTIPLE FILES

	eval_classify(SwDataBase, Families, DataBaseHit);

	// std::cout << "\nWriting results into file '" << dsoptions::OutDetect << ".swds' ... " << std::flush; 
	// Output.open(dsoptions::OutDetect+".swds");
	// Output << "#FamName\tSeqName\tORF\tMatch DbID\tPosition\tScore" << std::endl;
	// for(auto SwHit : DataBaseHit){
	// 	unsigned ORF = SwHit.orf();
	// 	std::string OrfPrint = "+";
	// 	if(ORF%2 == 1){
	// 		OrfPrint[0] = '-';
	// 	}
	// 	OrfPrint.push_back((ORF/2)+1+'0');
	// 	Output << Families[SwHit.fam_id()].id() << "\t" << Families[SwHit.fam_id()][SwHit.seq_num()].name() << "\t" << OrfPrint << "\t" << SwDataBase.id_fam_name(SwHit.db_fam_id()) << "\t" << SwHit.position() << "\t" << SwHit.score() << std::endl;
	// }
	// Output.close();
	// std::cout << "\n--> Done!" << std::endl;
}



std::vector<spacedhit> family_hit(sequence_family & Family, unsigned FamID, spacedword_db & SwDataBase){
	std::vector<spacedhit> PatternHits;
	#pragma omp parallel
	{
		std::vector<spacedhit> SequenceHits;
		#pragma omp for schedule(dynamic) nowait
		for(size_t i = 0; i < SwDataBase.size(); i++){
			pattern Pat = SwDataBase.pattern_families(i);
			std::vector<spacedword_family> DbSw = SwDataBase.spaced_families(i);
			std::sort(DbSw.begin(), DbSw.end());
			std::vector<spacedhit> SeqHits = sequence_hit(Family, FamID, Pat, DbSw);
			SequenceHits.insert(SequenceHits.end(), SeqHits.begin(), SeqHits.end());
		}
		if(omp_get_num_threads() == 1){
			std::swap(PatternHits, SequenceHits);
		}
		else{
			#pragma omp critical
			PatternHits.insert(PatternHits.end(), SequenceHits.begin(), SequenceHits.end());
		}
	}
	std::sort(PatternHits.begin(), PatternHits.end());
	return PatternHits;
}

std::vector<spacedhit> sequence_hit(sequence_family & Family, unsigned FamID, pattern & Pat, std::vector<spacedword_family> & DbSw){
	if(dsoptions::BitLength == alphabet::ProteinBitLength){
		protein_spacedword NOUSE;
		return sw_sequence_hit(NOUSE, Family, FamID, Pat, DbSw, FamSizeUsed[FamID]);
	}
	dna_spacedword NOUSE;
	return sw_sequence_hit(NOUSE, Family, FamID, Pat, DbSw, FamSizeUsed[FamID]);
}


void init_omp(){
	#ifdef _OPENMP
	omp_set_dynamic(0);
	omp_set_num_threads(dsoptions::ThreadNumber);
	#endif	
}


void eval_classify(spacedword_db & SwDataBase, std::vector<sequence_family> & Families, std::vector<spacedhit> & DataBaseHit){
	std::map<size_t, size_t> DbIdToFamId;
	for(size_t i = 0; i < SwDataBase.families_size(); i++){
		std::string FamName = SwDataBase.id_fam_name(i);
		auto Pos = std::find_if(Families.begin(), Families.end(), [FamName](const sequence_family & Fam){
			return FamName == Fam.id();
		});
		if(Pos != Families.end()){
			DbIdToFamId.insert(std::pair<size_t,size_t>(i,std::distance(Families.begin(), Pos)));
		}
	}

	std::vector<int64_t> FamTruePos(Families.size(), 0), FamFalsePos(Families.size(), 0), FamTrueNeg(Families.size(), 0), FamFalseNeg(Families.size(), 0);

	auto BucketStart = DataBaseHit.begin(), BucketIter = BucketStart;
	bool NotEnd = true;
	bool Border = false;
	size_t SwTotal = 0;
	while(NotEnd){
		if(BucketIter == DataBaseHit.end()){
			NotEnd = false;
			Border = true;
		}
		else{
			if(BucketIter->fam_id() != BucketStart->fam_id()){
				Border = true;
			}
			else if(BucketIter->seq_num() != BucketStart->seq_num()){
				Border = true;
			}
		}
		if(Border){
			auto MaxIter = BucketStart;
			for(auto BIter = BucketStart; BIter != BucketIter; BIter++){
				if(MaxIter->score() < BIter->score()){
					MaxIter = BIter;
				}
			}
			SwTotal++;
//			std::cout << SwDataBase.id_fam_name(MaxIter->db_fam_id()) << ": " << Families[MaxIter->fam_id()].id() << " == " << Families[DbIdToFamId[MaxIter->db_fam_id()]].id() << std::endl;
			if(Families[MaxIter->fam_id()].id() == SwDataBase.id_fam_name(MaxIter->db_fam_id())){
//				std::cout << Families[MaxIter->fam_id()].id() << " == " << SwDataBase.id_fam_name(MaxIter->db_fam_id()) << std::endl;
				FamTruePos[MaxIter->fam_id()]++;
			}
			else{
//				std::cout << Families[MaxIter->fam_id()].id() << " != " << SwDataBase.id_fam_name(MaxIter->db_fam_id()) << std::endl;
				size_t SwFamId = DbIdToFamId[MaxIter->db_fam_id()];
//				size_t FamId = MaxIter->fam_id();
				FamFalsePos[SwFamId]++;
//				FamFalseNeg[FamId]++;
			}
			Border = false;
			BucketStart = BucketIter;
		}
		BucketIter++;
	}

	for(size_t i = 0; i < FamTrueNeg.size(); i++){
//		FamFalseNeg[i] += FamSizeUsed[i] - FamTruePos[i];
		FamFalseNeg[i] = FamSizeUsed[i] - FamTruePos[i];
		FamTrueNeg[i] = SwTotal - FamTruePos[i] - FamFalsePos[i] - FamFalseNeg[i];
	}

	std::ofstream Output(dsoptions::OutDetect+".eval");
	Output << "#PfamID\tSensitivity\tSpecifity\tPrecision\tTP\t\tFP\t\tTN\t\tFN\t\tSeqSize" << std::endl;
	size_t TruePositiv = 0, FalsePositiv = 0, TrueNegative = 0, FalseNegative = 0;
	Output << std::fixed << std::setprecision(7);
	for(size_t i = 0; i < Families.size(); i++){
		Output << Families[i].id();
		if(FamTruePos[i] + FamFalseNeg[i] > 0){
			Output << "\t" << FamTruePos[i]/((double)FamTruePos[i]+FamFalseNeg[i]);
		}
		else{
			Output << "\t- "; 
		}
		if(FamTrueNeg[i] + FamFalsePos[i] > 0){
			Output << "\t" << FamTrueNeg[i]/((double)FamTrueNeg[i]+FamFalsePos[i]);
		}
		else{
			Output << "\t- ";
		}
		if(FamTruePos[i] + FamFalsePos[i] > 0){
			Output << "\t" << FamTruePos[i]/((double)FamTruePos[i]+FamFalsePos[i]);
		}
		else{
			Output << "\t- ";
		}
		Output << "\t" << FamTruePos[i] << "\t\t" << FamFalsePos[i] << "\t\t" << FamTrueNeg[i] << "\t\t" << FamFalseNeg[i] << "\t\t" << FamSizeUsed[i] << std::endl;
		TruePositiv += FamTruePos[i];
		FalsePositiv += FamFalsePos[i];
		FalseNegative += FamFalseNeg[i];
		TrueNegative += FamTrueNeg[i];
	}
	Output << "\n\nTotal\tSensitivity\tSpecifity\tPrecision" << std::endl;
	Output << "Pfam";
	if(TruePositiv + FalsePositiv > 0){
		Output << "\t" << TruePositiv/((double)TruePositiv+FalseNegative);
	}
	else{
		Output << "\t- "; 
	}
	if(TrueNegative + FalsePositiv > 0){
		Output << "\t" << TrueNegative/((double)TrueNegative+FalsePositiv);
	}
	else{
		Output << "\t- ";
	}
	if(TruePositiv + FalsePositiv > 0){
		Output << "\t" << TruePositiv/((double)TruePositiv+FalsePositiv);
	}
	else{
		Output << "\t- ";
	}
	Output << std::endl << std::endl;
	Output << "Sens = TP/(TP+FN)" << std::endl;
	Output << "Spec = TN/(TN+FP)" << std::endl;
	Output << "Prec = TP/(TP+FP)" << std::endl;

	Output.close();
}
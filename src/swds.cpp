#include <iostream>
#include <string>
#include "dsoptions.hpp"
#include "searchdatabase.hpp"


void SecureMessage(std::string && errmsg, std::string && opt);

int main(int argc, char* argv[]){
	bool OutSet = false;
	if(argc < 3){
		SecureMessage("param",std::string(argv[0]));
		return 0;
	}
	else{
		int ctr;
		for(ctr = 1; ctr < argc-2; ctr++){
			std::string tmp = argv[ctr];
			switch(tmp[1]){
				case 'B':
				case 'b':
					if (ctr < argc - 3) {
						ctr++;
						dsoptions::BlockSize = std::atoi(argv[ctr]);
					}
					break;
				case 'H':
				case 'h':
					if (ctr < argc - 3) {
						ctr++;
						dsoptions::HitThreshold = std::atof(argv[ctr]);
					}
					break;
				case 'R':
				case 'r':
					dsoptions::AlphaReduction = true;
					break;
				case 'S':
				case 's':
					if (ctr < argc - 3) {
						ctr++;
						dsoptions::WordPerSeq = std::atoi(argv[ctr]);
					}
					break;
				case 'T':
				case 't':
					if (ctr < argc - 3) {
						ctr++;
						dsoptions::ThreadNumber = std::atoi(argv[ctr]);
					}
					break;
				case '-':
					if(tmp == "--help"){
						SecureMessage("param",std::string(argv[0]));
						return 0;
					}
					else if(tmp == "--version"){
						SecureMessage("version","");
						return 0;
					}
					else if(tmp == "--outfile"){
						if (ctr < argc - 3) {
							OutSet = true;
							ctr++;
							dsoptions::OutDetect = std::string(argv[ctr]);
						}
					}
					break;
				default:
					SecureMessage("arg",std::string(argv[ctr]));
			}
		}
		if(ctr == argc-2){
			dsoptions::InSwDb = std::string(argv[ctr++]);
			dsoptions::InSeq = std::string(argv[ctr++]);
			if(OutSet == false){
				auto Pos = dsoptions::InSeq.find_last_of(".");
				if(Pos != std::string::npos){
					dsoptions::OutDetect = dsoptions::InSeq.substr(0,Pos);
				}
				else{
					dsoptions::OutDetect = dsoptions::InSeq;
				}
			}
			else{
				auto Pos = dsoptions::OutDetect.find_last_of(".");
				if(Pos != std::string::npos){
					dsoptions::OutDetect = dsoptions::OutDetect.substr(0,Pos);
				}
			}
			search_database();
		}
		else{
			std::cerr << "FATAL ERROR!\nINPUT FILES missing!" << std::endl;
		}
	}
	return 0;
}

/*===Functions===============================================================*/
void SecureMessage(std::string && errmsg, std::string && opt){
	if(errmsg == "arg"){
		std::cerr << "Argument '" << opt << "' is unknown!" << std::endl;
		return;
	}
	else if (errmsg == "param") {
		std::cout << "Usage:\n\t" << opt << " [options] <SWDB FILE> <SEQUENCE FILE>\n" << std::endl;
		std::cout << "Options:\n" << std::endl;
		std::cout << "\t\t -h <FLOAT>: \t\t Minimal threshold hit score for detection, default h = 0.\n" << std::endl;
		std::cout << "\t\t -b <INT>: \t\t Window size used for detection, default b = 3500 | 1167 (DNA|Protein).\n" << std::endl;
		std::cout << "\t\t -s <INT>: \t\t Number of minimal required spaced words per sequence, default s = 1.\n" << std::endl;
		//std::cout << "\t\t -r: \t\t\t Activate alphabet reduction for protein sequences.\n" << std::endl;
		std::cout << "\t\t -t <INT>: \t\t Number of CPU-Threads to be used, default: t = 1.\n" << std::endl;
		std::cout << "\t\t --outfile <FILE>: \t Save detection to <FILE-prefix>.swds instead of <SEQUENCE FILE-prefix>.swds.\n" << std::endl;
		std::cout << "\t=== Additional Parameters ====" << std::endl;
		std::cout << "\t\t --version: \t\t Print the program version.\n" << std::endl;
		std::cout << "\t\t --help: \t\t Print this help.\n" << std::endl;
		return;
	}
	else if (errmsg == "version") {
		std::cout << "Versionsnummer" << std::endl;
		return;
	}
}
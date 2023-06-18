#ifndef GOPTIONS_HPP_
#define GOPTIONS_HPP_

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <string>
#include <iostream>
#include <map>
#include <stdexcept>
#include <fstream>
using namespace std;

const std::string version("1.0");

class OptionsExitsProgram : public std::exception
{};

class Goptions {
public:
	// The constructor sets up all the various options that will be parsed
	Goptions() {
		SetOptions();
	}

	// Parse options runs through the heirarchy doing all the parsing
	void ParseOptions(int argc, char const *argv[]) {
		ParseCommandLine(argc, argv);
		NotDefined();
		CheckForHelp();
		CheckForVersion();
	}

	// Below is the interface to access the data, once ParseOptions has been run
	std::string Path() {
		return results["path"].as<std::string>();
	}

	std::string Verbosity() {
		return results["verbosity"].as<std::string>();
	}

	std::vector<std::string> IncludePath() {
		if (results.count("include-path")) {
			return results["include-path"].as<std::vector<std::string> >();
		}
		return std::vector<std::string>();
	}
/*
   std::string MasterFile() {
      if (results.count("master-file")) {
         return results["master-file"].as<std::string>();
      }
      return "";
   }
   std::vector<std::string> Files() {
      if (results.count("file")) {
         return results["file"].as<std::vector<std::string> >();
      }
      return std::vector<std::string>();
   }
*/

	bool GUI() {
		if (results["run-gui"].as<bool>()) {
			return true;
		}
		return false;
	}

//generic
	string GetGenericGenoFile() {
		return results["bfile"].as<string>();
	}

	string GetGenericOutFile() {
		return results["out"].as<string>();
	}

	string GetGenericExtractFile() {
		return generic_extract_file;
	}

	string GetGenericKeepFile() {
		return generic_keep_file;
	}

	string GetGenericForceRefAlleleFile() {
		return generic_ref_allele_file;
	}

	string GetGenericSNPTagFile() {
		return generic_snp_tag_file;
	}

	int GetGenericBlockSNPs() {
		return generic_block_snps;
	}

	float GetGenericBlockKB() {
		return generic_block_kb;
	}

	int GetGenericBlockNum() {
		return generic_block_num;
	}

	string GetGenericPhenoFile() {
//      if (results.count("pheno")) {
//         return results["pheno"].as<string>();
//      } else {
			return generic_pheno_file;
//      }
	}

	vector<int> GetGenericPhenoNum() {
		if (results.count("pheno-num")) {
			return results["pheno-num"].as<vector<int>>();
		} else {
			return generic_pheno_index;
		}
	}

	string GetGenericCovarFile() {
		if (results.count("covar")) {
			return results["covar"].as<string>();
		} else {
			return generic_covar_file;
		}
	}

	vector<int> GetGenericCovarNum() {
		if (results.count("covar-num")) {
			return results["covar-num"].as<vector<int>>();
		} else {
			return generic_covar_index;
		}
	}

	vector<string> GetGenericChr() {
		if (results.count("chr")) {
			return results["chr"].as<vector<string> >();
		} else {
			return generic_chr;
      }
   }

	bool IsGenericMemoryEfficient() {
		return results["mem"].as<bool>();
	}

	int GetGenericImput() {
		if (generic_imput != 0 && generic_imput != 1) {
			cerr << "unknow value for adjustment for imputation (0 for observed genotypes, or 1 for hwe if outbred or p vs (1-p) for inbred)." << endl;
			exit(1);      
		}
		return generic_imput;
	}

	bool IsGenericSubstractMean() {
		return true;
//      return !results["miss"].as<bool>();
//      return results["substract-mean"].as<bool>();
	}

	bool IsGenericNoMailman() {
		return results["no-mailman"].as<bool>();
	}

	bool IsGenericFastMode() {
		return !results["no-mailman"].as<bool>();
	}

	bool IsGenericInbred() {
		return inbred;
	}

	bool IsGenericMissing() {
		return results["miss"].as<bool>();
	}

	bool IsGenericVarNorm() {
		return results["var-norm"].as<bool>();
	}

	bool IsGenericQImput() {
		return results["q-imput"].as<bool>();
	}

	bool IsGenericDebug() {
		return results["debug"].as<bool>();
	}

	int GetGenericEigenvecNumber() {
		return results["evec"].as<int>();
	}

	int GetGenericMailmanBlockSize() {
		int k = (int) ceil(GetGenericEigenvecNumber()/10.0) * 10;

		if (CheckRandHEMasterOption() ) {
			k = GetRandHEB0() > GetRandHEB1() ? GetRandHEB0() : GetRandHEB1();
		} else if(CheckMeMasterOption() || CheckXLDMasterOption()) {
			k = GetGenericIteration();
		} else if (CheckEncMasterOption()) {
			k = GetEncK();
		}
		return k;
	}

	int GetGenericThreads() {
		return results["threads"].as<int>();
	}

	int GetGenericSeed() {
		return results["seed"].as<int>();
	}

	int GetGenericIteration() {
		return results["iter"].as<int>();
	}

	double GetGenericMAF() {
		return generic_maf;
	}

	double GetGenericLMISS() {
		return generic_lmiss;
	}

	double GetGenericHWE() {
		return generic_hwe;
	}

	bool GetGenericRemoveZeroGstd() {
		return generic_remove_zero_gstd;
	}

	int GetGenericAdjVar() {
		if (generic_adj_var != 0 && generic_adj_var != 1) {
			cerr << "unknow value for adjustment for locus standardization (0 or 1)" << endl;
			exit(1);
		}
		return generic_adj_var;
	}

//make-bed
	bool CheckMakebedMasterOption() {
		return results["make-bed"].as<bool>();
	}

//freq
	bool CheckFreqMasterOption() {
		return results["freq"].as<bool>();
	}

//propc
	bool CheckPropcMasterOption() {
		return results["propc"].as<bool>();
	}

	bool IsPropcAccuracy() {
		return results["accu"].as<bool>();
	}

	int GetPropcMaxIteration() {
		if(results["evec"].as<int>()!=2 && results["max-it"].as<int>()==4){
			return results["evec"].as<int>()+2;
		} else {
         return results["max-it"].as<int>();
		}
	}

	int GetPropcAcceleratedEM() {
		return results["accel-em"].as<int>();
	}

	double GetPropcConvergenceLimit() {
		return results["conv-limit"].as<double>();
	}

//eigengwas
	bool CheckEigenGWASMasterOption() {
		return results["eigengwas"].as<bool>();
	}
	string GetEigenGWASreadfile() {
		return eigengwas_read_v;
	}

//rand-he
	bool CheckRandHEMasterOption() {
		return results["rand-he"].as<bool>();
	}

	double GetRandHERoundErrOption() {
		if (round_err <= 0 || round_err >= 1) {
			cerr << "incorrect round error for rand-he: " << round_err << endl;
			exit(1);
		}
		return round_err;
	}

	int GetRandHEB0() {
		if (randhe_B0 < 0) {
			cerr << "incorrect randhe_B0 for rand-he: " << randhe_B0 << endl;
			exit(1);
		}
		return randhe_B0;
	}

	int GetRandHEB1() {
		if (randhe_B1 < 0 || randhe_B1 > randhe_B0) {
			cerr << "incorrect randhe_B1 for rand-he: " << randhe_B0 << endl;
			cerr << "randhe_B1 should 0< b1 < b0" << endl;
			exit(1);
		}
		return randhe_B1;
	}

	int GetRandHEMaxIt() {
		if (randhe_maxit < 0) {
			cerr << "incorect rand-he-max-it for rand-he " << randhe_maxit << endl;
			exit(1);
		}
		return randhe_maxit;
	}
//enc
	int GetEncK() {
		return results["enc-k"].as<int>();
	}

	bool CheckEncMasterOption() {
		return results["enc"].as<bool>();
	}

//encreg
	bool CheckEncRegMasterOption() {
		return enc_reg_switch;
	}

	string GetEncRegF1() {
		return enc_reg_f1;
	}

	string GetEncRegF2() {
		return enc_reg_f2;
	}

//cld
	bool CheckXLDMasterOption() {
		return results["xld"].as<bool>();
	}

	int GetXLDalg() {
		if (xld_alg != 0 && xld_alg != 1) {
			cerr << "Unknow value " << xld_alg << " for xld-alg (0 or 1 only)" << endl;
			exit(0);
		}
		return results["xld-alg"].as<int>();
	}

	bool GetXLDStage() {
		return xld_stage;
	}

	string GetXLDList() {
		return xld_list;
	}
//me
	bool CheckMeMasterOption() {
		return me_switch;
	}

	int GetMeLeap() {
		return me_leap;
	}

	double GetMeStop() {
		return me_stop;
	}
private:

	void NotDefined() {
		if (IsGenericMissing() && !IsGenericFastMode()) {
//   	if (missing && !fast_mode) {
			cerr << "Missing version works only with mailman i.e. fast mode\n EXITING..." << endl;
			exit(-1);
		}
		if (IsGenericFastMode() && IsGenericMemoryEfficient()) {
//	   if (fast_mode && memory_efficient) {
			cerr << "Memory effecient version for mailman EM not yet implemented" << endl;
			cerr << "Ignoring Memory effecient Flag" << endl;
		}
		if (IsGenericMissing() && IsGenericVarNorm()) {
//	   if (missing && var_normalize) {
			cout << "Missing version works only without variance normalization\n EXITING..." << endl;
			exit(-1);
		}
	}

	void SetOptions() {
		SetGenericOptions();
		SetMakebedOptions();
		SetFreqOptions();
		SetPropcOptions();
		SetEigenGWASOptions();
		SetEncOptions();
		SetEncRegOptions();
		SetRandHEOptions();
		SetXLDOptions();
		SetMeOptions();
		SetCommonOptions();
		SetEnvMapping();
	}

	void SetGenericOptions() {
		genericOpts.add_options()
			("help", "produce help message.")
			("version,v", "print version string.")

			("debug", po::bool_switch()->default_value(false), "debug mode.")
//			("txt", po::bool_switch()->default_value(false), "text pedigree files (default false).")
			("mem", po::bool_switch()->default_value(false), "the flag states whether to use a memory effecient version for the EM algorithm or not. The memory efficient version is a little slow than the not efficient version (default: false).")
			("var-norm", po::bool_switch()->default_value(true), "normalization for mailman.")
			("imput", po::value(&generic_imput), "imputation 0 for obs genotypes, and 1 for hwe if outbred and p vs (1-p) for inbred.")
			("no-mailman", po::bool_switch()->default_value(false), "no mailman (default, false).")
			("miss", po::bool_switch()->default_value(false), "no missing (default, false, when true the missing genotypes will be imputed.")
 //			("substract-mean", po::bool_switch()->default_value(true), "substract-mean for mailman (default, true).")

			("bfile", po::value<string>(&generic_bGeno_file), "root of plink binary pedigree files.")
			("out", po::value<string>()->default_value("out"), "root for output files.")
			("pheno", po::value<string>(&generic_pheno_file), "pheno file.")
			("pheno-num", po::value<vector<int> >(&generic_pheno_index)->multitoken(), "pheno index.")
			("covar", po::value<string>(&generic_covar_file), "covar file.")
			("covar-num", po::value<vector<int> >(&generic_covar_index)->multitoken(), "covar index.")
			("q-imput", po::bool_switch()->default_value(true), "quan-imputation.")

			("threads", po::value<int>(&generic_thread)->default_value(1), "thread number.")
			("seed", po::value<int>(&generic_seed)->default_value(2021), "seed for generating random numbers.")
			("evec", po::value<int>()->default_value(2), "eigenvector to estimate.")
			("iter", po::value<int>()->default_value(5), "iteration for randomization.")

			("config,c", po::value<string>(&config_file), "config files to parse (always parses default.cfg).")

			("extract", po::value<string>(&generic_extract_file), "extract snp file.")
			("chr", po::value<vector<string> >(&generic_chr)->multitoken(), "chromosome index.")
			("keep", po::value<string>(&generic_keep_file), "keep individual file.")
			("snp-tag", po::value<string>(&generic_snp_tag_file), "snp tag file.")
			("block-snps", po::value<int>(&generic_block_snps)->default_value(-1), "snp numbers for each block. 1000 by default.")
			("block-kb", po::value<float>(&generic_block_kb)->default_value(-1.0), "kb for each block. 50.0 by default.")
			("block-num", po::value<int>(&generic_block_num)->default_value(-1), "snp numbers for each block. 100 by default.")
			("force-ref-allele", po::value<string>(&generic_ref_allele_file), "force reference allele file")

			("inbred", po::bool_switch(&inbred)->default_value(false), "inbred population.")
			("maf", po::value<double>(&generic_maf)->default_value(-1.0), "maf cutoff.")
			("lmiss", po::value<double>(&generic_lmiss)->default_value(-1.0), "locus missing rate.")
			("hwe", po::value<double>(&generic_hwe)->default_value(-1.0), "hardy-weinberg cutoff.")
			("remove-zero-gstd", po::bool_switch(&generic_remove_zero_gstd)->default_value(false), "keep zero std loci.")
			("adj-var", po::value<int>(&generic_adj_var)->default_value(0), "0 for observed variance, 1 for expected variance 2pq (or 4pq if inbred).");
		;
	}

	void SetMakebedOptions() {
		makebedOpts.add_options()
			("make-bed", po::bool_switch(&makebed_switch)->default_value(false), "master option for make binary files.")
		;
	}

	void SetFreqOptions() {
		freqOpts.add_options()
			("freq", po::bool_switch()->default_value(false), "frequency option for make binary files.")
		;
	}

	void SetPropcOptions() {
		propcOpts.add_options()
			("propc", po::bool_switch(&propc_switch)->default_value(false), "master option for propc.")
			("accu", po::bool_switch()->default_value(false), "output the likelihood computation as a function of iterations.")
			("max-it", po::value<int>(&propc_max_it)->default_value(4), "maximun iteration for propc.")
			("conv-limit", po::value<double>()->default_value(0.05), "The value of the threshold telling the algorithm that it has converged (the value of -1, meaning no auto termination condition.")
			("accel-em", po::value<int>()->default_value(0), "The flag stating whether to use accelerated EM or not (default: 0, and can be 1 or 2).")
		;
	}

	void SetEigenGWASOptions() {
		eigengwasOpts.add_options()
			("eigengwas", po::bool_switch(&eigengwas_switch)->default_value(false), "master option for eigengwas.")
			("read-v", po::value<string>(&eigengwas_read_v), "read eigenvectors and eigenvalues.")
		;
	}

	void SetEncOptions() {
		encOpts.add_options()
			("enc", po::bool_switch(&enc_switch)->default_value(false), "master option for enc.")
			("enc-k", po::value<int>(&enc_k)->default_value(10), "iteration for generating tags.")
		;
	}

	void SetEncRegOptions() {
		encregOpts.add_options()
			("enc-reg", po::bool_switch(&enc_reg_switch)->default_value(false), "master option for enc-reg.")
			("enc-f1", po::value<string>(&enc_reg_f1), "encrypted file 1.")
			("enc-f2", po::value<string>(&enc_reg_f2), "encrypted file 2.")
		;
	}

	void SetRandHEOptions() {
		randheOpts.add_options()
			("rand-he", po::bool_switch(&randhe_switch)->default_value(false), "master option for randomized Haseman-Elston regression.")
			("round-err", po::value<double>(&round_err)->default_value(0.1), "rounding error for Haseman-Elston regression.")
			("rand-he-b0", po::value<int>(&randhe_B0)->default_value(10), "B0 for Haseman-Elston regression.")
			("rand-he-b1", po::value<int>(&randhe_B1)->default_value(5), "B1 for Haseman-Elston regression.")
			("rand-he-max-it", po::value<int>(&randhe_maxit)->default_value(100), "Max it for Haseman-Elston regression.")
		;
	}

	void SetXLDOptions() {
		xldOpts.add_options()
			("xld", po::bool_switch(&xld_switch)->default_value(false), "master option for cross chromosome linkage disequilibrium.")
			("xld-alg", po::value<int>(&xld_alg)->default_value(0), "iteration for generating tags. 0 for determistic, 1 for randomization.")
			("xld-stage", po::bool_switch(&xld_stage)->default_value(false), "print xxz under xld-alg 1.")
			("xld-list", po::value<string>(&xld_list), "xld-list for randomizated xld algorithm.")
		;
	}

	void SetMeOptions() {
		meOpts.add_options()
			("me", po::bool_switch(&me_switch)->default_value(false), "master option for me.")
			("me-leap", po::value<int>(&me_leap)->default_value(-1), "me step for iteration.")
			("me-stop", po::value<double>(&me_stop)->default_value(-1.0), "me stop threshold for iteration.")
		;
	}

	void SetCommonOptions() {
		common_options.add_options()
			("path", po::value<std::string>()->default_value(""),
				"the execution path to use (imports from environment if not specified)")
			("verbosity", po::value<std::string>()->default_value("INFO"),
				"set verbosity: DEBUG, INFO, WARN, ERROR, FATAL")
			("include-path,I", po::value<std::vector<std::string> >()->composing(),
				"paths to search for include files")
			("run-gui", po::bool_switch(), "start the GUI")
		;
	}

	void SetEnvMapping() {
		env_to_option["PATH"] = "path";
		env_to_option["EXAMPLE_VERBOSE"] = "verbosity";
	}

	void ParseCommandLine(int argc, char const *argv[]) {
		po::options_description cmd_opts;
		cmd_opts.add(genericOpts)
		.add(makebedOpts)
		.add(freqOpts)
		.add(propcOpts)
		.add(eigengwasOpts)
		.add(encOpts)
		.add(encregOpts)
		.add(randheOpts)
		.add(xldOpts)
		.add(meOpts)
		.add(common_options);

		store(po::command_line_parser(argc, argv).
			options(cmd_opts).run(), results);
		notify(results);
	}

	void CheckForHelp() {
		if (results.count("help")) {
			PrintHelp();
		}
	}

	void PrintHelp() {
		cout << "Program Options Example" << endl;
		cout << "Usage: example [OPTION]... MASTER-FILE [FILE]..." << endl;
		cout << "  or   example [OPTION] --run-gui" << endl;
		po::options_description help_opts;
		help_opts
		.add(genericOpts)
		.add(makebedOpts)
		.add(freqOpts)
		.add(propcOpts)
		.add(eigengwasOpts)
		.add(encOpts)
		.add(encregOpts)
		.add(randheOpts)
		.add(xldOpts)
		.add(meOpts)
		.add(common_options);
		cout << help_opts << endl;
		exit(1);
//		throw OptionsExitsProgram();
	}

	void CheckForVersion() {
		if (results.count("version")) {
		PrintVersion();
		}
	}

	void PrintVersion() {
		std::cout << "Program Options Example " << version << std::endl;
		throw OptionsExitsProgram();
	}

	std::string EnvironmentMapper(std::string env_var) {
		// ensure the env_var is all caps
		std::transform(env_var.begin(), env_var.end(), env_var.begin(), ::toupper);

		auto entry = env_to_option.find(env_var);
		if (entry != env_to_option.end()) {
			return entry->second;
		}
		return "";
	}

	po::options_description genericOpts;
	po::options_description makebedOpts;
	po::options_description freqOpts;
	
	po::options_description propcOpts;
	po::options_description eigengwasOpts;
	po::options_description encOpts;
	po::options_description encregOpts;
	po::options_description randheOpts;
	po::options_description xldOpts;
	po::options_description meOpts;

	std::map<std::string, std::string> env_to_option;
	po::options_description common_options;

	po::variables_map results;

	int opt;
	string config_file;
	string generic_bGeno_file;
	string generic_pheno_file;
	vector<int> generic_pheno_index;
	string generic_covar_file;
	vector<int> generic_covar_index;

	string generic_extract_file;
	vector<string> generic_chr;
	string generic_keep_file;
	string generic_snp_tag_file;
	string generic_ref_allele_file;

	int generic_thread;
	int generic_seed;

	int generic_block_snps;
	float generic_block_kb;
	int generic_block_num;

	double generic_maf;
	double generic_lmiss;
	double generic_hwe;
	bool generic_remove_zero_gstd;
	int generic_adj_var;
	int generic_imput;

	bool propc_switch;
	int propc_vec_num;
	int propc_max_it;

	bool eigengwas_switch;
	string eigengwas_read_v;
	int eigengwas_vec;

	bool enc_switch;
	int enc_k;

	bool randhe_switch;
	double round_err;
	int randhe_B0;
	int randhe_B1;
	int randhe_maxit;

	bool xld_switch;
	int xld_alg;
	bool xld_stage = false;
	string xld_list;

	bool makebed_switch;

	bool inbred;

	bool me_switch;
	int me_leap;
	double me_stop;

	bool enc_reg_switch;
	string enc_reg_f1;
	string enc_reg_f2;

};

#endif

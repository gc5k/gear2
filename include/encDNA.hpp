#ifndef ENCDNA_HPP_
#define ENCDNA_HPP_

#include "time.h"

#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <math.h>
#include <boost/math/distributions/students_t.hpp>
#include "genotype.h"
#include "Goptions.hpp"
#include "global.h"

using namespace std;

extern Goptions goptions;
extern genotype g;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;

class encDNA {
private:
	int K;
	int seed;
	MatrixXdr encGeno;
	MatrixXdr encGRM;
	vector<int> snpOrder;
	clock_t encDNA_begin;

public:
	encDNA(int k) {
		clock_t encDNA_begin = clock();

		K = k;
		seed = goptions.GetGenericSeed();
		checkOrder();
	}

	void checkOrder() {
		map<string, string> _saMap;
		map<string, int> _snpIdxMap;

		if (!goptions.GetGenericForceRefAlleleFile().empty()) {
			string saFile = goptions.GetGenericForceRefAlleleFile();
			ifstream inp(saFile);
			if (!inp.is_open()) {
				cerr << "Error reading file "<< saFile << endl;
				exit(1);
			} else {
				cout << "Reading SNP allele pairs from " << saFile << endl;
			}

			string line;
			string whitespaces (" \t\f\v\n\r");
			int LineCnt = 0;
			while (std::getline(inp, line)) {
				char c = line[0];
				if (c == '#') {
					continue;
				}
				if (line.empty()) continue;

				line.erase(0, line.find_first_not_of(" ")); //remove head whitespace;
				std::size_t found = line.find_last_not_of(whitespaces); //remove tail whitespace;
				if (found != std::string::npos) {
					line.erase(found + 1);
				}

				istringstream ss(line);
				string _s;
				string _a;
				string word;
				int cnt_word = 0;
				while (ss >> word && cnt_word < 2) {
					switch (cnt_word) {
						case 0:
							_s = word;
							break;
						default:
							_a = word;
							break;
					}
					cnt_word++;
				}
				if (cnt_word != 2) {
					cerr << "Line " << LineCnt + 1 << " invalid format" << endl;
					cerr << "Please make sure two columns in each line for --recode-ref-allele" << endl;
					exit(0);
				}
				auto it = _saMap.find(_s);
				if (it != _saMap.end()) {
					cerr << "Line: " << LineCnt + 1 << " snp " << _s << " is associated with ref allele " << _a << " and " << it->second << endl;
					cerr << "SNP " << _s << " has been mapped to more than one reference alleles" << endl;
					exit(0);
				}
				_saMap.insert(pair<string, string>(_s, _a));
				_snpIdxMap.insert(pair<string, int>(_s, LineCnt));
				LineCnt++;
			}

			if (g.act_snp.size() == _saMap.size()) {
				for (int i = 0; i < g.act_snp.size(); i++) {
					snpInfo _snp = g.get_snp_info(g.act_snp[i]);
					auto it = _saMap.find(_snp.snpID);
					if (it != _saMap.end()) {
						if (it->second == _snp.a1) {
							int _idx = (_snpIdxMap.find(_snp.snpID))->second;
							snpOrder.push_back(_idx);
						} else {
						cout << "SNP ref allele did not match: " << _snp.snpID << " " 
						<< _snp.a1 << " " << it->second << endl;
						cout << "Make sure all snps have reference allels matched." << endl;
						exit(0);
						}
					}
				}
			}
			cout << "Lined up " << snpOrder.size() << " unique snps for --enc" << endl;
		}
	}

	void encG() {

		int Nindv = g.get_active_sample_size();
		int Nsnp = g.get_active_snp_number();
		cout << "\nencG generates " << K << " tags for " << Nindv <<" samples." << endl;

		g.make_MailmanP();
		mailbox::setMem();

		srand(seed);
		std::default_random_engine generator(seed);
		std::normal_distribution<double> norm_dist(0, sqrt(1.0 / (1.0 * Nsnp)));
		cout<< "Random number is sampled from N[0, sqrt(1/" << Nsnp << ")]" << endl;
		MatrixXdr Bz(K, Nsnp);
		for (int i = 0; i < Bz.rows(); i++)
			for (int j = 0; j < Bz.cols(); j++)
				Bz(i, snpOrder[j]) = norm_dist(generator);

		MatrixXdr encGenoT(K, Nindv);

		mailbox::multiply_y_post(Bz, K, encGenoT, true);

		cout << "Transposing encGeno..." << endl;
		encGeno.resize(Nindv, K);
		encGeno = encGenoT.transpose();
		cout << "encGeno size " << encGeno.rows() << " " << encGeno.cols() << endl;
		cout << "done." <<endl;
	}

	void deCapKing() {
		cout << "Generating encGRM for " << encGeno.rows() << " samples using " << K << " tags." << endl;
		clock_t encGRM_begin = clock();

		encGRM.resize(encGeno.rows(), encGeno.rows());
		for (int i = 0; i < encGRM.rows(); i++) {
			for (int j = 0; j <= i; j++) {
				for (int p = 0; p < K; p++) {
					encGRM(i, j) += encGeno(i, p) * encGeno(j, p); //not very efficient 
				}
				encGRM(i, j) /= K;
			}
		}
		cout << "first couple of sample pairs" << endl;
		cout << encGRM(0, 0) << endl;
		cout << encGRM(1, 0) << " " << encGRM(1, 1) << endl;
		cout << encGRM(2, 0) << " " << encGRM(2, 1) << " " << encGRM(2, 2) << endl;
		clock_t encGRM_end = clock();
		double encGRM_time = double(encGRM_end - encGRM_begin) / CLOCKS_PER_SEC;
		cout << "encGRM time: " << encGRM_time << "s." << endl;

	}

	void print_encG() {
		cout << "Save encGeno to " << (goptions.GetGenericOutFile() + string(".enc")).c_str() <<endl;
		ofstream e_file;
		e_file.open((goptions.GetGenericOutFile() + string(".enc")).c_str(), ios::out);
		for (int i = 0; i < encGeno.rows(); i++) {
			indxInfo ind = g.get_fam_info(g.act_ind[i]);
			e_file << ind.fid << " " << ind.iid << " ";
			for (int j = 0; j < encGeno.cols(); j++) {
				e_file << encGeno(i, j);
				if (j != (encGeno.cols() - 1)) e_file << " ";
			}
			e_file << endl;
		}
		e_file.close();

		cout << "Saved encGRM into " << (goptions.GetGenericOutFile() + string(".encGRM")).c_str() <<endl;
		ofstream e_file_grm;
		e_file_grm.open((goptions.GetGenericOutFile() + string(".encGRM")).c_str(), ios::out);
		for (int i = 0; i < encGRM.rows(); i++) {
			for (int j = 0; j <= i; j++) {
				e_file_grm << i << " " << j << " " << encGRM(i, j) << endl;
			}
		}
		e_file_grm.close();
	}

	~encDNA() {
		cout << "Finishing encDNA ..." << endl;
		mailbox::cleanMem();
		clock_t encDNA_end = clock();
		double encDNA_time = double(encDNA_end - encDNA_begin) / CLOCKS_PER_SEC;
		cout << "encDNA total time " << encDNA_time << "s." << endl;
	}
};

#endif

/*
void ENC(int seed, int kval) {
	cout << "ENC generates " << goptions.GetEncK() << " tags for " << g.Nindv <<" samples" <<endl;
	clock_t ENC_begin = clock();

	srand(goptions.GetGenericSeed());
	std::default_random_engine generator(goptions.GetGenericSeed());
	std::normal_distribution<double> norm_dist(0, 1.0);

	MatrixXdr Bz(goptions.GetEncK(), g.Nsnp);
	for (int i = 0; i < Bz.rows(); i++)
		for (int j = 0; j < Bz.cols(); j++)
			Bz(i, j) = norm_dist(generator);

	MatrixXdr encG(goptions.GetEncK(), g.Nindv);
	multiply_y_post(Bz, goptions.GetEncK(), encG, true);

	ofstream e_file;
	e_file.open((goptions.GetGenericOutFile() + string(".enc.txt")).c_str());
	for (int i = 0; i < encG.cols(); i++) {
		for (int j = 0; j < encG.rows(); j++) {
			e_file << encG(j, i);
			if (j != (encG.rows() - 1)) e_file << " ";
		}
		e_file << endl;
	}

	e_file.close();
	clock_t ENC_end = clock();
	double ENC_time = double(ENC_end - ENC_begin) / CLOCKS_PER_SEC;
	cout<< "Save encG to " <<(goptions.GetGenericOutFile() + string(".enc.txt")).c_str() <<endl;
	cout << "ENC time " << ENC_time << endl;
}
*/
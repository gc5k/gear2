#ifndef ENCGREG_HPP_
#define ENCGREG_HPP_

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
#include "plinkNum.h"


using namespace std;
using namespace plinknum;


extern Goptions goptions;
extern genotype g;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;

class encGReg {
public:
	encGReg() {
		encReg_begin = clock();
		cout << "Start encG-regression " << endl;
		f1 = goptions.GetEncRegF1();
		f2 = goptions.GetEncRegF2();
		readEncG(f1, &fam1, &fmat1);
		readEncG(f2, &fam2, &fmat2);
		enc_k = fmat1.cols() * 1.0;
		//gReg();
	}

	void gReg() {

		using boost::math::students_t;
		students_t Tdist(enc_k - 1);
		using boost::math::complement;

		cout << "Running enc-reg regression ..." << endl;

		MatrixXdr encGRM = fmat1 * fmat2.transpose() / enc_k;

		ofstream e_file;
		cout << "Writing results to " << (goptions.GetGenericOutFile() + string(".encReg")).c_str() << endl;
		e_file.open((goptions.GetGenericOutFile() + string(".encReg")).c_str(), ios::out);
		e_file << "fid1 iid1" << "\t" << "fid2 iid2" << "\t" << "encGreg" <<"\t" << "se" << "\t" << "tstat" << "\t" << "pval" << endl;
		for (int i = 0; i < encGRM.rows(); i++) {
			for (int j = 0; j < encGRM.cols(); j++) {
				double _s = sqrt((1 - encGRM(i, j) * encGRM(i, j)) / (enc_k - 1));
				double _t = encGRM(i, j) / _s;

				string _p = boost::to_string(cdf(complement(Tdist, _t)));
				if (cdf(complement(Tdist, _t)) == 0) _p = plinknum::tstat(enc_k - 1, fabs(_t));
				e_file << fam1[i] << "\t" << fam2[j] << "\t" << encGRM(i, j) << "\t" << _s << "\t" << _t << "\t" << _p << endl;
			}
		}
		e_file.close();
	}

	void readEncG(string fname, vector<string> *fvec, MatrixXdr *fmat) {
		ifstream inp(fname.c_str());
		cout << "Reading " << fname.c_str() << " ..." << endl;

		vector<vector<double> > encGvec;
		int cntLine = 0;
		string line;
		string whitespaces (" \t\f\v\n\r");
		int _Nind = 0;
		while (std::getline(inp, line)) {
			cntLine++;
			line.erase(0, line.find_first_not_of(" ")); //remove head whitespace;
			std::size_t found = line.find_last_not_of(whitespaces); //remove tail whitespace;
			if (found != std::string::npos) line.erase(found + 1);

			char c = line[0];
			if (c == '#') continue;
			if (line.empty()) continue;
			istringstream ss(line);
			string word;
			int cnt_len = 0;
			vector<double> _encVec;
			double _sum = 0;
			double _sumSq = 0;
			string _id;
			while (ss >> word) {
				switch (cnt_len) {
					case 0:
						_id = word;
						break;
					case 1:
						_id = _id + " " + word;
						break;
					default:
						double _v = boost::lexical_cast<double>(word);
						_encVec.push_back(_v);
						_sum += _v;
						_sumSq += _v * _v;
				}
				cnt_len++;
			}
			double _len = _encVec.size();
			double _m = _sum / (_len * 1.0);
			double _sd = sqrt((_sumSq - _m * _m * cnt_len * 1.0)/ (_len - 1.0));

			//standardization
			for (int i = 0; i < _encVec.size(); i++) {
				_encVec[i] = (_encVec[i] - _m) / _sd;
			}
			encGvec.push_back(_encVec);
			(*fvec).push_back(_id);
		}

		(*fmat).resize(encGvec.size(), encGvec[0].size());
		for (int i = 0; i < encGvec.size(); i++) {
			vector<double> _encG = encGvec[i];
			for (int j = 0; j < _encG.size(); j++) {
				(*fmat)(i, j) = _encG[j];
			}
		}
		cout << "read " << encGvec.size() << " samples for " << encGvec[0].size() << " enc scores" << endl;
	}

    ~encGReg() {
		cout << "Finishing enc-reg ..." << endl;
		clock_t encReg_end = clock();
		double encReg_time = double(encReg_end - encReg_begin) / CLOCKS_PER_SEC;
		cout << "EigenGWAS total time " << encReg_time << "s." << endl;
    }

private:
    string f1;
	vector<string> fam1;
    string f2;
	vector<string> fam2;
    MatrixXdr fmat1;
    MatrixXdr fmat2;
	double enc_k;
	clock_t encReg_begin;

};
#endif

#ifndef XLD_HPP_
#define XLD_HPP_

#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
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

struct rSqUnit {
	long double psq;
	long double rsq;
	long double sd;
	int mi;
	int mj;
	string p1;

	// long double rsq_s;
	// long double sd_s;
	// string ps;
	rSqUnit() {
		psq = 0.0;
		rsq = 0.0;
		sd = 0.0;
		mi = 0;
		mj = 0;
		p1 = "0";

		// rsq_s = -1;
		// sd_s = -1;
		// ps = "NA";
	}
};
void estimate_rsq(int start, int end, int itB, vector<int>& snp_sizeVec, vector<MatrixXdr>& XXzVec,vector<rSqUnit>& rsqVec){
	int si = floor(sqrt(start * 2));
	if(si * (si + 1) / 2 > start) si--;
	int sj = start - (si * (si + 1) / 2);

	using boost::math::students_t;
	students_t Tdist(itB - 1);
	using boost::math::complement;
    
	int index = start, i = si, j = sj;
	while(index < end){

	long _Nind = XXzVec[i].cols();
	rSqUnit ru;
	VectorXd Lb = XXzVec[i].cwiseProduct(XXzVec[j]).rowwise().sum();

	double LbS;
	long double LbMean = 0;
	long double Lb_Sqsum = 0;

	// for (int k = 0; k < Lb.rows(); k++) {
    //         LbS=0;
	// 		for (int l = 0; l < Lb.cols(); l++) {
	// 			LbS += Lb(k, l);
	// 			//LbMean += Lb(k, l);
	// 		}
	// 		LbMean+=LbS;
	// 		Lb_Sqsum += LbS * LbS;
	// 	}
	Lb_Sqsum = Lb.array().square().sum()/ (1.0 * snp_sizeVec[i] * snp_sizeVec[j] * snp_sizeVec[i] * snp_sizeVec[j]);
	LbMean = Lb.sum() / (1.0 * snp_sizeVec[i] * snp_sizeVec[j] * itB);
	long double rho2 = (LbMean - _Nind) / (1.0 * _Nind * _Nind );
	long double sd_Lb = sqrt((Lb_Sqsum - itB * LbMean * LbMean) / (itB * 1.0));

	long double sd = sqrt(sd_Lb * sd_Lb / itB * 1.0) / (_Nind * _Nind  * 1.0) ;
	long double tstat = rho2 / sd;
	string p_tstat = boost::to_string(cdf(complement(Tdist, tstat)));
	if (cdf(complement(Tdist, tstat)) == 0) p_tstat = plinknum::tstat(itB - 1, tstat);
	ru.rsq = rho2;
	ru.sd = sd;
	ru.mi = snp_sizeVec[i];
	ru.mj = snp_sizeVec[j];
	ru.p1 = p_tstat;
	ru.psq = rho2 + 1.0 / _Nind;
	rsqVec[(i * (i+1) / 2) + j] = ru;

	j++;
	index++;
		if(j > i){
            i++;
			j = 0;
		}
	}
}

class xLD {
public:
	xLD() {
		clock_t xld_begin = clock();

		if (goptions.GetXLDStage()) {

			cout << "Starting XLD randomization algorithm (stage) ..." << endl;

			cout << "Examining tags ..." << endl;
			for (int c = 0; c < g.gMinfo.snpTagVec.size(); c++) {
				vector<string> _tag;
				_tag.push_back(g.gMinfo.snpTagVec[c]);
				g.generate_active_snp_sample(_tag);
				int _Nsnp = g.get_active_snp_number();
				if (_Nsnp > 0) {
					tag_nameVec.push_back(g.gMinfo.snpTagVec[c]);
				}
			}
			if (tag_nameVec.size() == 0) {
				cerr << "No valid tag associated with snps. Quit" << endl;
				exit(0);
			} else {
				cout << "Validated " << tag_nameVec.size() << " tags for xld randomization algorithm" << endl;
				cout << endl;
			}

			itB = goptions.GetGenericIteration();

			writeXXz();
			writeXXzList();
		}

		if (!goptions.GetXLDList().empty()) {
			readXXzList();
			generateXLD();
		}

		if (!goptions.GetXLDStage() && goptions.GetXLDList().empty()) {
			cout << "Starting XLD randomization algorithm ..." << endl;

			cout << "Examining tags ..." << endl;
			for (int c = 0; c < g.gMinfo.snpTagVec.size(); c++) {
				vector<string> _tag;
				_tag.push_back(g.gMinfo.snpTagVec[c]);
				g.generate_active_snp_sample(_tag);
				int _Nsnp = g.get_active_snp_number();
				if (_Nsnp > 0) {
					tag_nameVec.push_back(g.gMinfo.snpTagVec[c]);
				}
			}
			if (tag_nameVec.size() == 0) {
				cerr << "No valid tag associated with snps. Quit" << endl;
				exit(0);
			} else {
				cout << "Validated " << tag_nameVec.size() << " tags for xld randomization algorithm" << endl;
				cout << endl;
			}
			rsqVec = vector<rSqUnit>(tag_nameVec.size() * (tag_nameVec.size() + 1) / 2);
			generateXXz();
			generateXLD();
		}
	}

	void readXXzList() {
		string xxzFile;
		string xxzList = goptions.GetXLDList();
		ifstream inp(xxzList.c_str());
		cout << "Reading " << xxzList.c_str() << " ..." << endl;

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
			if (c == '#')	continue;
			if (line.empty()) continue;
			istringstream ss(line);
			string word;
			int cnt_len = 0;
			vector<string> strVec;
			while (ss >> word && cnt_len <5) {
				strVec.push_back(word);
				cnt_len++;
			}
			if (cnt_len == 5) {
				xxzFile = strVec[0];
				int m;
				try {
					m = boost::lexical_cast<int>(strVec[1]);
					_Nind = boost::lexical_cast<int>(strVec[2]);
					itB = boost::lexical_cast<int>(strVec[3]);
					tag_nameVec.push_back(strVec[4]);
				} catch (const boost::bad_lexical_cast &e) {
					cerr << e.what() << '\n';
					cerr << "Line " << cnt_len << " of " << xxzList << " the elements  are invalid number(s)" << endl; 
					exit(1);
				}
				snp_sizeVec.push_back(m);
			}

			MatrixXdr XXz(itB, _Nind);
			ifstream xxz_inp(xxzFile.c_str());
			cout << "Reading " << xxzFile.c_str() << ": iteration = " << itB << ", sample size = " << _Nind << endl;
			int _cntLine = 0;
			string _line;
			string _whitespaces (" \t\f\v\n\r");
			while (std::getline(xxz_inp, _line)) {
				vector<string> _strVec;
				string _word;
				istringstream _ss(_line);

				while (_ss >> _word) {
					_strVec.push_back(_word);
				}

				for (int j = 2; j < _strVec.size(); j++) {
					XXz(j - 2, _cntLine) = boost::lexical_cast<double>(_strVec[j]);
				}
				_cntLine++;
			}
			XXzVec.push_back(XXz);
		}
	}

	void writeXXzList() {
		string xxzList = goptions.GetGenericOutFile() + ".xxz.list";
		fstream xxzListOut(xxzList.c_str(), ios::out);
		for (int i = 0; i < tag_nameVec.size(); i++) {
			string xxzname = goptions.GetGenericOutFile() + "." + tag_nameVec[i] + ".xxz";
			xxzListOut << xxzname << " " << snp_sizeVec[i] << " " << g.act_ind.size() << " " << itB << " " << tag_nameVec[i] << endl;	
		}
		xxzListOut.close();
		cout << "Writing summary information into " << xxzList << endl;
	}

	void generateXXz() {

		cout << "--------------------Generating XXz---------------" << endl;
		MatrixXdr zRand; //(p,k)
		int Nsample = g.get_active_sample_size();
		itB = goptions.GetGenericIteration();
		zRand.resize(Nsample, itB);

		std::default_random_engine generator(goptions.GetGenericSeed());
		std::normal_distribution<double> norm_dist(0, 1.0);
		for (int i = 0; i < zRand.cols(); i++)
			for (int j = 0; j < zRand.rows(); j++)
				zRand(j, i) = norm_dist(generator);

		int _Nind = g.get_active_sample_size();
		for (int i = 0; i < tag_nameVec.size(); i++) {
			vector<string> _tag;
			_tag.push_back(tag_nameVec[i]);
			g.generate_active_snp_sample(_tag);
			int _Nsnp = g.get_active_snp_number();
			snp_sizeVec.push_back(_Nsnp);

			g.make_MailmanP();
			mailbox::setMem();
			MatrixXdr XXz(itB, _Nind);

			for (int j=0; j < itB; j++){

			//MatrixXdr zBeta(_Nsnp, itB);
			MatrixXdr zBeta(_Nsnp, 1);
			MatrixXdr zRandCol = zRand.col(j);
			MatrixXdr XXz_tmp(1, _Nind);
			//mailbox::multiply_y_pre(zRand, itB, zBeta, true);
			mailbox::multiply_y_pre(zRandCol, 1, zBeta, true);

			MatrixXdr zBetaT = zBeta.transpose();			
			mailbox::multiply_y_post(zBetaT, 1, XXz_tmp, true);
			XXz.row(j) = XXz_tmp;

			}

			mailbox::cleanMem();
			XXzVec.push_back(XXz);
		}
	}

	void writeXXz() {

		cout << "--------------------Generating XXz---------------" << endl;
		MatrixXdr zRand; //(p,k)
		int Nsample = g.get_active_sample_size();
		itB = goptions.GetGenericIteration();
		zRand.resize(Nsample, itB);

		std::default_random_engine generator(goptions.GetGenericSeed());
		std::normal_distribution<double> norm_dist(0, 1.0);
		for (int i = 0; i < zRand.cols(); i++)
			for (int j = 0; j < zRand.rows(); j++)
				zRand(j, i) = norm_dist(generator);

		int _Nind = g.get_active_sample_size();
		for (int i = 0; i < tag_nameVec.size(); i++) {
			vector<string> _tag;
			_tag.push_back(tag_nameVec[i]);
			g.generate_active_snp_sample(_tag);
			int _Nsnp = g.get_active_snp_number();
			snp_sizeVec.push_back(_Nsnp);

			g.make_MailmanP();
			mailbox::setMem();
			MatrixXdr zBeta(_Nsnp, itB);
			mailbox::multiply_y_pre(zRand, itB, zBeta, true);

			MatrixXdr zBetaT = zBeta.transpose();

			MatrixXdr XXz(itB, _Nind);
			mailbox::multiply_y_post(zBetaT, itB, XXz, true);
		
			string xxzname = goptions.GetGenericOutFile() + "." + tag_nameVec[i] + ".xxz";
			cout << "Writing XXz of tag " << tag_nameVec[i] << " into " << xxzname << endl;
			fstream xxzout(xxzname.c_str(), ios::out);
			for (int j = 0; j < g.act_ind.size(); j++) {
				indxInfo iIdx = g.get_fam_info(g.act_ind[j]);
				xxzout << iIdx.fid << " " << iIdx.iid << " ";
				for (int k = 0; k < itB; k++) {
					xxzout << XXz(k, j);
					if (k != (itB - 1)) { 
						xxzout << " ";
					}
				}
				xxzout << endl;
			}
			xxzout.close();
			mailbox::cleanMem();
		}
	}

	void generateXLD() {
		int nthreads = goptions.GetGenericThreads();
		int all_size = XXzVec.size() * (XXzVec.size()+1) / 2;
		nthreads = min(nthreads, all_size);
		std::thread th[nthreads];
		int perthread = all_size / nthreads;
		int t = 0;	
				
		//xld_stat.resize(XXzVec.size() * (XXzVec.size() + 1) / 2);
		for (; t < nthreads-1; t++) {
			th[t] = std::thread(estimate_rsq, t * perthread, (t+1) * perthread, itB, std::ref(snp_sizeVec),std::ref(XXzVec), std::ref(rsqVec));
		}
		th[t] = std::thread(estimate_rsq, t * perthread, all_size, itB, std::ref(snp_sizeVec),std::ref(XXzVec), std::ref(rsqVec));
		for (int t = 0; t < nthreads; t++) {
			th[t].join();
		}

		// for (int i = 0; i < XXzVec.size(); i++) {
		// 	int _Nind = XXzVec[i].cols();
		// 	vector<double> l1 = LbChr[(i + 1) * (i + 2)/2 - 1];
		// 	rSqUnit r1 = rsqVec[(i + 1) * (i + 2) / 2 - 1];
		// 	for (int j = 0; j <= i; j++) {
		// 		int idx_ij = i  * (i + 1) / 2 + j;
		// 		vector<double> l2 = LbChr[(j + 1) * (j + 2)/2 - 1];
		// 		vector<double> l12 = LbChr[idx_ij];
		// 		rSqUnit r2 = rsqVec[(j + 1) * (j + 2) / 2 - 1];
		// 		rSqUnit r12 = rsqVec[idx_ij];
		// 		double _r12sumSq = 0;
		// 		double _r12mean = 0;
		// 		int _effCnt = 0;
		// 		double _l12 = 0;
		// 		for (int k = 0; k < itB; k++) {
		// 			double _r12 = (l12[k] - _Nind) / (1.0 * _Nind * (_Nind + 1.0));
		// 			double _r1 = (l1[k] - _Nind) / (1.0 * _Nind * (_Nind + 1.0));
		// 			double _r2 = (l2[k] - _Nind) / (1.0 * _Nind * (_Nind + 1.0));
		// 			if ( ((_r1 * r1.mi - 1) / (r1.mi - 1)) > 0 && ((_r2 * r2.mj - 1) / (r2.mj - 1)) > 0) {
		// 				if (i != j) {
		// 					_l12 = _r12 / sqrt(((_r1 * r1.mi - 1) / (r1.mi - 1)) * ((_r2 * r2.mj - 1) / (r2.mj - 1)));
		// 				} else {
		// 					_l12 = 1;
		// 				}
		// 				_r12sumSq += _l12 * _l12;
		// 				_r12mean += _l12;
		// 				_effCnt++;
		// 			}
		// 		}

		// 		string p_tstat = "NA";
		// 		if (i != j) {
		// 			if (_effCnt > 1) {
		// 				_r12mean /= (_effCnt * 1.0);
		// 				double r12_s = -1;
		// 				double rsd = -1;
		// 				double sd_s = -1;
		// 				double tstat = 0;
		// 				bool _valid = true;
		// 				if ( ((r1.rsq * r1.mi - 1) / (r1.mi - 1)) > 0 && ((r2.rsq * r2.mj - 1) / (r2.mj - 1)) > 0) {
		// 					r12_s = r12.rsq / sqrt(((r1.rsq * r1.mi - 1) / (r1.mi - 1)) * ((r2.rsq * r2.mj - 1) / (r2.mj - 1)));
		// 				} else {
		// 					_valid = false;
		// 				}

		// 				if ((_r12sumSq - _effCnt * _r12mean * _r12mean) > 0 && r12_s > 0) {
		// 					rsd = sqrt((_r12sumSq - _effCnt * _r12mean * _r12mean) / (_effCnt * 1.0 - 1.0));
		// 					sd_s = rsd / sqrt(_effCnt * 1.0 - 1.0); 
		// 				} else {
		// 					_valid = false;
		// 				}

		// 				if (_valid) {
		// 					tstat = r12_s / sd_s;
		// 					rsqVec[idx_ij].rsq_s = r12_s;
		// 					rsqVec[idx_ij].sd_s = sd_s;
		// 					students_t _Tdist(_effCnt - 1);
		// 					p_tstat = boost::to_string(cdf(complement(_Tdist, tstat)));
		// 					if (cdf(complement(_Tdist, tstat)) == 0) p_tstat = plinknum::tstat(_effCnt - 1, tstat);
		// 					rsqVec[idx_ij].ps = p_tstat;
		// 				} else{
		// 					rsqVec[idx_ij].rsq_s = -1;
		// 					rsqVec[idx_ij].sd_s = -1;
		// 					rsqVec[idx_ij].ps = "NA";
		// 				}
		// 			}
		// 		} else {
		// 			if (_effCnt > 1) {
		// 				rsqVec[idx_ij].rsq_s = 1;
		// 				rsqVec[idx_ij].sd_s = 0;
		// 				rsqVec[idx_ij].ps = "0";
		// 			}
		// 		}
		// 	}
		// }

		ofstream e_file_r;
		string fname_r = goptions.GetGenericOutFile() + string(".xldr");
		e_file_r.open(fname_r.c_str());
		cout << "Writing randomized x-ld results to " << fname_r.c_str() << " ..." << endl;
		e_file_r << "Tagi\tSNPi\tTagj\tSNPj\trsq\tsd\tpval\tpsq" << endl;
		//cout << "Tagi\tSNPi\tTagj\tSNPj\trsq\tsd\tpval\trsq_s\tsd_s\tpval_s" << endl;

		for (int i = 0; i < XXzVec.size(); i++) {
			for (int j = 0; j <= i; j++) {
				rSqUnit rsq = rsqVec[i  * (i + 1) / 2 + j];
				e_file_r << tag_nameVec[i] << "\t" << rsq.mi << "\t" << tag_nameVec[j] << "\t" << rsq.mj << "\t" 
				<< rsq.rsq << "\t" << rsq.sd << "\t" << rsq.p1 << "\t" << rsq.psq << endl;
				//<< rsq.rsq_s << "\t" << rsq.sd_s << "\t" << rsq.ps << endl;
				// cout << tag_nameVec[i] << "\t" << rsq.mi << "\t" << tag_nameVec[j] << "\t" << rsq.mj << "\t" 
				// << rsq.rsq << "\t" << rsq.sd << "\t" << rsq.p1 << "\t"
				// << rsq.rsq_s << "\t" << rsq.sd_s << "\t" << rsq.ps << endl;
			}
		}
		e_file_r.close();
	}

	~xLD() {
		cout << "Finishing randomized xLD algorithm ..." << endl;
		// if (!goptions.GetXLDStage() && goptions.GetXLDList().empty()) {
		// 	mailbox::cleanMem();
		// }
		clock_t xld_end = clock();
		double xld_time = double(xld_end - xld_begin) / CLOCKS_PER_SEC;
		cout << "xld randomization total time " << xld_time << "s." << endl;
	}

private:
	vector<MatrixXdr> xld_stat;
	vector<rSqUnit> rsqVec;
	vector<MatrixXdr> XXzVec;
	vector<int> snp_sizeVec;
	vector<string> tag_nameVec;
	clock_t xld_begin;

	int itB;
};
#endif

#ifndef EIGENGWAS_HPP_
#define EIGENGWAS_HPP_

#include "time.h"
#include <limits>
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <math.h>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/exception/to_string.hpp>
#include "genotype.h"
#include "Goptions.hpp"
#include "global.h"
#include "mailbox.h"
#include "plinkNum.h"

using namespace std;
using namespace plinknum;
using namespace mailbox;

extern Goptions goptions;
extern genotype g;


typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;

class EigenGWAS {
	private:
		clock_t eg_begin;
		clock_t eg_end;

	public:
	EigenGWAS() {
		g.make_MailmanP();
		mailbox::setMem();
		proEigen proEigenCal;
		proEigenCal.RandPC();
	
		eg_begin = clock();
		Scan(proEigenCal.getEveP(), proEigenCal.getEigenValue());
		
	}

	EigenGWAS(const string & name) {
		g.make_MailmanP();
		mailbox::setMem();
		eg_begin = clock();
		Scan(read_eveP(), read_eval());
	}

	~EigenGWAS() {
		mailbox::cleanMem();
		eg_end = clock();
		cout << "Finishing up EigenGWAS..." << endl;
		double eg_time = double(eg_end - eg_begin) / CLOCKS_PER_SEC;
		cout << "EigenGWAS scan time: " << eg_time << "s" << endl;
	}

	void Scan(MatrixXdr IndEigenVec, MatrixXdr eigenVal) {

		if (IndEigenVec.rows() != g.act_ind.size()) {
			cerr << "The rows of eigenvector (" << IndEigenVec.rows() << ") does not match "
			<< "the active sample size (" << g.act_ind.size() << ")" << endl;
			exit(0);
		}

		if (IndEigenVec.cols() != eigenVal.rows()) {
			cerr << "The cols of eigenvector (" << IndEigenVec.cols() << ") does not match "
			<< "the eigen values " << endl;
			exit(0);
		}

//		double pmin = 4.940e-324;
		clock_t eg_begin = clock();
		int Nsnp = g.get_active_snp_number();
		int Nindv = g.get_active_sample_size();
		cout << "---------------EigenGWAS scan-------------" << endl;
		cout << "SNP number: " << Nsnp << endl;
		cout << "Sample size: " << Nindv << endl;

		cout << "Eigenvector number: " << IndEigenVec.cols() <<endl;

		cout << "Standardization eigenvec" << endl;
		MatrixXdr evePS(IndEigenVec.rows(), IndEigenVec.cols());
		Eigen::VectorXd _cm = IndEigenVec.colwise().sum() / (1.0 * IndEigenVec.rows());
		Eigen::VectorXd _css = IndEigenVec.colwise().squaredNorm().colwise().sum() / (1.0 * IndEigenVec.rows());
		Eigen::VectorXd _csq = (_css - _cm.cwiseProduct(_cm)).array().sqrt();
		Eigen::VectorXd _csq_inv = 1 / _csq.array();
		evePS = (IndEigenVec.rowwise() - _cm.transpose()) * _csq_inv.asDiagonal();

		//using mailman
		MatrixXdr EgBeta(Nsnp, IndEigenVec.cols());
		mailbox::multiply_y_pre(evePS, evePS.cols(), EgBeta, true);

		EgBeta = EgBeta / (1.0 * Nindv);

		using boost::math::students_t;
		students_t Tdist(Nindv - 1);

		using boost::math::chi_squared;
		using boost::math::non_central_chi_squared;
		using boost::math::quantile;
		using boost::math::complement;

		chi_squared Chisq_dist(1);
		for (int i = 0; i < EgBeta.cols(); i++) {
			double eval = eigenVal(i, 0);
			vector<int> indIdx1;
			vector<int> indIdx2;
			for (int j = 0; j < evePS.rows(); j++) {
				if (evePS(j, i) < 0) {
					indIdx1.push_back(g.act_ind[j]);
				} else {
					indIdx2.push_back(g.act_ind[j]);
				}
			}
			cout << "------------Scanning e-space " << i + 1 << ", the eigenvalue is " << eval << endl;

			vector<double> t_stat(EgBeta.rows());
			vector<double> seB(EgBeta.rows());
			vector<string> pt2tail(EgBeta.rows());

			for (int j = 0; j <EgBeta.rows(); j++) {
				seB[j] = sqrt((1 - pow(EgBeta(j, i), 2))/(Nindv - 1));
				t_stat[j] = EgBeta(j, i) / seB[j];
				pt2tail[j] = boost::to_string(cdf(complement(Tdist, fabs(t_stat[j]))) * 2);
				if (cdf(complement(Tdist, fabs(t_stat[j]))) == 0) pt2tail[j] = tstat(Nindv-1,fabs(t_stat[j]));

				if (goptions.IsGenericDebug()) {
					cout << "MK: " << j << " infor= " << g.get_bim_info(g.act_snp[j]) << " egB=" 
					<< EgBeta(j, i) << " seB=" << seB[j] << " t=" << t_stat[j] 
					<< " p=" << pt2tail[j] << endl;
				}
			}

			vector<double> t_temp(EgBeta.rows());
			for (int j = 0; j <EgBeta.rows(); j++) {
				t_temp[j]=fabs(t_stat[j]);
			}
			sort(t_temp.begin(), t_temp.end());

			double t_median;
			if (t_temp.size() % 2 != 0) {
				t_median = t_temp[t_temp.size()/2];
			} else {
				t_median = (t_temp[t_temp.size()/2] + t_temp[t_temp.size()/2+1])/2;
			}

			string p_median=boost::to_string(cdf(complement(Tdist, t_median)));
			if (cdf(complement(Tdist, t_median)) == 0) p_median = plinknum::tstat(Nindv-1, t_median);
			double gc = pow(t_median, 2)/0.455;
			cout << "Median p is " << p_median << ", and GC factor is " << gc << endl;
			////GGC
			vector<double> chisq;
			chisq.resize(EgBeta.rows());
			for (int i = 0; i < EgBeta.rows(); i++) {
				chisq[i] = t_stat[i]*t_stat[i];
			}
			//estimate ncp
			double ncp_est;
			int block = 21;
			int step = 0.9 / ((block-1) * 1.0) ;
			for (int i = 0; i < block; i++) {
				ncp_est += quantile2(chisq,0.05+i*step)/block;
			}
			//generate ggc for each snp
			non_central_chi_squared non_Chisq_dist(1, ncp_est);
			vector<double> chisq_ncp;
			vector<double> chisq_zero;
			chisq_ncp.resize(EgBeta.rows());
			chisq_zero.resize(EgBeta.rows());

			for (int i = 0; i < EgBeta.rows(); i++) {
				chisq_zero[i] = quantile(Chisq_dist,rand() / double(RAND_MAX));
			}
			//do seeds for two distribtions have to be different?
			//srand(goptions::GetGenericSeed()+1);           
			for (int i=0; i<EgBeta.rows(); i++){
				chisq_ncp[i] = quantile(non_Chisq_dist,rand() / double(RAND_MAX));
			}
			sort(chisq_zero.begin(), chisq_zero.end());
			sort(chisq_ncp.begin(), chisq_ncp.end());
			
			//vector<vector<double> > ggc(EgBeta.rows(), vector<double>(2));
			vector<double> ggc;
			ggc.resize(EgBeta.rows());
			for (int i = 0; i < EgBeta.rows(); i++) {
				ggc[i] = chisq_ncp[i]/chisq_zero[i];
			}

			//generate an index for original order and sort
			vector<vector<double> > chi_t_ggc(EgBeta.rows(), vector<double>(2));
			for (int i = 0; i < EgBeta.rows(); i++) {
				chi_t_ggc[i][1]=i;
				chi_t_ggc[i][0]=chisq[i];
			}
			sort(chi_t_ggc.begin(),chi_t_ggc.end());

			for (int i=0; i<EgBeta.rows(); i++) {
				chi_t_ggc[i][0]/=ggc[i];
			}

			sort(chi_t_ggc.begin(),chi_t_ggc.end(), [](const vector<double>& a, const vector<double>& b) {
				return  a[1]<b[1];
			});
			

			ofstream e_file;
			string fname = goptions.GetGenericOutFile() + string("."+to_string(i + 1) + ".egs");
			e_file.open(fname.c_str());

			cout << "Writing results to " << fname.c_str() << " ..." << endl;
			e_file << "CHR\tSNP\tPOS\tBP\tA1\tA2\tBeta\tSE\tT-stat\tP\tP|gc\tP|eval\tP|ggc\tFst" << endl;
			for (int j = 0; j < EgBeta.rows(); j++) {
				double p1 = g.cal_freq(g.act_snp[j], indIdx1);
				double p2 = g.cal_freq(g.act_snp[j], indIdx2);
				double pmean = (p1 * (indIdx1.size() * 1.0) + p2 * (indIdx2.size() * 1.0)) / (evePS.rows() * 1.0);
				double fst = (p1 - p2) * (p1 - p2) / (pmean * (1 - pmean));
				double chi_t_gc = t_stat[j] * t_stat[j]/gc;
				string pgc = boost::to_string(cdf(complement(Chisq_dist, chi_t_gc)));
				if (cdf(complement(Chisq_dist, chi_t_gc)) == 0) pgc = plinknum::Chisq(chi_t_gc);
				double chi_t_ev = t_stat[j] * t_stat[j]/eval;
				string pev = boost::to_string(cdf(complement(Chisq_dist, chi_t_ev)));
				if (cdf(complement(Chisq_dist, chi_t_ev)) == 0) pev = plinknum::Chisq(chi_t_ev);
				//GGC
				string pggc = boost::to_string(cdf(complement(Chisq_dist, chi_t_ggc[j][0])));
				if (cdf(complement(Chisq_dist, chi_t_ggc[j][0])) == 0) pggc = plinknum::Chisq(chi_t_ggc[j][0]);
				e_file << g.get_bim_info(g.act_snp[j]) << "\t" << std::setprecision(5) << EgBeta(j, i) << "\t" << seB[j] << "\t" << t_stat[j] << "\t" << pt2tail[j] 
				<< "\t" << pgc << "\t" << pev << "\t"<< pggc << "\t" << fst << endl;
			}
			e_file.close();
		}

		cout << "P, P|gc, P|eval, and P|ggc are original p-value, adjusted by gc, by eigenvalue, generalized gc, respectively" << endl;

		//cout << "minimal p-val was set " << pmin << endl;
		clock_t eg_end = clock();
		double eg_time = double(eg_end - eg_begin) / CLOCKS_PER_SEC;
		cout << "EigenGWAS total time " << eg_time << "s." << endl;
	}

	MatrixXdr read_eveP() {
		int Nind = g.get_active_sample_size();
		int k = goptions.GetGenericEigenvecNumber();
		MatrixXdr eveP_r(Nind, k);
	
		string pname = goptions.GetEigenGWASreadfile() + string(".projections.txt");

		ifstream eveP(pname);
		int cntLine = 0;
		string line;
		vector<string> strVec;
			//string _whitespaces (" \t\f\v\n\r");
		while (std::getline(eveP, line)) {
			string word;
			istringstream is(line);

			while(is >> word) {
				strVec.push_back(word);
			}

			for (int i = 0; i < k; i++) {
				eveP_r(cntLine, i) = boost::lexical_cast<double>(strVec[cntLine * (k + 2) + i + 2]);				
			}
			cntLine++;
		}
		return eveP_r;
	}

	MatrixXdr read_eval() {
		int k = goptions.GetGenericEigenvecNumber();
		MatrixXdr eval_r(k, 1);
		string vname = goptions.GetEigenGWASreadfile() + string(".evals.txt");

		ifstream eval(vname);
		int _cntLine = 0;
		string _line;
		//string _whitespaces (" \t\f\v\n\r");
		while (std::getline(eval, _line)) {
			eval_r(_cntLine, 0) = boost::lexical_cast<double>(_line);
			_cntLine++;
		}
		return eval_r;
	}

//quantile
	double quantile2(vector<double>&x, double q) {
		assert(q >= 0.0 && q <= 1.0);
		const int n = x.size();
		double id = (n-1) * q;
		int lo = floor(id);
		int hi = ceil(id);
		double qs = x[lo];
		double h = (id-lo);
		return (1.0 - h) * qs + h * x[hi];
	}
};

#endif
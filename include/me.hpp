#ifndef ME_HPP_
#define ME_HPP_

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
#include "mailman.h"
using namespace std;

extern Goptions goptions;
extern genotype g;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;


class me {
public:
	me() {
		clock_t heReg_begin = clock();

		cout << "--------------------Randomized estimation for effective number of markers---------------" << endl;
		cout << "Sample size: " << g.get_active_sample_size() << endl;
		cout << "SNP number: " << g.get_active_snp_number() << endl;
		cout << "Iteration: " << goptions.GetGenericIteration() << endl;
	}

	void meRandomization() {

		int Nsnp = g.get_active_snp_number();
		int Nindv = g.get_active_sample_size();
		int Iter = goptions.GetGenericIteration();
		int iter_leap = (goptions.GetMeLeap() >= 1 && goptions.GetMeLeap() < Iter) ? goptions.GetMeLeap() : Iter;
		cout << "Set iter_leap = " << iter_leap << endl; 
		int block = Iter / iter_leap;
		int shift = Iter % iter_leap;
		double me_stop = goptions.GetMeStop();

		vector<int> me_step;
		for (int i = 0; i < block; i++) me_step.push_back(iter_leap);
		if (shift != 0) me_step.push_back(shift);

		g.make_MailmanP();
		mailbox::setMem();

		srand(goptions.GetGenericSeed());
		std::default_random_engine generator(goptions.GetGenericSeed());
		std::normal_distribution<double> norm_dist(0, 1.0);

		bool stopFlag = false;
		for (int _it = 0; _it < me_step.size(); _it++) {
			cout << "-----Begining Leap " << (_it + 1) << " [";
			cout << _it * iter_leap << "-->";
			if (_it != me_step.size()) {
				cout << (_it + 1) * iter_leap;
			} else {
				cout << _it * iter_leap + me_step[_it];
			}
			cout << "]" << endl;
			int iter = me_step[_it];
			MatrixXdr Bz(Nindv, iter);
			for (int i = 0; i < Bz.rows(); i++) {
				for (int j = 0; j < Bz.cols(); j++) {
					Bz(i, j) = norm_dist(generator);
				}
			}

			//geno_matrix * Bz; //(p x n) * (n x iter) = p x iter
			MatrixXdr T1(Nsnp, iter);
			cout << "Performing T1 = G^T * z" << endl;
			mailbox::multiply_y_pre(Bz, iter, T1, true);

			cout << "Performing T2 = T1^T * G^T" << endl;
			MatrixXdr T1_T = T1.transpose(); //iter X p
			MatrixXdr T2(iter, Nindv);
			mailbox::multiply_y_post(T1_T, iter, T2, true);

			cout << "Performing T3 = T2^T * G^T" << endl;
			MatrixXdr T2_T = T2.transpose();
			MatrixXdr T3(Nsnp, iter);
			mailbox::multiply_y_pre(T2_T, iter, T3, true);

			cout << "Performing T4 = T3^T * G^T" << endl;
			MatrixXdr T3_T = T3.transpose();
			MatrixXdr T4(iter, Nindv);
			mailbox::multiply_y_post(T3_T, iter, T4, true);

			if (goptions.IsGenericDebug()) {
				cout << T2(0, 0) << " " << T2(0, 1) << " " << T2(1, 0) << " " << T2(1, 1) << endl;
			}

			for (int i = 0; i < T2.rows(); i++) {
				long double _LB2 = (T2.row(i).array() * T2.row(i).array()).sum() / (1.0 * Nsnp * Nsnp);
				long double _LB3 = (T3.col(i).array() * T3.col(i).array()).sum() / (1.0 * Nsnp * Nsnp * Nsnp);
				long double _LB4 = (T4.row(i).array() * T4.row(i).array()).sum() / (1.0 * Nsnp * Nsnp * Nsnp * Nsnp);
				LBVec2.push_back(_LB2);
				LBVec3.push_back(_LB3);
				LBVec4.push_back(_LB4);
			}

			double _ratio = 0;
			long double _lb_m2 = 0, _lb_sd2 = 0, _lb_sum2 = 0, _lb_ss2 = 0;
			long double _lb_m3 = 0, _lb_sd3 = 0, _lb_sum3 = 0, _lb_ss3 = 0;
			long double _lb_m4 = 0, _lb_sd4 = 0, _lb_sum4 = 0, _lb_ss4 = 0;
			long double _me = 0, _me_sd = 0;
			for (int k = 0; k < LBVec2.size(); k++) {
				_lb_sum2 += LBVec2[k];
				_lb_ss2 += LBVec2[k] * LBVec2[k];

				_lb_sum3 += LBVec3[k];
				_lb_ss3 += LBVec3[k] * LBVec3[k];

				_lb_sum4 += LBVec4[k];
				_lb_ss4 += LBVec4[k] * LBVec4[k];
			}
			_lb_m2 = _lb_sum2 / (1.0 * LBVec2.size());
			_lb_sd2 = LBVec2.size() > 1 ? (sqrt((_lb_ss2 - _lb_m2 * _lb_m2 * (1.0 * LBVec2.size())) / (LBVec2.size() - 1.0))) : 0;

			_lb_m3 = _lb_sum3 / (1.0 * LBVec3.size());
			_lb_sd3 = LBVec3.size() > 1 ? (sqrt((_lb_ss3 - _lb_m3 * _lb_m3 * (1.0 * LBVec3.size())) / (LBVec3.size() - 1.0))) : 0;

			_lb_m4 = _lb_sum4 / (1.0 * LBVec4.size());
			_lb_sd4 = LBVec4.size() > 1 ? (sqrt((_lb_ss4 - _lb_m4 * _lb_m4 * (1.0 * LBVec4.size())) / (LBVec4.size() - 1.0))) : 0;

			_me = 1 / ((_lb_m2 - Nindv) / (1.0 * Nindv * (Nindv + 1)));
			_me_sd = LBVec2.size() > 1 ? (_me / (Nindv * 1.0) * _me / (Nindv * 1.0) * _lb_sd2 / sqrt(LBVec2.size() * 1.0)) : 0;
			_ratio = LBVec2.size() > 1 ? (_me_sd / _me) : 1;

			cout <<"Iter " << (_it * iter_leap + me_step[_it]) << ", \tLB2: " << _lb_m2 << ", " << _lb_sd2 <<
			 ", \tLB3: " << _lb_m3 << ", " << _lb_sd3 << 
			 ", \tLB4: " << _lb_m4 << ", " << _lb_sd4 << 
			 "\tMe: " << _me << ", " << _me_sd << ", " << _ratio << endl;

			if (_ratio < me_stop) {
				cout << "Stopping rule applied: " << _ratio << " < " << me_stop <<", after " << (_it * iter_leap + me_step[_it]) << " iterations " << endl;
				cout << "mean LB2: " << _lb_m2 << endl;
				cout << "sd LB2: " <<_lb_sd2 << endl;
				cout << "mean LB3: " << _lb_m3 << endl;
				cout << "sd LB3: " <<_lb_sd3 << endl;
				cout << "mean LB4: " << _lb_m4 << endl;
				cout << "sd LB4: " <<_lb_sd4 << endl;
				cout << "mean Me: " << _me << endl;
				cout << "se Me: " << _me_sd << endl;
				cout << "ratio (se_Me/Me): " << _ratio << endl;
				stopFlag = true;
			} else {
				if (_it == me_step.size() - 1) {
					cout << "mean LB2: " << _lb_m2 << endl;
					cout << "sd LB2: " <<_lb_sd2 << endl;
					cout << "mean LB3: " << _lb_m3 << endl;
					cout << "sd LB3: " <<_lb_sd3 << endl;
					cout << "mean LB4: " << _lb_m4 << endl;
					cout << "sd LB4: " <<_lb_sd4 << endl;
					cout << "mean Me: " << _me << endl;
					cout << "se Me: " << _me_sd << endl;
					cout << "ratio (se_Me/Me): " << _ratio << endl;
				}
			}
			if (stopFlag) break;
		}
	}

	void writeResults() {
		int Nindv = g.get_active_sample_size();
		ofstream me_file;
		string fname = goptions.GetGenericOutFile() + string(".it.me");
		me_file.open(fname.c_str());
		me_file << "Iter\tMean_LB2\tSD_LB2\tMean_LB3\tSD_LB3\tMean_LB4\tSD_LB4\tMe\tSD_Me\tSD_Me/Me" << endl;

		double _ratio = 0;
		long double _lb_m2 = 0, _lb_sd2 = 0, _lb_sum2 = 0, _lb_ss2 = 0;
		long double _lb_m3 = 0, _lb_sd3 = 0, _lb_sum3 = 0, _lb_ss3 = 0;
		long double _lb_m4 = 0, _lb_sd4 = 0, _lb_sum4 = 0, _lb_ss4 = 0;
		long double _me = 0, _me_sd = 0;

		for (int k = 0; k < LBVec2.size(); k++) {
			_lb_sum2 += LBVec2[k];
			_lb_ss2 += LBVec2[k] * LBVec2[k];
			_lb_m2 = _lb_sum2 / (1.0 * (k + 1));
			_lb_sd2 = k > 0 ? sqrt((_lb_ss2 - _lb_m2 * _lb_m2 * (1.0 * (k + 1))) / ((k + 1) - 1.0)) : 0;

			_lb_sum3 += LBVec3[k];
			_lb_ss3 += LBVec3[k] * LBVec3[k];
			_lb_m3 = _lb_sum3 / (1.0 * (k + 1));
			_lb_sd3 = k > 0 ? sqrt((_lb_ss3 - _lb_m3 * _lb_m3 * (1.0 * (k + 1))) / ((k + 1) - 1.0)) : 0;

			_lb_sum4 += LBVec4[k];
			_lb_ss4 += LBVec4[k] * LBVec4[k];
			_lb_m4 = _lb_sum4 / (1.0 * (k + 1));
			_lb_sd4 = k > 0 ? sqrt((_lb_ss4 - _lb_m4 * _lb_m4 * (1.0 * (k + 1))) / ((k + 1) - 1.0)) : 0;

			_me = 1 / ((_lb_m2 - Nindv) / (1.0 * Nindv * (Nindv + 1)));
			_me_sd = k > 0 ? (_me / (Nindv * 1.0) * _me / (Nindv * 1.0) * _lb_sd2 / sqrt(((k + 1) - 1) * 1.0)) : 0;
			_ratio = k > 0 ? (_me_sd / _me) : 1;
			me_file << (k + 1) 
			<< "\t" << _lb_m2 << "\t" << _lb_sd2 
			<< "\t" << _lb_m3 << "\t" << _lb_sd3 
			<< "\t" << _lb_m4 << "\t" << _lb_sd4
			<< "\t" << _me << "\t" << _me_sd << "\t" << _ratio << endl;
		}
		me_file.close();
	}

	~me() {
		cout << "Finishing randomized estimation for Me ..." << endl;
		mailbox::cleanMem();
		clock_t me_end = clock();
		double me_time = double(me_end - me_begin) / CLOCKS_PER_SEC;
		cout << "Me time " << me_time << "s" << endl;
	}

private:
	vector<long double> LBVec2;
	vector<long double> LBVec3;
	vector<long double> LBVec4;
	clock_t me_begin;
};

#endif
#ifndef RHEREG_HPP_
#define RHEREG_HPP_

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

class RHEreg {

public:
	RHEreg(MatrixXdr yVal, int phe_idx) {

		Yval.resize(yVal.rows(), 1);
		long double _ym = yVal.array().sum() / (1.0 * yVal.rows());
		long double _yss = (yVal.cwiseProduct(yVal)).sum() / (1.0 * yVal.rows());
		long double _ysd = sqrt(_yss - _ym * _ym);
		for (int i = 0; i < yVal.rows(); i++) {
			Yval(i, 0) = (yVal(i, 0) - _ym) / _ysd; // it is needed to remove the mean
		}

		srand(goptions.GetGenericSeed());
		Nsnp = 1.0 * g.get_active_snp_number();
		Nindv = 1.0 * g.get_active_sample_size();
		iter = goptions.GetGenericIteration();
		cout << "Randomized HE (mailman), iteration for " << iter << " times" << endl;
		round_err = goptions.GetRandHERoundErrOption();
		cout << "rounding error is " << round_err << endl;
		heReg_begin = clock();

		g.make_MailmanP();
		mailbox::setMem();

		// Numerator
		RandHE();

		std::default_random_engine generator(goptions.GetGenericSeed());

		int B_step = goptions.GetGenericIteration();
		cout << "Burning in process, B0 = " << B_step << endl;
		It_total += B_step;
		cout << "Burn " << It_total << endl;
		generator.seed(goptions.GetGenericSeed() + It_total);
		iterLB(B_step, generator);
		updateStats();

		Stop_B = (trK4 * h2) / ((yK3y - 2 * yK2y + yK1y) + (yK2y - 2 * yK1y + yIy) * (v_E / h2)) / round_err;
		long double err = (1.0 / It_total * trK4 * h2) / ((yK3y - 2 * yK2y + yK1y) + (yK2y - 2 * yK1y + yIy) * (v_E/h2));
		cout << "err " << err << ", and B1 is " << Stop_B << endl;

		int B1 = Stop_B;
		while (B1 > 0) {
			It_total += B_step;
			cout << "Burn " << It_total << endl;
			generator.seed(goptions.GetGenericSeed() + It_total);
			iterLB(B_step, generator);
			updateStats();
			B1 -= B_step;
			err = (1.0 / It_total * trK4 * h2) / ((yK3y - 2 * yK2y + yK1y) + (yK2y - 2 * yK1y + yIy) * (v_E/h2));
			cout << "err " << err << endl;
		}
	}

	~RHEreg() {
		cout << "Finishing randomized estimation for h2 ..." << endl;
		mailbox::cleanMem();
		clock_t heReg_end = clock();

		double heReg_time = double(heReg_end - heReg_begin) / CLOCKS_PER_SEC;
		cout << "RHE time " << heReg_time << endl;
	}

private:
	void RandHE() {
		int _it = 1;

		MatrixXdr yT_X(static_cast<int> (Nsnp), _it);
		mailbox::multiply_y_pre(Yval, _it, yT_X, true); // y*X

		MatrixXdr yT_XT = yT_X.transpose();
		MatrixXdr yT_X_XT(_it, static_cast<int> (Nindv));
		mailbox::multiply_y_post(yT_XT, _it, yT_X_XT, true);

		MatrixXdr yT_X_XT_T = yT_X_XT.transpose();
		MatrixXdr yT_X_XT_X(static_cast<int> (Nsnp), _it);
		mailbox::multiply_y_pre(yT_X_XT_T, _it, yT_X_XT_X, true);

		yT_X_XT_X /= Nsnp;
		yK3y = (yT_X_XT_X.array() * yT_X_XT_X.array()).sum() / Nsnp;
		cout << "yK3y " << yK3y << endl;

		yK2y = (yT_X_XT.array() * yT_X_XT.array()).sum() / (Nsnp * Nsnp);
		cout << "yK2y " << yK2y << endl;

		yK1y = (yT_X.array() * yT_X.array()).sum() / (Nsnp);
		cout << "yK1y " << yK1y << endl;

		yIy = (Yval.array() * Yval.array()).sum();
		cout << "yIy " << yIy << endl;

		yvar = yIy / (Nindv - 1.0);
	}

	void iterLB(int it, std::default_random_engine generator) {
		MatrixXdr Bz(static_cast<int> (Nindv), it);
		std::normal_distribution<double> norm_dist(0, 1.0);
		for (int i = 0; i < Bz.rows(); i++) {
			for (int j = 0; j < Bz.cols(); j++) {
				Bz(i, j) = norm_dist(generator);
			}
		}

		//geno_matrix * Bz; //(p x n) * (n x iter) = p x iter
		MatrixXdr T1(static_cast<int> (Nsnp), it);
		mailbox::multiply_y_pre(Bz, it, T1, true);
		MatrixXdr T1_T = T1.transpose(); //iter X p
		MatrixXdr T2(it, static_cast<int> (Nindv));
		mailbox::multiply_y_post(T1_T, it, T2, true);

		for (int i = 0; i < T2.rows(); i++) {
			long double _lb = 0;
			for (int j = 0; j < T2.cols(); j++) {
				_lb += T2(i, j) * T2(i, j);
			}
			LB2.push_back(_lb);
		}

		MatrixXdr T2_T = T2.transpose();
		MatrixXdr T3(static_cast<int> (Nsnp), it);
		mailbox::multiply_y_pre(T2_T, it, T3, true);

		MatrixXdr T3_T = T3.transpose();
		MatrixXdr T4(it, static_cast<int> (Nindv));
		mailbox::multiply_y_post(T3_T, it, T4, true);

		for (int i = 0; i < T4.rows(); i++) {
			long double _lb4 = 0;
			for (int j = 0; j < T4.cols(); j++) {
				_lb4 += T4(i, j) * T4(i, j); 
			}
			LB4.push_back(_lb4);
		}

		long double lb2_sum = 0;
		long double lb2_ss = 0;
		long double lb4_sum = 0;
		long double lb4_ss = 0;
		for (int i = 0; i < LB4.size(); i++) {
			long double _lb2 = LB2[i] / (Nsnp * Nsnp);
			lb2_sum += _lb2;
			lb2_ss += _lb2 * _lb2;

			long double _lb4 = LB4[i] / (Nsnp * Nsnp * Nsnp * Nsnp);
			lb4_sum += _lb4;
			lb4_ss += _lb4 * _lb4;
		}
		trK2 = lb2_sum / (1.0 * LB2.size());
		trK4 = lb4_sum / (1.0 * LB4.size());

		trK2_var = (lb2_ss - trK2 * trK2 * LB2.size()) / (1.0 * LB2.size() - 1);
		trK2_var = 2 * trK4;
		long double _lb_sd2 = sqrt(trK2_var);

		trK4_var = (lb4_ss - trK4 * trK4 * LB4.size()) / (1.0 * LB4.size() - 1);
		cout << "trK4, trK4_var " << trK4 << " " << trK4_var << endl;

		me = 1 / ((trK2 - Nindv) / ((Nindv + 1) * Nindv));
		me_sd = (me / Nindv * me / Nindv * _lb_sd2 / sqrt(1.0 * LB2.size()));
	}

	void updateStats() {

		long double trK2_N = trK2 - Nindv;
		h2 = (yK1y - yIy) / trK2_N;
		v_E = (trK2 * yIy / Nindv - yK1y) / trK2_N;

		long double v_denom = trK2_N * trK2_N;
		long double v_num = 2 * (h2 * (yK3y - 2 * yK2y + yK1y) + v_E * (yK2y - 2 * yK1y + yIy));
		long double delta = trK2_var / v_denom;
		long double delta_B = delta / (It_total * 1.0);

		long double _s1 = (delta_B + delta_B * delta_B) * h2 * h2;
		long double h2_v = v_num / v_denom + _s1;
		cout << "------------------- " << It_total << endl;

		cout << "h2 h2se v_e " << h2 <<" " << sqrt(h2_v) << " " << v_E << endl;
		cout << "trK2, trK2_var " << trK2 << " " << trK2_var << endl;
		cout << "Me, Me_sd " << me << " " << me_sd << endl;
		cout << "trK4, trK4_var " << trK4 << " " << trK4_var << endl;
	}

	long double yK3y = 0;
	long double yK2y = 0;
	long double yK1y = 0;
	long double yIy = 0;

	vector<long double> LB2;
	vector<long double> LB4;
	double Nsnp;
	double Nindv;
	MatrixXdr Yval;
	long double yvar = 0;
	long double h2 = 0;
	long double v_E = 0;
	long double mse = 0;

	long double me = 0;
	long double me_sd = 0;
	long double trK2 = 0;
	long double trK2_var = 0;
	long double trK4 = 0;
	long double trK4_var = 0;

	int Stop_B = 0;
	int It_total = 0;
	int iter;

	double round_err = 0.1;
	clock_t heReg_begin;
};

#endif

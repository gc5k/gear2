#ifndef RHEREG2_HPP_
#define RHEREG2_HPP_

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

class RHEreg2 {

public:
	RHEreg2(MatrixXdr yVal, int phe_idx, MatrixXdr cVal) {

		Yval.resize(yVal.rows(), 1);
		Eigen::VectorXd _ym = yVal.colwise().sum() / (1.0 * yVal.rows());
		Eigen::VectorXd _yss = yVal.colwise().squaredNorm().colwise().sum() / (1.0 * yVal.rows());
		Eigen::VectorXd _ysq = (_yss - _ym.cwiseProduct(_ym)).array().sqrt();
		Eigen::VectorXd _ysq_inv = 1 / _ysq.array();
		Yval = (yVal.rowwise() - _ym.transpose()) * _ysq_inv.asDiagonal();

		Wval.resize(cVal.rows(), cVal.cols());
		Eigen::VectorXd _cm = cVal.colwise().sum() / (1.0 * cVal.rows());
		Eigen::VectorXd _css = cVal.colwise().squaredNorm().colwise().sum() / (1.0 * cVal.rows());
		Eigen::VectorXd _csq = (_css - _cm.cwiseProduct(_cm)).array().sqrt();
		Eigen::VectorXd _csq_inv = 1 / _csq.array();
		Wval = (cVal.rowwise() - _cm.transpose()) * _csq_inv.asDiagonal();

		srand(goptions.GetGenericSeed());
		Nsnp = 1.0 * g.get_active_snp_number();
		Nindv = 1.0 * g.get_active_sample_size();
		round_err = goptions.GetRandHERoundErrOption();
		BlockSize = goptions.GetGenericMailmanBlockSize();

		MatrixXdr _wx(static_cast<int> (Nindv), Wval.cols());
		_wx = MatrixXdr::Constant(static_cast<int> (Nindv), Wval.cols() + 1, 1.0);
		_wx.middleCols(1, Wval.cols()) = Wval.middleCols(0, Wval.cols());
		MatrixXdr _wwInv = (_wx.transpose() * _wx).inverse();
		MatrixXdr _wy = _wx.transpose() * Yval;

		beta = _wwInv * _wy;
		cout << "beta" << endl;
		cout << beta <<endl;
		yres = Yval - _wx * beta;
		cout << yres(0, 0) << " " << yres(1, 0) << " ~~~ " << Yval(0, 0) << " " << Yval(1, 0) << endl;
		cout << yres(yres.rows() - 2, 0) << " " << yres(yres.rows() -1, 0) << " ~~~ " << Yval(yres.rows() - 2, 0) << " " << Yval(yres.rows() - 1, 0) << endl;


		cout << "Randomized HE (mailman) parameters: " <<endl;
		cout << "Rounding error is " << round_err << endl;
		cout << "Max iteration is " << goptions.GetRandHEMaxIt() << endl;
		cout << "B0 and B1: " << goptions.GetRandHEB0() << " " << goptions.GetRandHEB1() << endl;

		heReg_begin = clock();

		g.make_MailmanP();
		mailbox::setMem();

		// Numerator
		RandHE();

		std::default_random_engine generator(goptions.GetGenericSeed());

		int B0 = goptions.GetRandHEB0();
		cout << "Burning in process, B0 = " << B0 << endl;
		It_total += B0;
		cout << "Burn " << It_total << endl;
		generator.seed(goptions.GetGenericSeed() + It_total);
		iterLB(B0, generator);
		updateStats();

		long double _v = (trK4 * h2) / (eta_h2 + eta_E * (v_E / h2));
		Stop_B = _v / round_err;
		long double err = _v / (1.0 * It_total);
		cout << "err " << err << " and estimated B1 is " << Stop_B << endl;

		long double _vw = (trV * trV * trK4 * h2w) / (eta_h2w + eta_Ew * (v_Ew / h2w));
		Stop_Bw = _vw / round_err;
		long double errw = _vw / (1.0 * It_total);
		cout << "errw " << errw << " and estimated B1w is " << Stop_Bw << endl;

		long double _vw_r = (trK4 * h2_r) / (eta_h2_r + eta_E_r * (v_E_r / h2_r));
		Stop_B_r = _vw_r / round_err;
		long double err_r = _vw_r / (1.0 * It_total);
		cout << "errw " << errw << " and estimated B1w is " << Stop_B_r << endl;

		if (Stop_Bw > goptions.GetRandHEMaxIt()) {
			Stop_Bw = goptions.GetRandHEMaxIt();
			cout << "Stop_B is set to " << Stop_Bw << endl;
		}

		int B_step = goptions.GetRandHEB1();
		int B1 = Stop_Bw;
		if (B1 > It_total) {
			while (B1 > 0) {
				It_total += B_step;
				cout << "Burn " << It_total << endl;
				generator.seed(goptions.GetGenericSeed() + It_total);
				iterLB(B_step, generator);
				updateStats();
				B1 -= B_step;
				err =  (trK4 * h2) / (eta_h2 + eta_E * (v_E / h2)) / (1.0 * It_total);
				cout << "Current estimated err is " << err << endl;
				_vw = (trV * trV * trK4 * h2w) / (eta_h2w + eta_Ew * (v_Ew / h2w));
				errw = _vw / (1.0 * It_total);
				cout << "Current estimated errw is " << errw << endl;

				err_r =  (trK4 * h2) / (eta_h2_r + eta_E_r * (v_E_r / h2_r)) / (1.0 * It_total);
			}
		}

		ofstream e_file_r;
		string fname_r = goptions.GetGenericOutFile() + string(".heReg");
		e_file_r.open(fname_r.c_str());
		cout << "Total iteration " << It_total << endl;
		cout << "N M Me Me_se " << Nindv << "\t" << Nsnp << "\t" << me << "\t" << me_sd << endl;
		cout << "Adjustment he_var\the_se (rounding error %) "  << h2w << "\t" << h2_sew << " (" << errw << ")" << endl;
		cout << "Adjustment y only, he_var\the_se (rounding error %) "  << h2_r << "\t" << h2_se << " (" << err_r << ")" << endl;
		cout << "No adjustment he_var\the_se (rounding error %) " << h2 << "\t" << h2_se << " (" << err << ")" << endl;

		e_file_r << "Total iteration " << It_total << endl;
		e_file_r << "N M Me Me_se " << Nindv << "\t" << Nsnp << "\t" << me << "\t" << me_sd << endl;
		e_file_r << "Adjustment he_var\the_se (rounding error %) "  << h2w << "\t" << h2_sew << " (" << errw << ")" << endl;
		e_file_r << "Adjustment y only, he_var\the_se (rounding error %) "  << h2_r << "\t" << h2_se << " (" << err_r << ")" << endl;
		e_file_r << "No adjustment he_var\the_se (rounding error %) " << h2 << "\t" << h2_se << " (" << err << ")" <<endl;
		e_file_r.close();
	}

	~RHEreg2() {
		cout << "Finishing randomized estimation for h2 ..." << endl;
		mailbox::cleanMem();
		clock_t heReg_end = clock();

		double heReg_time = double(heReg_end - heReg_begin) / CLOCKS_PER_SEC;
		cout << "RHE time " << heReg_time << endl;
	}

private:
	void RandHE() {
/*
		for (int i = 0; i < static_cast<int> (Nindv); i++) {
			long double _g = 0;
			for (int j = 0; j < static_cast<int> (Nsnp); j++ ) {
				double __g = g.get_geno_center(j, i);
				_g += __g * __g;
			}
			trK1 += _g / Nsnp;
		}
*/
		cout << "y1~~~~~~~~~~" << endl;
		int _it = 1;
		MatrixXdr yT_X(static_cast<int> (Nsnp), _it);
		mailbox::multiply_y_pre(Yval, _it, yT_X, true); // y*X
		MatrixXdr yT_Xr(static_cast<int> (Nsnp), _it);
		mailbox::multiply_y_pre(yres, _it, yT_Xr, true); // y*X

		cout << "y2~~~~~~~~~~" << endl;
		MatrixXdr yT_XT = yT_X.transpose();
		MatrixXdr yT_X_XT(_it, static_cast<int> (Nindv));
		mailbox::multiply_y_post(yT_XT, _it, yT_X_XT, true);
		MatrixXdr yT_XT_r = yT_Xr.transpose();
		MatrixXdr yT_X_XT_r(_it, static_cast<int> (Nindv));
		mailbox::multiply_y_post(yT_XT_r, _it, yT_X_XT_r, true);

		cout << "y3~~~~~~~~~~" << endl;
		MatrixXdr yT_X_XT_T = yT_X_XT.transpose();
		MatrixXdr yT_X_XT_X(static_cast<int> (Nsnp), _it);
		mailbox::multiply_y_pre(yT_X_XT_T, _it, yT_X_XT_X, true);
		MatrixXdr yT_X_XT_T_r = yT_X_XT_r.transpose();
		MatrixXdr yT_X_XT_X_r(static_cast<int> (Nsnp), _it);
		mailbox::multiply_y_pre(yT_X_XT_T_r, _it, yT_X_XT_X_r, true);

		cout << "y4~~~~~~~~~~" << endl;
		MatrixXdr yT_X_XT_X_T = yT_X_XT_X.transpose();
		MatrixXdr yT_X_XT_X_XT(_it, static_cast<int> (Nindv));
		mailbox::multiply_y_post(yT_X_XT_X_T, _it, yT_X_XT_X_XT, true);
		MatrixXdr yT_X_XT_X_T_r = yT_X_XT_X_r.transpose();
		MatrixXdr yT_X_XT_X_XT_r(_it, static_cast<int> (Nindv));
		mailbox::multiply_y_post(yT_X_XT_X_T_r, _it, yT_X_XT_X_XT_r, true);

		yIy = (Yval.array() * Yval.array()).sum();
		yvar = yIy / (Nindv - 1.0);
		yIy_r = (yres.array() * yres.array()).sum();
		yvar_r = yIy_r / (Nindv - 1.0);
		cout << "L_Iy " << yIy << " " << yIy_r << endl;

		yK1y = (yT_X.array() * yT_X.array()).sum() / (Nsnp);
		yK1y_r = (yT_Xr.array() * yT_Xr.array()).sum() / (Nsnp);
		cout << "L_1y " << yK1y << " " << yK1y_r << endl;

		yK2y = (yT_X_XT.array() * yT_X_XT.array()).sum() / (Nsnp * Nsnp);
		yK2y_r = (yT_X_XT_r.array() * yT_X_XT_r.array()).sum() / (Nsnp * Nsnp);
		cout << "L_2y " << yK2y << " " << yK2y_r << endl;;

		yK3y = (yT_X_XT_X.array() / Nsnp * yT_X_XT_X.array() / Nsnp).sum() / Nsnp;
		yK3y_r = (yT_X_XT_X_r.array() / Nsnp * yT_X_XT_X_r.array() / Nsnp).sum() / Nsnp;
		cout << "L_3y " << yK3y << " " << yK3y_r << endl;;

//////////
		int _itW = Wval.cols(); //it possibly crashes when _itW > iteration.
		wIw = (Wval.transpose()) * Wval;
		cout << "wIw " << wIw << endl;
		wIw_inv = wIw.inverse();
		cout << "wIw_inv " << wIw_inv << endl;

		MatrixXdr wT_X(static_cast<int> (Nsnp), _itW);
		for (int _c = 0; _c * BlockSize < _itW; _c++) {
			int _b = (BlockSize * (_c + 1)) <= _itW ? BlockSize : (_itW - BlockSize * _c);
			MatrixXdr _wT_X(static_cast<int> (Nsnp), _b);

			MatrixXdr _wval = Wval.middleCols(_c * BlockSize, _b);
			mailbox::multiply_y_pre(_wval, _b, _wT_X, true); // y*X
			wT_X.middleCols(_c * BlockSize, _b) = _wT_X.middleCols(0, _b);
		}
		wK1w = wT_X.transpose() * wT_X / Nsnp;
		cout << "wK1w " << wK1w << endl;

		MatrixXdr wT_X_T = wT_X.transpose();
		MatrixXdr wT_X_T_X(_itW, static_cast<int> (Nindv));
		for (int _c = 0; _c * BlockSize < _itW; _c++) {
			int _b = (BlockSize * (_c + 1)) <= _itW ? BlockSize : (_itW - BlockSize * _c);
			MatrixXdr _wT_X_T_X(_b, static_cast<int> (Nindv));

			MatrixXdr _wT_X_T = wT_X_T.middleRows(_c * BlockSize, _b);
			mailbox::multiply_y_post(wT_X_T, _itW, wT_X_T_X, true);
			wT_X_T.middleRows(_c * BlockSize, _b) = _wT_X_T.middleRows(0, _b);
		}
		wK2w = wT_X_T_X  * wT_X_T_X.transpose() / (Nsnp * Nsnp);
		cout << "wK2w " << wK2w << endl;

		yK = yT_X_XT / Nsnp; // 1 x n
		cout << "yK " << yK.rows() << " " << yK.cols() << endl;
		yP = ((Yval.transpose() * Wval) * wIw_inv) * Wval.transpose(); // 1 x n
		cout << "yP " << yP.rows() << " " << yP.cols() << endl;

		yKK = yT_X_XT_X_XT / (Nsnp * Nsnp); // 1 x n
		cout << "yKK " << yKK.rows() << " " << yKK.cols() << endl;
		yKP = ((yK * Wval) * wIw_inv) * Wval.transpose(); // 1 x n
		cout << "yKP " << yKP.rows() << " " << yKP.cols() << endl;

		MatrixXdr yKP_T = yKP.transpose();
		yKPX.resize(static_cast<int> (Nsnp), 1);
		mailbox::multiply_y_pre(yKP_T, 1, yKPX, true);
		yKPX = yKPX.transpose() / sqrt(Nsnp); // 1 x m
		cout << "yKPX " << yKPX.rows() << " " << yKPX.cols() << endl;

		MatrixXdr yP_T = yP.transpose();
		MatrixXdr yP_X(static_cast<int> (Nsnp), 1);
		mailbox::multiply_y_pre(yP_T, 1, yP_X, true);

		MatrixXdr yP_X_T = yP_X.transpose();
		MatrixXdr yP_X_T_X(1, static_cast<int> (Nindv));
		mailbox::multiply_y_post(yP_X_T, 1, yP_X_T_X, true);
		yPK = yP_X_T_X / Nsnp; // 1 x n
		cout << "yPK " << yPK.rows() << " " << yPK.cols() << endl;
	}

	void iterLB(int it, std::default_random_engine generator) {
		MatrixXdr Bz(static_cast<int> (Nindv), it);
		std::normal_distribution<double> norm_dist(0, 1.0);
		for (int i = 0; i < Bz.cols(); i++) {
		// Copy Yval
			std::vector<double> ShufYvalVector(Yval.data(), Yval.data() + Yval.size());
			std::shuffle(ShufYvalVector.begin(), ShufYvalVector.end(), generator);
			for (int j = 0; j < Bz.rows(); j++) {
				Bz(j, i) = ShufYvalVector[j];
//				Bz(i, j) = norm_dist(generator);
			}
			cout << "shuffledYval " << Bz(0, i) << " " << Bz(1, i) << " " << ShufYvalVector[0] << " " << ShufYvalVector[1] << " " << Yval(0, 0) << " " << Yval(1, 0) << endl;
		}

		//geno_matrix * Bz; //(p x n) * (n x iter) = p x iter
		MatrixXdr T1(static_cast<int> (Nsnp), it);
		mailbox::multiply_y_pre(Bz, it, T1, true);
		MatrixXdr T1_T = T1.transpose(); //iter X p
		MatrixXdr T2(it, static_cast<int> (Nindv));
		mailbox::multiply_y_post(T1_T, it, T2, true);

		MatrixXdr T2_T = T2.transpose();
		MatrixXdr T3(static_cast<int> (Nsnp), it);
		mailbox::multiply_y_pre(T2_T, it, T3, true);

		MatrixXdr T3_T = T3.transpose();
		MatrixXdr T4(it, static_cast<int> (Nindv));
		mailbox::multiply_y_post(T3_T, it, T4, true);

		for (int i = 0; i < T4.rows(); i++) {
			long double _lb1 = (T1_T.row(i).array() * T1_T.row(i).array()).sum();
			long double _lb2 = (T2.row(i).array() * T2.row(i).array()).sum();
			long double _lb4 = (T4.row(i).array() * T4.row(i).array()).sum();
			LB1.push_back(_lb1);
			LB2.push_back(_lb2);
			LB4.push_back(_lb4);
		}

		long double lb1_sum = 0;
		long double lb1_ss = 0;
		long double lb2_sum = 0;
		long double lb2_ss = 0;
		long double lb4_sum = 0;
		long double lb4_ss = 0;
		for (int i = 0; i < LB4.size(); i++) {
			long double _lb1 = LB1[i] / Nsnp;
			lb1_sum += _lb1;
			lb1_ss += _lb1 * _lb1;

			long double _lb2 = LB2[i] / (Nsnp * Nsnp);
			lb2_sum += _lb2;
			lb2_ss += _lb2 * _lb2;

			long double _lb4 = LB4[i] / (Nsnp * Nsnp * Nsnp * Nsnp);
			lb4_sum += _lb4;
			lb4_ss += _lb4 * _lb4;
		}
		trK1 = lb1_sum / (1.0 * LB1.size());
		trK2 = lb2_sum / (1.0 * LB2.size());
		trK4 = lb4_sum / (1.0 * LB4.size());
		trK1_var = (lb1_ss - trK1 * trK1 * LB1.size()) / (1.0 * LB1.size() - 1);
		trK2_var0 = (lb2_ss - trK2 * trK2 * LB2.size()) / (1.0 * LB2.size() - 1);
		trK2_var = 2 * trK4;
		long double _lb_sd2 = sqrt(trK2_var);

		trK4_var = (lb4_ss - trK4 * trK4 * LB4.size()) / (1.0 * LB4.size() - 1);
		me = 1 / ((trK2 - Nindv) / ((Nindv + 1) * Nindv));
		me_sd = (me / Nindv * me / Nindv * _lb_sd2 / sqrt(1.0 * LB2.size()));
		cout << "trK1, trK1_var " << trK1 << " " << trK1_var << endl;
		cout << "trK2, trK2_var (trK4*2), trK2_var0 " << trK2 << " " << trK2_var << " " << trK2_var0 << endl;
		cout << "trK4, trK4_var " << trK4 << " " << trK4_var << endl;
		cout << "me, me_sd " << me << " " << me_sd << endl;
		
	}

	void updateStats() {

		cout << endl << "Updating parameters --------------- after " << It_total << " iterations." << endl;

		long double trK2_N = trK2 - Nindv;

		trV = Nindv - (wIw_inv * wIw).trace();
		trKV = trK1 - (wIw_inv * wIw).trace();
		trKVKV = trK2 + (wK1w * wIw_inv * wK1w * wIw_inv).trace() - 2 * (wK2w * wIw_inv).trace();

		eta_h2 = yK3y - 2 * yK2y + yK1y;
		eta_E = yK2y - 2 * yK1y + yIy;

		long double vh_numw1 = trV * trV * (yK3y - (yKP * yKP.transpose()).sum() - (yKK * yPK.transpose()).sum() + (yKPX * yKPX.transpose()).sum());
		long double vh_numw2 = trKV * trV * (yK2y - (yK * yPK.transpose()).sum() - (yKP * yK.transpose()).sum() + (yPK * yPK.transpose()).sum());
		long double vh_numw4 = trKV * trKV * (yK1y - (yP * yK.transpose()).sum() - (yK * yP.transpose()).sum() + (yP * yPK.transpose()).sum());
		eta_h2w = vh_numw1 - 2 * vh_numw2 + vh_numw4;
		long double ve_numw1 = trV * trV * yK2y;
		long double ve_numw2 = trKV * trV * yK1y;
		long double ve_numw4 = trKV * trKV * yIy;
		eta_Ew = ve_numw1 - 2 * ve_numw2 + ve_numw4;

		cout << "trV trKV trKVKV " << trV << " " << trKV << " " << trKVKV << endl;

		h2 = (yK1y - yIy) / trK2_N;
		v_E = (trK2 * yIy / Nindv - yK1y) / trK2_N;
		long double v_denom = trK2_N * trK2_N;
		long double v_num = 2 * (eta_h2 * h2 + eta_E * v_E);
		long double v1 = v_num / v_denom;
		long double delta = 2 * trK4 / v_denom;
		long double delta_B = delta / (It_total * 1.0);
		long double _s1 = delta_B * h2 * h2;
		long double _s2 = delta_B * delta_B * h2 * h2;
		h2_se = sqrt(v1 + _s1 + _s2);
		mse = v1 + _s1;

		cout << "v_num v_denom v1 " << v_num << " " << v_denom << " " << v1 << endl;
		cout << "delta delta_B " << delta << " " << delta_B << endl;
		cout << "eta_h2 eta_E " << eta_h2 << " " << eta_E << endl;
		cout << "h2 h2se v_E " << h2 << " " << h2_se << " " << v_E << endl;
		cout << "mse " << mse << endl;

		Ezh2 = Nindv * Nindv / (sqrt(2) * me) * (h2) / sqrt(eta_h2 * h2 + eta_E * v_E);
		zh2 = Nindv * Nindv / (sqrt(2) * me) * (h2) / sqrt(eta_h2 * h2 + eta_E * v_E + trK4 * h2 * h2 / (1.0 * It_total));
		cout << "Ezh2 zh2 " << Ezh2 << " " << zh2 << endl << endl;

		h2w = (yK1y * trV - yIy * trKV) / (trKVKV * trV - trKV * trKV);
		v_Ew = (trKVKV * yIy - trKV * yK1y) / (trV * trKVKV - trKV * trKV);
		long double v_denomw = (trV * trKVKV - trKV * trKV) * (trV * trKVKV - trKV * trKV);
		long double v_numw = 2 * (eta_h2w * h2w + eta_Ew * v_Ew);
		long double v1w = v_numw / v_denomw;
		long double deltaw = trV * trV * 2 * trK4 / v_denomw;
		long double delta_Bw = deltaw / (It_total * 1.0);
		long double _s1w = delta_Bw * h2w * h2w;
		long double _s2w = delta_Bw * delta_Bw * h2w * h2w;
		h2_sew = sqrt(v1w + _s1w);
		msew = v1w + _s1w + _s2w;
		cout << "W:v_numw v_denomw v1w " << v_numw << " " << v_denomw << " " << v1w << endl;
		cout << "W:deltaw delta_Bw " << deltaw << " " << delta_Bw << endl;
		cout << "W:eta_h2w eta_Ew " << eta_h2w << " " << eta_Ew << endl;
		cout << "W:h2w h2wse v_Ew " << h2w << " " << h2_sew << " " << v_Ew << endl;
		cout << "msew " << msew << endl << endl;

		eta_h2_r = yK3y_r - 2 * yK2y_r + yK1y_r;
		eta_E_r = yK2y_r - 2 * yK1y_r + yIy_r;

		h2_r = (yK1y_r - yIy_r) / trK2_N;
		v_E_r = (trK2 * yIy_r / Nindv - yK1y_r) / trK2_N;
		long double v_denom_r = trK2_N * trK2_N;
		long double v_num_r = 2 * (eta_h2 * h2_r + eta_E * v_E_r);
		long double v1_r = v_num_r / v_denom_r;
		long double delta_r = 2 * trK4 / v_denom_r;
		long double delta_B_r = delta_r / (It_total * 1.0);
		long double _s1_r = delta_B_r * h2_r * h2_r;
		long double _s2_r = delta_B_r * delta_B_r * h2_r * h2_r;
		h2_se_r = sqrt(v1_r + _s1_r + _s2_r);
		mse_r = v1_r + _s1_r + _s2_r;
		cout << "W_1:v_num_r v_denom_r v1_r " << v_num_r << " " << v_denom_r << " " << v1_r << endl;
		cout << "W_1:delta_1 delta_B_1 " << delta_r << " " << delta_B_r << endl;
		cout << "W_1:eta_h2_1 eta_E_1 " << eta_h2_r << " " << eta_E_r << endl;
		cout << "W_1:h2w h2wse v_E_1 " << h2_r << " " << h2_se_r << " " << v_E_r << endl;
		cout << "W_1:mse_r " << mse_r << endl;
		Ezh2_w = Nindv * Nindv / (sqrt(2) * me) * (h2) / sqrt(eta_h2_r * h2_r + eta_E_r * v_E_r);
		zh2_w = Nindv * Nindv / (sqrt(2) * me) * (h2) / sqrt(eta_h2_r * h2_r + eta_E_r * v_E + trK4 * h2 * h2 / (1.0 * It_total));
		cout << "Ezh2_w zh2_w " << Ezh2_w << " " << zh2_w << endl << endl;
/*
		h2w_1 = (yK1y * trV - yIy * trKV) / Nindv / trK2_N;
		v_Ew_1 = (trKVKV * yIy - trKV * yK1y) / Nindv / trK2_N;
		long double v_denomw_1 = trK2_N * trK2_N;
		long double v_numw_1 = 2 * (eta_h2w * h2w_1 + eta_Ew * v_Ew_1) / (Nindv * Nindv);
		long double v1w_1 = v_numw_1 / v_denomw_1;
		long double deltaw_1 = ((trV * trKVKV - trKV * trKV) / Nindv / Nindv) * ((trV * trKVKV - trKV * trKV) / Nindv / Nindv) * trV * trV * 2 * trK4 / (v_denomw_1 * v_denomw_1);
		long double delta_Bw_1 = deltaw_1 / (It_total * 1.0);
		long double _s1w_1 = delta_Bw_1 * h2w_1 * h2w_1;
		long double _s2w_1 = trV * trV * 2 * trK4 / (trK2_N * trK2_N * trK2_N * trK2_N) * h2w_1 * h2w_1 / (It_total * It_total * 1.0);
		h2_sew_1 = sqrt(v1w_1 + _s1w_1 + _s2w_1);
		cout << "W_1:v_numw v_denomw v1w " << v_numw_1 << " " << v_denomw_1 << " " << v1w_1 << endl;
		cout << "W_1:deltaw delta_Bw " << deltaw_1 << " " << delta_Bw_1 << endl;
		cout << "W_1:eta_h2w eta_Ew " << eta_h2w << " " << eta_Ew << endl;
		cout << "W_1:h2w h2wse v_Ew " << h2w_1 << " " << h2_sew_1 << " " << v_Ew_1 << endl;
		*/
	}

	long double yvar = 0;
	long double h2 = 0, h2_se = 0;
	long double v_E = 0;
	long double mse = 0;

	long double h2w = 0, h2_sew = 0;
	long double v_Ew = 0;
	long double msew = 0;

	long double yvar_r = 0;
	long double h2_r = 0, h2_se_r = 0;
	long double v_E_r = 0;
	long double mse_r = 0;

	vector<long double> LB1;
	vector<long double> LB2;
	vector<long double> LB4;
	long double me = 0;
	long double me_sd = 0;
	long double trK1 = 0;
	long double trK1_var = 0;
	long double trK2 = 0;
	long double trK2_var = 0;
	long double trK2_var0 = 0;
	long double trK4 = 0;
	long double trK4_var = 0;

	long double yK3y = 0;
	long double yK2y = 0;
	long double yK1y = 0;
	long double yIy = 0;

	long double yK3y_r = 0;
	long double yK2y_r = 0;
	long double yK1y_r = 0;
	long double yIy_r = 0;

	long double trV = 0;
	long double trKV = 0;
	long double trKVKV = 0;

	MatrixXdr wK3w; 
	MatrixXdr wK2w;
	MatrixXdr wK1w;
	MatrixXdr wIw;
	MatrixXdr wIw_inv;

	MatrixXdr yK;
	MatrixXdr yKK;
	MatrixXdr yKP;
	MatrixXdr yKPX;
	MatrixXdr yP;
	MatrixXdr yPK;

	MatrixXdr beta;
	MatrixXdr yres;
	MatrixXdr wy;

	double Nsnp;
	double Nindv;
	MatrixXdr Yval;
    MatrixXdr Wval;

	int BlockSize = 0;
	int Stop_B = 0;
	int Stop_Bw = 0;
	int Stop_B_r = 0;
	long double eta_h2, eta_E;
	long double eta_h2w, eta_Ew;
	long double eta_h2_r, eta_E_r;

	double Ezh2 = 0;
	double zh2 = 0;

	double Ezh2_w = 0;
	double zh2_w = 0;
	int It_total = 0;

	double round_err = 0.1;
	clock_t heReg_begin;
};

#endif

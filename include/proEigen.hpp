#ifndef PROEIGEN_HPP_
#define PROEIGEN_HPP_

#include "time.h"
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <math.h>
#include <algorithm>
#include <boost/math/distributions/students_t.hpp>
#include "genotype.h"
#include "Goptions.hpp"
#include "global.h"
#include "mailman.h"

using namespace std;
using namespace mailbox;

extern Goptions goptions;
extern genotype g;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;

class proEigen {

public:
	proEigen() {

	}

	MatrixXdr getEveP() {
		return eveP;
	}

	void RandPC() {
		cout << "------------Estimating eigenvectors------------" << endl;
		int k = goptions.GetGenericMailmanBlockSize();
		int Nsnp = g.get_active_snp_number();
		int Nindv = g.get_active_sample_size();
		clock_t pc_begin = clock();

		pair<double,double> prev_error = make_pair(0.0, 0.0);
		bool toStop = false;
	
		// if (convergence_limit != -1)
		if (goptions.GetPropcConvergenceLimit() != -1)
			toStop = true;

		double prevnll = 0.0;

		c.resize(Nsnp, k);
		means.resize(Nsnp, 1);
		stds.resize(Nsnp, 1);

		for (int i = 0; i < Nsnp; i++) {
			means(i, 0) = g.get_col_mean(g.act_snp[i]);
			stds(i, 0) = g.get_col_std(g.act_snp[i]);
		}

		std::default_random_engine generator(goptions.GetGenericSeed());
		std::normal_distribution<double> norm_dist(0, 1.0);
		for (int i = 0; i < c.rows(); i++)
			for (int j = 0; j < c.cols(); j++)
				c(i, j) = norm_dist(generator);
		// Initial intermediate data structures
		// Operate in blocks to improve caching
		//

		ofstream c_file;
		if (goptions.IsGenericDebug()) {
			// c_file.open((string(command_line_opts.OUTPUT_PATH) + string("cvals_orig.txt")).c_str());
			c_file.open((goptions.GetGenericOutFile() + string(".cvals_orig.txt")).c_str());
			c_file << std::setprecision(15) << c << endl;
			c_file.close();
			printf("Read Matrix\n");
		}

		cout << "Running on dataset of " << Nsnp << " SNPs and " << Nindv << " samples" << endl;

		#if SSE_SUPPORT == 1
			if (fast_mode)
				cout<<"Using Optimized SSE FastMultiply"<<endl;
		#endif

		clock_t it_begin = clock();
		for (int i = 0; i < goptions.GetPropcMaxIteration(); i++) {
			MatrixXdr c1, c2, cint, r, v;
			double a, nll;
			clock_t propc_begin0 = clock();

			cout << "------------ Begin epoch " << (i+1) << "--------------" << endl;

	//		if (accelerated_em != 0) {
			if (goptions.GetPropcAcceleratedEM() != 0) {
				#if DEBUG == 1
				if (debug) {
					print_time();
					cout << "Before EM" << endl;
				}
				#endif
				c1 = run_EM(c);
				c2 = run_EM(c1);
				#if DEBUG == 1
				if (debug) {
					print_time();
					cout << "After EM but before acceleration" << endl;
				}
				#endif
				r = c1 - c;
				v = (c2 - c1) - r;
				a = -1.0 * r.norm() / (v.norm()) ;
				if (goptions.GetPropcAcceleratedEM() == 1) {
					if (a > -1) {
						a = -1;
						cint = c2;
					} else {
						cint = c - 2 * a * r + a * a * v;
						nll = get_error_norm(cint).second;
						if ( i > 0 ) {
							while ( nll>prevnll && a < -1) {
								a = 0.5 * (a-1);
								cint = c - 2 * a * r + (a * a * v);
								nll = get_error_norm(cint).second;
							}
						}
					}
					c = cint;
				} else if (goptions.GetPropcAcceleratedEM() == 2) {
						cint = c - 2 * a * r + a * a * v;
						c = cint;
						// c = run_EM(cint);
				}
			} else {
				c = run_EM(c);
			}

			if (goptions.GetPropcAcceleratedEM() == 1 || goptions.IsPropcAccuracy() || toStop) {
				pair<double, double> e = get_error_norm(c);
				prevnll = e.second;
				if (goptions.IsPropcAccuracy())
					cout << "Iteration " << i+1 << "  " << std::setprecision(15) << e.first << "  " << e.second << endl;
				if (abs(e.first - prev_error.first) <= goptions.GetPropcConvergenceLimit()) {
					cout << "Breaking after " << i+1 << " iterations" << endl;
					break;
				}
				prev_error = e;
			}

			clock_t propc_end0 = clock();
			double propc_time0 = double(propc_end0 - propc_begin0) / CLOCKS_PER_SEC;
			cout << "*********** End epoch " << (i+1) << " in " << propc_time0 << "s." << endl;
		}

		clock_t it_end = clock();

		print_vals();

		clock_t pc_end = clock();
		double avg_it_time = double(it_end - it_begin) / (goptions.GetPropcMaxIteration() * 1.0 * CLOCKS_PER_SEC);
		double total_time = double(pc_end - pc_begin) / CLOCKS_PER_SEC;
		cout << "\nAVG Iteration Time: " << avg_it_time << "\nTotal runtime: " << total_time << endl;
	}

	MatrixXdr getEigenValue() {
		return eigenVals;
	}
private:
	MatrixXdr c; //(p,k)
	MatrixXdr means; //(p,1)
	MatrixXdr stds; //(p,1)
	MatrixXdr eveP; //(n,k) //projection matrix for EigenGWAS
	MatrixXdr eigenVals;

	pair<double, double> get_error_norm(MatrixXdr &c) {

		int k = goptions.GetGenericMailmanBlockSize();
		int Nsnp = g.get_active_snp_number();
		int Nindv = g.get_active_sample_size();
		int k_orig = goptions.GetGenericEigenvecNumber();

		HouseholderQR<MatrixXdr> qr(c);
		MatrixXdr Q;
		Q = qr.householderQ() * MatrixXdr::Identity(Nsnp, k);
		MatrixXdr q_t(k, Nsnp);
		q_t = Q.transpose();
		MatrixXdr b(k, Nindv);
		// Need this for subtracting the correct mean in case of missing data
		if (goptions.IsGenericMissing()) {
			multiply_y_post(q_t, k, b, false);
			// Just calculating b from seen data
			MatrixXdr M_temp(k, 1);
			M_temp = q_t * means;
			for (int j = 0; j < Nindv; j++) {
				MatrixXdr M_to_remove(k, 1);
				M_to_remove = MatrixXdr::Zero(k, 1);
				for (int i = 0; i < g.act_na_ind[j].size(); i++) {
					int idx = g.act_na_ind[j][i];
					vector<int>::iterator it = find(g.act_snp.begin(), g.act_snp.end(), idx);
					int _idx = *it;
//					M_to_remove = M_to_remove + (Q.row(idx).transpose() * g.get_col_mean(idx));
					M_to_remove = M_to_remove + (Q.row(_idx).transpose() * g.get_col_mean(idx));
				}
				b.col(j) -= (M_temp - M_to_remove);
			}
		} else {
			multiply_y_post(q_t, k, b, true);
		}

		JacobiSVD<MatrixXdr> b_svd(b, ComputeThinU | ComputeThinV);
		MatrixXdr u_l, d_l, v_l;
		if (goptions.IsGenericFastMode())
			u_l = b_svd.matrixU();
		else
			u_l = Q * b_svd.matrixU();
		v_l = b_svd.matrixV();
		d_l = MatrixXdr::Zero(k, k);
		for (int kk = 0; kk < k; kk++)
			d_l(kk, kk) = (b_svd.singularValues())(kk);

		MatrixXdr u_k, v_k, d_k;
		u_k = u_l.leftCols(k_orig);
		v_k = v_l.leftCols(k_orig);
		d_k = MatrixXdr::Zero(k_orig, k_orig);
		for (int kk = 0; kk < k_orig; kk++)
			d_k(kk, kk) = (b_svd.singularValues())(kk);

		MatrixXdr b_l, b_k;
		b_l = u_l * d_l * (v_l.transpose());
		b_k = u_k * d_k * (v_k.transpose());

		if (goptions.IsGenericFastMode()) {
			double temp_k = b_k.cwiseProduct(b).sum();
			double temp_l = b_l.cwiseProduct(b).sum();
			double b_knorm = b_k.norm();
			double b_lnorm = b_l.norm();
			double norm_k = (b_knorm * b_knorm) - (2 * temp_k);
			double norm_l = (b_lnorm * b_lnorm) - (2 * temp_l);	
			return make_pair(norm_k, norm_l);
		} else {
			MatrixXdr e_l(Nsnp, Nindv);
			MatrixXdr e_k(Nsnp, Nindv);
 			for (int p = 0; p < Nsnp; p++) {
				for (int n = 0; n < Nindv; n++) {
					e_l(p, n) = g.get_geno_center(g.act_snp[p], g.act_ind[n]) - b_l(p, n);
					e_k(p, n) = g.get_geno_center(g.act_snp[p], g.act_ind[n]) - b_k(p, n);
				}
			}

			double ek_norm = e_k.norm();
			double el_norm = e_l.norm();
			return make_pair(ek_norm, el_norm);
		}
	}

/* Run one iteration of EM when genotypes are not missing
 * c_orig : p X k matrix
 * Output: c_new : p X k matrix
 */

	MatrixXdr run_EM_not_missing(MatrixXdr &c_orig) {

		int k = goptions.GetGenericMailmanBlockSize();
		int Nsnp = g.get_active_snp_number();
		int Nindv = g.get_active_sample_size();

		#if DEBUG==1
			if (debug) {
				print_time();
				cout << "Enter: run_EM_not_missing" << endl;
			}
		#endif

	 	// c_temp : k X p matrix: (C^T C)^{-1} C^{T}
		MatrixXdr c_temp(k, Nsnp);
		MatrixXdr c_new(Nsnp, k);
		c_temp = ((c_orig.transpose() * c_orig).inverse()) * (c_orig.transpose());

		#if DEBUG == 1
			if (debug) {
				print_timenl();
			}
		#endif

		/*E-step: Compute X = Z G
		* G: p X n genotype matrix
		* Z: k X p matrix: (C^T C)^{-1} C^{T}
		* X: k X n matrix
 		* x_fn: X
 		* c_temp: Z
 		*/
		MatrixXdr x_fn(k, Nindv);
		multiply_y_post(c_temp, k, x_fn, true);

		#if DEBUG == 1
			if (debug) {
				print_timenl();
			}
	    #endif

		//x_temp: n X k matrix X^{T} (XX^{T})^{-1}
		MatrixXdr x_temp(Nindv, k);
		x_temp = (x_fn.transpose()) * ((x_fn * (x_fn.transpose())).inverse());

		/* M-step: X = G Z
		* G: p X n genotype matrix
		* Z: n X k matrix: X^{T}(XX^{T})^{-1}
		* X = p X k matrix
		* c_new: X
		* x_temp: Z
		*/
		multiply_y_pre(x_temp, k, c_new, true);

		#if DEBUG==1
			if(debug) {
				print_time();
				cout << "Exiting: run_EM_not_missing" << endl;
			}
		#endif
		return c_new;
	}

	MatrixXdr run_EM_missing(MatrixXdr &c_orig) {
		int k = goptions.GetGenericMailmanBlockSize();
		int Nsnp = g.get_active_snp_number();
		int Nindv = g.get_active_sample_size();
		MatrixXdr c_new(Nsnp, k);
		MatrixXdr mu(k, Nindv);

		// E step
		MatrixXdr c_temp(k, k);
		c_temp = c_orig.transpose() * c_orig;

		MatrixXdr T(k, Nindv);
		MatrixXdr c_fn;
		c_fn = c_orig.transpose();
		multiply_y_post(c_fn, k, T, false);

		MatrixXdr M_temp(k, 1);
		M_temp = c_orig.transpose() * means;

		for (int j = 0; j < Nindv; j++) {
			MatrixXdr D(k, k);
			MatrixXdr M_to_remove(k, 1);
			D = MatrixXdr::Zero(k, k);
			M_to_remove = MatrixXdr::Zero(k, 1);
			for (int i = 0; i < g.act_na_ind[j].size(); i++) {
				int idx = g.act_na_ind[j][i];
				vector<int>::iterator it = find(g.act_ind.begin(), g.act_ind.end(), idx);
				int _idx = *it;
				D = D + (c_orig.row(idx).transpose() * c_orig.row(_idx));
				M_to_remove = M_to_remove + (c_orig.row(_idx).transpose() * g.get_col_mean(idx));
			}
			mu.col(j) = (c_temp - D).inverse() * (T.col(j) - M_temp + M_to_remove);
		}

		#if DEBUG == 1
			if (debug) {
				ofstream x_file;
//				x_file.open((string(command_line_opts.OUTPUT_PATH)+string("x_in_fn_vals.txt")).c_str());
				x_file.open((goptions.GetGenericOutFile() + string("x_in_fn_vals.txt")).c_str());
				x_file<<std::setprecision(15)<<mu<<endl;
				x_file.close();
			}
		#endif

		// M step
		MatrixXdr mu_temp(k, k);
		mu_temp = mu * mu.transpose();
		MatrixXdr T1(Nsnp, k);
		MatrixXdr mu_fn;
		mu_fn = mu.transpose();
		multiply_y_pre(mu_fn, k, T1, false);
		MatrixXdr mu_sum(k, 1);
		mu_sum = MatrixXdr::Zero(k, 1);
		mu_sum = mu.rowwise().sum();

		for (int i = 0; i < Nsnp; i++) {
			MatrixXdr D(k, k);
			MatrixXdr mu_to_remove(k, 1);
			D = MatrixXdr::Zero(k, k);
			mu_to_remove = MatrixXdr::Zero(k, 1);
			for (int j = 0; j < g.act_na_snp[i].size(); j++) {
				int idx = g.act_na_snp[i][j];
				vector<int>::iterator it = find(g.act_ind.begin(), g.act_ind.end(), idx);
				int _idx = *it;
				D = D + (mu.col(_idx) * mu.col(_idx).transpose());
				mu_to_remove = mu_to_remove + (mu.col(_idx));
			}
			c_new.row(i) = (((mu_temp - D).inverse()) * (T1.row(i).transpose() - (g.get_col_mean(g.act_snp[i]) * (mu_sum - mu_to_remove)))).transpose();
			double mean;
			mean = g.get_col_sum(g.act_snp[i]);
			mean = mean -  (c_orig.row(i) * (mu_sum - mu_to_remove))(0, 0);
			mean = mean * 1.0 / (Nindv - g.act_na_snp[i].size());
			g.update_col_mean(g.act_snp[i], mean);
		}

		// IMPORTANT: Update the value of means variable present locally, so that for next iteration, updated value of means is used.
		for (int i = 0; i < Nsnp; i++) {
			means(i, 0) = g.get_col_mean(g.act_snp[i]);
			// Also updating std, just for consistency, though, it is not used presently.
			stds(i, 0) = g.get_col_std(g.act_snp[i]);
		}

		return c_new;
	}

	MatrixXdr run_EM(MatrixXdr &c_orig) {
		if (goptions.IsGenericMissing()) {
			return run_EM_missing(c_orig);
		} else {
			return run_EM_not_missing(c_orig);
		}
	}

	void print_vals() {
		int k = goptions.GetGenericMailmanBlockSize();
		int Nsnp = g.get_active_snp_number();
		int Nindv = g.get_active_sample_size();
		int k_orig = goptions.GetGenericEigenvecNumber();

		HouseholderQR<MatrixXdr> qr(c);
		MatrixXdr Q;
		Q = qr.householderQ() * MatrixXdr::Identity(Nsnp, k);
		MatrixXdr q_t(k, Nsnp);
		q_t = Q.transpose();
		MatrixXdr b(k, Nindv);

		// Need this for subtracting the correct mean in case of missing data
		if (goptions.IsGenericMissing()) {
			multiply_y_post(q_t, k, b, false);
			// Just calculating b from seen data
			MatrixXdr M_temp(k, 1);
			M_temp = q_t * means;
			for (int j = 0; j < Nindv; j++) {
				MatrixXdr M_to_remove(k, 1);
				M_to_remove = MatrixXdr::Zero(k, 1);
				for (int i = 0; i < g.act_na_ind[j].size(); i++) {
					int idx = g.act_na_ind[j][i];
					vector<int>::iterator it = find(g.act_snp.begin(), g.act_snp.end(), idx);
					int _idx = *it;
					M_to_remove = M_to_remove + (Q.row(_idx).transpose() * g.get_col_mean(idx));
				}
				b.col(j) -= (M_temp - M_to_remove);
			}
		} else {
			multiply_y_post(q_t, k, b, true);
		}

		JacobiSVD<MatrixXdr> b_svd(b, ComputeThinU | ComputeThinV);
		MatrixXdr u_l;
		u_l = b_svd.matrixU();
		MatrixXdr v_l;
		v_l = b_svd.matrixV();
		MatrixXdr u_k;
		MatrixXdr v_k, d_k;
		u_k = u_l.leftCols(k_orig);
		v_k = v_l.leftCols(k_orig);

		ofstream evec_file;
		evec_file.open((goptions.GetGenericOutFile() + string(".evecs.txt")).c_str());
		MatrixXd SNPQ = Q * u_k;
		for (int i = 0; i < g.get_active_snp_number(); i++) {
			snpInfo _snp = g.get_snp_info(g.act_snp[i]);
			evec_file << _snp.snpID << " ";
			for (int j = 0; j < k_orig; j++) {
				evec_file << std::setprecision(8) << SNPQ(i, j);
				if (j != (k_orig - 1)) {
					evec_file << " ";
				}
			}
			evec_file << endl;
		}
//		evec_file << std::setprecision(15) << Q * u_k << endl;
		evec_file.close();

		ofstream eval_file;
		eval_file.open((goptions.GetGenericOutFile() + string(".evals.txt")).c_str());
		eigenVals.resize(k_orig, 1);
		for (int kk = 0; kk < k_orig; kk++) {
			eigenVals(kk, 0) = (b_svd.singularValues())(kk) * (b_svd.singularValues())(kk)/ (1.0 * Nsnp);
			eval_file << std::setprecision(15) << eigenVals(kk, 0) <<endl;
		}
		eval_file.close();

		eveP = v_l.leftCols(k_orig);
		ofstream proj_file;
		proj_file.open((goptions.GetGenericOutFile() + string(".projections.txt")).c_str());
		for (int i = 0; i < g.get_active_sample_size(); i++) {
			indxInfo _ind = g.get_fam_info(g.act_ind[i]);
			proj_file << _ind.fid << " " << _ind.iid << " ";
			for (int j = 0; j < k_orig; j++) {
				proj_file << std::setprecision(8) << v_k(i, j);
				if (j != (k_orig - 1)) {
					proj_file << " ";
				}
			}
			proj_file << endl;
		}
//		proj_file << std::setprecision(15)<< v_k << endl;
		proj_file.close();

		if (goptions.IsGenericDebug()) {
			ofstream c_file;
			c_file.open((goptions.GetGenericOutFile()+string(".cvals.txt")).c_str());
			c_file << std::setprecision(15) << c << endl;
			c_file.close();

			ofstream means_file;
			means_file.open((goptions.GetGenericOutFile()+string(".means.txt")).c_str());
			means_file << std::setprecision(15) << means << endl;
			means_file.close();

			d_k = MatrixXdr::Zero(k_orig, k_orig);
			for (int kk =0; kk < k_orig; kk++)
				d_k(kk, kk)  =(b_svd.singularValues())(kk);
			MatrixXdr x_k;
			x_k = d_k * (v_k.transpose());
			ofstream x_file;
			x_file.open((goptions.GetGenericOutFile() + string(".xvals.txt")).c_str());
			x_file << std::setprecision(15) << x_k.transpose() << endl;
			x_file.close();
		}
	}
};

#endif


/*
MatrixXdr c; //(p,k)
MatrixXdr means; //(p,1)
MatrixXdr stds; //(p,1)
MatrixXdr eveP; //(n,k) //projection matrix for EigenGWAS

pair<double, double> get_error_norm(MatrixXdr &c) {

	int k = goptions.GetGenericMailmanBlockSize();
	int Nsnp = g.Nsnp;
	int Nindv = g.Nindv;
	int k_orig = goptions.GetGenericEigenvecNumber();

	HouseholderQR<MatrixXdr> qr(c);
	MatrixXdr Q;
	Q = qr.householderQ() * MatrixXdr::Identity(Nsnp, k);
	MatrixXdr q_t(k, Nsnp);
	q_t = Q.transpose();
	MatrixXdr b(k, Nindv);
	// Need this for subtracting the correct mean in case of missing data
	if (goptions.IsGenericMissing()) {
		multiply_y_post(q_t, k, b, false);
		// Just calculating b from seen data
		MatrixXdr M_temp(k, 1);
		M_temp = q_t * means;
		for (int j = 0; j < Nindv; j++) {
			MatrixXdr M_to_remove(k, 1);
			M_to_remove = MatrixXdr::Zero(k, 1);
			for (int i = 0; i < g.not_O_j[j].size(); i++) {
				int idx = g.not_O_j[j][i];
				M_to_remove = M_to_remove + (Q.row(idx).transpose() * g.get_col_mean(idx));
			}
			b.col(j) -= (M_temp - M_to_remove);
		}
	} else {
		multiply_y_post(q_t, k, b, true);
	}

	JacobiSVD<MatrixXdr> b_svd(b, ComputeThinU | ComputeThinV);
	MatrixXdr u_l, d_l, v_l;
	if (goptions.IsGenericFastMode())
        u_l = b_svd.matrixU();
    else
        u_l = Q * b_svd.matrixU();
	v_l = b_svd.matrixV();
	d_l = MatrixXdr::Zero(k, k);
	for (int kk = 0; kk < k; kk++)
		d_l(kk, kk) = (b_svd.singularValues())(kk);

	MatrixXdr u_k, v_k, d_k;
	u_k = u_l.leftCols(k_orig);
	v_k = v_l.leftCols(k_orig);
	d_k = MatrixXdr::Zero(k_orig, k_orig);
	for (int kk = 0; kk < k_orig; kk++)
		d_k(kk, kk) = (b_svd.singularValues())(kk);

	MatrixXdr b_l, b_k;
    b_l = u_l * d_l * (v_l.transpose());
    b_k = u_k * d_k * (v_k.transpose());

    if (goptions.IsGenericFastMode()) {
        double temp_k = b_k.cwiseProduct(b).sum();
        double temp_l = b_l.cwiseProduct(b).sum();
        double b_knorm = b_k.norm();
        double b_lnorm = b_l.norm();
        double norm_k = (b_knorm * b_knorm) - (2 * temp_k);
        double norm_l = (b_lnorm * b_lnorm) - (2 * temp_l);	
        return make_pair(norm_k, norm_l);
    } else {
        MatrixXdr e_l(Nsnp, Nindv);
        MatrixXdr e_k(Nsnp, Nindv);
        for (int p_iter = 0; p_iter < Nsnp; p_iter++) {
            for (int n_iter = 0; n_iter < Nindv; n_iter++) {
                e_l(p_iter, n_iter) = g.get_geno(p_iter, n_iter, goptions.IsGenericVarNorm()) - b_l(p_iter, n_iter);
                e_k(p_iter, n_iter) = g.get_geno(p_iter, n_iter, goptions.IsGenericVarNorm()) - b_k(p_iter, n_iter);
            }
        }

        double ek_norm = e_k.norm();
        double el_norm = e_l.norm();
        return make_pair(ek_norm, el_norm);
    }
}
*/

/* 
Run one iteration of EM when genotypes are not missing
c_orig : p X k matrix
Output: c_new : p X k matrix
MatrixXdr run_EM_not_missing(MatrixXdr &c_orig) {

	int k = goptions.GetGenericMailmanBlockSize();
	int Nsnp = g.Nsnp;
	int Nindv = g.Nindv;

	#if DEBUG==1
		if (debug) {
			print_time();
			cout << "Enter: run_EM_not_missing" << endl;
		}
	#endif

 	// c_temp : k X p matrix: (C^T C)^{-1} C^{T}
	MatrixXdr c_temp(k, Nsnp);
	MatrixXdr c_new(Nsnp, k);
	c_temp = ((c_orig.transpose() * c_orig).inverse()) * (c_orig.transpose());

	#if DEBUG == 1
		if (debug) {
			print_timenl();
		}
	#endif

 	//E-step: Compute X = Z G
 	//G: p X n genotype matrix
 	//Z: k X p matrix: (C^T C)^{-1} C^{T}
 	//X: k X n matrix
 	//x_fn: X
 	//c_temp: Z
 	//
	MatrixXdr x_fn(k, Nindv);
	multiply_y_post(c_temp, k, x_fn, true);

	#if DEBUG == 1
		if (debug) {
			print_timenl();
		}
	#endif

	//x_temp: n X k matrix X^{T} (XX^{T})^{-1}
	MatrixXdr x_temp(Nindv, k);
	x_temp = (x_fn.transpose()) * ((x_fn*(x_fn.transpose())).inverse());

	//M-step: X = G Z
 	//G: p X n genotype matrix
 	//Z: n X k matrix: X^{T}(XX^{T})^{-1}
 	//X = p X k matrix
 	//c_new: X
 	//x_temp: Z

	multiply_y_pre(x_temp, k, c_new, true);

	#if DEBUG==1
		if(debug) {
			print_time();
			cout << "Exiting: run_EM_not_missing" << endl;
		}
	#endif
	return c_new;
}
*/

/*
MatrixXdr run_EM_missing(MatrixXdr &c_orig) {
	int k = goptions.GetGenericMailmanBlockSize();
	int Nsnp = g.Nsnp;
	int Nindv = g.Nindv;
	MatrixXdr c_new(Nsnp, k);
	MatrixXdr mu(k, Nindv);

	// E step
	MatrixXdr c_temp(k, k);
	c_temp = c_orig.transpose() * c_orig;

	MatrixXdr T(k, Nindv);
	MatrixXdr c_fn;
	c_fn = c_orig.transpose();
	multiply_y_post(c_fn, k, T, false);

	MatrixXdr M_temp(k, 1);
	M_temp = c_orig.transpose() * means;

	for (int j = 0; j < Nindv; j++) {
		MatrixXdr D(k, k);
		MatrixXdr M_to_remove(k, 1);
		D = MatrixXdr::Zero(k, k);
		M_to_remove = MatrixXdr::Zero(k, 1);
		for (int i = 0; i < g.not_O_j[j].size(); i++) {
			int idx = g.not_O_j[j][i];
			D = D + (c_orig.row(idx).transpose() * c_orig.row(idx));
			M_to_remove = M_to_remove + (c_orig.row(idx).transpose() * g.get_col_mean(idx));
		}
		mu.col(j) = (c_temp - D).inverse() * (T.col(j) - M_temp + M_to_remove);
	}

	#if DEBUG == 1
		if (debug) {
			ofstream x_file;
//			x_file.open((string(command_line_opts.OUTPUT_PATH)+string("x_in_fn_vals.txt")).c_str());
			x_file.open((goptions.GetGenericOutFile() + string("x_in_fn_vals.txt")).c_str());
			x_file<<std::setprecision(15)<<mu<<endl;
			x_file.close();
		}
	#endif

	// M step
	MatrixXdr mu_temp(k, k);
	mu_temp = mu * mu.transpose();
	MatrixXdr T1(Nsnp, k);
	MatrixXdr mu_fn;
	mu_fn = mu.transpose();
	multiply_y_pre(mu_fn, k, T1, false);
	MatrixXdr mu_sum(k, 1);
	mu_sum = MatrixXdr::Zero(k, 1);
	mu_sum = mu.rowwise().sum();

	for (int i = 0; i < Nsnp; i++) {
		MatrixXdr D(k, k);
		MatrixXdr mu_to_remove(k, 1);
		D = MatrixXdr::Zero(k, k);
		mu_to_remove = MatrixXdr::Zero(k, 1);
		for (int j = 0; j < g.not_O_i[i].size(); j++) {
			int idx = g.not_O_i[i][j];
			D = D + (mu.col(idx) * mu.col(idx).transpose());
			mu_to_remove = mu_to_remove + (mu.col(idx));
		}
		c_new.row(i) = (((mu_temp-D).inverse()) * (T1.row(i).transpose() - (g.get_col_mean(i) * (mu_sum-mu_to_remove)))).transpose();
		double mean;
		mean = g.get_col_sum(i);
		mean = mean -  (c_orig.row(i) * (mu_sum-mu_to_remove))(0, 0);
		mean = mean * 1.0 / (Nindv - g.not_O_i[i].size());
		g.update_col_mean(i, mean);
	}

	// IMPORTANT: Update the value of means variable present locally, so that for next iteration, updated value of means is used.
	for (int i = 0; i < Nsnp; i++) {
		means(i, 0) = g.get_col_mean(i);
		// Also updating std, just for consistency, though, it is not used presently.
		stds(i, 0) = g.get_col_std(i);
	}

	return c_new;
}
*/

/*
MatrixXdr run_EM(MatrixXdr &c_orig) {
	if (goptions.IsGenericMissing()) {
		return run_EM_missing(c_orig);
	} else {
		return run_EM_not_missing(c_orig);
	}
}
*/

/*
void print_vals() {
	int k = goptions.GetGenericMailmanBlockSize();
	int Nsnp = g.Nsnp;
	int Nindv = g.Nindv;
	int k_orig = goptions.GetGenericEigenvecNumber();

	HouseholderQR<MatrixXdr> qr(c);
	MatrixXdr Q;
	Q = qr.householderQ() * MatrixXdr::Identity(Nsnp, k);
	MatrixXdr q_t(k, Nsnp);
	q_t = Q.transpose();
	MatrixXdr b(k, Nindv);

	// Need this for subtracting the correct mean in case of missing data
	if (goptions.IsGenericMissing()) {
		multiply_y_post(q_t, k, b, false);
		// Just calculating b from seen data
		MatrixXdr M_temp(k, 1);
		M_temp = q_t * means;
		for (int j = 0; j < Nindv; j++) {
			MatrixXdr M_to_remove(k, 1);
			M_to_remove = MatrixXdr::Zero(k, 1);
			for (int i = 0; i < g.not_O_j[j].size(); i++) {
				int idx = g.not_O_j[j][i];
				M_to_remove = M_to_remove + (Q.row(idx).transpose() * g.get_col_mean(idx));
			}
			b.col(j) -= (M_temp - M_to_remove);
		}
	} else {
		multiply_y_post(q_t, k, b, true);
	}

	JacobiSVD<MatrixXdr> b_svd(b, ComputeThinU | ComputeThinV);
	MatrixXdr u_l;
	u_l = b_svd.matrixU();
	MatrixXdr v_l;
	v_l = b_svd.matrixV();
	MatrixXdr u_k;
	MatrixXdr v_k, d_k;
	u_k = u_l.leftCols(k_orig);
	v_k = v_l.leftCols(k_orig);

	ofstream evec_file;
	evec_file.open((goptions.GetGenericOutFile() + string("evecs.txt")).c_str());
	evec_file << std::setprecision(15) << Q * u_k << endl;
	evec_file.close();
	ofstream eval_file;
	eval_file.open((goptions.GetGenericOutFile() + string("evals.txt")).c_str());
	for(int kk = 0; kk < k_orig; kk++)
		eval_file << std::setprecision(15) << (b_svd.singularValues())(kk) * (b_svd.singularValues())(kk)/g.Nsnp<<endl;
	eval_file.close();

	eveP = v_l.leftCols(k_orig);
	ofstream proj_file;
	proj_file.open((goptions.GetGenericOutFile() + string("projections.txt")).c_str());
	proj_file << std::setprecision(15)<< v_k << endl;
	proj_file.close();

	if (goptions.IsGenericDebug()) {
		ofstream c_file;
		c_file.open((goptions.GetGenericOutFile()+string("cvals.txt")).c_str());
		c_file << std::setprecision(15) << c << endl;
		c_file.close();

		ofstream means_file;
		means_file.open((goptions.GetGenericOutFile()+string("means.txt")).c_str());
		means_file << std::setprecision(15) << means << endl;
		means_file.close();

		d_k = MatrixXdr::Zero(k_orig, k_orig);
		for (int kk =0; kk < k_orig; kk++)
			d_k(kk,kk)  =(b_svd.singularValues())(kk);
		MatrixXdr x_k;
		x_k = d_k * (v_k.transpose());
		ofstream x_file;
		x_file.open((goptions.GetGenericOutFile() + string("xvals.txt")).c_str());
		x_file << std::setprecision(15) << x_k.transpose() << endl;
		x_file.close();
	}
}
*/

/*
void ProPC() {

	int k = goptions.GetGenericMailmanBlockSize();
	int Nsnp = g.Nsnp;
	clock_t pc_begin = clock();

	pair<double,double> prev_error = make_pair(0.0, 0.0);
	bool toStop = false;
	
	//	if (convergence_limit != -1)
	if (goptions.GetPropcConvergenceLimit() != -1)
		toStop = true;

	double prevnll = 0.0;

	c.resize(Nsnp, k);
	means.resize(Nsnp, 1);
	stds.resize(Nsnp, 1);
	for (int i = 0; i < Nsnp; i++) {
		means(i, 0) = g.get_col_mean(i);
		stds(i, 0) = g.get_col_std(i);
	}

	std::default_random_engine generator(goptions.GetGenericSeed());
	std::normal_distribution<double> norm_dist(0, 1.0);
	for (int i = 0; i < c.rows(); i++)
		for (int j = 0; j < c.cols(); j++)
			c(i, j) = norm_dist(generator);
	// Initial intermediate data structures
	// Operate in blocks to improve caching
	//

	ofstream c_file;
	if (goptions.IsGenericDebug()) {
		// c_file.open((string(command_line_opts.OUTPUT_PATH) + string("cvals_orig.txt")).c_str());
		c_file.open((goptions.GetGenericOutFile() + string("cvals_orig.txt")).c_str());
		c_file << std::setprecision(15) << c << endl;
		c_file.close();
		printf("Read Matrix\n");
	}

	cout << "Running on Dataset of " << g.Nsnp << " SNPs and " << g.Nindv << " Individuals" << endl;

	#if SSE_SUPPORT == 1
		if (fast_mode)
			cout<<"Using Optimized SSE FastMultiply"<<endl;
	#endif

	clock_t it_begin = clock();
	for (int i = 0; i < goptions.GetPropcMaxIteration(); i++) {
		MatrixXdr c1, c2, cint, r, v;
		double a, nll;
		if (goptions.IsGenericDebug()) {
			print_time ();
			cout << "*********** Begin epoch " << i << "***********" << endl;
		}

//		if (accelerated_em != 0) {
		if (goptions.GetPropcAcceleratedEM() != 0) {
			#if DEBUG == 1
			if (debug) {
				print_time();
				cout << "Before EM" << endl;
			}
			#endif
			c1 = run_EM(c);
			c2 = run_EM(c1);
			#if DEBUG == 1
			if (debug) {
				print_time();
				cout << "After EM but before acceleration" << endl;
			}
			#endif
			r = c1 - c;
			v = (c2 - c1) - r;
			a = -1.0 * r.norm() / (v.norm()) ;
			if (goptions.GetPropcAcceleratedEM() == 1) {
				if (a > -1) {
					a = -1;
					cint = c2;
				} else {
					cint = c - 2 * a * r + a * a * v;
					nll = get_error_norm(cint).second;
					if ( i > 0 ) {
						while ( nll>prevnll && a < -1) {
							a = 0.5 * (a-1);
							cint = c - 2 * a * r + (a * a * v);
							nll = get_error_norm(cint).second;
						}
					}
				}
				c = cint;
			} else if (goptions.GetPropcAcceleratedEM() == 2) {
				cint = c - 2 * a * r + a * a * v;
				c = cint;
				// c = run_EM(cint);
			}
		} else {
			c = run_EM(c);
		}

		if (goptions.GetPropcAcceleratedEM() == 1 || goptions.IsPropcAccuracy() || toStop) {
			pair<double, double> e = get_error_norm(c);
				prevnll = e.second;
				if (goptions.IsPropcAccuracy())
					cout << "Iteration " << i+1 << "  " << std::setprecision(15) << e.first << "  " << e.second << endl;
				if (abs(e.first - prev_error.first) <= goptions.GetPropcConvergenceLimit()) {
					cout << "Breaking after " << i+1 << " iterations" << endl;
					break;
				}
				prev_error = e;
		}

		if (goptions.IsGenericDebug()) {
				print_time();
				cout << "*********** End epoch " << i << "***********" << endl;
		}
	}

	clock_t it_end = clock();

	print_vals();

	clock_t pc_end = clock();
	double avg_it_time = double(it_end - it_begin) / (goptions.GetPropcMaxIteration() * 1.0 * CLOCKS_PER_SEC);
	double total_time = double(pc_end - pc_begin) / CLOCKS_PER_SEC;
	cout << "\nAVG Iteration Time: " << avg_it_time << "\nTotal runtime: " << total_time << endl;
}
*/

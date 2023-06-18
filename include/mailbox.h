#ifndef MAILBOX_H
#define MAILBOX_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include "mailman.h"

#if SSE_SUPPORT == 1
	#define fastmultiply fastmultiply_sse
	#define fastmultiply_pre fastmultiply_pre_sse
#else
	#define fastmultiply fastmultiply_normal
	#define fastmultiply_pre fastmultiply_pre_normal
#endif

using namespace Eigen;
using namespace std;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;

extern Goptions goptions;
extern genotype g;
extern double **partialsums;
extern double *sum_op;

extern double **yint_e;
extern double ***y_e;

extern double **yint_m;
extern double ***y_m;

/*
void multiply_y_post_naive(MatrixXdr &op, int Nrows_op, MatrixXdr &res) {
	res = op * g.geno_matrix;
}
*/

namespace mailbox {

void multiply_y_post_naive_mem(MatrixXdr &op, int Nrows_op, MatrixXdr &res) {
	int Nsnp = g.get_active_snp_number();
	int Nindv = g.get_active_sample_size();
	for (int n = 0; n < Nindv; n++) {
		for (int k = 0; k < Nrows_op; k++) {
			double temp = 0;
			for (int p = 0; p < Nsnp; p++)
				temp += op(k, p) * g.get_geno_center(g.act_snp[p], g.act_ind[n]);
			res(k, n) = temp;
		}
	}
}

void multiply_y_post_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_e, double **y_e, double *partialsums) {
	int blocksize = goptions.GetGenericMailmanBlockSize();
	for (int i = 0; i < g.get_active_sample_size(); i++) {
		memset(y_e[i], 0, blocksize * sizeof(double));
	}

	for (int seg = begin; seg < end; seg++) {
		mailman::fastmultiply_pre(g.segment_size_hori, g.get_active_sample_size(), 
			Ncol_op, seg * g.segment_size_hori, 
			g.mailman_p[seg], op, yint_e, partialsums, y_e);
	}
}

/*
 * E-step: Compute X = D Y 
 * Y : p X n genotype matrix
 * D : k X p matrix: (C^T C)^{-1} C^{T}
 * X : k X n matrix
 *
 * op_orig : D
 * Nrows_op : k
 * res : X
 * subtract_means :
 */
void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res, bool subtract_means) {

	int Nsnp = g.get_active_snp_number();
	int Nindv = g.get_active_sample_size();
	MatrixXdr op;
	op = op_orig.transpose();

	if (goptions.IsGenericVarNorm() && goptions.IsGenericSubstractMean()) {
		for (int p = 0; p < Nsnp; p++) {
			for (int k = 0; k < Nrows_op; k++)
				op(p, k) = op(p, k) / (g.get_col_std(g.act_snp[p]));
		}
	}

	#if DEBUG == 1
		if (debug) {
			print_time ();
			cout << "Starting mailman on postmultiply" << endl;
		}
	#endif

	int Ncol_op = Nrows_op;
	int nthreads = goptions.GetGenericThreads();

	nthreads = (nthreads > g.Nsegments_hori) ? g.Nsegments_hori : nthreads;

	std::thread th[nthreads];
	int perthread = g.Nsegments_hori/nthreads;
	int t = 0;
	for (; t < nthreads - 1; t++) {
		th[t] = std::thread(multiply_y_post_fast_thread, t * perthread, 
			(t+1) * perthread, std::ref(op), Ncol_op, yint_e[t], y_e[t], partialsums[t]);
	}

	th[t] = std::thread(multiply_y_post_fast_thread, t * perthread, 
		g.Nsegments_hori - 1, std::ref(op), Ncol_op, yint_e[t], y_e[t], partialsums[t]);

	for (int t = 0; t < nthreads; t++) {
		th[t].join();
	}

/*
	int seg_iter;
	for(seg_iter = 0; seg_iter < g.Nsegments_hori-1; seg_iter++){
		mailman::fastmultiply_pre (g.segment_size_hori, g.Nindv, Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter], op, yint_e, partialsums[0], y_e);
	}
*/

	for (int t = 1; t < nthreads; t++) {
		for (int n = 0; n < Nindv; n++)
			for (int k = 0; k < Ncol_op; k++)
				y_e[0][n][k] += y_e[t][n][k];
	}

	int last_seg_size = (g.get_active_snp_number() % g.segment_size_hori !=0) ? g.get_active_snp_number() % g.segment_size_hori : g.segment_size_hori;
	mailman::fastmultiply_pre(last_seg_size, g.get_active_sample_size(), Ncol_op, (g.Nsegments_hori - 1) * g.segment_size_hori, g.mailman_p[g.Nsegments_hori - 1], op, yint_e[0], partialsums[0], y_e[0]);

	for (int n = 0; n < Nindv; n++) {
		for (int k = 0; k < Ncol_op; k++) {
			res(k, n) = y_e[0][n][k];
			y_e[0][n][k] = 0;
		}
	}

	#if DEBUG == 1
		if (debug) {
			print_time (); 
			cout << "Ending mailman on postmultiply" << endl;
		}
	#endif

	if (!goptions.IsGenericSubstractMean())
		return;

	double *sums_elements = new double[Ncol_op];
 	memset(sums_elements, 0, Nrows_op * sizeof(int));

 	for (int k = 0; k < Ncol_op; k++) {
 		double sum_to_calc = 0.0;
 		for (int p = 0; p < Nsnp; p++)
 			sum_to_calc += g.get_col_mean(g.act_snp[p]) * op(p, k);
 		sums_elements[k] = sum_to_calc;
 	}

 	for (int k = 0; k < Ncol_op; k++) {
 		for (int n = 0; n < Nindv; n++)
 			res(k, n) = res(k, n) - sums_elements[k];
 	}
}

void multiply_y_post(MatrixXdr &op, int Nrows_op, MatrixXdr &res, bool subtract_means) {
	if (goptions.IsGenericFastMode())
		multiply_y_post_fast(op, Nrows_op, res, subtract_means);
	else {
		multiply_y_post_naive_mem(op, Nrows_op, res);
/*
		if (goptions.IsGenericMemoryEfficient())
			multiply_y_post_naive_mem(op, Nrows_op, res);
		else
			multiply_y_post_naive(op, Nrows_op, res);
*/
	}
}

/*
void multiply_y_pre_naive(MatrixXdr &op, int Ncol_op, MatrixXdr &res) {
	res = g.geno_matrix * op;
}
*/

void multiply_y_pre_naive_mem(MatrixXdr &op, int Ncol_op, MatrixXdr &res) {
	int Nsnp = g.get_active_snp_number();
	int Nindv = g.get_active_sample_size();
	for (int p = 0; p < Nsnp; p++) {
		for (int k = 0; k < Ncol_op; k++) {
			double temp = 0;
			for (int n = 0; n < Nindv; n++)
				temp += g.get_geno_center(g.act_snp[p], g.act_ind[n]) * op(n, k);
			res(p, k) = temp;
		}
	}
}

void multiply_y_pre_fast_thread(int begin, int end, MatrixXdr &op, int Ncol_op, double *yint_m, double **y_m, double *partialsums, MatrixXdr &res) {
	for (int seg = begin; seg < end; seg++) {
		mailman::fastmultiply(g.segment_size_hori, g.get_active_sample_size(), Ncol_op, g.mailman_p[seg], op, yint_m, partialsums, y_m);
		int p_base = seg * g.segment_size_hori;
		for (int p = p_base; (p < p_base + g.segment_size_hori) && (p < g.get_active_snp_number()); p++) {
			for (int k = 0; k < Ncol_op; k++)
				res(p, k) = y_m[p - p_base][k];
		}
	}
}

/*
 * M-step: Compute C = Y E 
 * Y: p X n genotype matrix
 * E: n K k matrix: X^{T} (XX^{T})^{-1}
 * C = p X k matrix
 *
 * op: E
 * Ncol_op: k
 * res: C
 * subtract_means:
 */
void multiply_y_pre_fast(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means) {

	int Nsnp = g.get_active_snp_number();

	for (int k = 0; k < Ncol_op; k++) {
		sum_op[k] = op.col(k).sum();
	}

	#if DEBUG == 1
		if (debug) {
			print_time();
			cout << "Starting mailman on premultiply" << endl;
			cout << "Nops = " << Ncol_op << "\t" << g.Nsegments_hori << endl;
			cout << "Segment size = " << g.segment_size_hori << endl;
			cout << "Matrix size = " << g.segment_size_hori << "\t" << g.Nindv << endl;
			cout << "op = " << op.rows () << "\t" << op.cols () << endl;
		}
	#endif

	//TODO: Memory Effecient SSE FastMultipy
	int nthreads = goptions.GetGenericThreads();
	nthreads = (nthreads > g.Nsegments_hori) ? g.Nsegments_hori : nthreads;

	std::thread th[nthreads];
	int perthread = g.Nsegments_hori / nthreads;

	int t = 0;
	for (; t < nthreads - 1; t++) {
//		cout << "Launching thread " << t << endl;
		th[t] = std::thread(multiply_y_pre_fast_thread, t * perthread, (t+1) * perthread, 
			std::ref(op), Ncol_op, yint_m[t], y_m[t], partialsums[t], std::ref(res));
	}
	
	th[t] = std::thread(multiply_y_pre_fast_thread, t * perthread, g.Nsegments_hori - 1, 
		std::ref(op), Ncol_op, yint_m[t], y_m[t], partialsums[t], std::ref(res));

	for (int t = 0; t < nthreads; t++) {
		th[t].join();
	}

/*
	for(int seg_iter = 0; seg_iter < g.Nsegments_hori - 1; seg_iter++){
		mailman::fastmultiply ( g.segment_size_hori, g.Nindv, Ncol_op, g.p[seg_iter], op, yint_m, partialsums, y_m);
		int p_base = seg_iter * g.segment_size_hori; 
		for(int p_iter=p_base; (p_iter < p_base + g.segment_size_hori) && (p_iter < g.Nsnp) ; p_iter++ ){
			for(int k_iter = 0; k_iter < Ncol_op; k_iter++) 
				res(p_iter, k_iter) = y_m [p_iter - p_base][k_iter];
		}
	}
*/

	int last_seg_size = (g.get_active_snp_number() % g.segment_size_hori !=0) ? g.get_active_snp_number() % g.segment_size_hori : g.segment_size_hori;
	mailman::fastmultiply(last_seg_size, g.get_active_sample_size(), Ncol_op, 
		g.mailman_p[g.Nsegments_hori - 1], op, yint_m[0], partialsums[0], y_m[0]);		
	int p_base = (g.Nsegments_hori - 1) * g.segment_size_hori;
	for (int p = p_base; (p < p_base + g.segment_size_hori) && (p < g.get_active_snp_number()); p++) {
		for (int k = 0; k < Ncol_op; k++)
			res(p, k) = y_m[0][p - p_base][k];
	}

	#if DEBUG == 1
		if (debug) {
			print_time (); 
			cout <<"Ending mailman on premultiply"<<endl;
		}
	#endif

	if (!goptions.IsGenericSubstractMean())
		return;

	for (int p = 0; p < Nsnp; p++) {
 		for (int k = 0; k < Ncol_op; k++) {
			res(p, k) = res(p, k) - (g.get_col_mean(g.act_snp[p]) * sum_op[k]);
			if (goptions.IsGenericVarNorm())
				res(p, k) = res(p, k) / (g.get_col_std(g.act_snp[p]));
 		}
 	}
}

//y*X
void multiply_y_pre(MatrixXdr &op, int Ncol_op, MatrixXdr &res, bool subtract_means) {
    if (goptions.IsGenericFastMode()) {
        multiply_y_pre_fast(op, Ncol_op, res, subtract_means);
	} else {
		multiply_y_pre_naive_mem(op, Ncol_op, res);
/*
		if (goptions.IsGenericMemoryEfficient())
			multiply_y_pre_naive_mem(op, Ncol_op, res);
		else
			multiply_y_pre_naive(op, Ncol_op, res);
*/
	}
}

void setMem() {
	int blocksize = goptions.GetGenericMailmanBlockSize();
	int hsegsize = g.segment_size_hori;	// = log_3(n)
	int hsize = pow(3, hsegsize);
//	int vsegsize = g.segment_size_ver;	// = log_3(p)
//	int vsize = pow(3, vsegsize);

	sum_op = new double[blocksize];
	partialsums = new double*[goptions.GetGenericThreads()];
	yint_m = new double*[goptions.GetGenericThreads()];
	yint_e = new double*[goptions.GetGenericThreads()];

	for (int t = 0; t < goptions.GetGenericThreads(); t++) {
		partialsums[t] = new double [blocksize];
		yint_m[t] = new double [hsize * blocksize];
		memset(yint_m[t], 0, hsize * blocksize * sizeof(double));
		yint_e[t] = new double [hsize * blocksize];
		memset(yint_e[t], 0, hsize * blocksize * sizeof(double));
	}

	y_e = new double**[goptions.GetGenericThreads()];
	y_m = new double**[goptions.GetGenericThreads()];
	for (int t = 0; t < goptions.GetGenericThreads(); t++) {
		y_e[t] = new double*[g.get_active_sample_size()];
		for (int i = 0; i < g.get_active_sample_size(); i++) {
			y_e[t][i] = new double[blocksize];
			memset(y_e[t][i], 0, blocksize * sizeof(double));
		}
		y_m[t] = new double*[hsegsize];
		for (int i = 0; i < hsegsize; i++) {
			y_m[t][i] = new double[blocksize];
			memset(y_m[t][i], 0, blocksize * sizeof(double));
		}
	}
}

void cleanMem() {
	int nthreads = goptions.GetGenericThreads();
	int hsegsize = g.segment_size_hori;	// = log_3(n)
	delete[] sum_op;
	for (int t = 0; t < nthreads; t++) {
		delete[] yint_e[t];
	}
	delete[] yint_e;

	for (int t = 0; t < nthreads; t++) {
		delete[] yint_m[t];
		delete[] partialsums[t];
	}
	delete[] yint_m;
	delete[] partialsums;

	for (int t = 0; t < nthreads; t++) {
		for (int i  = 0; i < hsegsize; i++)
			delete[] y_m[t][i];
		delete[] y_m[t];
	}
	delete[] y_m;

	for (int t = 0; t < nthreads; t++) {
		for (int i  = 0; i < g.get_active_sample_size(); i++)
			delete[] y_e[t][i]; 
		delete[] y_e[t];
	}
	delete[] y_e;
}

} //namespace
#endif

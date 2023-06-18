#ifndef XLDGRM_HPP_
#define XLDGRM_HPP_

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

struct rSqUnit_grm {
	double rsq;
	double sd;
	int mi;
	int mj;
	string p1;

	double rsq_s;
	double sd_s;
	string ps;
	rSqUnit_grm() {
		rsq = 0.0;
		sd = 0.0;
		mi = 0;
		mj = 0;
		p1 = "0";

		rsq_s = -1;
		sd_s = -1;
		ps = "NA";
	}
};

void compute_rsq(int start, int end, int itB, vector<int>& snp_sizeVec, vector<double>& grmM, vector<double>& grmSq, vector<VectorXd>& grmVec, vector<rSqUnit_grm>& rsqVec){
	int si = floor(sqrt(start * 2));
	if(si * (si + 1) / 2 > start) si--;
	int sj = start - (si * (si + 1) / 2);
	double Nind = g.act_ind.size();
	double f = (Nind * (Nind - 1.0) / 2.0);

	using boost::math::students_t;
	students_t Tdist(itB - 1);
	using boost::math::complement;
    
	int index = start, i = si, j = sj;
	while(index < end){
        rSqUnit_grm ru;
		double grmSqSum = i == j ? grmSq[i] : (grmVec[i].cwiseProduct(grmVec[j]).sum() / f);
		double grmM1 = grmM[i];
		double grmSq1 = grmSq[i];
		double grmM2 = grmM[j];
		double grmSq2 = grmSq[j];
		double v1 = grmSq1 - grmM1 * grmM1;
		double v2 = grmSq2 - grmM2 * grmM2;

		double cv = grmSqSum - grmM1 * grmM2;
		double v3 = 2 / (Nind * (Nind - 1)) * (v1 * v2 + cv * cv);
		double rsq = grmSqSum;
		double tstat = rsq / sqrt(v3);
		//cout << " tstat: "<<tstat << " rsq: "<< rsq<< " v3:"<< v3<<" v1:"<<v1<<" v2:"<<v2 <<" cv:"<<cv<<" grmSqSum:"<<grmSqSum<< endl;
		string p_tstat = boost::to_string(cdf(complement(Tdist, tstat)));
		if (cdf(complement(Tdist, tstat)) == 0) p_tstat = plinknum::tstat(itB - 1, tstat);

		ru.rsq = rsq;
		ru.sd = sqrt(v3);
		ru.mi = snp_sizeVec[i];
		ru.mj = snp_sizeVec[j];
		ru.p1 = p_tstat;
		rsqVec[(i * (i+1) / 2) + j] = ru;

		j++;
		index++;
		if(j > i){
            i++;
			j = 0;
		}
	}
    
	index = start;
	i = si;
	j = sj;
	double v1 = grmSq[i] - grmM[i] * grmM[i];
	int r1_mi = snp_sizeVec[i];
	double r1_rsq = grmSq[i];
	while(index < end){
		double v2 = grmSq[j] - grmM[j] * grmM[j];
		int r2_mj = snp_sizeVec[j];
		double r2_rsq = grmSq[j];
		int idx_ij = i * (i + 1) / 2 + j;
		rSqUnit_grm r12 = rsqVec[idx_ij];
		double v12 = (r12.sd * r12.sd) * (Nind * (Nind - 1) / 2) - v1 * v2;
		string p_tstat = "NA";
		if (i != j) {
			double r12_s = -1;
			double v12_s = 0;
			double sd_s = -1;
			double tstat = 0;
			bool _valid = true;
			if (((r1_mi * r1_rsq - 1) / (r1_mi - 1)) > 0 && ((r2_mj * r2_rsq - 1) / (r2_mj - 1)) > 0) {
				r12_s =  r12.rsq / sqrt(((r1_rsq * r1_mi - 1) / (r1_mi - 1)) * ((r2_rsq * r2_mj - 1) / (r2_mj - 1)));
			} else {
				_valid = false;
			}

			v12_s = 2 * r12_s * r12_s / (Nind * (Nind - 1)) * (v1 * v2 / v12 + v12 / (v1 * v2) - 2);
			if (v12_s < 0 || r12_s < 0) {
				_valid = false;
			}
			
			if (_valid) {
				rsqVec[idx_ij].rsq_s = r12_s;
				rsqVec[idx_ij].sd_s = sqrt(v12_s);
				tstat = r12_s / sqrt(v12_s);
				p_tstat = boost::to_string(cdf(complement(Tdist, tstat)));
				if (cdf(complement(Tdist, tstat)) == 0) p_tstat = plinknum::tstat(itB - 1, tstat);
				rsqVec[idx_ij].ps = p_tstat;
			} else {
				rsqVec[idx_ij].rsq_s = -1;
				rsqVec[idx_ij].sd_s = -1;
				rsqVec[idx_ij].ps = "NA";
			}
		} else {
			rsqVec[idx_ij].rsq_s = 1;
			rsqVec[idx_ij].sd_s = 0;
			rsqVec[idx_ij].ps = "0";
		}

		j++;
		index++;
		if(j > i){
            i++;
			j = 0;
			v1 = grmSq[i] - grmM[i] * grmM[i];
	        r1_mi = snp_sizeVec[i];
	        r1_rsq = grmSq[i];
		}
	}
}

class xLDgrm {
public:
	xLDgrm() {
		xldgrm_begin = clock();

		cout << "Starting determistic XLD algorithm ..." << endl;
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
			cout << "Validated " << tag_nameVec.size() << " tags for xld" << endl;
			cout << endl;
		}
        rsqVec = vector<rSqUnit_grm>(tag_nameVec.size() * (tag_nameVec.size() + 1) / 2);
		itB = goptions.GetGenericIteration();
		generateGRM();
		calXld();
	}

	void calXld() {
		using boost::math::students_t;
		students_t Tdist(itB - 1);
		using boost::math::complement;

		double Nind = g.act_ind.size();
		double f = (Nind * (Nind - 1.0) / 2.0);
		vector<double> grmM(tag_nameVec.size());
		vector<double> grmSq(tag_nameVec.size());
        for (int i = 0; i < tag_nameVec.size(); i++){
			grmM[i] = grmVec[i].sum() / f;
			grmSq[i] = grmVec[i].cwiseProduct(grmVec[i]).sum() / f;
		}

		int nthreads = goptions.GetGenericThreads();
		int all_size = tag_nameVec.size() * (tag_nameVec.size()+1) / 2;
		nthreads = min(nthreads, all_size);
		std::thread th[nthreads];
		int perthread = all_size / nthreads;
		int t = 0;
		for (; t < nthreads-1; t++) {
			th[t] = std::thread(compute_rsq, t * perthread, (t+1) * perthread, itB, std::ref(snp_sizeVec), std::ref(grmM), std::ref(grmSq), std::ref(grmVec), std::ref(rsqVec));
		}
		th[t] = std::thread(compute_rsq, t * perthread, all_size, itB, std::ref(snp_sizeVec), std::ref(grmM), std::ref(grmSq), std::ref(grmVec), std::ref(rsqVec));
		for (int t = 0; t < nthreads; t++) {
			th[t].join();
		}

		ofstream e_file_r;
		string fname_r = goptions.GetGenericOutFile() + string(".xld");
		e_file_r.open(fname_r.c_str());
		cout << "Writing grm-based x-ld results to " << fname_r.c_str() << " ..." << endl;
		e_file_r << "Tagi\tSNPi\tTagj\tSNPj\trsq\tsd\tpval\trsq_s\tsd_s\tpval_s" << endl;
		cout << "Tagi\tSNPi\tTagj\tSNPj\trsq\tsd\tpval\trsq_s\tsd_s\tpval_s" << endl;

		for (int i = 0; i < tag_nameVec.size(); i++) {
			for (int j = 0; j <= i; j++) {
				rSqUnit_grm rsq = rsqVec[i  * (i + 1) / 2 + j];
				e_file_r << tag_nameVec[i] << "\t" << rsq.mi << "\t" << tag_nameVec[j] << "\t" << rsq.mj << "\t" 
				<< rsq.rsq << "\t" << rsq.sd << "\t" << rsq.p1 << "\t"
				<< rsq.rsq_s << "\t" << rsq.sd_s << "\t" << rsq.ps << endl;

				cout << tag_nameVec[i] << "\t" << rsq.mi << "\t" << tag_nameVec[j] << "\t" << rsq.mj << "\t" 
				<< rsq.rsq << "\t" << rsq.sd << "\t" << rsq.p1 << "\t"
				<< rsq.rsq_s << "\t" << rsq.sd_s << "\t" << rsq.ps << endl;
			}
		}
		e_file_r.close();
	}

	~xLDgrm() {
		cout << "Finishing determistic xLD algorithm ..." << endl;
		mailbox::cleanMem();
		clock_t xldgrm_end = clock();
		double xldgrm_time = double(xldgrm_end - xldgrm_begin) / CLOCKS_PER_SEC;
		cout << "xld total time " << xldgrm_time << "s." << endl;
	}

private:
	int itB;
	vector<VectorXd> grmVec;
	vector<int> snp_sizeVec;
	vector<string> tag_nameVec;
	vector<rSqUnit_grm> rsqVec;
	clock_t xldgrm_begin;

	void generateGRM() {
		int _Nind = g.act_ind.size();

		VectorXd grm(((_Nind-1) * _Nind) / 2);
		int grm_index = 0;
		MatrixXdr gSet;
		for (int c = 0; c < tag_nameVec.size(); c++) {
			vector<string> _tag;
			_tag.push_back(tag_nameVec[c]);
			g.generate_active_snp_sample(_tag);
			int _Nsnp = g.get_active_snp_number();
			snp_sizeVec.push_back(_Nsnp);

			g.make_MailmanP();
			mailbox::setMem();
			grm_index = 0;

			int i = 0;
			for (; i+itB < _Nind; i+=itB) {
				gSet.resize(itB, _Nsnp);
				for (int _i = 0; _i < itB; _i++) {
					for (int _j = 0; _j < _Nsnp; _j++) {
						gSet(_i, _j) = g.get_geno_center(g.act_snp[_j], g.act_ind[i + _i]);
					}
				}

				MatrixXdr gXset(itB, _Nind);
				mailbox::multiply_y_post(gSet, itB, gXset, true);
				gXset /= (1.0 * _Nsnp);

				for (int _i = 0; _i < itB; _i++) {
					for (int _j = 0; _j < i+_i; _j++) {
						grm(grm_index) = gXset(_i, _j);
						grm_index++;
					}
				}
			}

			int _itB = _Nind - i;
			gSet.resize(_itB, _Nsnp);
			for (int _i = 0; _i < _itB; _i++) {
				for (int _j = 0; _j < _Nsnp; _j++) {
					gSet(_i, _j) = g.get_geno_center(g.act_snp[_j], g.act_ind[_i + i]);
				}
			}

			MatrixXdr gXset(_itB, _Nind);
			mailbox::multiply_y_post(gSet, _itB, gXset, true);
			gXset /= (1.0 * _Nsnp);

			for (int _i = 0; _i < _itB; _i++) {
				for (int _j = 0; _j < i+_i; _j++) {
					grm(grm_index) = gXset(_i, _j);
					grm_index++;
				}
			}
			grmVec.push_back(grm);
		}
	}
};

#endif

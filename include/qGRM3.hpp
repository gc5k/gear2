#ifndef QGRM3_HPP_
#define QGRM3_HPP_

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

class qGRM3 {

public:
	qGRM3() {
/*
		clock_t qGRM_begin0 = clock();

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
		itB = goptions.GetGenericIteration();
		generateGRM();
		clock_t qGRM_end0 = clock();
		double qGRM_time0 = double(qGRM_end0 - qGRM_begin0) / CLOCKS_PER_SEC;
		cout << "qGRM total time0 " << qGRM_time0 << "s." << endl;
		mailbox::cleanMem();
*/
		bed_transpose();

		bed_preset();
		int gByteCol = ceil((g.act_snp.size() * 1.0) / 4);

        GRM = new double *[g.act_ind.size()];
        for (int i = 0; i < g.act_ind.size(); i++) {
            GRM[i] = new double [g.act_ind.size()]();
        }

		cout << "Calculating GRM" << endl;
		clock_t qGRM_begin1 = clock();
		for (int i = 0; i < gByteCol; i++) {
			update_grm_unit(i);
			for (int l1 = 0; l1 < g.act_ind.size(); l1++) {
				for (int l2 = 0; l2 <= l1; l2++) {
					GRM[l1][l2] += grmUnit[gBedT[l1][i]][gBedT[l2][i]];
				}
			}
		}
/*
		for (int i = 0; i < g.act_ind.size(); i++) {
			cout << i << " " << GRM[i][i]/g.act_snp.size() << endl;
		}
*/
		clock_t qGRM_end1 = clock();
		double qGRM_time1 = double(qGRM_end1 - qGRM_begin1) / CLOCKS_PER_SEC;
		cout << "qGRM total time " << qGRM_time1 << "s." << endl;

	}

	~qGRM3() {
		int gByteCol = ceil((g.act_snp.size() * 1.0) / 4);

		for (int i = 0; i < g.act_ind.size(); i++) {
			delete[] gBedT[i];
			delete[] GRM[i];
		}
		delete[] gBedT;
		delete[] GRM;

		for (int i = 0; i < 256; i++) {
			delete[] grmUnit[i];
		}
		delete[] grmUnit;
	}

private:
	void bed_transpose() {
		cout << "Transposing genotypes ..." << endl;
		int _n = g.act_ind.size();
		int _m = g.act_snp.size();
		int gByteCol = ceil((g.act_snp.size() * 1.0 / 4));
		gBedT = new unsigned char *[_n];
		for (int i = 0; i < _n; i++) {
			gBedT[i] = new unsigned char [gByteCol]();
		}

		for (int i = 0; i < _n; i++) {
			unsigned char _byte = 0;
			int _shift = 0;

			for (int j = 0; j < _m; j++) {
				unsigned char _v = 0;
				if (g.is_geno_NA(g.act_snp[j], g.act_ind[i])) {
					_v = 1;
				} else {
					_v = g.get_geno_add(g.act_snp[j], g.act_ind[i]);
					_v = _v == 0 ? _v : (_v + 1);
				}
				_byte += _v << (_shift * 2);
				_shift++;

				if (_shift == 4) {
					int byte_index = j / 4;
					gBedT[i][byte_index] = _byte;
					_shift = 0;
					_byte = 0;
				}
			}
			if (_shift != 0) gBedT[i][gByteCol - 1] = _byte;
		}
		cout << "Finishing transposing ..." << endl;
	}

	void bed_preset() {
		for (int i = 0; i < 256; i++) {
			g_cnt[i][3] = i / pow(4, 3);
			g_cnt[i][2] = (i % (int) pow(4, 3)) / pow(4, 2);
			g_cnt[i][1] = (i % (int) pow(4, 2)) / pow(4, 1);
			g_cnt[i][0] = i % (int) pow(4, 1);
			if (g_cnt[i][3] != 1 && g_cnt[i][2] != 1
			&& g_cnt[i][1] != 1 && g_cnt[i][0] != 1) g_cnt_idx.push_back(i);
		}
	}

	void update_grm_unit(int idx) {
		int i_end = ((idx + 1) << 2) > g.act_snp.size() ? g.act_snp.size() : ((idx + 1) << 2);
		double _g[4][4] = { {0.0, 0.0, 0.0, 0.0},
							{0.0, 0.0, 0.0, 0.0},
							{0.0, 0.0, 0.0, 0.0},
							{0.0, 0.0, 0.0, 0.0}};
		for (int i = (idx << 2), j = 0; i < i_end; i++, j++) {
			double _f = g.get_locus_freq(g.act_snp[i]);
			double _st = g.get_col_std(g.act_snp[i]);
			_g[j][0] = (0 - 2 * _f) / _st;
			_g[j][1] = 0;
			_g[j][2] = (1 - 2 * _f) / _st;
			_g[j][3] = (2 - 2 * _f) / _st;
//			cout << _f << " " << _st << " " << _g[j][0] << " " << _g[j][1] << " " << _g[j][2] << " " << _g[j][3] << endl;
		}

		grmUnit = new double* [256];
		for (int i = 0; i < 256; i++) {
			grmUnit[i] = new double[256]();
		}

		if (goptions.GetGenericImput() > -1) {
			/*
			for (int i = 0; i < g_cnt_idx.size(); i++) {
				int _i = g_cnt_idx[i];
				for (int j = 0; j < g_cnt_idx.size(); j++) {
					int _j = g_cnt_idx[j];
					grmUnit[_i][_j] = _g[0][g_cnt[_i][0]] * _g[0][g_cnt[_j][0]]
						+ _g[1][g_cnt[_i][1]] * _g[1][g_cnt[_j][1]]
						+ _g[2][g_cnt[_i][2]] * _g[2][g_cnt[_j][2]]
						+ _g[3][g_cnt[_i][3]] * _g[3][g_cnt[_j][3]]; 
				}
			}
			*/
			for (int i = 0; i < g_cnt_idx.size(); i++) {
				int _i = g_cnt_idx[i];
				for (int j = 0; j <= i; j++) {
					int _j = g_cnt_idx[j];
					grmUnit[_i][_j] = _g[0][g_cnt[_i][0]] * _g[0][g_cnt[_j][0]]
						+ _g[1][g_cnt[_i][1]] * _g[1][g_cnt[_j][1]]
						+ _g[2][g_cnt[_i][2]] * _g[2][g_cnt[_j][2]]
						+ _g[3][g_cnt[_i][3]] * _g[3][g_cnt[_j][3]];
					grmUnit[_j][_i] = grmUnit[_i][_j];
				}
			}
		} else {
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < 256; j++) {
					grmUnit[i][j] = _g[0][g_cnt[i][0]] * _g[0][g_cnt[j][0]]
						+ _g[1][g_cnt[i][1]] * _g[1][g_cnt[j][1]]
						+ _g[2][g_cnt[i][2]] * _g[2][g_cnt[j][2]]
						+ _g[3][g_cnt[i][3]] * _g[3][g_cnt[j][3]]; 
				}
			}
		}
	}

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

	int itB;
	vector<string> tag_nameVec;
	vector<VectorXd> grmVec;
	vector<int> snp_sizeVec;

	unsigned char g_cnt[256][4];
	vector<int> g_cnt_idx;
	double** grmUnit;
	unsigned char** gBedT;
	double** GRM;
};
#endif

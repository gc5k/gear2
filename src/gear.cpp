/** 
Haha, not all of this code is written by Aman Agrawal 
 (Indian Institute of Technology, Delhi)
*/

#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <math.h>
#include <boost/math/distributions/students_t.hpp>

#include "time.h"
#include <thread>
#include <chrono>

#include "global.h"
#include "genotype.h"
#include "mailman.h"
#include "mailbox.h"

#include "Goptions.hpp"
#include "helper.h"
#include "proEigen.hpp"
#include "EigenGWAS.hpp"
#include "encDNA.hpp"
#include "encGReg.hpp"
#include "xLD.hpp"
#include "xLDgrm.hpp"
#include "me.hpp"
#include "RHEreg.hpp"
#include "RHEreg2.hpp"
#include "lineup.hpp"
#include "qGRM3.hpp"
#include "pheno.hpp"

using namespace Eigen;
using namespace std;
using namespace mailbox;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;

//options command_line_opts;
extern Goptions goptions;
extern genotype g;

//memory set
extern double **partialsums;
extern double *sum_op;

extern double **yint_e;
extern double ***y_e;
extern double **yint_m;
extern double ***y_m;

struct timespec t0;

int main(int argc, char const *argv[]) {

	pheno phe;
	pheno cov;

	try {
		goptions.ParseOptions(argc, argv);
	}
	catch (OptionsExitsProgram){}

	auto start = std::chrono::system_clock::now();
	clock_t io_begin = clock();
	clock_gettime(CLOCK_REALTIME, &t0);

	srand(goptions.GetGenericSeed());
	g.revvup(goptions.IsGenericFastMode(), goptions.IsGenericMissing(),
		goptions.IsGenericMemoryEfficient(), goptions.IsGenericVarNorm());
	g.read_plinkTrio(goptions.GetGenericGenoFile());
	
	g.QC();
	g.generate_active_snp_sample();

	if (!goptions.IsGenericFastMode() && !goptions.IsGenericMemoryEfficient()) {
		cout << "Genotype standardization..." << endl;
	}
//	qGRM3 q;

	//TODO: Implement these codes.
/*	
	if (missing && !fast_mode) {
		cout << "Missing version works only with mailman i.e. fast mode\n EXITING..." << endl;
		exit(-1);
	}
	if (fast_mode && memory_efficient) {
		cout << "Memory effecient version for mailman EM not yet implemented" << endl;
		cout << "Ignoring Memory effecient Flag" << endl;
	}
	if (missing && var_normalize) {
		cout << "Missing version works only without variance normalization\n EXITING..." << endl;
		exit(-1);
	}
*/


	clock_t io_end = clock();
	double io_time = double(io_end - io_begin) / CLOCKS_PER_SEC;
	cout << "IO Time: " << io_time << "s" << endl;

	if (goptions.CheckMakebedMasterOption()) {
		g.writeBed();
	} else if (goptions.CheckFreqMasterOption()) {
		g.writeFreq();
	} else if (goptions.CheckXLDMasterOption()) {
		if (g.act_ind.size() != g.gMinfo.fam_N) g.post_update_freq();

		if (goptions.GetXLDalg() == 0) {
			xLDgrm xld0;
		} else if (goptions.GetXLDalg() == 1 || goptions.GetXLDStage() || !goptions.GetXLDList().empty()) {
			xLD xld1;
		}
	} else if (goptions.CheckPropcMasterOption()) {
		if (g.act_ind.size() != g.gMinfo.fam_N) g.post_update_freq();
		g.make_MailmanP();
		mailbox::setMem();
		proEigen proEigenCal;
		proEigenCal.RandPC();
		mailbox::cleanMem();
	} else if (goptions.CheckEigenGWASMasterOption()) {
		if (!goptions.GetEigenGWASreadfile().empty()) {
			EigenGWAS eg(goptions.GetEigenGWASreadfile());
		} else {
			if (g.act_ind.size() != g.gMinfo.fam_N) g.post_update_freq();
			EigenGWAS eg;
		}
	} else if (goptions.CheckRandHEMasterOption()) {
		vector<string> gSubVec;
		g.get_subGInfoVec(gSubVec);
		lineup lup(gSubVec);
		vector<string> qSubVec;
		if (!goptions.GetGenericPhenoFile().empty()) {
			phe.read_q_file(goptions.GetGenericPhenoFile(), goptions.GetGenericPhenoNum());
			phe.get_NotMissSub(qSubVec);
			lup.lineup2(qSubVec);
		} else {
			cout << "no phenotype file specified." << endl;
		}
		vector<string> cSubVec;
		if (!goptions.GetGenericCovarFile().empty()) {
			cov.read_q_file(goptions.GetGenericCovarFile(), goptions.GetGenericCovarNum());
			cov.get_NotMissSub(cSubVec);
			lup.lineup2(cSubVec);
		} else {
			cout << "no covar file specified." << endl;
		}
		vector<string> chrVec;
		g.update_active(lup.get_subinterset(), chrVec);
		MatrixXdr qvec = phe.getQpheno(lup.get_subinterset(), 0);
		cout << qvec.rows() << " " << qvec.cols() << " " << qvec(0, 0) << " " << qvec(1, 0) << endl;
		cout << "active samples " << g.act_ind.size() << endl;
		if (g.act_ind.size() != g.gMinfo.fam_N) g.post_update_freq();

		if (goptions.GetGenericCovarFile().empty()) {
			RHEreg rhe(qvec, 0);
		} else {
			MatrixXdr cvec = cov.getQpheno(lup.get_subinterset(), cov.get_q_index());
			RHEreg2 rhe(qvec, 0, cvec);
		}

	} else if (goptions.CheckEncMasterOption()) {
		encDNA eDNA(goptions.GetEncK());
		eDNA.encG();
		eDNA.deCapKing();
		eDNA.print_encG();
	} else if (goptions.CheckEncRegMasterOption()) {
		encGReg Greg;
		Greg.gReg();
	} else if (goptions.CheckMeMasterOption()) {
//		qGRM3 qgrm3;
//		exit(0);
		if (g.act_ind.size() != g.gMinfo.fam_N) g.post_update_freq();
		me me1;
		me1.meRandomization();
		me1.writeResults();
	}

	std::chrono::duration<double> wctduration = std::chrono::system_clock::now() - start;
	cout << "Wall clock time " << wctduration.count() << "s" << endl;
	return 0;
}

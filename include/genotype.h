#ifndef GENOTYPE_H
#define GENOTYPE_H
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <map>
#include "Goptions.hpp"

using namespace Eigen;
using namespace std;

extern Goptions goptions;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;

struct snpInfo {
	string chr;
	string snpID;
	float pos;
	int bp;
	string a1;
	string a2;
	string snpTag;
	bool keepIt;

	double Freq;
	double FreqImputed;
	double std;
	double adj_std;
	int Not_NA_cnt;

	vector<int> locusCnt; //AA, Aa, aa, NA
	bool mafFlag;
	bool hweFlag;
	bool lmissFlag;
	bool zero_stdFlag;
	vector<int> na_snp;
};

struct indxInfo {
	string sid;
	string fid;
	string iid;
	string ppid;
	string mmid;
	string sex;
	string qval;
	int cnt;
	int line_index;
	bool keepIt;
};

struct quickUnit {
	double frq;
	double std;
	int sum;
	int SqSum;
	int Not_NA_cnt;
	int NA_cnt;
	vector<int> naVec;
	vector<int> geno_cnt; //AA Aa aa NA
	quickUnit(vector<int> _geno_cnt, vector<int>_naVec) {
		copy(_geno_cnt.begin(), _geno_cnt.end(), std::back_inserter(geno_cnt));
		copy(_naVec.begin(), _naVec.end(), std::back_inserter(naVec));
		sum = geno_cnt[1] + geno_cnt[2] * 2;
		SqSum = geno_cnt[1] + geno_cnt[2] * 4;
		Not_NA_cnt = geno_cnt[0] + geno_cnt[1] + geno_cnt[2];
		NA_cnt = geno_cnt[3];
		if (Not_NA_cnt >= 2) {
			frq = sum * 0.5 / (Not_NA_cnt * 1.0);
			std = sqrt((SqSum * 1.0 - frq * frq * 4 * Not_NA_cnt) / (1.0 * Not_NA_cnt - 1));
		} else {
			frq = 0.0;
			std = 0.0;
		}
	}
};

struct genoParam {
	bool isMailman;
	bool isMemEfficient;
	bool allowMissing;
	bool isGenoNormalization;
};

struct gMatrix {
	vector<string> snpTagVec;
	vector<string> chrInfoVec;
	int bim_N;
	int fam_N;
};

class genotype {

public:
	genoParam genoPt;
	bool isRevvup;

	gMatrix gMinfo;
	int Nsegments_hori, segment_size_hori; //segment_size_ver, Nsegments_ver;

	vector<vector<int> > act_na_snp; //not_O_i (snp) in active snp matrix
	vector<vector<int> > act_na_ind; //not_O_j (individual) in active snp matrix

	vector<int> act_snp; //not_O_i
	vector<int> act_ind; //not_O_j

	vector<vector<int> > mailman_p;

	void revvup(bool isFast, bool isMissing, bool isMemEff, bool isVarNorm);
	void read_plinkTrio(string geno_file);
	void make_MailmanP();

	int get_geno_add(int snpindex, int indvindex);
	bool is_geno_NA(int snpindex, int indvindex);
	double get_geno_center(int snpindex, int indvindex);
	double get_col_mean(int snpindex);
	double get_col_sum(int snpindex);
	double get_col_std(int snpindex);
	double get_locus_freq_orig(int snpindex);
	double get_locus_freq(int snpindex);
	string get_bim_info(int snpIdx);
	snpInfo get_snp_info(int snpIdx);
	indxInfo get_fam_info(int indIdx);

	void update_col_mean(int snpindex, double value);

	int get_sample_size();
	int get_active_sample_size();
	int get_snp_number();
	int get_active_snp_number();

	void get_subGInfoVec(vector<string>& indxInfoVec);
	void update_active(vector<string> subVec, vector<string> snpTagVec);
	void writeBed();
	void writeFreq();

	void QC();
	void generate_active_snp_sample();
	void generate_active_snp_sample(vector<string> chrVec);
	double cal_freq(int snpIdx, const vector<int> indIdx);
	quickUnit cal_geno(int snpIdx, vector<int> indIdx);
	void quickFreq(string bedfilename);
	void post_update_freq();
	~genotype();

private:
	
	char PLINK_magic[3];
	bool isImputed;

	unsigned char mask;
//	vector<unsigned char> maskG;
	vector<indxInfo> indVec; //information for all samples in fam file
	vector<snpInfo> snpVec; //information for all snps in bim file
	string gfile;

	int read_fam(string famfilename);
	int read_bim(string bimfilename);
	void read_bed(string bedfilename);

	float cal_observed_freq(const unsigned char* line);
	quickUnit cal_observed_freq_pro(const unsigned char* gBlock);
	void generate_bed_mailmanP2();
	int simulate_geno_hwe(double p_j);
	int simulate_geno_obs(int snpIdx);

	void imputation();
	void setImputed();
	bool hasImputed();

	void snpFilterQC();
	void sampleFilterQC();
	void MetricQC();
	void bed_preset();
	void bed_flip_preset();
	void recode_ref_allele();

	vector<vector<int> > naVec; //bed missing genotype index
	int **geno_cnt; //genotype counts
	unsigned char *bed_flip;
	int nbyte;

	unsigned char PLINK_BED_BYTE1 = 0x6c;
	unsigned char PLINK_BED_BYTE2 = 0x1b;
	unsigned char PLINK_BED_BYTE3 = 0x1;

	int genoCnt; //genosize = 3 for outbred, =2 for inbred
	unsigned int unitPerGenos;
	int gBitSize;

	unsigned char **gBedBlock;
};

#endif

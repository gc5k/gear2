#include <bits/stdc++.h>
#include <boost/lexical_cast.hpp>
#include <map>
#include <algorithm>
#include "genotype.h"
#include "time.h"
#include "Goptions.hpp"

using namespace std;
using namespace boost;

// Functions to read binary files
template<typename T>
static std::istream & binary_read(std::istream& stream, T& value) {
	return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}

template <class T>
inline void printvector(vector<T> &t, string delim = " ", bool newline = false) {
		for (int i = 0; i < t.size(); i++)
			cout << t[i] << delim;
		if (newline)
			cout << endl;
}

template <class T>
inline void printvectornl(vector<T> &t, string delim = " ") {
	printvector(t, delim, true);
}

struct snpHelp {
	int index;
	int bp;
};

bool compareSNP (snpHelp a, snpHelp b) {
	return a.bp < b.bp;
}
genotype::~genotype() {
	cout << "Finishing genotype..." << endl;
	for (int i = 0; i < gMinfo.bim_N; i++) {
		delete[] gBedBlock[i];
	} 
	delete[] gBedBlock;

	int _len = pow(256, nbyte);
	for (int i = 0; i < _len; i++) {
		delete[] geno_cnt[i];
	}
	delete[] geno_cnt;

	if (!goptions.GetGenericForceRefAlleleFile().empty()) delete[] bed_flip;

}

void genotype::revvup(bool isFast, bool isMissing, bool isMemEff, bool isVarNorm) {
	genoPt.isMailman = isFast;
	genoPt.allowMissing = isMissing;
	genoPt.isMemEfficient = isMemEff;
	genoPt.isGenoNormalization = isVarNorm;

	isRevvup = true;
	isImputed = false;

	gBitSize = 2;
	genoCnt = 3;
	unitPerGenos = 4; //wordsize / gBitSize;
	mask = 0;
	for (int i = 0; i < gBitSize; i++)
		mask = mask | (0x1 << i);
/*
	maskG.resize(unitPerGenos);
	maskG[0] = 255 - pow(2, 0) - pow(2, 1);
	maskG[1] = 255 - pow(2, 2) - pow(2, 3);
	maskG[2] = 255 - pow(2, 4) - pow(2, 5);
	maskG[3] = 255 - pow(2, 6) - pow(2, 7);
*/
	bed_preset(); //generate two arrays
	bed_flip_preset(); //--force-ref-allele
}

void genotype::read_plinkTrio(string geno_file) {

	gfile = geno_file;
	if (!isRevvup) {
		cerr << "Haven't properly setup parameters.\n EXITING..." << endl;
		exit(-1);
	}

	clock_t plink_read_begin = clock();
	cout << "--------------------Read plink triplet files---------------" << endl;

	std::stringstream fbim;
	fbim << geno_file << ".bim";

	gMinfo.bim_N = read_bim(fbim.str());

	std::stringstream ffam;
	ffam << geno_file << ".fam";
	gMinfo.fam_N = read_fam(ffam.str());

	std::stringstream fbed;
	fbed << geno_file << ".bed";
	read_bed(fbed.str());

	clock_t plink_read_end = clock();
	double plink_io_time = double(plink_read_end - plink_read_begin) / CLOCKS_PER_SEC;
	cout << "Reading plink data in " << plink_io_time << "s." <<endl;
	//genotype::quickFreq(fbed.str());//test

	imputation();
	isImputed = true;
}

void genotype::make_MailmanP() {
	generate_bed_mailmanP2();
}

void genotype::update_active(vector<string> subVec, vector<string> snpTagVec) {
	clock_t mailman_begin0 = clock();
	vector<int>().swap(act_ind);
	if (!subVec.empty()) {
		set<string> _indSet;
		for (int i = 0; i < subVec.size(); i++) _indSet.insert(subVec[i]);
		for (int i = 0; i < indVec.size(); i++)
			if (_indSet.count(indVec[i].sid) > 0 && indVec[i].keepIt) act_ind.push_back(i);
	} else {
		for (int i = 0; i < indVec.size(); i++) if (indVec[i].keepIt) act_ind.push_back(i);
	}

	vector<int>().swap(act_snp);
	if (!snpTagVec.empty()) {
		set<string> _tagSet;
		for (int i = 0; i < snpTagVec.size(); i++) _tagSet.insert(snpTagVec[i]);
		for (int i = 0; i < snpVec.size(); i++)
			if (_tagSet.count(snpVec[i].snpTag) > 0 && snpVec[i].keepIt) act_snp.push_back(i);
	} else {
		for (int i = 0; i < snpVec.size(); i++) if (snpVec[i].keepIt) act_snp.push_back(i);
	}

	act_na_ind.resize(act_ind.size());
	act_na_snp.resize(act_snp.size());
	if (!hasImputed()) {
		for (int i = 0; i < act_snp.size(); i++) {
			for (int j = 0; j < act_ind.size(); j++) {
				bool isNA = is_geno_NA(act_snp[i], act_ind[j]);
				if (isNA) {
					act_na_snp[i].push_back(act_ind[j]);
					act_na_snp[j].push_back(act_snp[i]);
				}
			}
		}
	}

	clock_t mailman_end0 = clock();
	double mailman_time0 = double(mailman_end0 - mailman_begin0) / CLOCKS_PER_SEC;
	cout << "Updating active " << act_snp.size() <<" SNPs and " << act_ind.size() << " samples in " << mailman_time0 << "s." <<endl;
}

void genotype::generate_active_snp_sample() {
	clock_t mailman_begin0 = clock();

	vector<int>().swap(act_snp);
	for (int i = 0; i < snpVec.size(); i++) if (snpVec[i].keepIt) act_snp.push_back(i);
	cout << "Generate the list of " << act_snp.size() << " SNPs." << endl;

	vector<int>().swap(act_ind);
	for (int i = 0; i < indVec.size(); i++) if (indVec[i].keepIt) act_ind.push_back(i);
	cout << "Generate the list of " << act_ind.size() << " samples." << endl;

	act_na_ind.resize(act_ind.size());
	act_na_snp.resize(act_snp.size());
	if (!hasImputed()) {
		for (int i = 0; i < act_snp.size(); i++) {
			for (int j = 0; j < act_ind.size(); j++) {
				bool isNA = is_geno_NA(act_snp[i], act_ind[j]);
				if (isNA) {
					act_na_snp[i].push_back(act_ind[j]);
					act_na_snp[j].push_back(act_snp[i]);
				}
			}
		}
	}

	clock_t mailman_end0 = clock();
	double mailman_time0 = double(mailman_end0 - mailman_begin0) / CLOCKS_PER_SEC;
	cout << "Generating active SNPs and samples in " << mailman_time0 << "s." <<endl;
}

void genotype::generate_active_snp_sample(vector<string> chrVec) {
	clock_t mailman_begin0 = clock();

	vector<int>().swap(act_snp);
	for (int i = 0; i < snpVec.size(); i++) {
		if (snpVec[i].keepIt) {
			vector<string>::iterator result = find(chrVec.begin( ), chrVec.end( ), snpVec[i].snpTag);
			if (result != chrVec.end())	act_snp.push_back(i);
		}
	}

	cout << "Generate the list of " << act_snp.size() << " SNPs for chromosome (";
	for (int i = 0; i < chrVec.size(); i++) {
		cout << chrVec[i];
		if (i < (chrVec.size() - 1)) {
			cout << ", ";
		} else {
			cout << ")" << endl;
		}
	}

	vector<int>().swap(act_ind);
	for (int i = 0; i < indVec.size(); i++) if (indVec[i].keepIt) act_ind.push_back(i);
	cout << "Generate the list for " << act_ind.size() << " samples." << endl;

	act_na_ind.resize(act_ind.size());
	act_na_snp.resize(act_snp.size());
	if (!hasImputed()) {
		for (int i = 0; i < act_snp.size(); i++) {
			for (int j = 0; j < act_ind.size(); j++) {
				bool isNA = is_geno_NA(act_snp[i], act_ind[j]);
				if (isNA) {
					act_na_snp[i].push_back(act_ind[j]);
					act_na_snp[j].push_back(act_snp[i]);
				}
			}
		}
	}

	clock_t mailman_end0 = clock();
	double mailman_time0 = double(mailman_end0 - mailman_begin0) / CLOCKS_PER_SEC;
	cout << "Generating active SNPs and samples in " << mailman_time0 << "s." <<endl;
}

void genotype::generate_bed_mailmanP2() {
	clock_t mailman_begin0 = clock();
	int _N_act_snp = get_active_snp_number();
	int _N_act_indv = get_active_sample_size();

	//cgb note: log_3 (Nindv) = [log_a (Nindv)] / [log_a (3)]
	segment_size_hori = ceil(log(_N_act_indv * 1.0) / log(genoCnt * 1.0));
	Nsegments_hori = ceil((_N_act_snp * 1.0) / (segment_size_hori * 1.0));

	cout << "--------------------Generating Mailman A=U*P---------------" << endl;
	cout << _N_act_snp << " active SNPs, " << _N_act_indv << " active samples" << endl;
	cout << "It generates " << Nsegments_hori << " blocks, each of " << segment_size_hori << " SNPs (" << pow(genoCnt, segment_size_hori) << ")" << endl;

	if (mailman_p.size() > 0) {//for recurrent calling purpose
		mailman_p.clear();
	}
	mailman_p.resize(Nsegments_hori, vector<int>(_N_act_indv));

	for (int i = 0; i < _N_act_snp; i++) {
		int _snp_idx = act_snp[i];
		int horiz_seg_no = i / segment_size_hori;
		for (int j = 0; j < _N_act_indv; j++) {
			int _ind_idx = act_ind[j];
			int val = get_geno_add(_snp_idx, _ind_idx);	
			mailman_p[horiz_seg_no][j] = mailman_p[horiz_seg_no][j] * genoCnt + val;
		}
	}

	if (goptions.IsGenericDebug()) {
		for (int i = 0; i < 1; i++)
			for (int j = 0; j < mailman_p[i].size(); j++)
				cout << mailman_p[i][j] << " ";
			cout << endl;
	}

	clock_t mailman_end0 = clock();
	double mailman_time0 = double(mailman_end0 - mailman_begin0) / CLOCKS_PER_SEC;
	cout << "Constructed Mailman P2 matrix in " << mailman_time0 << "s." <<endl;
}

void genotype::read_bed(string bedfilename) {

	ifstream ifs(bedfilename.c_str(), ios::in|ios::binary);

	if (!ifs.is_open()) {
		cerr << "Error reading file " << bedfilename << endl;
		exit(1);
	} else {
		cout << "Reading " << bedfilename << " ..." << endl;
	}

	binary_read(ifs, PLINK_magic);
	int ByteCol = ceil((gMinfo.fam_N * 1.0) / unitPerGenos);

	gBedBlock = new unsigned char* [gMinfo.bim_N];
	for (int i = 0; i < gMinfo.bim_N; i++) {
		gBedBlock[i] = new unsigned char [ByteCol]();
	}

//	na_snp.resize(gMinfo.bim_N);
//	na_ind.resize(gMinfo.fam_N);
	vector<double> frqVec;
	double _F = 0;
	for (int i = 0; i < gMinfo.bim_N; i++) {
		ifs.read(reinterpret_cast<char*> (gBedBlock[i]), ByteCol * sizeof(unsigned char));
	}

	if (!goptions.GetGenericForceRefAlleleFile().empty()) {
		recode_ref_allele();
	}

	long missCnt = 0;
	for (int i = 0; i < gMinfo.bim_N; i++) {
		quickUnit qU = cal_observed_freq_pro(gBedBlock[i]);

		snpVec[i].Not_NA_cnt = qU.Not_NA_cnt;
		snpVec[i].Freq = qU.frq;
		snpVec[i].FreqImputed = qU.frq;

		snpVec[i].std = qU.std;
		if (goptions.GetGenericAdjVar() == 0) {//default observed variance
			snpVec[i].adj_std = snpVec[i].std;
		} else {
			if (!goptions.IsGenericInbred()) { // outbred expected var
				snpVec[i].adj_std = sqrt(2 * qU.frq * (1 - qU.frq));
			} else { //inbred expected var
				snpVec[i].adj_std = sqrt(4 * qU.frq * (1 - qU.frq));
			}
		}

		copy(qU.geno_cnt.begin(), qU.geno_cnt.end(), std::back_inserter(snpVec[i].locusCnt));
		copy(qU.naVec.begin(), qU.naVec.end(), std::back_inserter(snpVec[i].na_snp));

		missCnt += snpVec[i].na_snp.size();
		frqVec.push_back(qU.frq);
		_F += qU.frq;
	}
	ifs.close();

	cout << "Read " << gMinfo.bim_N << " SNPs and " 
	<< gMinfo.fam_N << " samples" << endl;
	cout << "Read " << missCnt << " missing genotypes (the overall missing rate is " << ((1.0 * missCnt)/(1.0 * gMinfo.fam_N * gMinfo.bim_N)) << ")" << endl;
	sort(frqVec.begin(), frqVec.end());
	_F /= (1.0 * gMinfo.bim_N);
	cout << "Frequency range (before imputation): (" << frqVec[0] << ", " << frqVec[frqVec.size()-1] << "), and mean freq is " << _F << endl;
}

int genotype::read_bim(string bimfilename) {

	ifstream inp(bimfilename);
	if (!inp.is_open()) {
		cerr << "Error reading file " << bimfilename << endl;
		exit(1);
	} else {
		cout << "Reading " << bimfilename << " ..." << endl;
	}

	int cnt_anno = 0;
	int cntLine = 0;
	map<string, int> _snpMap;
	map<string, int> _chrMap;
	vector<string> tagVec;

	string line;
	string whitespaces (" \t\f\v\n\r");
	while (std::getline(inp, line)) {
		char c = line[0];
		if (c=='#') {
			cnt_anno++;
			continue;
		}

		if (line.empty()) continue;

		cntLine++;
		line.erase(0, line.find_first_not_of(" ")); //remove head whitespace;

		std::size_t found = line.find_last_not_of(whitespaces); //remove tail whitespace;

		if (found != std::string::npos) {
			line.erase(found+1);
		}

		istringstream ss(line);

		string word;
		int cnt_len = 0;
		vector<string> strVec;

		while (ss >> word && cnt_len <6) {
			strVec.push_back(word);
			cnt_len++;
		}

		snpInfo snpIn;
		if (cnt_len == 6) {
			try {
				snpIn.chr = boost::lexical_cast<string>(strVec[0]);
			} catch (const boost::bad_lexical_cast &e) {
				cerr << "Line " << cntLine << ": " << e.what() << '\n';
			}
			string _chr = snpIn.chr;

			if (_chrMap.count(_chr) == 0) {
				_chrMap.insert(pair<string, int> (_chr, 1));
				tagVec.push_back(_chr);
			} else {
				map<string, int>::iterator iter;
				iter = _chrMap.find(_chr);
 				int val = (int) iter->second;
				val++;
				_chrMap.erase(_chr);
				_chrMap.insert(pair<string, int> (_chr, val));
			}

			snpIn.snpID = strVec[1];
			if (_snpMap.count(strVec[1]) == 0) {
				_snpMap.insert(pair<string, int> (strVec[1], 1));
			} else {
				map<string, int>::iterator iter;
				iter = _snpMap.find(strVec[1]);
				int val = (int) iter->second;
				val++;
				_snpMap.erase(strVec[1]);
				_snpMap.insert(pair<string, int> (strVec[1], val));
			}

			try {
				snpIn.pos = boost::lexical_cast<float>(strVec[2]);
			} catch (const boost::bad_lexical_cast &e) {
				cerr << e.what() << '\n';
			}
			try {
				snpIn.bp = boost::lexical_cast<int>(strVec[3]);
			} catch (const boost::bad_lexical_cast &e) {
				cerr << e.what() << '\n';
			}
			snpIn.a1 = strVec[4];
			snpIn.a2 = strVec[5];
			snpIn.snpTag = _chr;
			snpIn.keepIt = true;
			snpVec.push_back(snpIn);
		}
	}

	inp.close();

	int dup_cnt = 0;

	for (map<string, int>::iterator it = _snpMap.begin(); it != _snpMap.end(); it++) {
		if (((int) it->second) > 1) {
			cerr << "Duplicated SNPs: " << (string) it->first << ", (" << (int) it->second << " counts)." <<endl;
			dup_cnt++;
		}
	}

    if (dup_cnt > 0) {
        cerr << dup_cnt << " duplicated SNP(s)." << endl;
        exit(1);
    }
	cout << "Read " << _snpMap.size() << " SNPs from " << bimfilename << endl;
	cout << "By default snp tag is associated with chromosome id " << endl; 

	for (auto it = tagVec.begin(); it != tagVec.end(); it++) {
		gMinfo.snpTagVec.push_back((string) (*it));
		gMinfo.chrInfoVec.push_back((string) (*it));
		cout << "Chr" << (string) (*it) << " has " << (int) (_chrMap.find((string) *it))->second << " SNPs" << endl;
	}

	return _snpMap.size();
}

string genotype::get_bim_info(int snpIdx) {
	if (snpIdx > snpVec.size()) {
		cerr << "Error getting snp information:" << snpIdx << "out of range (" << snpVec.size() << ")" <<endl;
	}
	snpInfo snpIn = snpVec[snpIdx];
	string snp = boost::lexical_cast<string>(snpIn.chr) + "\t" + snpIn.snpID + "\t" + boost::lexical_cast<string> (snpIn.pos) + "\t"
	+ boost::lexical_cast<string> (snpIn.bp) + "\t" + snpIn.a1 + "\t" + snpIn.a2;
	return snp;
}

snpInfo genotype::get_snp_info(int snpIdx) {
	if (snpIdx > snpVec.size()) {
		cerr << "Error getting snp information:" << snpIdx << "out of range (" << snpVec.size() << ")" <<endl;
	}
	return snpVec[snpIdx];
}

indxInfo genotype::get_fam_info(int indIdx) {
	return indVec[indIdx];
}

int genotype::read_fam(string famfilename) {

	map<string, int> _subMap;
	ifstream inp(famfilename.c_str());
	if (!inp.is_open()) {
		cerr << "Error reading file " << famfilename.c_str() << endl;
		exit(1);
	} else {
		cout << "Reading " << famfilename.c_str() << " ..." << endl;
	}

	string line;
	string whitespaces (" \t\f\v\n\r");
	int eff_cnt = 0;
	int famFile_line_cnt = 0;
	int cnt_anno = 0;

	while (std::getline(inp, line)) {
		famFile_line_cnt++;

		char c = line[0];
		if (c == '#') {
			cnt_anno++;
			continue;
		}
		if (line.empty()) continue;

		line.erase(0, line.find_first_not_of(" ")); //remove head whitespace;
		std::size_t found = line.find_last_not_of(whitespaces); //remove tail whitespace;
		if(found != std::string::npos) {
			line.erase(found+1);
		}

		eff_cnt++;
		istringstream ss(line);

		string word;
		int cnt_len = 0;
		vector<string> strVec;
		indxInfo iInfo;
		iInfo.cnt = 0;
		while (ss >> word && cnt_len < 6) {
			strVec.push_back(word);
			cnt_len++;
			switch (cnt_len) {
				case (1):
					iInfo.fid = word;
					break;
				case (2):
					iInfo.iid = word;
					break;
				case (3):
					iInfo.ppid = word;
					break;
				case (4):
					iInfo.mmid = word;
					break;
				case (5):
					iInfo.sex = word;
					break;
				case (6):
					iInfo.qval = word;
			}
		}

		iInfo.sid = iInfo.fid + " " + iInfo.iid;
		iInfo.keepIt = true;
		if (strVec.size() == 6) {
			string sid = iInfo.fid + " " + iInfo.iid;
			int cnt = _subMap.count(sid);
			if (cnt == 0) {
				iInfo.line_index = eff_cnt;
				_subMap.insert(pair<string, int>(sid, cnt));
				indVec.push_back(iInfo);
			} else {
				map<string, int>::iterator iter = _subMap.find(sid);
				int val = iter->second;
				val++;
				_subMap.erase(sid);
				_subMap.insert(pair<string, int>(sid, val));
				cerr << "line " << eff_cnt << " duplicated id: " << sid << endl;
			}
		}
	}

	inp.close();

	if (cnt_anno > 0) {
		cout << "Read " << cnt_anno << " annotation lines ('#')." <<endl;
	}

	int dup_cnt = 0;
	for (map<string, int>::iterator iter = _subMap.begin(); iter != _subMap.end(); iter++) {
		if ((iter->second) > 1) {
			cerr << "Duplicated SNPs: " << iter->first << ", (" <<iter->second << " counts)." <<endl;
			dup_cnt++;
		}
	}

	if (dup_cnt > 0) {
		cerr << dup_cnt << " duplicated id(s)." << endl;
		exit(1);
	}
	cout << "Read " << _subMap.size() << " samples from " << famfilename << endl;
	return _subMap.size();
}

void genotype::imputation() {

	cout << "--------------------Naive Imputation---------------" << endl;
	cout << "Imputation after ";
	if (goptions.GetGenericImput() == 0) {
		cout << "observed genotypic proportions ";
	} else {
		if (!goptions.IsGenericInbred()) {
			cout << " HWE for outbred population ";
		} else {
			cout << " p vs (1-p) for inbred population " ;
		}
	}
	cout << "according to selected options." << endl;
	clock_t imputation_begin0 = clock();

	int total_cnt = 0;
	vector<double> frqVec;
	double _F = 0; 

	long int _gCnt[4] = {0, 0, 0, 0};
	for (int i = 0; i < snpVec.size(); i++) {
		for (int j = 0; j < snpVec[i].na_snp.size(); j++) {
			unsigned char val;
			if (goptions.GetGenericImput() == 0) {
				val = static_cast<unsigned char> (simulate_geno_obs(i));
			} else {
				val = static_cast<unsigned char> (simulate_geno_hwe(snpVec[i].FreqImputed));
			}

			_gCnt[val]++;

			if (val != 0) val++;
			int byte_index = snpVec[i].na_snp[j] / unitPerGenos;
			int shift_index = snpVec[i].na_snp[j] % unitPerGenos;
			gBedBlock[i][byte_index] += ((val - 1) << (shift_index << 1));
			snpVec[i].locusCnt[val == 0 ? val : (val - 1)]++;
			snpVec[i].locusCnt[3]--;
			total_cnt++;
		}
		if (snpVec[i].na_snp.size() > 0) {
			snpVec[i].Not_NA_cnt = snpVec[i].locusCnt[0] + snpVec[i].locusCnt[1] + snpVec[i].locusCnt[2];
			int _locusSum = snpVec[i].locusCnt[1] + snpVec[i].locusCnt[2] * 2;
			snpVec[i].FreqImputed = (_locusSum * 0.5) / (snpVec[i].Not_NA_cnt * 1.0);
			double _n =  snpVec[i].Not_NA_cnt * 1.0;
			double _frq = snpVec[i].FreqImputed;
			double _locusSqSum = snpVec[i].locusCnt[1] * 1.0 * 1.0 + snpVec[i].locusCnt[2] * 2.0 * 2.0;
			snpVec[i].std = sqrt( (_locusSqSum - _n * 4 * _frq * _frq) / (_n - 1.0) );
		}
		frqVec.push_back(snpVec[i].FreqImputed);
		_F += snpVec[i].FreqImputed;

		if (goptions.GetGenericAdjVar() == 0) {//default observed variance
			snpVec[i].adj_std = snpVec[i].std;
		} else {
			if (!goptions.IsGenericInbred()) { // outbred expected var
				snpVec[i].adj_std = sqrt(snpVec[i].FreqImputed * (1 - snpVec[i].FreqImputed) * 2);
			} else { //inbred expected var
				snpVec[i].adj_std = sqrt(snpVec[i].FreqImputed * (1 - snpVec[i].FreqImputed) * 4);
			}
		}

	}
	sort(frqVec.begin(), frqVec.end());
	_F /= (frqVec.size() * 1.0);
	clock_t imputation_end0 = clock();
	double imputation_time0 = double(imputation_end0 - imputation_begin0) / CLOCKS_PER_SEC;
	cout << "Frequency range (after imputation): (" << frqVec[0] << ", " << frqVec[frqVec.size()-1] << "), and mean freq is " << _F << endl;
	cout << "Imputed " << total_cnt << " genotypes (aa:Aa:AA = " << _gCnt[0] << ":" << _gCnt[1] << ":" << _gCnt[2] << ") in " 
	<< imputation_time0 << "s" << endl;
}

// Accessor Functions
int genotype::get_geno_add(int snpindex, int indvindex) {
	int byte_index = indvindex / unitPerGenos;
	int shift_index = indvindex % unitPerGenos;
	unsigned char gbyte = gBedBlock[snpindex][byte_index];
	int geno = (gbyte >> (shift_index << 1)) & mask;

	// PLINK coding: 00->0, 01->missing, 10->1, 11->2
	if (geno == 1) { //missing value, what to do?
		return 0;
	} else {
		if (geno > 1) --geno;
		return geno;
	}
}

bool genotype::is_geno_NA(int snpindex, int indvindex) {
	int byte_index = indvindex / unitPerGenos;
	int shift_index = indvindex % unitPerGenos;
	unsigned char gbyte = gBedBlock[snpindex][byte_index];
	int geno = (gbyte >> (shift_index << 1)) & mask;

	return geno == 1 ? true : false; 
}

double genotype::get_geno_center(int snpindex, int indvindex) {
	int byte_index = indvindex / unitPerGenos;
	int shift_index = indvindex % unitPerGenos;
	unsigned char gbyte = gBedBlock[snpindex][byte_index];
	unsigned char geno0 = (gbyte >> (shift_index << 1)) & mask;
	// PLINK coding: 00->0, 01->missing, 10->1, 11->2
	if (geno0 == 1) { //missing value, what to do?
		return 0;
	} else {
		double geno = (geno0 == 0? geno0 : (geno0 - 1)) * 1.0;
		geno = (geno - 2 * snpVec[snpindex].FreqImputed) / snpVec[snpindex].adj_std;
		return geno;
	}
}

double genotype::get_col_mean(int snpindex) {
	double temp = snpVec[snpindex].FreqImputed * 2;
	return temp;
}

double genotype::get_locus_freq_orig(int snpindex) {
	return snpVec[snpindex].Freq;
}

double genotype::get_locus_freq(int snpindex) {
	return snpVec[snpindex].FreqImputed;
}

double genotype::get_col_sum(int snpindex) {
	double cnt = snpVec[snpindex].locusCnt[1] * 1.0 + snpVec[snpindex].locusCnt[2] * 2.0;
	return cnt;
}

double genotype::get_col_std(int snpindex) {
	return snpVec[snpindex].adj_std;
}

// Modifier Functions
void genotype::update_col_mean(int snpindex, double value) {
	snpVec[snpindex].FreqImputed = value / 2;
}

int genotype::simulate_geno_hwe(double p) {
	double rval = rand() / double(RAND_MAX);
	double _p = 1 - p;
	if (!goptions.IsGenericInbred()) {
		double P[3] = {_p * _p, 2 * p * _p, p * p};
		if (rval < P[0])
			return 0;
		else if (rval >= P[0] && rval < (P[0] + P[1])) {
			return 1;
		} else {
			return 2;
		}
	} else {
		if (rval < _p) {
			return 0;
		} else {
			return 2;
		}
	}
}

int genotype::simulate_geno_obs(int snpIdx) {
	double _N =  (snpVec[snpIdx].locusCnt[0] + snpVec[snpIdx].locusCnt[1] + snpVec[snpIdx].locusCnt[2]) * 1.0;
	double _P[3] = { (1.0 * snpVec[snpIdx].locusCnt[0])/_N, 
	((snpVec[snpIdx].locusCnt[0] + snpVec[snpIdx].locusCnt[1]) * 1.0)/_N, 
		1.0};

	double rval = rand() / double(RAND_MAX);
	if (rval <= _P[0]) {
		return 0;
	} else if (rval > _P[0] && rval < _P[1]) {
		return 1;
	} else {
		return 2;
	}
}

void genotype::get_subGInfoVec(vector<string>& indxVec) {
	for (int i = 0; i < indVec.size(); i++) {
		indxVec.push_back(indVec[i].sid);
	}
}

void genotype::snpFilterQC() {
	cout << "Applying SNP filters (if any)..." <<endl;
	map<string, string> _tagMap;
	bool _tagFlag = false;
	set<string> _snpSet;
	bool _snpFlag = false;
	set<string> _chrSet;
	bool _chrFlag = false;

	if (!goptions.GetGenericSNPTagFile().empty()) {
		vector<string> _tagVec;
		string TagFile = goptions.GetGenericSNPTagFile();
		ifstream inp(TagFile);
		if (!inp.is_open()) {
			cerr << "Error reading file "<< TagFile << endl;
			exit(1);
		} else {
			cout << "Reading SNP tag from " << TagFile << endl;
		}

		string line;
		string whitespaces (" \t\f\v\n\r");
		int LineCnt = 0;
		while (std::getline(inp, line)) {
			char c = line[0];
			if (c == '#') {
				continue;
			}
			if (line.empty()) continue;

			line.erase(0, line.find_first_not_of(" ")); //remove head whitespace;
			std::size_t found = line.find_last_not_of(whitespaces); //remove tail whitespace;
			if (found != std::string::npos) {
				line.erase(found + 1);
			}
			LineCnt++;

			istringstream ss(line);
			string _s;
			string _t;
			string word;
			int cnt_word = 0;
			while (ss >> word && cnt_word < 2) {
				switch (cnt_word) {
					case 0:
						_s = word;
						break;
					default:
						_t = word;
						break;
				}
				cnt_word++;
			}
			if (cnt_word != 2) {
				cerr << "Line " << LineCnt << " invalid format " << endl;
				cerr << "Please make sure two columns in each line for --snp-tag " << endl;
				exit(0);
			}
			auto it = _tagMap.find(_s);
			if (it != _tagMap.end()) {
				cerr << "Line: " << LineCnt << " snp " << _s << " is associated with tag " << _t << " and " << it->second << endl;
				cerr << "SNP " << _s << " has been mapped to more than one tags" << endl;
				exit(0);
			}

			_tagMap.insert(pair<string, string>(_s, _t));
			if (find(_tagVec.begin(), _tagVec.end(), _t) == _tagVec.end()) {
				_tagVec.push_back(_t);
			}
		}
		if (LineCnt != _tagMap.size()) {
			cerr << "Read " << LineCnt << " lines but mapped to " << _tagMap.size() << " unique snp-tag pairs" << endl;
			cerr << "Make sure each SNP is assigned to one tag only." << endl;
			exit(0);
		}
		cout << "read " << _tagMap.size () << " SNP tags from " << TagFile.c_str() << endl;
		gMinfo.snpTagVec.clear(); //removed the original tag, which is chromosome id.
		for (auto it = _tagVec.begin(); it != _tagVec.end(); it++) {
			gMinfo.snpTagVec.push_back(*it);
		}
		int _tagCnt = 0;
		for (int i = 0; i < snpVec.size(); i++) {
			snpInfo sInfo = snpVec[i];
			auto it = _tagMap.find(sInfo.snpID);
			if (it != _tagMap.end()) {
				snpVec[i].snpTag = it->second;
				_tagCnt++;
			}
		}
		cout << "Updated tags for " << _tagCnt << " matched SNPs" << endl;
		_tagFlag = true;
		inp.close();
	}

	if (!goptions.GetGenericExtractFile().empty()) { // exact snp
		cout << "Exacting SNPs from " << goptions.GetGenericExtractFile() << endl;
		string ExFile = goptions.GetGenericExtractFile();
		ifstream inp(ExFile.c_str());
		if (!inp.is_open()) {
			cerr << "Error reading file "<< ExFile <<endl;
			exit(1);
		} else {
			cout << "Reading " << ExFile <<endl;
		}

		string line;
		string whitespaces (" \t\f\v\n\r");
		while (std::getline(inp, line)) {
			line.erase(0, line.find_first_not_of(" ")); //remove head whitespace;
    		std::size_t found = line.find_last_not_of(whitespaces); //remove tail whitespace;
			if (found != std::string::npos) {
				line.erase(found+1);
			}

			char c = line[0];
			if (c == '#') {
				continue;
			}
			if (line.empty()) continue;

			istringstream ss(line);

			string word;
			ss >> word;
			_snpSet.insert(word);
		}
		cout << "read " << _snpSet.size () << " unique SNPs from " << ExFile.c_str() <<endl;
		if (_snpSet.size() > 0) _snpFlag = true;
		inp.close();
	}

	if (!goptions.GetGenericChr().empty()) {
		cout << "Chromosome extracting ..." << endl;
		vector<string> chrV = goptions.GetGenericChr();
		for (int i = 0; i < chrV.size(); i++) {
			_chrSet.insert(chrV[i]);
		}
		if (_chrSet.size() > 0) {
			gMinfo.chrInfoVec.clear();
			copy(chrV.begin(), chrV.end(), std::back_inserter(gMinfo.chrInfoVec));
			_chrFlag = true;
			if (!_tagFlag) {
				gMinfo.snpTagVec.clear();
				copy(chrV.begin(), chrV.end(), std::back_inserter(gMinfo.snpTagVec));
			}
		}
	}

	if (_snpFlag || _chrFlag || _tagFlag) {
		int filterCnt = 0;
		for (int i = 0; i < snpVec.size(); i++) {
			bool flagSNP = true;
			bool flagChr = true;
			bool flagTag = true;
			snpInfo sInfo = snpVec[i];
			if (_snpFlag) flagSNP = _snpSet.count(sInfo.snpID) > 0 ? true : false;
			if (_chrFlag) flagChr = _chrSet.count(sInfo.chr) > 0 ? true : false;
			if (_tagFlag) flagTag = _tagMap.find(sInfo.snpID) != _tagMap.end() ? true : false;
			snpVec[i].keepIt = flagSNP & flagChr & flagTag;
			if (snpVec[i].keepIt) filterCnt++;
		}
		cout << filterCnt << " snps were remained for analysis" << endl;
	}

	if (goptions.GetGenericBlockNum() > 0 || goptions.GetGenericBlockSNPs() > 0 || goptions.GetGenericBlockKB() > 0) {
		set<string> _tagSet;
		vector<string> _tagVec;
		map<string, vector<snpHelp>> _chrSNPmap;
		//putting snps within chromosome
		for (int i = 0; i < snpVec.size(); i++) {
			if (snpVec[i].keepIt) {
				map<string, vector<snpHelp>>::iterator _it = _chrSNPmap.find(snpVec[i].chr);
				if (_it != _chrSNPmap.end()) {
					vector<snpHelp>* _snpHelpVec = &(_it->second);
					snpHelp _snpHelp;
					_snpHelp.index = i;
					_snpHelp.bp = snpVec[i].bp;
					_snpHelpVec->push_back(_snpHelp);
				} else {
					vector<snpHelp> _snpHelpVec;
					snpHelp _snpHelp;
					_snpHelp.index = i;
					_snpHelp.bp = snpVec[i].bp;
					_snpHelpVec.push_back(_snpHelp);
					_chrSNPmap.insert(pair<string, vector<snpHelp>>(snpVec[i].chr, _snpHelpVec));
				}
			}
		}

		map<string, vector<snpHelp>>::iterator _it = _chrSNPmap.begin();
		for (; _it != _chrSNPmap.end(); _it++) {
			vector<snpHelp> _snpHelpVec = _it->second;
			string _chr = _it->first;
			sort(_snpHelpVec.begin(), _snpHelpVec.end(), compareSNP); //sorting snps according to pos

			if (goptions.GetGenericBlockNum() > 0) {
				int _block_num = goptions.GetGenericBlockNum();
				int _interval = ceil(_snpHelpVec.size() / (1.0 * _block_num));
				cout << "Assigning " << _snpHelpVec.size() << " snps on chr " << _chr << " according to block numbers: " << _block_num << ", and " << _interval << " SNPs in each block." << endl;
				for (int k = 0; k < _snpHelpVec.size(); k++) {
					int _intv = ceil((k + 1) / (1.0 * _interval));
					string _tag = _chr + "_" + boost::lexical_cast<string> (_intv);
					snpVec[_snpHelpVec[k].index].snpTag = _tag;
					if (_tagSet.count(_tag) == 0) {
						_tagSet.insert(_tag);
						_tagVec.push_back(_tag);
					}
				}
			} else if (goptions.GetGenericBlockKB() > 0) {
				int _bp_start = _snpHelpVec[0].bp;
				int _bp_end = _snpHelpVec[_snpHelpVec.size() -1].bp;
				int _bp_len = _bp_end - _bp_start;
				float _kb = goptions.GetGenericBlockKB() * 1000.0;
				int _bp_block = ceil(_bp_len / _kb);
				int _cnt_block = 0;
				for (int k = 0; k < _snpHelpVec.size(); k++) {
					int _intv = ceil((_snpHelpVec[k].bp - _bp_start) / _kb);
					if (k == 0) _intv++;
					string _tag = _chr + "_" + boost::lexical_cast<string> (_intv);
					snpVec[_snpHelpVec[k].index].snpTag = _tag;
					if (_tagSet.count(_tag) == 0) {
						_tagSet.insert(_tag);
						_tagVec.push_back(_tag);
						_cnt_block++;
					}
				}
				cout << "Assigning " << _snpHelpVec.size() << " snps on chr " 
				<< _chr <<" [bp:" << _bp_start << "-" << _bp_end << "] according to block kb: " 
				<< goptions.GetGenericBlockKB() << ", and generated " << _cnt_block << " tags" << endl;
			} else if (goptions.GetGenericBlockSNPs() > 0) {
					int _block_snp = goptions.GetGenericBlockSNPs();
					cout << "Assigning " << _snpHelpVec.size() << " snps on chr " << _chr << " according to block snps: " << _block_snp << endl;
					for (int k = 0; k < _snpHelpVec.size(); k++) {
						int _intv = ceil((k + 1) / (1.0 * _block_snp));
						string _tag = _chr + "_" + boost::lexical_cast<string> (_intv);
						snpVec[_snpHelpVec[k].index].snpTag = _tag;
						if (_tagSet.count(_tag) == 0) {
							_tagSet.insert(_tag);
							_tagVec.push_back(_tag);
						}
					}
			}
		}
		gMinfo.snpTagVec.clear();
		copy(_tagVec.begin(), _tagVec.end(), std::back_inserter(gMinfo.snpTagVec));
		cout << _tagVec.size() << " tags are generated\n" << endl;
	}
}

void genotype::sampleFilterQC() {

	cout << "Appling sample filters (if any)..." << endl;
	set<string> _indSet;
	bool _subFlag = false;
	if (!goptions.GetGenericKeepFile().empty()) { // exact snp
		cout << "Keeping samples..." << endl;
		string keepFile = goptions.GetGenericKeepFile();
		ifstream inp(keepFile.c_str());
		if (!inp.is_open()) {
			cerr << "Error reading file "<< keepFile <<endl;
			exit(1);
		} else {
			cout << "Reading " << keepFile <<endl;
		}

		string line;
		string whitespaces (" \t\f\v\n\r");
		while (std::getline(inp, line)) {
			line.erase(0, line.find_first_not_of(" ")); //remove head whitespace;
    		std::size_t found = line.find_last_not_of(whitespaces); //remove tail whitespace;
			if (found !=std::string::npos) {
				line.erase(found+1);
			}

			char c = line[0];
			if (c=='#') {
				continue;
			}
			if (line.empty()) continue;

			istringstream ss(line);
			string word;
			string sid;
			int cnt_len = 0;
			while(ss >> word && cnt_len < 2) {
				if (cnt_len == 0) sid = word;
				if (cnt_len == 1) sid = sid + " " + word;
				cnt_len++;
			}
			_indSet.insert(sid);
		}
		cout << "Read " << _indSet.size () << " unique samples from " << keepFile.c_str() <<endl;
		if (_indSet.size() > 0) {
			for (int i = 0; i < indVec.size(); i++) {
				if (_indSet.count(indVec[i].sid) == 0) indVec[i].keepIt = false;
			}
		}
		inp.close();
	}
}

void genotype::setImputed() {
	isImputed = true;
}

bool genotype::hasImputed() {
	return isImputed;
}

void genotype::QC() {
	cout << "--------------------QC---------------" << endl;
	cout << "1) SNP filter qc" << endl;
	snpFilterQC();
	cout << "2) Sample filter qc" << endl;
	sampleFilterQC();
	cout << "3) Metric qc" << endl;
	MetricQC();
}

void genotype::MetricQC() {

	cout << "Conducting SNP metric QC ..." << endl;
	double _maf = goptions.GetGenericMAF();
	double _lmiss = goptions.GetGenericHWE();
	double _hwe = goptions.GetGenericHWE();
	bool _remove_zero_gstd = goptions.GetGenericRemoveZeroGstd();

	if (goptions.CheckEigenGWASMasterOption() || goptions.CheckPropcMasterOption()) {
		if (_maf < 0.0) { // default maf was -1 in goption
			cout << "Since maf has not been specified, maf threshold is set to 0.01 for EigenGWAS scan" << endl;
 			_maf = 0.01;
		}
		if (!_remove_zero_gstd) {
			cout << "Loci of no variation are removed for EigenGWAS." << endl;
			_remove_zero_gstd = true;
		}
	} else if (goptions.CheckXLDMasterOption()) {
		if (!_remove_zero_gstd) {
			cout << "Loci of no variation are removed for XLD." << endl;
			_remove_zero_gstd = true;
		}
	} else if (goptions.CheckMeMasterOption()) {
		if (!_remove_zero_gstd) {
			cout << "Loci of no variation are removed for ME." << endl;
			_remove_zero_gstd = true;
		}
		if (_maf < 0.0) { // default maf was -1 in goption
			cout << "Since maf has not been specified, maf threshold is set to 0.01 for ME." << endl;
 			_maf = 0.01;
		}
	}

	int _lmafCnt = 0, _lmissCnt = 0, _zeroCnt = 0;
	for (int i = 0; i < snpVec.size(); i++) {
		bool flag = snpVec[i].keepIt;
		if (_maf >= 0.0) {
			double lmaf = snpVec[i].Freq < 0.5 ? snpVec[i].Freq : (1.0 - snpVec[i].Freq);
			if (lmaf < _maf) {
				snpVec[i].mafFlag = false;
				flag = false;
				_lmafCnt++;
			}
		}
		if (_lmiss >= 0) {
			double lmiss = 1 - (1.0 * snpVec[i].Not_NA_cnt) / (1.0 * (snpVec[i].locusCnt[0] + snpVec[i].locusCnt[1] + snpVec[i].locusCnt[2] + snpVec[i].locusCnt[3]));
			if (lmiss < _lmiss) {
				snpVec[i].lmissFlag = false;
				flag = false;
				_lmissCnt++;
			}
		}
		if (_remove_zero_gstd) {
			if (snpVec[i].std == 0) {
				snpVec[i].mafFlag = false;
				flag = false;
				_zeroCnt++;
			}
		}
		snpVec[i].keepIt = flag;
	}
	cout <<"Removed SNP numbers: " << _lmafCnt << " (MAF qc), " << _lmissCnt << " (locus missing qc), " << _zeroCnt << " (zero s.d. qc)" <<endl;
}

void genotype::writeBed() {

	string gfamname = goptions.GetGenericOutFile() + ".fam";
	fstream famout(gfamname.c_str(), ios::out);

	string gbimname = goptions.GetGenericOutFile() + ".bim";
	fstream bimout(gbimname.c_str(), ios::out);

	string gbedname = goptions.GetGenericOutFile() + ".bed";
	fstream bedout(gbedname.c_str(), ios::out|ios::binary);
	bedout.write((char*) PLINK_magic, 3);

	cout << "Writing " << act_snp.size() << " snps, " << act_ind.size() << " samples." << endl;

	if (act_ind.size() == gMinfo.fam_N) {
		for (int i = 0; i < indVec.size(); i++) {
			indxInfo iInfo = indVec[i];
			famout << iInfo.fid << "  " << iInfo.iid << "  " << iInfo.ppid << "  " << iInfo.mmid 
				<< "  " << iInfo.sex << "  " << iInfo.qval << endl;
		}
		if (act_snp.size() == gMinfo.bim_N) {
			for (int i = 0; i < snpVec.size(); i++) {
				snpInfo sInfo = snpVec[i];
				bimout << sInfo.chr << "  " << sInfo.snpID << "  " << sInfo.pos << "  " << sInfo.bp 
					<< "  " << sInfo.a1 << "  " << sInfo.a2 << endl; 
			}
			int _ByteCol = ceil((gMinfo.fam_N * 1.0) / unitPerGenos);
			for (int i = 0; i < gMinfo.bim_N; i++) {
				for (int j = 0; j < _ByteCol; j++) {
					bedout.put(gBedBlock[i][j]); //can be faster?
				}
			}
		} else {
			for (int i = 0; i < act_snp.size(); i++) {
				int _i = act_snp[i];
				snpInfo sInfo = snpVec[_i];
				bimout << sInfo.chr << "  " << sInfo.snpID << "  " << sInfo.pos << "  " << sInfo.bp 
				<< "  " << sInfo.a1 << "  " << sInfo.a2 << endl;
			}
			int _ByteCol = ceil((gMinfo.fam_N * 1.0) / unitPerGenos);
			for (int i = 0; i < act_snp.size(); i++) {
				int _i = act_snp[i];
				for (int j = 0; j < _ByteCol; j++) {
					bedout.put(gBedBlock[_i][j]); //can be faster?
				}
			}
		}
	} else {
		for (int i = 0; i < act_ind.size(); i++) {
			int _idx = act_ind[i];
			indxInfo iInfo = indVec[_idx];
			famout << iInfo.fid << "  " << iInfo.iid << "  " << iInfo.ppid << "  " << iInfo.mmid 
				<< "  " << iInfo.sex << "  " << iInfo.qval << endl;
		}
		for (int i = 0; i < act_snp.size(); i++) {
			int _i = act_snp[i];
			snpInfo sInfo = snpVec[_i];
			bimout << sInfo.chr << "  " << sInfo.snpID << "  " << sInfo.pos << "  " << sInfo.bp 
			<< "  " << sInfo.a1 << "  " << sInfo.a2 << endl;
		}

		for (int i = 0; i < act_snp.size(); i++) {
			int _i = act_snp[i];
			unsigned char _byte = 0;
			int _shift = 0;
			for (int j = 0; j < act_ind.size(); j++) {
				int _j = act_ind[j];
				unsigned char _v = 0;
				if (is_geno_NA(_i, _j)) {
					_v = 1;
				} else {
					_v = get_geno_add(_i, _j);
					_v = _v == 0 ? _v : (_v + 1);
				}
				_byte += _v << (_shift * 2);
				_shift++;

				if (_shift == 4) {
					bedout.put(_byte);
					_shift = 0;
					_byte = 0;
				}
			}
			if (_shift != 0) bedout.put(_byte);
		}
	}
	famout.close();
	bimout.close();
	bedout.close();
}

void genotype::writeFreq() {

	cout << "Writing freqiencies for " << act_snp.size() << " SNPs and " << act_ind.size() << " samples ..." << endl;
	string frqname = goptions.GetGenericOutFile() + ".frq";
	fstream frqout(frqname.c_str(), ios::out);
	frqout << "chr\tsnp\tA\ta\tfreq(A)\tfreq(A-imputation)\taa\tAa\tAA\tna\tsd" << endl;

	if (gMinfo.fam_N == act_ind.size()) {
		for (int i = 0; i < act_snp.size(); i++) {
			int _i = act_snp[i];
			snpInfo sInfo = snpVec[_i];
			frqout << sInfo.chr << "\t" << sInfo.snpID << "\t"
				<< sInfo.a1 << "\t" << sInfo.a2 << "\t"
				<< snpVec[_i].Freq << "\t"
				<< snpVec[_i].FreqImputed << "\t"
				<< snpVec[_i].locusCnt[0] << "\t"
				<< snpVec[_i].locusCnt[1] << "\t"
				<< snpVec[_i].locusCnt[2] << "\t"
				<< snpVec[_i].locusCnt[3] << "\t"
				<< snpVec[_i].std << endl;
		}
	} else {
		for (int i = 0; i < act_snp.size(); i++) {
			int _i = act_snp[i];
			snpInfo sInfo = snpVec[_i];
			quickUnit qu = cal_geno(_i, act_ind);
			frqout << sInfo.chr << "\t" << sInfo.snpID << "\t"
				<< sInfo.a1 << "\t" << sInfo.a2 << "\t"
				<< qu.frq << "\t"
				<< qu.frq << "\t"
				<< qu.geno_cnt[0] << "\t"
				<< qu.geno_cnt[1] << "\t"
				<< qu.geno_cnt[2] << "\t"
				<< qu.geno_cnt[3] << "\t"
				<< qu.std << endl;
		}
	}
	frqout.close();
	cout << "Writing frequencies into " << frqname.c_str() << endl;
}

int genotype::get_sample_size() {
	return indVec.size();
}

int genotype::get_active_sample_size() {
	return act_ind.size();
}

int genotype::get_snp_number() {
	return snpVec.size();
}

int genotype::get_active_snp_number() {
	return act_snp.size();
}

double genotype::cal_freq(int snpIdx, const vector<int> indIdx) {
	int misscnt = 0;
	double cnt = 0;
	for (int i = 0; i < indIdx.size(); i++) {
		int g = get_geno_add(snpIdx, indIdx[i]);
//		msb[snpIdx][indIdx[i]] * 2 + lsb[snpIdx][indIdx[i]];
		if (g != 1) {
			cnt += 0.5 * g;
		} else {
			misscnt++;
		}
	}

	if (indIdx.size() == misscnt) {
		return 0.0;
	} else {
		return cnt/(1.0 * indIdx.size() - misscnt);
	}
}

quickUnit genotype::cal_geno(int snpIdx, vector<int> indIdx) {
	vector<int> _gCnt(4);
	vector<int> _naVec;
	for (int i = 0; i < indIdx.size(); i++) {
		int _ind = indIdx[i];
		if (is_geno_NA(snpIdx, _ind)) {
			_gCnt[3]++;
			_naVec.push_back(_ind);
		} else {
			_gCnt[get_geno_add(snpIdx, _ind)]++;
		}
	}
	quickUnit qU(_gCnt, _naVec);
	return qU;
}

float genotype::cal_observed_freq(const unsigned char* line) {
	int y[4];
	int observed_sum = 0;
	int observed_ct = 0;
	int gByteCol = ceil((gMinfo.fam_N * 1.0) / unitPerGenos);

	for (int k = 0; k < gByteCol; k++) {
		unsigned char c = line[k];
		y[0] = (c) & mask;
		y[1] = (c>>2) & mask;
		y[2] = (c>>4) & mask;
		y[3] = (c>>6) & mask;
		int j0 = k * unitPerGenos;
		int lmax = 4;

		if (k == gByteCol - 1) {
			lmax = gMinfo.fam_N % 4;
			lmax = (lmax==0) ? 4 : lmax;
		}

		for (int l = 0; l < lmax; l++) {
			int j = j0 + l;
			// int ver_seg_no = j/segment_size_ver;
			// Extract PLINK coded genotype and convert into 0/1/2
			// PLINK coding: 00->0, 01->missing, 10->1, 11->2
			int val = y[l];
			if (val != 1) {
				val = (val == 0) ? val : (--val);
				observed_sum += val;
				observed_ct++;
			}
		}
	}
	return observed_sum * 0.5 / observed_ct;
}

void genotype::bed_flip_preset() {
	bed_flip = new unsigned char [256]();
	for (int i = 0; i < 256; i++) {
//		bed_flip[i] = 0;
		int _v[4];
		for (int j = 0; j < 4; j++) {
			_v[j] = (i % (int) pow(4, j + 1))/ pow(4, j);
			switch (_v[j]) {
				case 0:
					bed_flip[i] += 3 * (int) pow(4, j);
					break;
				case 3:
					break;
				default:
					bed_flip[i] += _v[j] * (int) pow(4, j);
					break;
			}
		}
	}
}

void genotype::bed_preset() {

//preset
	nbyte = 2; //if necessary,change the length of index and nbyte

	int _len = pow(256, nbyte);
	naVec.resize(_len); //missing genotype location 
	geno_cnt = new int *[_len];
	for (int i = 0; i < _len; i++) {
		geno_cnt[i] = new int[4]();
//		for (int j = 0; j < 4; j++) {
//			geno_cnt[i][j] = 0;
//		}
	}

	//PLINK bed: two bytes as a unit, from right to left,so the naIdx maybe like [7,6,5,4,3,2,1,0] for two bytes
	//or [15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0] for four bytes.
	int naIdx[8 / gBitSize * nbyte];
	int _b = (8 / gBitSize) * nbyte; 
	int idx = 0;
	for (int i = 0; i < _b; i++) naIdx[i] = _b - i - 1;

	int pw = 8 / gBitSize * nbyte - 1; // pw = 7 if we use int;

	for (int i = 0; i < _len; i++) {
		int remain = i;
		vector<int> _naV;
		int P = (int) pow(4, pw);
		for (int j = 0; j < 4 * nbyte; j++) {
			int _v = remain / P;
			switch (_v) {
				case 0:
					geno_cnt[i][0]++;
					break;
				case 1:
					geno_cnt[i][3]++;
					_naV.push_back(naIdx[j]);
					break;
				case 2:
					geno_cnt[i][1]++;
					break;
				default:
					geno_cnt[i][2]++;
			}
			remain = remain % P;
			P /= 4;
		}
		copy(_naV.begin(), _naV.end(), std::back_inserter(naVec[i]));
	}
	if (goptions.IsGenericDebug()) {
		for (int i = 0; i < _len; i++) {
			for (int j = 0; j < naVec[i].size(); j++) {
				cout << naVec[i][j] << " ";
			}
			cout << endl;
		}
	}
}

quickUnit genotype::cal_observed_freq_pro(const unsigned char* gBlock) {
	vector<int> _gCnt(4);
	vector<int> _naVec;
	int gConst = unitPerGenos * nbyte;
	int gByteCol = ceil((gMinfo.fam_N * 1.0) / (gConst));

	for (int k = 0; k < gByteCol; k++) {
		unsigned long int c = 0;
		for (int j = 0; j < nbyte; j++)	c += (gBlock[k * nbyte + j] << (gConst * j));

		if (k == gByteCol - 1) {
			c = 0;
			int _nbyte = (gMinfo.fam_N % (gConst) == 0) ? nbyte : ceil((gMinfo.fam_N % (gConst) * 1.0) / unitPerGenos);
			for (int _j = 0; _j < _nbyte; _j++) c += (gBlock[k * nbyte + _j] << (gConst * _j));
			_gCnt[0] -= gByteCol * (gConst) - gMinfo.fam_N;
		}

		for (int l = 0; l < _gCnt.size(); l++) _gCnt[l] += geno_cnt[c][l];
		if (geno_cnt[c][3] > 0) {
			for (int j = 0; j < naVec[c].size(); j++) {
				_naVec.push_back(naVec[c][j] + k * (gConst));
			}
		}
	}
	quickUnit qU(_gCnt, _naVec);
	return qU;
}

void genotype::quickFreq(string bedfilename) {

//old method
	clock_t freq_begin0 = clock();
 	cout << "Reading " << bedfilename << " for " << gMinfo.bim_N << " SNPs of " << gMinfo.fam_N << " samples" << endl;

 	ifstream ifs(bedfilename.c_str(), ios::in|ios::binary);
 	int gByteCol = ceil((gMinfo.fam_N * 1.0) / unitPerGenos);

 	binary_read(ifs, PLINK_magic);////???
	int _ByteCol = ceil((gMinfo.fam_N * 1.0) / unitPerGenos);
 	unsigned char* gBlock = new unsigned char[_ByteCol]();
 	for (int i = 0; i < gMinfo.bim_N; i++) {
 		ifs.read(reinterpret_cast<char*> (gBlock), _ByteCol * sizeof(unsigned char));
 		float p_j = cal_observed_freq(gBlock);
 	}
 	clock_t freq_end0 = clock();
 	double freq_time0 = double(freq_end0 - freq_begin0) / CLOCKS_PER_SEC;
 	cout << "Freq method 1 " << freq_time0 << "s." <<endl;
 	ifs.close();

 //new method
	clock_t freq_begin1 = clock();
 	cout << "Reading " << bedfilename << " for " << gMinfo.bim_N << " SNPs of " << gMinfo.fam_N << " samples" << endl;

  	ifstream ifs1(bedfilename.c_str(), ios::in|ios::binary);

  	binary_read(ifs1, PLINK_magic);
 	unsigned char* gBlock1 = new unsigned char[_ByteCol]();
  	for (int i = 0; i < gMinfo.bim_N; i++) {
  		ifs1.read(reinterpret_cast<char*> (gBlock1), _ByteCol * sizeof(unsigned char));
  		quickUnit new_p_j = cal_observed_freq_pro(gBlock1);
  	}
  	clock_t freq_end1 = clock();
  	double freq_time1 = double(freq_end1 - freq_begin1) / CLOCKS_PER_SEC;
  	cout << "Freq method 2 " << freq_time1 << "s." <<endl;
  	ifs1.close();
}

void genotype::post_update_freq() {
	for (int i = 0; i < act_snp.size(); i++) {
		int _i = act_snp[i];
		quickUnit qu = cal_geno(_i, act_ind);

		snpVec[_i].Not_NA_cnt = qu.Not_NA_cnt;
		snpVec[_i].Freq = qu.frq;
		snpVec[_i].FreqImputed = qu.frq;
		snpVec[_i].std = qu.std;

		if (goptions.GetGenericAdjVar() == 0) {//default observed variance
			snpVec[_i].adj_std = snpVec[_i].std;
		} else {
			if (!goptions.IsGenericInbred()) { // outbred expected var
				snpVec[_i].adj_std = sqrt(snpVec[_i].FreqImputed * (1 - snpVec[_i].FreqImputed) * 2);
			} else { //inbred expected var
				snpVec[i].adj_std = sqrt(snpVec[_i].FreqImputed * (1 - snpVec[_i].FreqImputed) * 4);
			}
		}

		snpVec[_i].locusCnt.clear();
		copy(qu.geno_cnt.begin(), qu.geno_cnt.end(), std::back_inserter(snpVec[_i].locusCnt));
		snpVec[_i].na_snp.clear();
		copy(qu.naVec.begin(), qu.naVec.end(), std::back_inserter(snpVec[_i].na_snp));
	}
}

void genotype::recode_ref_allele() {
	cout << "Updating SNP reference alleles ..." <<endl;
	map<string, string> _saMap;
	vector<int> _recodeIndexVec;

	if (!goptions.GetGenericForceRefAlleleFile().empty()) {
		string saFile = goptions.GetGenericForceRefAlleleFile();
		ifstream inp(saFile);
		if (!inp.is_open()) {
			cerr << "Error reading file "<< saFile << endl;
			exit(1);
		} else {
			cout << "Reading SNP allele pairs from " << saFile << endl;
		}

		string line;
		string whitespaces (" \t\f\v\n\r");
		int LineCnt = 0;
		while (std::getline(inp, line)) {
			char c = line[0];
			if (c == '#') {
				continue;
			}
			if (line.empty()) continue;

			line.erase(0, line.find_first_not_of(" ")); //remove head whitespace;
			std::size_t found = line.find_last_not_of(whitespaces); //remove tail whitespace;
			if (found != std::string::npos) {
				line.erase(found + 1);
			}
			LineCnt++;

			istringstream ss(line);
			string _s;
			string _a;
			string word;
			int cnt_word = 0;
			while (ss >> word && cnt_word < 2) {
				switch (cnt_word) {
					case 0:
						_s = word;
						break;
					default:
						_a = word;
						break;
				}
				cnt_word++;
			}
			if (cnt_word != 2) {
				cerr << "Line " << LineCnt << " invalid format" << endl;
				cerr << "Please make sure two columns in each line for --recode-ref-allele" << endl;
				exit(0);
			}
			auto it = _saMap.find(_s);
			if (it != _saMap.end()) {
				cerr << "Line: " << LineCnt << " snp " << _s << " is associated with ref allele " << _a << " and " << it->second << endl;
				cerr << "SNP " << _s << " has been mapped to more than one reference alleles" << endl;
				exit(0);
			}
			_saMap.insert(pair<string, string>(_s, _a));
		}

		if (LineCnt != _saMap.size()) {
			cerr << "Read " << LineCnt << " lines but mapped to " << _saMap.size() << " unique snp-allele pairs" << endl;
			cerr << "Make sure each SNP is assigned to one allele only." << endl;
			exit(0);
		}

		cout << "Read " << _saMap.size () << " SNP-allele pairs from " << saFile.c_str() << endl;
		for (int i = 0; i < snpVec.size(); i++) {
			auto it = _saMap.find(snpVec[i].snpID);
			if (it != _saMap.end()) {
				if (snpVec[i].a1 == it->second) {
					continue;
				} else if (snpVec[i].a2 == it->second) {
					string __a = snpVec[i].a1;
					snpVec[i].a1 = snpVec[i].a2;
					snpVec[i].a2 = __a;
					_recodeIndexVec.push_back(i);
				} else {
					cerr << "Matched neither alleles for SNP " << snpVec[i].snpID << endl;
					cerr << "Make sure the reference allele in " << saFile.c_str() << endl;
					exit(0);
				}
			}
		}
		cout << "Matched " << _recodeIndexVec.size() << " SNP for flipping ref alleles" << endl;
		inp.close();
	}

	if (_recodeIndexVec.size() != 0) {
		int gByteCol = ceil((gMinfo.fam_N * 1.0) / unitPerGenos);
		int gShift = 4 - gMinfo.fam_N % unitPerGenos;
		for (int i = 0; i < _recodeIndexVec.size(); i++) {
			int snpIdx = _recodeIndexVec[i];
			for (int j = 0; j < gByteCol; j++) {
				if (j != (gByteCol - 1)) {
					gBedBlock[snpIdx][j] = (unsigned char) bed_flip[gBedBlock[snpIdx][j]];
				} else { //last byte
					unsigned char m = 255 >> (gShift >> 1);
					unsigned char idx = bed_flip[gBedBlock[snpIdx][j]];
					idx &= m;
					gBedBlock[snpIdx][j] = idx;
				}
			}
		}
	}
}

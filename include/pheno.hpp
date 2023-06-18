#ifndef PHENO_HPP_
#define PHENO_HPP_

#include <bits/stdc++.h>
#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <map>
#include <algorithm>
#include <vector>
#include "Goptions.hpp"

using namespace Eigen;
using namespace std;
using namespace boost;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;

extern Goptions goptions;

struct subQInfo {
	string fid;
	string iid;
	string sid;
	int line_index;
	int naCnt;
	vector<double> qval;
	vector<bool> isNA;
};

class pheno {

public:
	map<string, int> sub_index;

	int Ntraits;
	int Nindv;
	vector<double> tSum;
	vector<double> tSqSum;
	vector<int> tObs;
	vector<double> tMean;
	vector<double> tSd;

	vector<subQInfo> subQInfoVec;
    MatrixXdr qvec;
	string pfile;
	vector<int> q_index;

    pheno() {

    }

	void read_q_file(string filename, vector<int> q_num) {

		pfile = filename;

		ifstream inp(filename.c_str());
		if (!inp.is_open()) {
			cerr << "Error reading file "<< filename << endl;
			exit(1);
		} else {
			cout << "Reading " << filename << endl;
		}

		string line;
		string whitespaces(" \t\f\v\n\r");

		int line_cnt = 0;
		int eff_line_cnt = 0;
		int cnt_anno = 0;

		map<string, int> _subMap;
		while (std::getline(inp, line)) {
			line_cnt++;
			line.erase(0, line.find_first_not_of(" ")); //remove head whitespace;
			std::size_t found = line.find_last_not_of(whitespaces); //remove tail whitespace;

			if (found != std::string::npos) {
				line.erase(found+1);
			}

			char c = line[0];
			if (c == '#') {
				cnt_anno++;
				continue;
			}

			if (line.empty()) continue;

			eff_line_cnt++;
			istringstream ss(line);

			string word;

			string sid;
			string fid;
			string iid;
			vector<int> na_t;
			vector<double> valVec;

			subQInfo _sQInfo;
			int naCnt = 0;
			int word_cnt = 0;
			while (ss >> word) {
				if (word_cnt == 0) {
					sid = word;
					fid = word;
				} else if (word_cnt == 1) {
					sid = sid + " " + word;
					iid = word;
				} else {
					if (word == "na" || word == "NA" || word =="nan" 
					    || word == "-9" || word == "-Inf" || word == "Inf" 
						|| word == "." || word == "-" || word == "?") {
						valVec.push_back(0.0);
						na_t.push_back(true);
						naCnt++;
					} else {
						try {
							valVec.push_back((boost::lexical_cast<double>(word)));
							na_t.push_back(false);
						} catch (const boost::bad_lexical_cast &e) {
							cerr << e.what() << '\n';
							exit(0);
						}
					}
				}
				word_cnt++;
			}

			if (eff_line_cnt == 1) {
				Ntraits = valVec.size();
				if (Ntraits == 0) {
					cerr << "No values observed in the first line." << endl;
					exit(1);
				}
			} else {
				if (valVec.size() != Ntraits) {
					cerr << "line " << line_cnt << " does not have " << Ntraits << " elements." << endl;
				}
			}

			_sQInfo.fid = fid;
			_sQInfo.iid = iid;
			_sQInfo.sid = fid + " " + iid;
			_sQInfo.line_index = eff_line_cnt;
			_sQInfo.naCnt = naCnt;
			_sQInfo.qval.resize(valVec.size());
			copy(valVec.begin(), valVec.end(), _sQInfo.qval.begin());
			_sQInfo.isNA.resize(na_t.size());
			copy(na_t.begin(), na_t.end(), _sQInfo.isNA.begin());
			subQInfoVec.push_back(_sQInfo);

			int cnt = _subMap.count(sid);
			if (cnt == 0) {
				_subMap.insert(pair<string, int>(sid, 1));
			} else {
				map<string, int>::iterator iter = _subMap.find(sid);
				int val = (int) iter->second;
				val++;
				_subMap.erase(sid);
				_subMap.insert(pair<string, int>(sid, val));
				cerr << "line " << eff_line_cnt << " duplicated id: " << sid << endl;
			}
		}

		bool flag = true;
		for (map<string, int>::iterator it = _subMap.begin(); it != _subMap.end(); it++) {
			if (it->second > 1) {
				cout << it->first << " duplicated " << it->second << " times." << endl;
				flag = false;
			}
		}

		if (!flag) {
			exit(1);
		}

		if (goptions.IsGenericDebug()) {
			for (int i = 0; i < subQInfoVec.size(); i++) {
				cout << subQInfoVec[i].fid << " " << subQInfoVec[i].iid << endl;
				for (int j = 0; j < subQInfoVec[i].qval.size(); j++) {
					cout << subQInfoVec[i].qval[j] << " ";
				}
				cout << endl;
			}

			for (int i = 0; i < subQInfoVec.size(); i++) {
				subQInfo sq = subQInfoVec[i];
				if (sq.naCnt == 0) {
					cout << " no missing in " << i <<endl;
				} else {
					cout << " missed " << sq.naCnt << " values for individual " << sq.sid << endl;
				}
			}
		}

		if (q_num.empty()) {
			q_index.resize(Ntraits);
			for (int i = 0; i < q_index.size(); i++) {
				q_index[i] = i;
			}
		} else {
			q_index.assign(q_num.begin(), q_num.end());
			for (int i = 0; i < q_index.size(); i++) {
				q_index[i]--;
			}
		}

		for (int i = 0; i < subQInfoVec.size(); i++) {
			sub_index.insert(pair<string, int>(subQInfoVec[i].sid, i));
		}
/*
		if (goptions.IsGenericQImput()) {
			Qimputation();
		}
*/
	}

	MatrixXdr getQpheno(vector<string> subVec, int q_i) {
		MatrixXdr qvec(subVec.size(), 1);
		for (int i = 0; i < subVec.size(); i++) {
			map<string, int>::iterator it = sub_index.find(subVec[i]);
//			cout << (it->first) << " " << (int) (it->second) << endl;
			qvec(i, 0) = subQInfoVec[(int) it->second].qval[q_i];
		}
		return qvec;
	}

	MatrixXdr getQpheno(vector<string> subVec, vector<int> q_vec) {
		MatrixXdr qvec(subVec.size(), q_vec.size());
		for (int i = 0; i < subVec.size(); i++) {
			map<string, int>::iterator it = sub_index.find(subVec[i]);
//			cout << (it->first) << " " << (int) (it->second) << endl;
			for (int j = 0; j < q_vec.size(); j++) {
				qvec(i, j) = subQInfoVec[(int) it->second].qval[q_vec[j]];
			}
		}
		return qvec;
	}


	void get_NotMissSub(vector<string>& qSubVec) {
		for (int i = 0; i < subQInfoVec.size(); i++) {
			subQInfo sq = subQInfoVec[i];
			bool flag = true;
			for (int j = 0; j < q_index.size(); j++) {
				if (sq.isNA[q_index[j]]) {
					flag = false;
					break;
				}
			}
			if (flag) qSubVec.push_back(sq.sid);
		}
	}

	void get_NotMissSub(vector<string>& qSubVec, int q_i) {
		for (int i = 0; i < subQInfoVec.size(); i++) {
			subQInfo sq = subQInfoVec[i];
			bool flag = true;
			if (sq.isNA[q_index[q_i]]) {
				flag = false;
				break;
			}
			if (flag) qSubVec.push_back(sq.sid);
		}
	}

	vector<int> get_q_index() {
		return q_index;
	}

	void Qimputation() {

		tSum.resize(Ntraits);
		tSqSum.resize(Ntraits);
		tObs.resize(Ntraits);
		tMean.resize(Ntraits);
		tSd.resize(Ntraits);

		for (int i = 0; i < subQInfoVec.size(); i++) {
			subQInfo sq = subQInfoVec[i];
			for (int j = 0; j < sq.isNA.size(); j++) {
				if (!subQInfoVec[i].isNA[j]) {
					tSum[j] += sq.qval[j];
					tSqSum[j] += sq.qval[j] * sq.qval[j];
					tObs[j]++;
				}
			}
		}

		for (int i = 0; i < tMean.size(); i++) {
			tMean[i] = tSum[i] / (1.0 * tObs[i]);
			tSd[i] = sqrt((tSqSum[i] - (tMean[i] * tMean[i]) * (1.0 * tObs[i]))/(1.0 * tObs[i] - 1));
		}

		srand(goptions.GetGenericSeed());
		std::default_random_engine generator(goptions.GetGenericSeed());

		for (int i = 0; i < subQInfoVec.size(); i++) {
			for (int j = 0; j < subQInfoVec[i].isNA.size(); j++) {
				if (subQInfoVec[i].isNA[j]) {
					std::normal_distribution<double> norm_dist(tMean[j], tSd[j]);
					subQInfoVec[i].qval[j] = norm_dist(generator);
				}
			}
		}
	}

};

#endif

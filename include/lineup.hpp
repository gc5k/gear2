#ifndef LINEUP_HPP_
#define LINEUP_HPP_

#include <vector>
#include <algorithm>
#include <iterator>
#include <genotype.h>

using namespace std;

class lineup {
public:
    vector<string> subSet;

    lineup (vector<string> vG) {
        subSet.resize(vG.size());
        copy(vG.begin(), vG.end(), subSet.begin());
        sort(subSet.begin(), subSet.end());
    }

    void lineup2 (vector<string> v3) {
        sort(v3.begin(), v3.end());
        vector<string> v;
        v.resize(min(subSet.size(), v3.size()));
        vector<string>::iterator v_end = set_intersection(subSet.begin(), subSet.end(), v3.begin(), v3.end(), v.begin());
        subSet.clear();
        for (vector<string>::iterator it = v.begin(); it != v_end; it++) {
            subSet.push_back(*it);
        }
        sort(subSet.begin(), subSet.end());
        cout << "line up " << subSet.size() << " samples. " << endl;
    }

    vector<string> get_subinterset() {
        return subSet;
    }
};

#endif

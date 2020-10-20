#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <iostream>
#include "omp.h"
using namespace std;

int main(int argc, char const *argv[])
{
    vector<int> values {};
    int max_iter {1};

    for (int i=0;i<3;i++) 
    {
        values.push_back(i);
        max_iter *= i+1;
    }

    for (auto n: values) {
            cout << n << " ";
    }
    cout << endl;

    omp_set_num_threads(5);

    #pragma omp parallel for
    for (int iter=1; iter<max_iter; iter++) {
        next_permutation(values.begin(), values.end());
        for (auto n: values) {
            cout << n << " ";
        }

        cout << endl;

        int code = iter;
		// the index correspond to the queen's number and the queen's column
		// we only generate configurations where there's only one queen per column
		for (int i = 0; i < 5; i++)
		{
			cout << code % 5 << " ";
			code /= 5;
		}
    }

    return 0;
}

#include<omp.h>
#include<iostream>
#include<array>

using std::cin;
using std::cout;
using std::endl;
using std::array;

int main() {
    constexpr int PROBLEMSIZE {12};
    constexpr int MAX_THREAD {6};
    
    array<int, PROBLEMSIZE> a {0};
    array<int, PROBLEMSIZE> t {0};

    // for (int &i: a) {
    //     cout << i << endl;
    // }

    // cout << "----------------------------" << endl;

    #pragma omp parallel for num_threads(6)
    for (int i = 0; i < PROBLEMSIZE; i++) {
        a[i] = i;
        // t[i] = omp_get_thread_num();
        cout << " " << a[i] << " ";
    }

    cout << "\n----------------------------" << endl;

    for (int i = 0; i < PROBLEMSIZE; i++) {
        cout << "a[" << i << "] = " << a[i] << " done by thread " << t[i] << endl;
    }

    return 0;
}

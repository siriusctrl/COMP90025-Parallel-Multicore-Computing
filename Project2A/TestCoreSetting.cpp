#include <omp.h>
#include <iostream>
#include <memory>
#include <vector>

using namespace std;

// typedef struct seqalign_parallel {
//     int number = 0;
//     int padding[15];
// } content;

class content {
public:
    int number = 0;
content(int i) {number = i;}
};

int main(int argc, char const *argv[]) {
    auto dp = make_unique<vector<unique_ptr<vector<unique_ptr<content>>>>>();
    auto b = make_unique<vector<unique_ptr<content>>>();

    
    for (size_t i = 0; i < 10; i++) {
        auto a = make_unique<content>(i);
        // auto b = make_unique<vector<unique_ptr<content>>>();
        // b->push_back(a);
        // dp->push_back(b);
        // cout << b->size() << endl;
        cout << a->number << endl;
        b->push_back(a);
    }

    // for (auto const &a:*dp) {
    //     cout << a->size() << endl;
    // }
    
    

    return 0;
}

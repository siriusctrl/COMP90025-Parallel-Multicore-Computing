#include<omp.h>
#include<iostream>
#include<array>
#include<sstream>
#include<string>
#include "openmp.h"

using std::cin;
using std::cout;
using std::endl;
using std::array;
using std::string;

int main(int argc, char* argv[]) {

    int n_thread {0};

    std::stringstream input {argv[1]};
    input >> n_thread;

    cout << "doing job with " << n_thread << " threads" << endl;

    cout << "pi = " << cal_pi_reduction(n_thread) << endl;

    return 0;
}


double cal_pi() {
    constexpr int N_THREAD {6};

    array<double, N_THREAD> res {0.0};

    const int NUM_STEPS {100'000};
    double step {1.0/(double) NUM_STEPS};

    double x {0};

    #pragma omp parallel for num_threads(N_THREAD)
    for (size_t i = 0; i < NUM_STEPS; i++) {
        x = (i+0.5)*step;
        res[omp_get_thread_num()] += 4.0/(1.0+x*x);
    }

    double sum {0.0};

    for (auto &i: res) {
        sum += i;
    }
    
    return step * sum;
}

double cal_pi_non_for() {
    constexpr int N_THREAD {6};
    array<double, N_THREAD> res {0.0};

    const int NUM_STEPS {100'000'000};
    double step {1.0/(double) NUM_STEPS};

    int nthreads;

    omp_set_num_threads(N_THREAD);

    #pragma omp parallel
    {
        int id {omp_get_thread_num()};
        double x {0.0};
        
        for (size_t i = id; i < NUM_STEPS; i+=omp_get_num_threads()) {
            x = (i+0.5)*step;
            res[id] += 4.0/(1.0+x*x);
        }
        
    }

    double sum {0.0};

    for (auto &i: res) {
        sum += i;
    }
    
    return step * sum;

}

/**
 * synchronization
 */
double cal_pi_atom(int n_thread) {

    const int NUM_STEPS {100'000'000};
    double step {1.0/(double) NUM_STEPS};

    double x {0};
    double pi {0.0};

    omp_set_num_threads(n_thread);

    #pragma omp parallel
    {
        int id {omp_get_thread_num()};
        double x {0.0}, sum {0.0};
        
        for (size_t i = id; i < NUM_STEPS; i+=omp_get_num_threads()) {
            x = (i+0.5)*step;
            sum += 4.0/(1.0+x*x);
        }

        // it is critical important to only sum the value here instead of in the for loop
        #pragma omp atomic
        pi += sum*step;
    }
    
    return pi;
}

/**
 * loop construct
 */ 
double cal_pi_reduction(int n_thread) {

    const int NUM_STEPS {100'000'000};
    double step {1.0/(double) NUM_STEPS};

    double x {0}, sum {0};

    omp_set_num_threads(n_thread);

    // private x is a must
    #pragma omp parallel for private(x) reduction(+:sum)
    for (size_t i = 0; i < NUM_STEPS; i++) {
        x = (i+0.5)*step;
        sum += 4.0/(1.0+x*x);
    }
    
    return step * sum;
}

/**
 * random number (Monte Carlo Calculations)
 */
double cal_pi_monte_carlo(int n_thread) {
    
}

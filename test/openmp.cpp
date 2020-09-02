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

    if (n_thread == -1) {
        // disable dynamic adjustment of # of threads
        omp_set_dynamic(0);
        n_thread = omp_get_max_threads();
        omp_set_num_threads(n_thread);
    } else {
        omp_set_num_threads(n_thread); 
    }

    cout << "doing job with " << n_thread << " threads" << endl;

    cout << "pi = " << cal_pi_reduction() << endl;

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
double cal_pi_atom() {

    const int NUM_STEPS {100'000'000};
    double step {1.0/(double) NUM_STEPS};

    double x {0};
    double pi {0.0};

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
double cal_pi_reduction() {

    const int NUM_STEPS {100'000'000};
    double step {1.0/(double) NUM_STEPS};

    double x {0}, sum {0};

    // private x is a must, can use firstprivate(x) to init the value, but its not
    // necessary here, since we do not do things like +=
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
double cal_pi_monte_carlo() {
    // constexpr long num_trails {100'00};
    // long N_circle {0};
    // double pi, x, y;
    // // radius of circle
    // double r {1.0};
    return 0.0;
}

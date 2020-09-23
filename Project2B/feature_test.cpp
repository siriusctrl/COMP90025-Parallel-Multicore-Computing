#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include <omp.h>

using namespace std;

int n_threads = 5;

void gg() {
    #pragma omp parallel for 
    for (int i=0; i<10; i++) {
        printf("I am %d i= %d\n", omp_get_thread_num(), i);
    }    
}


int main(int argc, char **argv) {

    // omp_set_nested(1);   /* make sure nested parallism is on */
    // int nprocs = 4;
    // int group1 = 1;
    // int group2 = nprocs - group1;

    // cout << nprocs << " " << group1 << " " << group2 << endl;

    // #pragma omp parallel default(none) shared(group1, group2) num_threads(2)
    // #pragma omp single
    // {
    //     #pragma omp task
    //     #pragma omp parallel num_threads(group1)
    //     {
    //         printf("I am %d from group %d\n",omp_get_thread_num(), omp_get_ancestor_thread_num(1));
    //     }

    //     #pragma omp task
    //     #pragma omp parallel num_threads(group2)
    //     {
    //         printf("I am %d from group %d\n",omp_get_thread_num(), omp_get_ancestor_thread_num(1));
    //     }
    // }
    omp_set_num_threads(3);

    // #pragma omp parallel 
    // {
    //     switch (omp_get_thread_num())
    //     {
    //         case 0:
    //         { 
    //             cout << "1" << endl; 
    //             n_threads--;
    //             break;
    //         }

    //         default: gg();
    //     }
    // }

    omp_set_nested(1);

    long a = 0;

    #pragma omp parallel num_threads(2)
    {
        if (omp_get_thread_num() == 0){
            cout << "1" << endl; 
            
            for (long i=0; i<3000000000; i++) {
                a++;
            }

            cout << "end" << endl;
        }else{
            gg();
        }
    }

    // omp_set_nested(1);   /* make sure nested parallism is on */
    // int nprocs = 4;
    // int group1 = 1;
    // int group2 = nprocs - group1;

    // bool job_done {false};

    // #pragma omp parallel default(shared) shared(group1, group2) num_threads(1)
    // {
    //     #pragma omp task
    //     #pragma omp parallel num_threads(group1)
    //     {   
    //         if (!job_done) {
    //             printf("I am %d from group %d\n",omp_get_thread_num(), omp_get_ancestor_thread_num(1));
    //         } else {
    //             printf("job done!\n");
    //         }
            
    //     }

    //     #pragma omp task
    //     #pragma omp parallel num_threads(group2)
    //     {
    //         printf("I am %d from group %d\n",omp_get_thread_num(), omp_get_ancestor_thread_num(1));
    //     }

    //     printf("go\n");
    // }

    return 0;
}
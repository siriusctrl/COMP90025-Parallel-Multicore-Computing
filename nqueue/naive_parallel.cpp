/* solving the N-Queens problem using OpenMP
   usage with gcc (version 4.2 or higher required):
     g++ -fopenmp -o naive_parallel naive_parallel.c
     ./naive_parallel n numWorkers
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <vector>

#include <omp.h>

using namespace std;

int check_acceptable(int queen_rows[], int n)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = i+1; j < n; j++)
		{
			// two queens in the same row => not a solution!
			if (queen_rows[i] == queen_rows[j]) return 0;
			
			// two queens in the same diagonal => not a solution!
			if (queen_rows[i] - queen_rows[j] == i - j ||
			    queen_rows[i] - queen_rows[j] == j - i)
			    return 0;
		}
	}

	return 1;
}

int main(int argc, char* argv[])
{
    int n;
    int max_iter = 1;
    
    double start_time, end_time;
    int number_solutions = 0;

    std::vector<int> values;

	// config
	{
	    int num_workers;
        int i;
	
        n = (argc > 1) ? atoi(argv[1]) : 8;
        num_workers = (argc > 2) ? atoi(argv[2]) : 12;
        
        omp_set_num_threads(num_workers);
	    
        // maximum of n factorial possible combination
        for (i = 0; i < n; i++)
        {
            max_iter *= i;
            values.push_back(i);
        }
    }
  
    start_time = omp_get_wtime();
    
	#pragma omp parallel for
	for (int iter = 0; iter < max_iter; iter++)
	{
	    int queen_rows[n];
		// the index correspond to the queen's number and the queen's column
		// we only generate configurations where there's only one queen per column
            // for (i = 0; i < n; i++)
            // {
            //     queen_rows[i] = code % n;
                
            //     code /= n;
            // }
		
		if (check_acceptable(queen_rows, n))
		{
			#pragma omp atomic
		    number_solutions++;
		}
	}

    // get end time
    end_time = omp_get_wtime();
    // print results
    printf("The execution time is %g sec\n", end_time - start_time);
    printf("Number of found solutions is %d\n", number_solutions);
    
	return 0;
}

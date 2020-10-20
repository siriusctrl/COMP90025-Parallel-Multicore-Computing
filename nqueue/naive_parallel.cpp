#include <iostream>
#include <vector>
#include <algorithm>
#include <mpi.h>

using std::vector;
using std::cin;
using std::cout;
using std::endl;

void find_solution(int n, int max_solution);
bool check_acceptable(const vector<int> &queen_rows);

const MPI_Comm comm = MPI_COMM_WORLD;

int main(int argc, char *argv[])
{
    int n, max_solutions;

    n = (argc > 1) ? atoi(argv[1]) : 8;
    max_solutions = (argc > 2) ? atoi(argv[2]) : 92;

    int rank, prov;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov);
    MPI_Comm_rank(comm, &rank);

    cout << rank << endl;

    find_solution(n, max_solutions);
    MPI::Finalize();

    return 0;
}


void find_solution(int n, int max_solution) 
{
    int solution[n] = {0};
    int n_solutions {0};
    vector<int> permutation{};

    for (int i=0; i<n; ++i) {
        permutation.push_back(i);
    }

    do {
        if (check_acceptable(permutation)) {
            n_solutions++;

            for (auto p: permutation) {
                cout << p << " ";
            }

            cout << endl;
        }
    } 
    while(n_solutions < max_solution && std::next_permutation(permutation.begin(), permutation.end()));

    cout << max_solution << " solutions found!" << endl;
}


bool check_acceptable(const vector<int> &queen_rows)
{
	for (int i = 0; i < queen_rows.size(); ++i)
	{
		for (int j = i+1; j < queen_rows.size(); ++j)
		{
			// two queens in the same diagonal => not a solution!
			if (queen_rows[i] - queen_rows[j] == i - j ||
			    queen_rows[i] - queen_rows[j] == j - i)
			    return false;
		}
	}

	return true;
}
// mpicxx -std=c++14 naive_parallel.cpp -o naive_parallel -O3
#include <iostream>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <random>
#include <set>

using std::vector;
using std::cin;
using std::cout;
using std::endl;
using std::set;

void find_solution(int n, int max_solution, int rank);
bool check_acceptable(const vector<int> &queen_rows);

const MPI_Comm comm = MPI_COMM_WORLD;
constexpr int TERMINATION_SIGN = -1;
constexpr int CHECK_STATE_SIGN = -2;

int main(int argc, char *argv[])
{
    int n, max_solutions;

    n = (argc > 1) ? atoi(argv[1]) : 8;
    max_solutions = (argc > 2) ? atoi(argv[2]) : 92;

    int rank, prov;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov);
    MPI_Comm_rank(comm, &rank);

    int size {0};
    MPI_Comm_size(comm, &size);

    std::set<vector<int>> solution_set {};

    if (rank == 0)
    {
        MPI_Status status;
        int finished {0};
        int n_solutions {0};

        while (true)
        {
            if (finished == size-1) {
                cout << "found " << n_solutions << " solutions!" << endl;
                break;
            }

            // cout << "finished " << finished << endl;

            vector<int> res(n, 0);

            MPI_Recv(&res.front(), n, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

            if (res[0] == TERMINATION_SIGN) {
                finished++;
            } else if (res[0] == CHECK_STATE_SIGN) {
                // cout << "sync message received" << endl;
            } else {
                solution_set.insert(res);
                n_solutions = solution_set.size();

                if (n_solutions <= max_solutions) {
                    for (auto i:res) {
                        cout << i << " ";
                    }

                    cout << endl;
                }
            }

            if (n_solutions > max_solutions)
            {
                bool temp {false};
                MPI_Send(&temp, 1, MPI_CXX_BOOL, status.MPI_SOURCE, 0, comm);
                finished++;
            } else {
                bool temp {true};
                MPI_Send(&temp, 1, MPI_CXX_BOOL, status.MPI_SOURCE, 0, comm);
            }
        }

    } else {
        find_solution(n, max_solutions, rank);
    }
    
    MPI::Finalize();

    return 0;
}


void find_solution(int n, int max_solution, int rank) 
{
    int solution[n] = {0};
    int n_solutions {0};
    vector<int> permutation{};

    auto rng = std::default_random_engine {(unsigned) rank};

    for (int i=0; i<n; ++i) {
        permutation.push_back(i);
    }

    std::shuffle(permutation.begin(), permutation.end(), rng);

    cout << rank << " has permutation of ";
    for (auto i:permutation) {
        cout << i << " ";
    }

    cout << endl;

    int size;
    MPI_Comm_size(comm, &size);

    int counter {0};
    bool feedback;

    do {
        // frequent sync the state with master node to prevent working for
        // too long
        if (counter % 10000000 == 0) {
            vector<int> check_state(n, 0);
            check_state[0] = CHECK_STATE_SIGN;
            MPI::Send(&check_state.front(), n, MPI_INT, 0, 0, comm);
            MPI::Recv(&feedback, 1, MPI_CXX_BOOL, 0, MPI_ANY_TAG, comm, nullptr);

            if (!feedback)
            {
                cout << "node " << rank << " stopped" << endl;
                return;
            }
        }

        if (check_acceptable(permutation)) 
        {
            MPI_Send(&permutation.front(), n, MPI_INT, 0, 0, comm);
            MPI_Status status;
            MPI_Recv(&feedback, 1, MPI_CXX_BOOL, 0, MPI_ANY_TAG, comm, &status);

            if (!feedback) 
            {
                cout << "node " << rank << " stopped" << endl;
                return;
            }
        }

        counter++;
    }
    while(n_solutions < max_solution && std::next_permutation(permutation.begin(), permutation.end()));

    cout << "node " << rank << " has no more solution to try, and tried " << counter << " solutions" << endl;

    // no more candidates to try
    permutation[0] = TERMINATION_SIGN;
    MPI_Send(&permutation.front(), n, MPI_INT, 0, 0, comm);
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
// mpicxx -std=c++14 naive_shuffle.cpp -o naive_shuffle -O3
// mpirun -np 4 naive_shuffle 8 1
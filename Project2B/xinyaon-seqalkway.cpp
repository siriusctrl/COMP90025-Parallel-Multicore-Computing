// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/
// with many fixes and changes for multiple sequence alignment and to include an MPI driver
#include <mpi.h>
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include "sha512.hh"

using namespace std;

std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap, int *penalties);
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);
void do_MPI_task(int rank);

/*
Examples of sha512 which returns a std::string
sw::sha512::calculate("SHA512 of std::string") // hash of a string, or
sw::sha512::file(path) // hash of a file specified by its path, or
sw::sha512::calculate(&data, sizeof(data)) // hash of any block of data
*/

// Return current time, for performance measurement
uint64_t GetTimeStamp()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
}

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;

// Driver code
int main(int argc, char **argv)
{
    int rank;

    int prov;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov);
    // cout << "threads: " << prov << endl;

    MPI_Comm_rank(comm, &rank);
    if (rank == root)
    {
        int misMatchPenalty;
        int gapPenalty;
        int k;
        std::cin >> misMatchPenalty;
        std::cin >> gapPenalty;
        std::cin >> k;
        std::string genes[k];
        for (int i = 0; i < k; i++)
            std::cin >> genes[i];

        int numPairs = k * (k - 1) / 2;

        int penalties[numPairs];

        uint64_t start = GetTimeStamp();

        // return all the penalties and the hash of all allignments
        std::string alignmentHash = getMinimumPenalties(genes,
                                                        k, misMatchPenalty, gapPenalty,
                                                        penalties);

        // print the time taken to do the computation
        printf("Time: %ld us\n", (uint64_t)(GetTimeStamp() - start));

        // print the alginment hash
        std::cout << alignmentHash << std::endl;

        for (int i = 0; i < numPairs; i++)
        {
            std::cout << penalties[i] << " ";
        }
        std::cout << std::endl;
    }
    else
    {
        // do stuff for MPI tasks that are not rank==root
        do_MPI_task(rank);
    }
    MPI_Finalize();
    return 0;
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/
#include "omp.h"
#include "math.h"

#include <vector>
// #include <tuple>
#include <algorithm>
#include <queue>

constexpr int TASK_DISTRIBUTION_TAG {0};
constexpr int RESULT_COLLECTION_TAG = 1;
constexpr int STOP_SYMBOL = -9;
int n_threads = 16;
MPI_Datatype MPI_TASK;
MPI_Datatype MPI_RESULT;

// define the type for a task
typedef struct
{
    int i, j, id;
} TASK_t;

// define the type for a result
typedef struct
{
    int penalty, id;
    char problem_hash[129];
} RES_t;

inline MPI_Datatype create_MPI_TASK() {
    // define the task type for MPI
    MPI_Datatype MPI_TASK;
    MPI_Type_contiguous(3, MPI_INT, &MPI_TASK);
    MPI_Type_commit(&MPI_TASK);
    return MPI_TASK;
}

inline MPI_Datatype create_MPI_RESULT() {
    MPI_Datatype MPI_RESULT;
    int block_len_arry[3];
    MPI_Aint displacements[3];
    MPI_Datatype old[3];
    block_len_arry[0] = 1;
    displacements[0] = offsetof(RES_t, penalty);
    old[0] = MPI_INT;
    block_len_arry[1] = 1;
    displacements[1] = offsetof(RES_t, id);
    old[1] = MPI_INT;
    block_len_arry[2] = 129;
    displacements[2] = offsetof(RES_t, problem_hash);
    old[2] = MPI_CHAR;
    MPI_Type_create_struct(3, block_len_arry, displacements, old, &MPI_RESULT);
    MPI_Type_commit(&MPI_RESULT);
    return MPI_RESULT;
}

inline RES_t do_task(std::string x, std::string y, int task_id, int misMatchPenalty, int gapPenalty);

inline bool task_CMP(const RES_t &a, const RES_t &b) {
    return a.id < b.id;
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
inline int **new2d(int width, int height) {
    int **dp = new int *[width];
    size_t size = width;
    size *= height;
    int *dp0 = new int[size];
    if (!dp || !dp0)
    {
        std::cerr << "getMinimumPenalty: new failed" << std::endl;
        exit(1);
    }
    dp[0] = dp0;
    for (int i = 1; i < width; i++)
        dp[i] = dp[i - 1] + height;

    return dp;
}

// called by the root MPI task only
// this procedure should distribute work to other MPI tasks
// and put together results, etc.
std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap,
                                int *penalties) {

    std::string alignmentHash {};

    int size;
    MPI_Comm_size(comm, &size);

    int meta[3] = {k, pxy, pgap}; // k, pxy, pgap
    MPI_Bcast(meta, 3, MPI_INT, root, comm);

    // Broadcast the sequence length list
    int seq_length[k];

    for (int i = 0; i < k; ++i) {
        seq_length[i] = genes[i].length();
    }

    MPI_Bcast(seq_length, k, MPI_INT, root, comm);

    MPI_TASK = create_MPI_TASK();
    MPI_RESULT = create_MPI_RESULT();

    // Broadcast the sequences
    for (int i = 0; i < k; i++) {
        char buffer[seq_length[i]];
        memcpy(buffer, genes[i].c_str(), seq_length[i]);
        MPI_Bcast(buffer, seq_length[i], MPI_CHAR, root, comm);
    }

    omp_set_nested(1);
    n_threads--;

    #pragma omp parallel default(shared) num_threads(2)
    {
        if (omp_get_thread_num() == 0) {
            queue<TASK_t> task_queue;
            int task_id = 0;
            for (int i = 1; i < k; i++) {
                for (int j = 0; j < i; j++) {
                    task_queue.push({i, j, task_id++});
                }
            }

            MPI_Datatype MPI_TASK = create_MPI_TASK();
            MPI_Datatype MPI_RESULT = create_MPI_RESULT();
            
            
            // keep a list of whether each worker is done
            int finished = 0;
            int workers = size;

            for (int i = 0; i < size; i++) {
                if (!task_queue.empty()) {
                    // get the task
                    TASK_t task = task_queue.front();
                    task_queue.pop();
                    // send to worker i
                    MPI_Send(&task, 1, MPI_TASK, i, TASK_DISTRIBUTION_TAG, comm);
                } else {
                    // ask the worker to stop
                    TASK_t task = {STOP_SYMBOL, STOP_SYMBOL, STOP_SYMBOL};
                    MPI_Send(&task, 1, MPI_TASK, i, TASK_DISTRIBUTION_TAG, comm);
                    finished++;
                }
            }

            MPI_Status status;
            vector<RES_t> results = {};

            while (finished < workers) {   
                RES_t temp;
                MPI_Recv(&temp, 1, MPI_RESULT, MPI_ANY_SOURCE, RESULT_COLLECTION_TAG, comm, &status);
                results.push_back(temp);

                // if there are still task_queue to do
                if (!task_queue.empty()) {
                    // get the task
                    TASK_t task = task_queue.front();
                    task_queue.pop();
                    // send to worker i
                    MPI_Send(&task, 1, MPI_TASK, status.MPI_SOURCE, TASK_DISTRIBUTION_TAG, comm);
                } else {
                    // tell them to stop if no more works
                    TASK_t task = {STOP_SYMBOL, STOP_SYMBOL, STOP_SYMBOL};
                    MPI_Send(&task, 1, MPI_TASK, status.MPI_SOURCE, TASK_DISTRIBUTION_TAG, comm);
                    finished++;
                }
            }

            std::sort(results.begin(), results.end(), task_CMP);

            for (int i = 0; i < results.size(); ++i) {
                alignmentHash = sw::sha512::calculate(alignmentHash.append(results[i].problem_hash));
                penalties[i] = results[i].penalty;
            }
        } else {
            MPI_Status status;
            // receive my initial task
            TASK_t my_task;
            MPI_Recv(&my_task, 1, MPI_TASK, root, TASK_DISTRIBUTION_TAG, comm, &status);
            int STOP = my_task.i;
            // uint64_t start, end, start1, end1;
            // start = GetTimeStamp();
            while (STOP != STOP_SYMBOL) {
                RES_t result = do_task(genes[my_task.i], genes[my_task.j], my_task.id, pxy, pgap);
                MPI_Send(&result, 1, MPI_RESULT, root, RESULT_COLLECTION_TAG, comm);
                MPI_Recv(&my_task, 1, MPI_TASK, root, TASK_DISTRIBUTION_TAG, comm, &status);
                STOP = my_task.i;
            }
            // end = GetTimeStamp();
            // cout << "rank[" << 0 << "] computes: " << end - start << endl;
        }
    }

    MPI_Type_free(&MPI_TASK);
    MPI_Type_free(&MPI_RESULT);

    return alignmentHash;
}

// called for all tasks with rank!=root
// do stuff for each MPI task based on rank
void do_MPI_task(int rank) {
    MPI_Status status;
    int meta[3]; // k, pxy, pgap
    MPI_Bcast(meta, 3, MPI_INT, root, comm);

    int misMatchPenalty = meta[1];
    int gapPenalty = meta[2];
    int k = meta[0];

    // receive the sequneces length list
    int seq_length[k];
    MPI_Bcast(seq_length, k, MPI_INT, root, comm);

    MPI_TASK = create_MPI_TASK();
    MPI_RESULT = create_MPI_RESULT();

    // receive the gene sequences
    string genes[k];
    for (int i = 0; i < k; i++) {
        char buffer[seq_length[i] + 1];
        MPI_Bcast(buffer, seq_length[i], MPI_CHAR, root, comm);
        buffer[seq_length[i]] = '\0';
        genes[i] = string(buffer, seq_length[i]);
    }

    // receive my initial task
    TASK_t my_task;
    MPI_Recv(&my_task, 1, MPI_TASK, root, TASK_DISTRIBUTION_TAG, comm, &status);
    int STOP = my_task.i;

    // uint64_t start, end, start1, end1;
    // start = GetTimeStamp();

    while (STOP != STOP_SYMBOL) {
        RES_t result = do_task(genes[my_task.i], genes[my_task.j], my_task.id, misMatchPenalty, gapPenalty);
        MPI_Send(&result, 1, MPI_RESULT, root, RESULT_COLLECTION_TAG, comm);
        MPI_Recv(&my_task, 1, MPI_TASK, root, TASK_DISTRIBUTION_TAG, comm, &status);
        // cout << "rank-" << rank << ": i=" << my_task.i << ", j=" << my_task.j << ", task-id=" << my_task.id << endl;
        STOP = my_task.i;
    }

    // end = GetTimeStamp();
    // cout << "rank[" << rank << "] computes: " << end - start << endl;

    MPI_Type_free(&MPI_TASK);
    MPI_Type_free(&MPI_RESULT);
}

// function to find out the minimum penalty
// return the minimum penalty and put the aligned sequences in xans and yans
inline int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans) {

    int i, j;

    int m = x.length();
    int n = y.length();

    const int TILE_WIDTH {(int)ceil((1.0 * m) / thread_weight)};
    const int TILE_HEIGHT {(int)ceil((1.0 * n) / thread_weight)};
    const int TILE_M {(int)ceil((1.0 * m) / TILE_WIDTH)};
    const int TILE_N {(int)ceil((1.0 * n) / TILE_HEIGHT)};

    // bool pr = true;

    omp_set_num_threads(n_threads);
    int **dp = new2d(m + 1, n + 1);

    #pragma omp parallel
    {
        #pragma omp for nowait
        for (i = 0; i <= m; ++i) {
            dp[i][0] = i * pgap;
        }

        #pragma omp for nowait
        for (i = 0; i <= n; ++i) {
            dp[0][i] = i * pgap;
        }
    }

    int thread_weight = n_threads * 4;


    int total_diagonal = TILE_M + TILE_N - 1;
    int row_min, row_max, k;

    for (int line = 1; line <= total_diagonal; ++line) {
        row_min = max(1, line - TILE_N + 1);
        row_max = min(line, TILE_M);

        #pragma omp parallel for
        for (k = row_min; k <= row_max; ++k) {
            int i_lb = 1 + (k - 1) * TILE_WIDTH;
            int i_ub = min(i_lb + TILE_WIDTH, m + 1);
            int j_lb = 1 + (line - k) * TILE_HEIGHT;
            int j_ub = min(j_lb + TILE_HEIGHT, n + 1);

            for (int ii = i_lb; ii < i_ub; ++ii) {
                for (int jj = j_lb; jj < j_ub; ++jj) {
                    if (x[ii - 1] == y[jj - 1]) {
                        dp[ii][jj] = dp[ii - 1][jj - 1];
                    } else {
                        dp[ii][jj] = min(dp[ii - 1][jj - 1] + pxy,
                                        min(dp[ii - 1][jj] + pgap,
                                          dp[ii][jj - 1] + pgap));
                    }
                }
            }
        }
    }

    // Reconstructing the solution
    int l = n + m;

    i = m;
    j = n;

    int xpos = l;
    int ypos = l;

    while (!(i == 0 || j == 0)) {
        if (x[i - 1] == y[j - 1]) {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)y[j - 1];
            i--;
            j--;
        } else if (dp[i - 1][j - 1] + pxy == dp[i][j]) {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)y[j - 1];
            i--;
            j--;
        } else if (dp[i - 1][j] + pgap == dp[i][j]) {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)'_';
            i--;
        } else {
            xans[xpos--] = (int)'_';
            yans[ypos--] = (int)y[j - 1];
            j--;
        }
    }

    int x_diff = xpos - i, y_diff = ypos - j;
    #pragma omp parallel for
    for (int ii = i; ii > 0; --ii) {
        xans[ii + x_diff] = (int)x[ii - 1];
    }

    #pragma omp parallel for
    for (int x_pos2 = xpos - i; x_pos2 > 0; --x_pos2) {
        xans[x_pos2] = (int)'_';
    }

    #pragma omp parallel for
    for (int jj = j; jj > 0; --jj) {
        yans[jj + y_diff] = (int)y[jj - 1];
    }

    #pragma omp parallel for
    for (int y_pos2 = ypos - j; y_pos2 > 0; --y_pos2) {
        yans[y_pos2] = (int)'_';
    }

    int ret = dp[m][n];

    delete[] dp[0];
    delete[] dp;

    return ret;
}

inline RES_t do_task(std::string gene1, std::string gene2, int task_id, int misMatchPenalty, int gapPenalty)
{

    int m = gene1.length(); // length of gene1
    int n = gene2.length(); // length of gene2
    int l = m + n;
    int xans[l + 1], yans[l + 1];
    int penalty = getMinimumPenalty(gene1, gene2, misMatchPenalty, gapPenalty, xans, yans);

    int id = 1;
    int i;

    // find the start of the extra gap
    for (i = l; i >= 1; i--) {
        if ((char)yans[i] == '_' && (char)xans[i] == '_') {
            id = i + 1;
            break;
        }
    }

    // extract the exact alignment for both string
    std::string align1 = "";
    std::string align2 = "";
    for (i = id; i <= l; i++) {
        align1.append(1, (char)xans[i]);
    }
    
    for (i = id; i <= l; i++) {
        align2.append(1, (char)yans[i]);
    }

    // alignmentHash = hash(alignmentHash ++ hash(hash(align1)++hash(align2)))
    std::string align1hash = sw::sha512::calculate(align1);
    std::string align2hash = sw::sha512::calculate(align2);
    std::string problem_hash = sw::sha512::calculate(align1hash.append(align2hash));

    RES_t result;
    result.penalty = penalty;
    result.id = task_id;
    strcpy(result.problem_hash, problem_hash.c_str());

    return result;
}

// mpicxx -std=c++14 -fopenmp -o xinyaon-seqalkway xinyaon-seqalkway.cpp -O3
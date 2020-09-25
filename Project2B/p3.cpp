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

constexpr int JOB_DISTRIBUTION_TAG {0};
constexpr int RESULT_COLLECTION_TAG = 1;
constexpr int STOP_SYMBOL = -9;
int n_threads = 16;

// MPI_Datatype MPI_JOB = create_MPI_JOB();
// MPI_Datatype MPI_RESULT = create_MPI_RESULT();
MPI_Datatype MPI_JOB;
MPI_Datatype MPI_RESULT;

// define the type for a job
typedef struct {
    int i, j, id;
} JOB_t;

// define the type for a result
typedef struct {
    int penalty, id;
    char hash[129];
} RES_t;

inline MPI_Datatype create_MPI_JOB() {
    // define the job type for MPI
    MPI_Datatype MPI_JOB;
    MPI_Type_contiguous(3, MPI_INT, &MPI_JOB);
    MPI_Type_commit(&MPI_JOB);
    return MPI_JOB;
}

inline MPI_Datatype create_MPI_RESULT() {
    MPI_Datatype MPI_RESULT;
    int block_len[3];
    MPI_Aint displacement[3];
    MPI_Datatype oldtypes[3];
    block_len[0] = 1;
    displacement[0] = offsetof(RES_t, penalty);
    oldtypes[0] = MPI_INT;
    block_len[1] = 1;
    displacement[1] = offsetof(RES_t, id);
    oldtypes[1] = MPI_INT;
    block_len[2] = 129;
    displacement[2] = offsetof(RES_t, hash);
    oldtypes[2] = MPI_CHAR;
    MPI_Type_create_struct(3, block_len, displacement, oldtypes, &MPI_RESULT);
    MPI_Type_commit(&MPI_RESULT);
    return MPI_RESULT;
}

inline RES_t do_job(std::string x, std::string y, int job_id, int misMatchPenalty, int gapPenalty);

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
                                int *penalties)
{

    MPI_JOB = create_MPI_JOB();            
    MPI_RESULT = create_MPI_RESULT();

    std::string alignmentHash = "";

    int size;
    MPI_Comm_size(comm, &size);

    int meta[3] = {k, pxy, pgap}; // k, pxy, pgap
    MPI_Bcast(meta, 3, MPI_INT, root, comm);

    int seq_length[k];

    for (int i = 0; i < k; ++i) {
        seq_length[i] = genes[i].length();
    }
    MPI_Bcast(seq_length, k, MPI_INT, root, comm);

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
            queue<JOB_t> jobs;
            int job_id = 0;
            for (int i = 1; i < k; i++) {
                for (int j = 0; j < i; j++) {
                    jobs.push({i, j, job_id++});
                }
            }
            
            // keep a list of whether each worker is done
            int finished = 0;
            int workers = size;

            for (int i = 0; i < size; i++) {
                if (!jobs.empty()) {
                    // get the job
                    JOB_t job = jobs.front();
                    jobs.pop();
                    // send to worker i
                    MPI_Send(&job, 1, MPI_JOB, i, JOB_DISTRIBUTION_TAG, comm);
                } else {
                    // ask the worker to stop
                    JOB_t job = {STOP_SYMBOL, STOP_SYMBOL, STOP_SYMBOL};
                    MPI_Send(&job, 1, MPI_JOB, i, JOB_DISTRIBUTION_TAG, comm);
                    finished++;
                }
            }

            MPI_Status status;
            vector<RES_t> results = {};

            while (finished < workers) {
                RES_t temp;
                MPI_Recv(&temp, 1, MPI_RESULT, MPI_ANY_SOURCE, RESULT_COLLECTION_TAG, comm, &status);

                results.push_back(temp);
                // if there are still jobs to do
                if (!jobs.empty()) {
                    // get the job
                    JOB_t job = jobs.front();
                    jobs.pop();
                    // send to worker i
                    MPI_Send(&job, 1, MPI_JOB, status.MPI_SOURCE, JOB_DISTRIBUTION_TAG, comm);
                    // if there are nothing to do
                } else {
                    // ask the worker to stop
                    JOB_t job = {STOP_SYMBOL, STOP_SYMBOL, STOP_SYMBOL};
                    MPI_Send(&job, 1, MPI_JOB, status.MPI_SOURCE, JOB_DISTRIBUTION_TAG, comm);
                    finished++;
                }
            }

            // re-order the job into correct order
            std::sort(results.begin(), results.end(), [](auto const &a, auto const &b) {
                return a.id < b.id;
            });

            for (int i = 0; i < results.size(); ++i) {
                alignmentHash = sw::sha512::calculate(alignmentHash.append(results[i].hash));
                penalties[i] = results[i].penalty;
            }
        } else {
            MPI_Status status;
            // receive my initial job
            JOB_t my_job;
            MPI_Recv(&my_job, 1, MPI_JOB, root, JOB_DISTRIBUTION_TAG, comm, &status);

            int STOP = my_job.i;

            uint64_t start, end, start1, end1;
            // start = GetTimeStamp();
            while (STOP != STOP_SYMBOL) {
                RES_t result = do_job(genes[my_job.i], genes[my_job.j], my_job.id, pxy, pgap);
                MPI_Send(&result, 1, MPI_RESULT, root, RESULT_COLLECTION_TAG, comm);
                MPI_Recv(&my_job, 1, MPI_JOB, root, JOB_DISTRIBUTION_TAG, comm, &status);
                STOP = my_job.i;
            }
            // end = GetTimeStamp();
            // cout << "rank[" << 0 << "] computes: " << end - start << endl;
        }
    }

    MPI_Type_free(&MPI_JOB);
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

    MPI_JOB = create_MPI_JOB();            
    MPI_RESULT = create_MPI_RESULT();

    // receive the gene sequences
    string genes[k];
    for (int i = 0; i < k; i++) {
        char buffer[seq_length[i] + 1];
        MPI_Bcast(buffer, seq_length[i], MPI_CHAR, root, comm);
        buffer[seq_length[i]] = '\0';
        genes[i] = string(buffer, seq_length[i]);
    }

    // receive my initial job
    JOB_t my_job;
    MPI_Recv(&my_job, 1, MPI_JOB, root, JOB_DISTRIBUTION_TAG, comm, &status);
    int STOP = my_job.i;

    uint64_t start, end, start1, end1;
    // start = GetTimeStamp();
    while (STOP != STOP_SYMBOL) {
        RES_t result = do_job(genes[my_job.i], genes[my_job.j], my_job.id, misMatchPenalty, gapPenalty);
        MPI_Send(&result, 1, MPI_RESULT, root, RESULT_COLLECTION_TAG, comm);
        MPI_Recv(&my_job, 1, MPI_JOB, root, JOB_DISTRIBUTION_TAG, comm, &status);
        STOP = my_job.i;
    }
    // end = GetTimeStamp();
    // cout << "rank[" << rank << "] computes: " << end - start << endl;

    MPI_Type_free(&MPI_JOB);
    MPI_Type_free(&MPI_RESULT);
}

// function to find out the minimum penalty
// return the minimum penalty and put the aligned sequences in xans and yans
inline int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans) {

    int i, j; // intialising variables

    int m = x.length(); // length of gene1
    int n = y.length(); // length of gene2

    // table for storing optimal substructure answers
    omp_set_num_threads(n_threads);
    int **dp = new2d(m + 1, n + 1);

    // intialising the table
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (i = 0; i <= m; ++i)
        {
            dp[i][0] = i * pgap;
        }
        
        #pragma omp for nowait
        for (i = 0; i <= n; ++i)
        {
            dp[0][i] = i * pgap;
        }
    }


    const int TILE_WIDTH {(int)ceil((1.0 * m) / (n_threads*2))};
    const int TILE_LENGTH {(int)ceil((1.0 * n) / (n_threads*2))};

    const int W_TILES {(int)ceil((1.0 * m) / TILE_WIDTH)};
    const int H_TILES {(int)ceil((1.0 * n) / TILE_LENGTH)};

    const int TOTAL_LINE = W_TILES + H_TILES - 1;
    int k;

    for (int line = 1; line <= TOTAL_LINE; ++line) {
        const int start_col = max(1, line - H_TILES + 1);
        const int count_of_line = min(line, W_TILES);

        #pragma omp parallel for
        for (k = start_col; k <= count_of_line; ++k) {
            int i_lb = 1 + (k - 1) * TILE_WIDTH;
            int i_ub = min(i_lb + TILE_WIDTH, m + 1);
            int j_lb = 1 + (line - k) * TILE_LENGTH;
            int j_ub = min(j_lb + TILE_LENGTH, n + 1);

            for (i = i_lb; i < i_ub; ++i) {
                for (j = j_lb; j < j_ub; ++j) {
                    if (x[i - 1] == y[j - 1]) {
                        dp[i][j] = dp[i - 1][j - 1];
                    } else {
                        dp[i][j] = min(dp[i - 1][j - 1] + pxy,
                                        min(dp[i - 1][j] + pgap,
                                          dp[i][j - 1] + pgap));
                    }
                }
            }
        }
    }

    int l = n + m;

    i = m;
    j = n;

    int xpos = l;
    int ypos = l;

    while (!(i == 0 || j == 0))
    {
        if (x[i - 1] == y[j - 1])
        {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)y[j - 1];
            i--;
            j--;
        }
        else if (dp[i - 1][j - 1] + pxy == dp[i][j])
        {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)y[j - 1];
            i--;
            j--;
        }
        else if (dp[i - 1][j] + pgap == dp[i][j])
        {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)'_';
            i--;
        }
        else
        {
            xans[xpos--] = (int)'_';
            yans[ypos--] = (int)y[j - 1];
            j--;
        }
    }

    int x_offset { xpos - i }, y_offset {ypos - j};

    // try to parallel the reconstruction step

    #pragma omp parallel for
    for (int ii = i; ii > 0; --ii)
    {
        xans[ii + x_offset] = (int)x[ii - 1];
    }

    #pragma omp parallel for
    for (int x_pos2 = xpos - i; x_pos2 > 0; --x_pos2) {
        xans[x_pos2] = (int)'_';
    }

    #pragma omp parallel for
    for (int jj = j; jj > 0; --jj) {
        yans[jj + y_offset] = (int)y[jj - 1];
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

inline RES_t do_job(std::string gene1, std::string gene2, int job_id, int misMatchPenalty, int gapPenalty) {

    int m = gene1.length(); // length of gene1
    int n = gene2.length(); // length of gene2
    int l = m + n;
    int xans[l + 1], yans[l + 1];
    int penalty = getMinimumPenalty(gene1, gene2, misMatchPenalty, gapPenalty, xans, yans);

    int id = 1;

    // find the start of the extra gap
    for (int i = l; i >= 1; i--) {
        if ((char)yans[i] == '_' && (char)xans[i] == '_') {
            id = i + 1;
            break;
        }
    }

    // extract the exact alignment for both string
    std::string align1 = "";
    std::string align2 = "";
    for (int i = id; i <= l; i++) {
        align1.append(1, (char)xans[i]);
    }

    for (int i = id; i <= l; i++) {
        align2.append(1, (char)yans[i]);
    }

    // alignmentHash = hash(alignmentHash ++ hash(hash(align1)++hash(align2)))
    std::string align1hash = sw::sha512::calculate(align1);
    std::string align2hash = sw::sha512::calculate(align2);
    std::string hash = sw::sha512::calculate(align1hash.append(align2hash));

    RES_t result;
    result.penalty = penalty;
    result.id = job_id;
    strcpy(result.hash, hash.c_str());

    return result;
}

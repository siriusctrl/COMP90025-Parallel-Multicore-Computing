// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ and
// fixed an error when initializing the dp array :-)
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <iomanip>

using namespace std;

int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);

// Return current time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t) 1000000 + tv.tv_usec;
}


// Driver code
int main() {
    int misMatchPenalty;
    int gapPenalty;
    std::string gene1;
    std::string gene2;
    std::cin >> misMatchPenalty;
    std::cin >> gapPenalty;
    std::cin >> gene1;
    std::cin >> gene2;
    std::cout << "misMatchPenalty=" << misMatchPenalty << std::endl;
    std::cout << "gapPenalty=" << gapPenalty << std::endl;

    int m = gene1.length(); // length of gene1
    int n = gene2.length(); // length of gene2
    int l = m + n;
    int xans[l + 1], yans[l + 1];

    uint64_t start = GetTimeStamp();

    // calling the function to calculate the result
    int penalty = getMinimumPenalty(gene1, gene2,
                                    misMatchPenalty, gapPenalty,
                                    xans, yans);

    // print the time taken to do the computation
    printf("Time: %ld us\n", (uint64_t)(GetTimeStamp() - start));

    // postprocessing of the answer, for printing results

    // Since we have assumed the answer to be n+m long,
    // we need to remove the extra gaps in the starting
    // id represents the index from which the arrays
    // xans, yans are useful
    int id = 1;
    int i;
    for (i = l; i >= 1; i--) {
        if ((char) yans[i] == '_' && (char) xans[i] == '_') {
            id = i + 1;
            break;
        }
    }

    // Printing the final answer
    std::cout << "Minimum Penalty in aligning the genes = ";
    std::cout << penalty << std::endl;
	std::cout << "The aligned genes are :" << std::endl;
	for (i = id; i <= l; i++)
	{
		std::cout<<(char)xans[i];
	}
	std::cout << "\n";
	for (i = id; i <= l; i++)
	{
		std::cout << (char)yans[i];
	}
	std::cout << "\n";

    return 0;
}

/******************************************************************************/
/* Do not change any lines above here.            */
/* All of your changes should be below this line. */
/******************************************************************************/

inline int min3(int a, int b, int c) {
    if (a <= b && a <= c) {
        return a;
    } else if (b <= a && b <= c) {
        return b;
    } else {
        return c;
    }
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
inline int **new2d(int width, int height) {
    int **dp = new int *[width];
    size_t size = width;
    size *= height;
    int *dp0 = new int[size];
    if (!dp || !dp0) {
        std::cerr << "getMinimumPenalty: new failed" << std::endl;
        exit(1);
    }
    dp[0] = dp0;
    for (int i = 1; i < width; i++)
        dp[i] = dp[i - 1] + height;

    return dp;
}

// function to find out the minimum penalty
// return the maximum penalty and put the aligned sequences in xans and yans
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
                      int *xans, int *yans) {
    
    constexpr int n_threads {6};
    omp_set_num_threads(n_threads);

    const int ROW {(int) (x.length() + 1)};
    const int COL {(int) (y.length() + 1)};

    // table for storing optimal substructure answers
    int **dp = new2d(ROW, COL);
//	size_t size = m + 1;
//	size *= n + 1;
//	memset (dp[0], 0, size);

    #pragma omp parallel
    {

        // intialising the table
        #pragma omp for nowait
        for (int i=0; i <= ROW-1; i++) {
            dp[i][0] = i * pgap;
        }

        #pragma omp for nowait
        for (int i=0; i <= COL-1; i++) {
            dp[0][i] = i * pgap;
        }
    }

    // Tile parallel
    

    const int tile_width { (int)ceil( ((float) x.length()) / n_threads )};
    const int tile_length { (int)ceil( ((float) y.length()) / n_threads )};

    const int width_tiles { (int)ceil( ((float) x.length()) / tile_width )};
    const int hight_tiles { (int)ceil( ((float) y.length()) / tile_length )};

    // cout << "tile_width :" << tile_width << " "<< "tile_length :" << tile_length<< endl;
    // cout << "width_tiles :" << width_tiles << " hight_tiles :" << hight_tiles<<endl;

    // There will be tile_width+COL-1 lines in the output
    for (int line = 1; line <= (width_tiles + hight_tiles - 1); line++) {

        /* Get column index of the first element in this line of output.
           The index is 0 for first tile_width lines and line - tile_width for remaining
           lines  */
        int start_col = max(0, line - width_tiles);

        /* Get count of elements in this line. The count of elements is
           equal to minimum of line number, tile_length-start_col and width_tiles */
        int count = min3(line, hight_tiles - start_col, width_tiles);

        /* Print elements of this line */
        #pragma omp parallel for
        for (int z = 0; z < count; z++) {

            // define bound for both i and j
            const int i_lb { (min(width_tiles, line) - z - 1) * tile_width + 1 };
            const int j_lb { (start_col + z) * tile_length + 1 };
            const int i_ub { min(i_lb + tile_width, ROW) };
            const int j_ub { min(j_lb + tile_length, COL) };

            for (int i = i_lb; i < i_ub; i++) {
                for (int j = j_lb; j < j_ub; j++) {

                    if (x[i - 1] == y[j - 1]) {
                        dp[i][j] = dp[i - 1][j - 1];
                    } else {
                        dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
                                dp[i - 1][j] + pgap ,
                                dp[i][j - 1] + pgap);
                    }
                }
            }
        }
    }

    // exit(1);

    // Reconstructing the solution
    int l = ROW + COL - 2; // maximum possible length

    int i = ROW - 1;
    int j = COL - 1;

    int xpos = l;
    int ypos = l;

    while (!(i == 0 || j == 0)) {
        // cout << "(i, j) ("<< i << ", " << j << ")" << endl;

        if (x[i - 1] == y[j - 1]) {
            xans[xpos--] = (int) x[i - 1];
            yans[ypos--] = (int) y[j - 1];
            i--;
            j--;
        } else if (dp[i - 1][j - 1] + pxy == dp[i][j]) {
            xans[xpos--] = (int) x[i - 1];
            yans[ypos--] = (int) y[j - 1];
            i--;
            j--;
        } else if (dp[i - 1][j] + pgap == dp[i][j]) {
            xans[xpos--] = (int) x[i - 1];
            yans[ypos--] = (int) '_';
            i--;
        } else if (dp[i][j - 1] + pgap == dp[i][j]) {
            xans[xpos--] = (int) '_';
            yans[ypos--] = (int) y[j - 1];
            j--;
        }
    }
    while (xpos > 0) {
        if (i > 0) xans[xpos--] = (int) x[--i];
        else xans[xpos--] = (int) '_';
    }
    while (ypos > 0) {
        if (j > 0) yans[ypos--] = (int) y[--j];
        else yans[ypos--] = (int) '_';
    }

    int ret = dp[ROW-1][COL-1];

    delete[] dp[0];
    delete[] dp;

    return ret;
}

//g++ -std=c++14 -fopenmp xinyaon-seqalignomp.cpp -o xinyaon-seqalignomp -O3
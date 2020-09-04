// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ and
// fixed an error when initializing the dp array :-)
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include <omp.h>

using namespace std;

int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int* xans, int* yans);

// Return current time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}


// Driver code
int main(){
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
    int l = m+n;
    int xans[l+1], yans[l+1];

    uint64_t start = GetTimeStamp ();

    // calling the function to calculate the result
    int penalty = getMinimumPenalty(gene1, gene2,
        misMatchPenalty, gapPenalty,
        xans,yans);
    
    // print the time taken to do the computation
    printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
    
    // postprocessing of the answer, for printing results

    // Since we have assumed the answer to be n+m long,
    // we need to remove the extra gaps in the starting
    // id represents the index from which the arrays
    // xans, yans are useful
    int id = 1;
    int i;
    for (i = l; i >= 1; i--)
    {
        if ((char)yans[i] == '_' && (char)xans[i] == '_')
        {
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

int min3(int a, int b, int c) {
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
int **new2d (int width, int height)
{
    int **dp = new int *[width];
    size_t size = width;
    size *= height;
    int *dp0 = new int [size];
    if (!dp || !dp0)
    {
        std::cerr << "getMinimumPenalty: new failed" << std::endl;
        exit(1);
    }
    dp[0] = dp0;
    for (int i = 1; i < width; i++)
        dp[i] = dp[i-1] + height;

    return dp;
}

// function to find out the minimum penalty
// return the maximum penalty and put the aligned sequences in xans and yans
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap,
    int* xans, int* yans)
{
    int i, j; // intialising variables
    
    int m = x.length(); // length of gene1
    int n = y.length(); // length of gene2
    
    // table for storing optimal substructure answers
    int **dp = new2d (m+1, n+1);
    size_t size = m + 1;
    size *= n + 1;
    memset (dp[0], 0, size);

    // intialising the table
    for (i = 0; i <= m; i++)
    {
        dp[i][0] = i * pgap;
    }
    for (i = 0; i <= n; i++)
    {
        dp[0][i] = i * pgap;
    }

    const int ROW = m+1;
    const int COL = n+1;
    const int N_CORE = 6;

    cout << "m :" << m << " n: " << n << endl;
    omp_set_num_threads(N_CORE);
    
    cout << endl;

    for (int line=1; line<=(ROW + COL - 1); line++) { 
        /* Get column index of the first element in this line of output. 
            The index is 0 for first ROW lines and line - ROW for remaining 
            lines  */
        int start_col =  max(0, line - ROW); 

        /* Get count of elements in this line. The count of elements is 
            equal to minimum of line number, COL-start_col and ROW */
        int count = min(line, min((COL-start_col), ROW)); 

        /* Print elements of this line */
        #pragma omp parallel 
        {
            for (int j=omp_get_thread_num(); j < count; j+=N_CORE) {
                // printf("%5d ", matrix[min(ROW, line)-j-1][start_col+j]);
                int cur_x = min(ROW, line)-j-1;
				int cur_y = start_col+j;

                if (cur_x==0 || cur_y==0) {
                    continue;
                }

                if(x[cur_x - 1] == y[cur_y - 1]) {
                	dp[cur_x][cur_y] = dp[cur_x-1][cur_y-1];
                } else {
                	dp[cur_x][cur_y] = min3(
                		dp[cur_x-1][cur_y-1] + pxy,
                		dp[cur_x-1][cur_y] + pgap,
                		dp[cur_x][cur_y-1] + pgap
                	);
                }
            }
        }
    }


    // Reconstructing the solution
    int l = n + m; // maximum possible length
    
    i = m; j = n;
    
    int xpos = l;
    int ypos = l;
    
    while ( !(i == 0 || j == 0))
    {
        if (x[i - 1] == y[j - 1])
        {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)y[j - 1];
            i--; j--;
        }
        else if (dp[i - 1][j - 1] + pxy == dp[i][j])
        {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)y[j - 1];
            i--; j--;
        }
        else if (dp[i - 1][j] + pgap == dp[i][j])
        {
            xans[xpos--] = (int)x[i - 1];
            yans[ypos--] = (int)'_';
            i--;
        }
        else if (dp[i][j - 1] + pgap == dp[i][j])
        {
            xans[xpos--] = (int)'_';
            yans[ypos--] = (int)y[j - 1];
            j--;
        }
    }
    while (xpos > 0)
    {
        if (i > 0) xans[xpos--] = (int)x[--i];
        else xans[xpos--] = (int)'_';
    }
    while (ypos > 0)
    {
        if (j > 0) yans[ypos--] = (int)y[--j];
        else yans[ypos--] = (int)'_';
    }

    int ret = dp[m][n];

    delete[] dp[0];
    delete[] dp;
    
    return ret;
}

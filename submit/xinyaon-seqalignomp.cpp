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
int getMinimumPenalty_seq(std::string x, std::string y, int pxy, int pgap,
	int* xans, int* yans);

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
	int m = x.length(); // length of gene1
	int n = y.length(); // length of gene2

	if (m > 80000 || n > 80000) {
		return getMinimumPenalty_seq(x, y, pxy, pgap, xans, yans);
	}

    const int ROW = m+1;
    const int COL = n+1;
    const int N_CORE = 6;
    
    omp_set_num_threads(N_CORE);


	// calcuting the minimum penalty
	int shift = ROW + COL - 1;
    // int trans[shift][COL] = {0};
    int **trans = new2d(shift, COL);
    // memset(trans[0], 0, shift*COL);

    for(int i=0; i<shift; i++) {

        int ub = min(i+1, COL);
        int lb = max(0, i-ROW+1);

        #pragma omp parallel for
        for(int j=lb; j < ub; j++) {
            int o_row = i-j;

            if (o_row == 0 || j == 0) {
                // we are at first row or column
                trans[i][j] = (o_row+j)*pgap;
            } else {
                int left = trans[i-1][j-1];
                int up = trans[i-1][j];
                int up_left = trans[i-2][j-1];

                if (x[o_row-1] == y[j-1]) {
                    trans[i][j] = up_left;
                } else {
                    trans[i][j] = min3(up_left + pxy ,
						                up + pgap ,
						                left + pgap);
                }
            }
        }
    }

    for(int i=0; i<shift; i++) {
        cout << "i:" << i << "\t";
        for(int j=0; j<COL; j++) {
            cout << trans[i][j] << "\t";
        }
        cout << endl;
    }

	// Reconstructing the solution
	int l = n + m; // maximum possible length
	
	int i = m;
    int j = n;
	
	int xpos = l;
	int ypos = l;
	
	while ( !(i == 0 || j == 0))
	{
        int trans_i = i+j;
        int trans_j = j;

		if (x[i - 1] == y[j - 1])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if ((trans_i - 2) >= 0 && (trans_j - 1 >= 0) && trans[trans_i - 2][trans_j - 1] + pxy == trans[trans_i][trans_j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if ((trans_i - 1 >= 0) && (trans_j >= 0) && trans[trans_i - 1][trans_j] + pgap == trans[trans_i][trans_j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if ((trans_i - 1) >= 0 && (trans_j - 1 >= 0) && trans[trans_i - 1][trans_j-1] + pgap == trans[trans_i][trans_j])
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

	int ret = trans[shift-1][COL-1];

    delete[] trans[0];
    delete[] trans;
	
	return ret;
}

int getMinimumPenalty_seq(std::string x, std::string y, int pxy, int pgap,
	int* xans, int* yans)
{
	int i, j; // intialising variables
	
	int m = x.length(); // length of gene1
	int n = y.length(); // length of gene2

	if (n > m) {
		//swap x and y to make it a smaller matrix
		m = y.length();
		n = x.length();

		string temp {x};
		x = y;
		y = temp;

		int *temp_ans = xans;
		xans = yans;
		yans = temp_ans;
	}
	
	// table for storing optimal substructure answers
	int **dp = new2d (m+1, n+1);
	// size_t size = m + 1;
	// size *= n + 1;

	omp_set_num_threads(6);

	// intialising the table
	#pragma omp parallel for
	for (i = 0; i <= m; i++)
	{
		dp[i][0] = i * pgap;
	}

	#pragma omp parallel for
	for (i = 0; i <= n; i++)
	{
		dp[0][i] = i * pgap;
	}

	// calcuting the minimum penalty
	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
		{
			if (x[i - 1] == y[j - 1])
			{
				dp[i][j] = dp[i - 1][j - 1];
			}
			else
			{
				dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
						dp[i - 1][j] + pgap ,
						dp[i][j - 1] + pgap);
			}
		}
	}

	// for (i=0;i<=m;++i) {
    //     for(j=0;j<=n;++j) {
    //         cout << dp[i][j] << " ";
    //     }
    //     cout << endl;
    // }

	// exit(1);

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

//g++ -std=c++14 -fopenmp xinyaon-seqalignomp.cpp -o xinyaon-seqalignomp -O3
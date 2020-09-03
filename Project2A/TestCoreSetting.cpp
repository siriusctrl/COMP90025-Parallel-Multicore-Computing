#include <omp.h>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char const *argv[]) {
    const int ROW = 5;
    const int COL = 6;
    const int N_CORE = 1;

    omp_set_num_threads(N_CORE);
    
    int matrix[ROW][COL] = {1, 2, 3, 4, 5, 6,  
                            7, 8, 9, 10, 11, 12, 
                            13, 14, 15, 16, 17, 18,  
                            19, 20, 21, 22, 23, 24,  
                            25, 26, 27, 28, 29, 30};

    // There will be ROW+COL-1 lines in the output 
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

				// printf("(%d,%d)\t", cur_x, cur_y);
                cout << matrix[cur_x][cur_y] << "\t";
            }
        }

        cout << endl;
    }

    return 0;
}

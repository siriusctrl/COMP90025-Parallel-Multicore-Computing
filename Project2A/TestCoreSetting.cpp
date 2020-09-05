#include <omp.h>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>

using namespace std;

void test_diagonal() {
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
}

void matrix_transform_vertical() {
    const int ROW = 5;
    const int COL = 6;
    const int N_CORE = 1;

    omp_set_num_threads(N_CORE);
    
    int matrix[ROW][COL] = {1,  2,  3,  4,  5,  6,  
                            7,  8,  9,  10, 11, 12, 
                            13, 14, 15, 16, 17, 18,  
                            19, 20, 21, 22, 23, 24,  
                            25, 26, 27, 28, 29, 30};

    int shift = COL+min(ROW, COL)-1;
    int trans[ROW][COL+min(ROW, COL)] = {0};

    for (int i=1; i < ROW; i++) {

        int upper_bound = min(shift, i+COL);
        
        for (int j=i+1; j < upper_bound; j++) {

            int o_col = j-i;

            printf("i am (%d, %d) back to (%d, %d)\n", i, j, i, o_col);

            trans[i][j] = matrix[i][o_col];
        }
    }


    for (int i=0;i<ROW;i++) {
        for (int j=0;j<shift;j++) {
            cout << trans[i][j] << "\t";
        }
        cout << endl;
    }

}

void matrix_transform_horizontal() {
    const int ROW = 15;
    const int COL = 2;
    const int N_CORE = 1;

    omp_set_num_threads(N_CORE);
    
    int matrix[ROW][COL] = {1,  2,  3,  4,  5,  6,  
                            7,  8,  9,  10, 11, 12, 
                            13, 14, 15, 16, 17, 18,  
                            19, 20, 21, 22, 23, 24,  
                            25, 26, 27, 28, 29, 30};

    // int matrix[ROW][COL] = {1,  2,  3,  4,  5,  
    //                         6,  7,  8,  9,  10, 
    //                         11, 12, 13, 14, 15, 
    //                         16, 17, 18, 19, 20, 
    //                         21, 22, 23, 24, 25, 
    //                         26, 27, 28, 29, 30};


    // int matrix[ROW][COL] = {1,  2,   
    //                         3,  4,  
    //                         5,  6,  
    //                         7,  8,  
    //                         9,  10, 
    //                         11, 12, 
    //                         13, 14, 
    //                         15, 16, 
    //                         17, 18, 
    //                         19, 20, 
    //                         21, 22, 
    //                         23, 24, 
    //                         25, 26, 
    //                         27, 28, 
    //                         29, 30};

    int shift = ROW + COL - 1;
    int trans[shift][COL] = {0};

    for(int a=0; a<shift; a++) {
        // cout << "i:" << i << "\t";
        for(int b=0; b<COL; b++) {
            trans[a][b] = 0;
        }
    }


    for(int i=0; i<shift; i++) {

        // int ub = min(i+1, COL);
        // int lb = max(0, i - min(ROW, COL) + (ROW < COL));

        int ub = min(i+1, COL);
        int lb = max(0, i-ROW+1);
        // int lb;
        // if (i<=ROW) {
        //     lb = 0;
        // } else {
        //     lb = i - ROW;
        // }
        

        // printf("lb:%d, ub:%d\n", lb, ub);


        // cout << "i:" << i << endl;
        // for(int a=0; a<shift; a++) {
        //     cout << "a:" << a << "\t";
        //     for(int b=0; b<COL; b++) {
        //         cout << trans[a][b] << "\t";
        //     }
        //     cout << endl;
        // }
        
        // cout << endl;

        #pragma omp parallel for
        for(int j=lb; j < ub; j++) {
            int o_row = i-j;
            // int left = trans[i-1][j-1];
            // int up = trans[i-1][j];
            // int up_left = trans[i-2][j-1];

            // printf("I am (%d, %d), back to (%d, %d) ", i, j, o_row, j);
            // printf("left:%d, up:%d, up_left:%d", left, up, up_left);
            // cout << endl;

            // if(o_row < 0 || o_row >= ROW) {
            //     cout << "wrong ";
            //     printf("I am (%d, %d), back to (%d, %d) \n", i, j, o_row, j);
            //     continue;
            // }

            // if (i == 14) {
            //     cout << "j:" << j << endl;
            // }

            trans[i][j] = matrix[o_row][j];
        }
    }

    for(int i=0; i<shift; i++) {
        cout << "i:" << i << "\t";
        for(int j=0; j<COL; j++) {
            cout << trans[i][j] << "\t";
        }
        cout << endl;
    }
}

int min3(int a, int b, int c) {
    return min(min(a,b),c);
}

void h_test() {
    string x = "AGGGCT";
    string y = "AGGGA";

    const int pxy = 3;
    const int pgap = 2;

    const int ROW = x.length()+1;
    const int COL = y.length()+1;
    const int N_CORE = 1;
    
    
    omp_set_num_threads(N_CORE);

	// calcuting the minimum penalty
    int shift = COL + min(ROW-1, COL);
    int trans[shift][COL] = {0};

    for(int a=0; a<shift; a++) {
        for(int b=0; b<COL; b++) {
            trans[a][b] = 0;
        }
    }

    for(int i=0; i<shift; i++) {

        int ub = min(i+1, COL);
        int lb = max(0, i - min(ROW, COL) + (ROW < COL));

        printf("lb:%d, ub:%d\n", lb, ub);
        #pragma omp parallel for
        for(int j=lb; j < ub; j++) {
            int o_row = i-j;

            printf("accessing (%d, %d)\n", i, j);

            for(int a=0; a<shift; a++) {
                cout << "i:" << a << "\t";
                for(int b=0; b<COL; b++) {
                    cout << trans[a][b] << "\t";
                }
                cout << endl;
            }

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

            cout << endl;
        }
    }

    for(int a=0; a<shift; a++) {
        cout << "i:" << a << "\t";
        for(int b=0; b<COL; b++) {
            cout << trans[a][b] << "\t";
        }
        cout << endl;
    }
}

int main(int argc, char const *argv[]) {
    matrix_transform_horizontal();
    // h_test();

    return 0;
}

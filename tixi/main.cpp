#include <iostream>
#include<ctime>
#include <vector>
#include <armadillo>
#include <thread>
#include <immintrin.h>
using namespace std;
using namespace arma;
//randnum from 1 to 1000
int my_rand() {
    int ran = rand();
    return ran % 2;
}
void initial(vector<vector<float>>& p, int n) {
   
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            p[i][j] = my_rand();
        }
    }
}


void solution2(vector<vector<short>>& a, vector<vector<short>>& b, int n) {
    
    vector<vector<double>> res(n, vector<double>(n));
        for(int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j){
                 double sum = 0;
                for(int k = 0; k < n; ++k)
                         sum += a[i][k] * b[k][j];
                 res[i][j] += sum;
            }

       
}
void solution3(vector < vector < float >> &a, vector<vector<float>>& b, int n, vector<vector<float>>& res) {
    
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k) {
            double r = a[i][k];
            for (int j = 0; j < n; ++j)
                res[i][j] += r * b[k][j];
        }

}

void solution1(vector<vector<int>>& a, vector<vector<int>>& b, int n) {
    vector<vector<double>> res(n, vector<double>(n));
    for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k)
        {
            double r = b[k][j];
            for (int i = 0; i < n; ++i)
                res[i][j] += a[i][k] * r;
        }

    

   

}
//多线程12
void subs4(vector<vector<float>>& res, vector<vector<float>>& a, vector<vector<float>>& b, int n, int start, int subn) {
    for (int i = start; i < subn; ++i)
        for (int k = 0; k < n; ++k) {
            double r = a[i][k];
            for (int j = 0; j < n; ++j)
                res[i][j] += r * b[k][j];
        }

}

void solution4(vector<vector<float>>& a, vector<vector<float>>& b, int n) {
    vector<vector<float>> res(n, vector<float>(n));
    int subn = n / 12;
    std::thread th1(subs4, ref(res),ref(a), ref(b), n, 0, subn * 1);
    std::thread th2(subs4, ref(res), ref(a), ref(b), n, subn * 1, subn * 2);
    std::thread th3(subs4, ref(res), ref(a), ref(b), n, subn * 2, subn * 3);
    std::thread th4(subs4, ref(res), ref(a), ref(b), n, subn * 3, subn * 4);
    std::thread th5(subs4, ref(res), ref(a), ref(b), n, subn * 4, subn * 5);
    std::thread th6(subs4, ref(res), ref(a), ref(b), n, subn * 5, subn * 6);
    std::thread th7(subs4, ref(res), ref(a), ref(b), n, subn * 6, subn * 7);
    std::thread th8(subs4, ref(res), ref(a), ref(b), n, subn * 7, subn * 8);
    std::thread th9(subs4, ref(res), ref(a), ref(b), n, subn * 8, subn * 9);
    std::thread th10(subs4, ref(res), ref(a), ref(b), n, subn * 9, subn * 10);
    std::thread th11(subs4, ref(res), ref(a), ref(b), n, subn * 10, subn * 11);
    std::thread th12(subs4, ref(res), ref(a), ref(b), n, subn * 11, n);
    th1.join();
    th2.join();
    th3.join();
    th4.join();
    th5.join();
    th6.join();
    th7.join();
    th8.join();
    th9.join();
    th10.join();
    th11.join();
    th12.join();

   
   

}




//分块矩阵,向下为x+ 1 2
//                  3 4
void add(vector<vector<double>>&res, vector<vector<double>>&temp1, vector<vector<double>>&temp2,int n,int startx,int starty) {
   // cout << "resin:" << &res << endl;
    for (int i = 0; i < n ; i++) {
        for (int j = 0; j < n; j++) {
            res[i + startx][j + starty] = temp1[i][j] + temp2[i][j];
        }
    }
}

void subs5_1(vector<vector<short>>& a, vector<vector<short>>& b, int n, int startx, int starty, int startx_1, int starty_1, vector<vector<double>> temp) {
    vector<vector<double>> res(n, vector<double>(n));

    if (n == 1) {
        temp[0][0] = a[startx][starty] * b[startx_1][starty_1];

    }
    else {
        vector<vector<double>> temp1(n/2, vector<double>(n/2));
        vector<vector<double>> temp2(n/2, vector<double>(n/2));
        subs5_1(a, b, n / 2, startx, starty, startx_1, starty_1,temp1);
        
        subs5_1(a, b, n / 2, startx, starty + n / 2, startx_1 + n / 2, starty_1,temp2);
        add(res, temp1, temp2, n / 2, 0, 0);//11 23

        subs5_1(a, b, n / 2, startx, starty, startx_1, starty_1 + n / 2, temp1);
        subs5_1(a, b, n / 2, startx, starty + n / 2, startx_1 + n / 2, starty_1 + n / 2, temp2);
        add(res, temp1, temp2, n / 2, 0, n / 2);

        subs5_1(a, b, n / 2, startx + n / 2, starty, startx_1, starty_1, temp1);
        subs5_1(a, b, n / 2, startx + n / 2, starty + n / 2, startx_1 + n / 2, starty_1, temp2);
        add(res, temp1, temp2, n / 2, n / 2, 0);

        subs5_1(a, b, n / 2, startx + n / 2, starty, startx_1, starty_1 + n / 2, temp1);
        subs5_1(a, b, n / 2, startx + n / 2, starty + n / 2, startx_1 + n / 2, starty_1 + n / 2, temp2);
        add(res, temp1, temp2, n / 2, n / 2, n / 2);
    }
  
  


}

void sub6(vector<vector<float>>& a, vector<vector<float>>& b, vector<vector<float>>& c, int startx, int endx, int starty, int endy,int n,int type) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n ; j++) {
            c[i][j] = a[startx+i][endx+j] + type * b[starty + i][endy+j];
        }
    }
}
void sub6_11(vector<vector<float>>& p5, vector<vector<float>>& p4, vector<vector<float>>& p2, vector<vector<float>>& p6,vector<vector<float>>& res,int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j] = p5[i][j] + p4[i][j] - p2[i][j] + p6[i][j];
        }
    }
}
void sub6_22(vector<vector<float>>& p5, vector<vector<float>>& p1, vector<vector<float>>& p3, vector<vector<float>>& p7, vector<vector<float>>& res, int n) {
    for (int i = 0; i < n; i++) {
        for (int j =0; j < n; j++) {
            res[i+n][j+n] = p5[i][j] + p1[i][j] - p3[i][j] + p7[i][j];
        }
    }
}
void sub6_12(vector<vector<float>>& p1, vector<vector<float>>& p2, vector<vector<float>>& res, int n) {
    for (int i = 0; i <  n; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j+n] = p2[i][j] + p1[i][j];
        }
    }
}
void sub6_21(vector<vector<float>>& p3, vector<vector<float>>& p4, vector<vector<float>>& res, int n) {
    for (int i = 0; i <n; i++) {
        for (int j = 0; j < n; j++) {
            res[i+n][j] = p3[i][j] + p4[i][j];
        }
    }
}
void sub6_3(vector<vector<float>>& a, vector<vector<float>>& b, int n, vector<vector<float>>& res, int startx, int endx, int starty, int endy) {

    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k) {
            float r = a[i+startx][k+endx];
            for (int j = 0; j < n; ++j)
                res[i][j] += r * b[k+starty][j+endy];
        }

}


void solution10(vector<vector<float>>& a, vector<vector<float>>& b, int n, vector<vector<float>>& res, int startx, int endx, int starty, int endy) {

    __m256 c1, c2, bm[8], c3, r, c4, re, cload;
    vector<vector<float>> bT(n, vector<float>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            bT[j][i] = b[starty + i][endy + j];
        }
    }
    for (int i = startx; i < n + startx; i++)
        for (int k = endx; k < n + endx; k += 8) {
            r = _mm256_loadu_ps(&a[i][k]);
            for (int j = endy; j < n + endy; j += 8) {

                for (int h = 0; h < 8; h++) {
                    bm[h] = _mm256_loadu_ps(&bT[j + h-endy][k-endx]);
                }

                c2 = _mm256_setzero_ps();
                c1 = _mm256_dp_ps(r, bm[0], 0b11110001);
                c2 = _mm256_blend_ps(c2, c1, 0b00010001);
                c1 = _mm256_dp_ps(r, bm[1], 0b11110010);
                c2 = _mm256_blend_ps(c2, c1, 0b00100010);
                c1 = _mm256_dp_ps(r, bm[2], 0b11110100);
                c2 = _mm256_blend_ps(c2, c1, 0b01000100);
                c1 = _mm256_dp_ps(r, bm[3], 0b11111000);
                c2 = _mm256_blend_ps(c2, c1, 0b10001000);
                c3 = _mm256_insertf128_ps(c2, _mm256_extractf128_ps(c2, 0b00000001), 0b00000000);
                c3 = _mm256_add_ps(c3, c2);//低位有效


                c1 = _mm256_dp_ps(r, bm[4], 0b11110001);
                c2 = _mm256_blend_ps(c2, c1, 0b00010001);
                c1 = _mm256_dp_ps(r, bm[5], 0b11110010);
                c2 = _mm256_blend_ps(c2, c1, 0b00100010);
                c1 = _mm256_dp_ps(r, bm[6], 0b11110100);
                c2 = _mm256_blend_ps(c2, c1, 0b01000100);
                c1 = _mm256_dp_ps(r, bm[7], 0b11111000);
                c2 = _mm256_blend_ps(c2, c1, 0b10001000);
                c4 = _mm256_insertf128_ps(c2, _mm256_extractf128_ps(c2, 0b00000000), 0b00000001);
                c4 = _mm256_add_ps(c4, c2);//高位有效
                re = _mm256_blend_ps(c3, c4, 0b11110000);
                cload = _mm256_loadu_ps(&res[i - startx][j - endy]);
                re = _mm256_add_ps(re, cload);
                _mm256_store_ps(&res[i - startx][j - endy], re);
            }
        }


}



void solution8(vector<vector<float>>& a, vector<vector<float>>& b, int n, vector<vector<float>>& res, int startx, int endx, int starty, int endy) {

    __m256 c1, c2, bm[8], c3, r, c4, re, cload;
    for (int i = startx; i < n + startx; i++)
        for (int k = endx; k < n + endx; k += 8) {
            r = _mm256_loadu_ps(&a[i][k]);
            for (int j = endy; j < n + endy; j += 8) {

                for (int h = 0; h < 8; h++) {
                    bm[h] = _mm256_set_ps(b[k + 7 - startx + starty][h + j], b[k + 6 - startx + starty][h + j], b[k + 5 - startx + starty][h + j], b[k + 4 - startx + starty][h + j], b[k + 3 - startx + starty][h + j], b[k + 2 - startx + starty][h + j], b[k + 1 - startx + starty][h + j], b[k - startx + starty][h + j]);
                }

                c2 = _mm256_setzero_ps();
                c1 = _mm256_dp_ps(r, bm[0], 0b11110001);
                c2 = _mm256_blend_ps(c2, c1, 0b00010001);
                c1 = _mm256_dp_ps(r, bm[1], 0b11110010);
                c2 = _mm256_blend_ps(c2, c1, 0b00100010);
                c1 = _mm256_dp_ps(r, bm[2], 0b11110100);
                c2 = _mm256_blend_ps(c2, c1, 0b01000100);
                c1 = _mm256_dp_ps(r, bm[3], 0b11111000);
                c2 = _mm256_blend_ps(c2, c1, 0b10001000);
                c3 = _mm256_insertf128_ps(c2, _mm256_extractf128_ps(c2, 0b00000001), 0b00000000);
                c3 = _mm256_add_ps(c3, c2);//低位有效


                c1 = _mm256_dp_ps(r, bm[4], 0b11110001);
                c2 = _mm256_blend_ps(c2, c1, 0b00010001);
                c1 = _mm256_dp_ps(r, bm[5], 0b11110010);
                c2 = _mm256_blend_ps(c2, c1, 0b00100010);
                c1 = _mm256_dp_ps(r, bm[6], 0b11110100);
                c2 = _mm256_blend_ps(c2, c1, 0b01000100);
                c1 = _mm256_dp_ps(r, bm[7], 0b11111000);
                c2 = _mm256_blend_ps(c2, c1, 0b10001000);
                c4 = _mm256_insertf128_ps(c2, _mm256_extractf128_ps(c2, 0b00000000), 0b00000001);
                c4 = _mm256_add_ps(c4, c2);//高位有效
                re = _mm256_blend_ps(c3, c4, 0b11110000);
                cload = _mm256_loadu_ps(&res[i - startx][j - endy]);
                re = _mm256_add_ps(re, cload);
                _mm256_store_ps(&res[i - startx][j - endy], re);
            }
        }


}


//startx为x的起始纵坐标，endx为起始横坐标
void solution6(int n, vector<vector<float>>& a, vector<vector<float>>& b, int startx, int endx, int starty, int endy, vector<vector<float>>& res) {
    /*if (n == 1) {
        res[0][0] = a[startx][endx] * b[starty][endy];
        return;
    }*/
    if (n <= 64) {
        sub6_3(a, b, n, res, startx, endx, starty, endy);
        return;
    }
    vector<vector<float>> s1(n / 2, vector<float>(n / 2));
    vector<vector<float>> s2(n / 2, vector<float>(n / 2));
    vector<vector<float>> s3(n / 2, vector<float>(n / 2));
    vector<vector<float>> s4(n / 2, vector<float>(n / 2));
    vector<vector<float>> s5(n / 2, vector<float>(n / 2));
    vector<vector<float>> s6(n / 2, vector<float>(n / 2));
    vector<vector<float>> s7(n / 2, vector<float>(n / 2));
    vector<vector<float>> s8(n / 2, vector<float>(n / 2));
    vector<vector<float>> s9(n / 2, vector<float>(n / 2));
    vector<vector<float>> s10(n / 2, vector<float>(n / 2));
    sub6(b, b, s1, starty, endy + n / 2, starty + n / 2, endy + n / 2, n / 2, -1);
    sub6(a, a, s2, startx, endx, startx, endx + n / 2, n / 2, 1);
    sub6(a, a, s3, startx + n / 2, endx, startx + n / 2, endx + n / 2, n / 2, 1);
    sub6(b, b, s4, starty + n / 2, endy, starty, endy, n / 2, -1);
    sub6(a, a, s5, startx, endx, startx + n / 2, endx + n / 2, n / 2, 1);
    sub6(b, b, s6, starty, endy, starty + n / 2, endy + n / 2, n / 2, 1);
    sub6(a, a, s7, startx, endx + n / 2, startx + n / 2, endx + n / 2, n / 2, -1);
    sub6(b, b, s8, starty + n / 2, endy, starty + n / 2, endy + n / 2, n / 2, 1);
    sub6(a, a, s9, startx, endx, startx + n / 2, endx, n / 2, -1);
    sub6(b, b, s10, starty, endy, starty, endy + n / 2, n / 2, 1);
    vector<vector<float>> p1(n / 2, vector<float>(n / 2));
    vector<vector<float>> p2(n / 2, vector<float>(n / 2));
    vector<vector<float>> p3(n / 2, vector<float>(n / 2));
    vector<vector<float>> p4(n / 2, vector<float>(n / 2));
    vector<vector<float>> p5(n / 2, vector<float>(n / 2));
    vector<vector<float>> p6(n / 2, vector<float>(n / 2));
    vector<vector<float>> p7(n / 2, vector<float>(n / 2));
    solution6(n / 2, a, s1, startx, endx, 0, 0, p1);
    solution6(n / 2, s2, b, 0, 0, starty + n / 2, endy + n / 2, p2);
    solution6(n / 2, s3, b, 0, 0, starty, endy, p3);
    solution6(n / 2, a, s4, startx + n / 2, endx + n / 2, 0, 0, p4);
    solution6(n / 2, s5, s6, 0, 0, 0, 0, p5);
    solution6(n / 2, s7, s8, 0, 0, 0, 0, p6);
    solution6(n / 2, s9, s10, 0, 0, 0, 0, p7);
    sub6_11(p5, p4, p2, p6, res, n / 2);
    sub6_12(p1, p2, res, n / 2);
    sub6_21(p3, p4, res, n / 2);
    sub6_22(p5, p1, p3, p7, res, n / 2);

    return;


}

//_mm256_dp_ps
// _mm256_set_ps
//_mm256_setzero_ps  __m256 avx_sum0 = _mm256_setzero_ps();
//_mm256_loadu_ps    r0 = _mm256_loadu_ps(&a[i][k]);
//_mm256_permute_ps
// _mm256_blend_ps
// _mm256_insertf128_ps
// _mm256_extractf128_ps
// _mm256_add_ps
// _mm256_store_ps
//只会处理8*8的
void avx_matrix(vector<vector<float>>& a, vector<vector<float>>& b,vector<vector<float>>& res,int startx, int endx, int starty, int endy) {
    //加载a的数组的寄存器，行row加载，连续存储
    __m256 am[8],bm[8], c1,c2,c3,c4,re;
   
    for (int i = 0; i < 8; i++) {
        am[i] = _mm256_loadu_ps(&a[i+startx][0+endx]);
        
        bm[i] = _mm256_set_ps(b[7 + starty][i+endy], b[6 + starty][i + endy], b[5 + starty][i + endy], b[4 + starty][i + endy], b[3 + starty][i + endy], b[2 + starty][i + endy], b[1 + starty][i + endy], b[0+starty][i + endy]);
    }
    for (int i = 0; i < 8; i++) {
       c2 = _mm256_setzero_ps();
       c1 = _mm256_dp_ps(am[i], bm[0], 0b11110001);
       c2 = _mm256_blend_ps(c2,c1, 0b00010001);
       c1 = _mm256_dp_ps(am[i], bm[1], 0b11110010);
       c2 = _mm256_blend_ps(c2, c1, 0b00100010);
       c1 = _mm256_dp_ps(am[i], bm[2], 0b11110100);
       c2 = _mm256_blend_ps(c2, c1, 0b01000100);
       c1 = _mm256_dp_ps(am[i], bm[3], 0b11111000);
       c2 = _mm256_blend_ps(c2, c1, 0b10001000);
       c3 = _mm256_insertf128_ps(c2,_mm256_extractf128_ps(c2, 0b00000001),0b00000000);
       c3 = _mm256_add_ps(c3, c2);//低位有效


       c1 = _mm256_dp_ps(am[i], bm[4], 0b11110001);
       c2 = _mm256_blend_ps(c2, c1, 0b00010001);
       c1 = _mm256_dp_ps(am[i], bm[5], 0b11110010);
       c2 = _mm256_blend_ps(c2, c1, 0b00100010);
       c1 = _mm256_dp_ps(am[i], bm[6], 0b11110100);
       c2 = _mm256_blend_ps(c2, c1, 0b01000100);
       c1 = _mm256_dp_ps(am[i], bm[7], 0b11111000);
       c2 = _mm256_blend_ps(c2, c1, 0b10001000);
       c4 = _mm256_insertf128_ps(c2, _mm256_extractf128_ps(c2, 0b00000000), 0b00000001);
       c4 = _mm256_add_ps(c4, c2);
       re = _mm256_blend_ps(c3, c4, 0b11110000);
       _mm256_store_ps(&res[i][0], re);
    }
    

}

void solution7(int n, vector<vector<float>>& a, vector<vector<float>>& b, int startx, int endx, int starty, int endy, vector<vector<float>>& res) {
    /*if (n == 1) {
        res[0][0] = a[startx][endx] * b[starty][endy];
        return;
    }*/
    if (n <= 512) {
        solution10(a, b, n, res, startx, endx, starty, endy);
        return;
    }
    vector<vector<float>> s1(n / 2, vector<float>(n / 2));
    vector<vector<float>> s2(n / 2, vector<float>(n / 2));
    vector<vector<float>> s3(n / 2, vector<float>(n / 2));
    vector<vector<float>> s4(n / 2, vector<float>(n / 2));
    vector<vector<float>> s5(n / 2, vector<float>(n / 2));
    vector<vector<float>> s6(n / 2, vector<float>(n / 2));
    vector<vector<float>> s7(n / 2, vector<float>(n / 2));
    vector<vector<float>> s8(n / 2, vector<float>(n / 2));
    vector<vector<float>> s9(n / 2, vector<float>(n / 2));
    vector<vector<float>> s10(n / 2, vector<float>(n / 2));
    sub6(b, b, s1, starty, endy + n / 2, starty + n / 2, endy + n / 2, n / 2, -1);
    sub6(a, a, s2, startx, endx, startx, endx + n / 2, n / 2, 1);
    sub6(a, a, s3, startx + n / 2, endx, startx + n / 2, endx + n / 2, n / 2, 1);
    sub6(b, b, s4, starty + n / 2, endy, starty, endy, n / 2, -1);
    sub6(a, a, s5, startx, endx, startx + n / 2, endx + n / 2, n / 2, 1);
    sub6(b, b, s6, starty, endy, starty + n / 2, endy + n / 2, n / 2, 1);
    sub6(a, a, s7, startx, endx + n / 2, startx + n / 2, endx + n / 2, n / 2, -1);
    sub6(b, b, s8, starty + n / 2, endy, starty + n / 2, endy + n / 2, n / 2, 1);
    sub6(a, a, s9, startx, endx, startx + n / 2, endx, n / 2, -1);
    sub6(b, b, s10, starty, endy, starty, endy + n / 2, n / 2, 1);
    vector<vector<float>> p1(n / 2, vector<float>(n / 2));
    vector<vector<float>> p2(n / 2, vector<float>(n / 2));
    vector<vector<float>> p3(n / 2, vector<float>(n / 2));
    vector<vector<float>> p4(n / 2, vector<float>(n / 2));
    vector<vector<float>> p5(n / 2, vector<float>(n / 2));
    vector<vector<float>> p6(n / 2, vector<float>(n / 2));
    vector<vector<float>> p7(n / 2, vector<float>(n / 2));
    solution7(n / 2, a, s1, startx, endx, 0, 0, p1);
    solution7(n / 2, s2, b, 0, 0, starty + n / 2, endy + n / 2, p2);
    solution7(n / 2, s3, b, 0, 0, starty, endy, p3);
    solution7(n / 2, a, s4, startx + n / 2, endx + n / 2, 0, 0, p4);
    solution7(n / 2, s5, s6, 0, 0, 0, 0, p5);
    solution7(n / 2, s7, s8, 0, 0, 0, 0, p6);
    solution7(n / 2, s9, s10, 0, 0, 0, 0, p7);
    sub6_11(p5, p4, p2, p6, res, n / 2);
    sub6_12(p1, p2, res, n / 2);
    sub6_21(p3, p4, res, n / 2);
    sub6_22(p5, p1, p3, p7, res, n / 2);

    return;


}
void solution9(int n, vector<vector<float>>& a, vector<vector<float>>& b, int startx, int endx, int starty, int endy, vector<vector<float>>& res) {
    /*if (n == 1) {
        res[0][0] = a[startx][endx] * b[starty][endy];
        return;
    }*/
    if (n <= 64) {
        solution8(a, b, n,res, startx, endx, starty, endy);
        return;
    }
    vector<vector<float>> s1(n / 2, vector<float>(n / 2));
    vector<vector<float>> s2(n / 2, vector<float>(n / 2));
    vector<vector<float>> s3(n / 2, vector<float>(n / 2));
    vector<vector<float>> s4(n / 2, vector<float>(n / 2));
    vector<vector<float>> s5(n / 2, vector<float>(n / 2));
    vector<vector<float>> s6(n / 2, vector<float>(n / 2));
    vector<vector<float>> s7(n / 2, vector<float>(n / 2));
    vector<vector<float>> s8(n / 2, vector<float>(n / 2));
    vector<vector<float>> s9(n / 2, vector<float>(n / 2));
    vector<vector<float>> s10(n / 2, vector<float>(n / 2));
    sub6(b, b, s1, starty, endy + n / 2, starty + n / 2, endy + n / 2, n / 2, -1);
    sub6(a, a, s2, startx, endx, startx, endx + n / 2, n / 2, 1);
    sub6(a, a, s3, startx + n / 2, endx, startx + n / 2, endx + n / 2, n / 2, 1);
    sub6(b, b, s4, starty + n / 2, endy, starty, endy, n / 2, -1);
    sub6(a, a, s5, startx, endx, startx + n / 2, endx + n / 2, n / 2, 1);
    sub6(b, b, s6, starty, endy, starty + n / 2, endy + n / 2, n / 2, 1);
    sub6(a, a, s7, startx, endx + n / 2, startx + n / 2, endx + n / 2, n / 2, -1);
    sub6(b, b, s8, starty + n / 2, endy, starty + n / 2, endy + n / 2, n / 2, 1);
    sub6(a, a, s9, startx, endx, startx + n / 2, endx, n / 2, -1);
    sub6(b, b, s10, starty, endy, starty, endy + n / 2, n / 2, 1);
    vector<vector<float>> p1(n / 2, vector<float>(n / 2));
    vector<vector<float>> p2(n / 2, vector<float>(n / 2));
    vector<vector<float>> p3(n / 2, vector<float>(n / 2));
    vector<vector<float>> p4(n / 2, vector<float>(n / 2));
    vector<vector<float>> p5(n / 2, vector<float>(n / 2));
    vector<vector<float>> p6(n / 2, vector<float>(n / 2));
    vector<vector<float>> p7(n / 2, vector<float>(n / 2));
    solution9(n / 2, a, s1, startx, endx, 0, 0, p1);
    solution9(n / 2, s2, b, 0, 0, starty + n / 2, endy + n / 2, p2);
    solution9(n / 2, s3, b, 0, 0, starty, endy, p3);
    solution9(n / 2, a, s4, startx + n / 2, endx + n / 2, 0, 0, p4);
    solution9(n / 2, s5, s6, 0, 0, 0, 0, p5);
    solution9(n / 2, s7, s8, 0, 0, 0, 0, p6);
    solution9(n / 2, s9, s10, 0, 0, 0, 0, p7);
    sub6_11(p5, p4, p2, p6, res, n / 2);
    sub6_12(p1, p2, res, n / 2);
    sub6_21(p3, p4, res, n / 2);
    sub6_22(p5, p1, p3, p7, res, n / 2);

    return;


}


int main() {
    clock_t startTime, endTime;

    int n1 = 1024;

    vector<vector<float>> matrix1(n1, vector<float>(n1));
    vector<vector<float>> matrix2(n1, vector<float>(n1));
    vector<vector<float>> res(n1, vector<float>(n1));
    vector<vector<float>> res3(n1, vector<float>(n1));
    vector<vector<float>> res6(n1, vector<float>(n1));
    vector<vector<float>> res7(n1, vector<float>(n1));
    vector<vector<float>> res8(n1, vector<float>(n1));
    vector<vector<float>> res9(n1, vector<float>(n1));
    initial(matrix1, n1);
    initial(matrix2, n1);
   // cout << &matrix1[0][0] << " " << &matrix1[0][1] << endl;
    /*
    startTime = clock();//计时开始
    solution1(matrix1,matrix2, n1);
    endTime = clock();//计时结束
    cout << endTime << endl;
    cout << startTime << endl;
    cout << "不要局部性The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
   

    /*
    startTime = clock();//计时开始
    solution2(matrix1, matrix2, n1);
    endTime = clock();//计时结束
    cout << endTime << endl;
    cout << startTime << endl;
    cout << "局部性一般The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

        */

    startTime = clock();//计时开始
    solution3(matrix1, matrix2, n1,res3);
    endTime = clock();//计时结束
    cout << endTime << endl;
    cout << startTime << endl;
    cout << "局部性较好The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    
    /*for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1;j++) {
            cout << res3[i][j] << " ";
        }
        cout << endl;
    }
    cout << "------------------------------------------------------------------------------------------" << endl;*/

 

    startTime = clock();//计时开始
    solution4(matrix1, matrix2, n1);
    endTime = clock();//计时结束
    cout << "12线程The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    
   
   /* startTime = clock();//计时开始
    subs5_1(matrix1, matrix2, n1,0,0,0,0,res5);
    endTime = clock();//计时结束
    cout << endTime << endl;
    cout << startTime << endl;
    cout << "strassen " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;*/

    mat a = randu<mat>(n1, n1);
    mat b = randu<mat>(n1, n1);

    startTime = clock();//计时开始
    mat c =a* b;
    endTime = clock();//计时结束
    cout << endTime << endl;
    cout << startTime << endl;
    cout << "strassen " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

   // startTime = clock();//计时开始
   // solution6(n1,matrix1,matrix2,0,0,0,0,res6);
   // endTime = clock();//计时结束
   // cout << "局部性+strassen " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
   /* for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            cout << res6[i][j]<<" ";
        }
        cout << endl;
    }*/
    /*cout << "------------------------------------------------------------------------------------------" << endl;
    
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            cout << matrix1[i][j] << " ";
        }
        cout << endl;
    }
    cout << "------------------------------------------------------------------------------------------" << endl;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            cout << matrix2[i][j] << " ";
        }
        cout << endl;
    }
    cout << "------------------------------------------------------------------------------------------" << endl;*/
    startTime = clock();//计时开始
    solution7(n1, matrix1, matrix2, 0, 0, 0, 0, res7);
    endTime = clock();//计时结束
    cout << "AVX3 +strassem" << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
   /* for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1; j++) {
            cout << res7[i][j]<<" ";
        }
        cout << endl;
    }
    cout << "------------------------------------------------------------------------------------------" << endl;-*/

    startTime = clock();//计时开始
   // solution8( matrix1, matrix2,n1, res7,0,0,0,0);
    endTime = clock();//计时结束
    cout << "AVX2 " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
   /* for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            cout << res8[i][j] << " ";
        }
        cout << endl;
    }*/

    startTime = clock();//计时开始
    //solution9(n1, matrix1, matrix2, 0, 0, 0, 0, res9);
    endTime = clock();//计时结束
    cout << "AVX2+strassen " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
   /* for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1; j++) {
            cout << res9[i][j] << " ";
        }
        cout << endl;
    }*/

    startTime = clock();//计时开始
    solution10(matrix1, matrix2, n1, res8, 0, 0, 0, 0);
    endTime = clock();//计时结束
    cout << "AVX3 " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    /*for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1; j++) {
            cout << res8[i][j] << " ";
        }
        cout << endl;
    }*/
    
    return 0;
}

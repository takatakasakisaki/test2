// memcpy.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <intrin.h>
#include <immintrin.h>
#include <ammintrin.h>
#include <emmintrin.h>
#include <xmmintrin.h>


#include <iostream>
#include <windows.h>
#undef max
#undef min
#include <stdlib.h>   // For _MAX_PATH definition
#include <stdio.h>
#include <malloc.h>
#include <valarray>
#include <algorithm>
#include <execution>
#include <omp.h>
#if 0
#define PITCH 2496
#define WIDTH 2440 
#else
#define PITCH 2240
#define WIDTH 2240 
#endif
#define HEIGHT 220
//#define HEIGHT 221
//#define PITCH 64
//#define WIDTH 48
//#define HEIGHT 220
//#define NUMLOOP 10000
#define NUMLOOP 1
using namespace std;
void satsub(uint8_t* src0, uint8_t* src1, uint8_t* dst, int size)
{
    LARGE_INTEGER PerformanceCount00;
    QueryPerformanceCounter(&PerformanceCount00);
    for (int i = 0; i < size; i++) {
        uint8_t c0 = *src0++;
        uint8_t c1 = *src1++;
        uint8_t diff = c0 > c1 ? c0-c1: 0;
        *dst++ = diff;
    }
    LARGE_INTEGER PerformanceCount01;
    QueryPerformanceCounter(&PerformanceCount01);
    printf("satsub= %d\n", (uint32_t)(PerformanceCount01.QuadPart - PerformanceCount00.QuadPart));;
}
static void avxsatsub256(uint8_t* src0, uint8_t* src1, uint8_t *dst)
{
        __m256i c0, c1, c2, c3, c4, c5, c6, c7;
		c0 = _mm256_loadu_si256(((const __m256i*)src0) + 0);
		c1 = _mm256_loadu_si256(((const __m256i*)src1) + 0);
		__m256i s0 = _mm256_subs_epu8(c0, c1);
		_mm256_storeu_si256(((__m256i*)dst) + 0, s0);

		c2 = _mm256_loadu_si256(((const __m256i*)src0) + 1);
		c3 = _mm256_loadu_si256(((const __m256i*)src0) + 1);
		__m256i s1 = _mm256_subs_epu8(c2, c3);
		_mm256_storeu_si256(((__m256i*)dst) + 1, s1);

		c4 = _mm256_loadu_si256(((const __m256i*)src1) + 2);
		c5 = _mm256_loadu_si256(((const __m256i*)src1) + 2);
		__m256i s2 = _mm256_subs_epu8(c4, c5);
		_mm256_storeu_si256(((__m256i*)dst) + 2, s2);

		c6 = _mm256_loadu_si256(((const __m256i*)src1) + 3);
		c7 = _mm256_loadu_si256(((const __m256i*)src1) + 3);
		__m256i s3 = _mm256_subs_epu8(c6, c7);
		_mm256_storeu_si256(((__m256i*)dst) + 3, s3);

}
void avxsatsub(uint8_t* src0, uint8_t* src1, uint8_t* dst, int size)
{
    LARGE_INTEGER PerformanceCount00;
    QueryPerformanceCounter(&PerformanceCount00);
    _mm_prefetch((const char*)(src0), _MM_HINT_NTA);
    _mm_prefetch((const char*)(src1), _MM_HINT_NTA);

    long long i;
    int n = size / 256;
//#pragma omp parallel for num_threads(4)
    for (i = 0; i < n; i++) {// per 256 bytes
        avxsatsub256( &src0[ i * 256], &src1[i * 256], &dst[i * 256]);
        _mm_prefetch(((const char*)src0)+i* 512, _MM_HINT_NTA);
        _mm_prefetch(((const char*)src1)+i* 512, _MM_HINT_NTA);
    }
    size -= 256 * n;
    for (; size > 32; size -= 32) {
		__m256i  c0 = _mm256_loadu_si256(((const __m256i*)src0) + 0);
		__m256i  c1 = _mm256_loadu_si256(((const __m256i*)src1) + 0);
        src0 += 32;
        src1 += 32;
		__m256i s0 = _mm256_subs_epu8(c0, c1);
		_mm256_storeu_si256(((__m256i*)dst) + 0, s0);
        dst += 32;

    }
    if (size > 0) {
        satsub(src0, src1, dst, size);
    }
    LARGE_INTEGER PerformanceCount01;
    QueryPerformanceCounter(&PerformanceCount01);
    printf("avxsatsub= %d\n", (uint32_t)(PerformanceCount01.QuadPart - PerformanceCount00.QuadPart));;

}
char  avx2max(uint8_t* src, int size)
{
    LARGE_INTEGER PerformanceCount00;
    QueryPerformanceCounter(&PerformanceCount00);
    _mm_prefetch((const char*)(src), _MM_HINT_NTA);

    __m256i max_c0 = { 0 };
    __m256i max_c1 = { 0 };
    __m256i max_c2 = { 0 };
    __m256i max_c3 = { 0 };
    __m256i max_c4 = { 0 };
    __m256i max_c5 = { 0 };
    __m256i max_c6 = { 0 };
    __m256i max_c7 = { 0 };
    __m256i max_256;
    for (; size >= 256; size -= 256) {// per 256 bytes
        __m256i c0, c1, c2, c3, c4, c5, c6, c7;
        c0 = _mm256_loadu_si256(((const __m256i*)src) + 0); //load 256 BITS = 32bytes
        c1 = _mm256_loadu_si256(((const __m256i*)src) + 1);
        c2 = _mm256_loadu_si256(((const __m256i*)src) + 2);
        c3 = _mm256_loadu_si256(((const __m256i*)src) + 3);
        c4 = _mm256_loadu_si256(((const __m256i*)src) + 4);
        c5 = _mm256_loadu_si256(((const __m256i*)src) + 5);
        c6 = _mm256_loadu_si256(((const __m256i*)src) + 6);
        c7 = _mm256_loadu_si256(((const __m256i*)src) + 7);
        max_c0 = _mm256_max_epu8(max_c0, c0);
        max_c1 = _mm256_max_epu8(max_c1, c1);
        max_c2 = _mm256_max_epu8(max_c2, c2);
        max_c3 = _mm256_max_epu8(max_c3, c3);
        max_c4 = _mm256_max_epu8(max_c4, c4);
        max_c5 = _mm256_max_epu8(max_c5, c5);
        max_c6 = _mm256_max_epu8(max_c6, c6);
        max_c7 = _mm256_max_epu8(max_c7, c7);
        _mm_prefetch((const char*)(src+512), _MM_HINT_NTA);
        src += 256;
    }
    for (; size > 128; size -= 128) {// per 128 bytes
        __m256i c0, c1, c2, c3;
        c0 = _mm256_loadu_si256(((const __m256i*)src) + 0); //load 256 BITS = 32bytes
        c1 = _mm256_loadu_si256(((const __m256i*)src) + 1);
        c2 = _mm256_loadu_si256(((const __m256i*)src) + 2);
        c3 = _mm256_loadu_si256(((const __m256i*)src) + 3);
        max_c0 = _mm256_max_epu8(max_c0, c0);
        max_c1 = _mm256_max_epu8(max_c1, c0);
        max_c2 = _mm256_max_epu8(max_c2, c0);
        max_c3 = _mm256_max_epu8(max_c3, c0);
        _mm_prefetch((const char*)(src+256), _MM_HINT_NTA);
        src += 128;
    }
    if(size >= 64){
        __m256i c0, c1;
        c0 = _mm256_loadu_si256(((const __m256i*)src) + 0); //load 256 BITS = 32bytes
        c1 = _mm256_loadu_si256(((const __m256i*)src) + 1);
        max_c0 = _mm256_max_epu8(max_c0, c0);
        max_c1 = _mm256_max_epu8(max_c1, c1);
        size -= 64;
        src += 64;
    }
    if(size >= 32){
        __m256i c0, c1, c2, c3, c4, c5, c6, c7;
        c0 = _mm256_loadu_si256(((const __m256i*)src) + 0); //load 256 BITS = 32bytes
        max_c0 = _mm256_max_epu8(max_c0, c0);
        src += 32;
        size -= 32;//  // per 32 bytes
    }
    max_c0 = _mm256_max_epu8(max_c0, max_c1);
    max_c2 = _mm256_max_epu8(max_c2, max_c3);
    max_c4 = _mm256_max_epu8(max_c4, max_c5);
    max_c6 = _mm256_max_epu8(max_c6, max_c7);

    max_c0 = _mm256_max_epu8(max_c0, max_c2);
    max_c4 = _mm256_max_epu8(max_c4, max_c6);

    max_c0 = _mm256_max_epu8(max_c0, max_c4);

    uint8_t ret = 0;
    for(int i = 0; i < 32; i++){
        if (max_c0.m256i_u8[i] > ret) {
            ret = max_c0.m256i_u8[i];
        }
    }
    if (size > 0) {

        for(int i = 0; i < size; i++){
            uint8_t v = *src++;
            if (v > ret) {
                ret = v;
            }
        }
    }
#if 0
    __m128i i128_0, i128_1;
    i128_0 = _mm256_extracti128_si256(max_c0, 0);
    i128_1 = _mm256_extracti128_si256(max_c0, 1);
    __m128i max128;
    max128 = _mm_max_epu8(i128_0, i128_1);
    uint8_t* maxp = std::max_element(std::execution::par_unseq, max128.m128i_u8, max128.m128i_u8 + 16);
    uint8_t ret = *maxp;
#endif
    LARGE_INTEGER PerformanceCount01;
    QueryPerformanceCounter(&PerformanceCount01);
    printf("max=%d,%d\n", ret, (uint32_t)(PerformanceCount01.QuadPart - PerformanceCount00.QuadPart));;

    return ret;

}
uint8_t  maxvalarray(uint8_t* src, int size)
#if 1
{
    uint8_t maxval = 0;
    for (int i = 0; i < size; i++) {
        uint8_t d = *src++;
        if (d > maxval) {
            maxval = d;
        }
    }
    return maxval;
}
#endif
int main()
{
    std::cout << "Hello World!\n";
    for (volatile int kk = 0; kk < 100000; kk++) {
    }
    valarray<double> elappsed(NUMLOOP);
    valarray<double> elappsed1(NUMLOOP);
    for (int ii = 0; ii < NUMLOOP; ii++) {
        LARGE_INTEGER PerformanceCount01;
        QueryPerformanceCounter(&PerformanceCount01);
        size_t memsz = (((WIDTH * HEIGHT) + 0xfff) & ~0xfff);
        size_t memszp = (((PITCH * HEIGHT) + 0xfff) & ~0xfff);
#if 0
        char* memp = (char*)malloc(PITCH * HEIGHT);
        char* dst = (char*)malloc(WIDTH * HEIGHT);
#else
        uint8_t *dst = (uint8_t*)VirtualAlloc(nullptr, memsz, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
        uint8_t *memp = (uint8_t*)VirtualAlloc(nullptr, memszp, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
        uint8_t *src0 = (uint8_t*)VirtualAlloc(nullptr, memszp, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
        uint8_t *src1 = (uint8_t*)VirtualAlloc(nullptr, memszp, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
        uint8_t *dst1 = (uint8_t*)VirtualAlloc(nullptr, memszp, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
#endif
        //VirtualLock(dst, memsz);
        //VirtualLock(memp, memsz);
        //cout << (void*)dst << ","<<(void*)memp << endl;
        if (memp == nullptr) {
            return 1;
        }
        if (dst == nullptr) {
            return 1;
        }
#if 1
        //for (int i = 0; i < PITCH * HEIGHT; i+=4096) {
        for (int i = 0; i < PITCH * HEIGHT; i+=1)  {
            memp[i] = (i+1)&0x3f;
        }
        memp[WIDTH * HEIGHT - 1] = 0xfe;
#endif
        int size = WIDTH*HEIGHT;
#if 0
        uint8_t src0[64];
        uint8_t src1[64];
        uint8_t dst1[64];
        for (int i = 0; i < 32; i++) {
            src0[i] = i;
            src1[i] = i+32;
            dst[i] = 0xff;
            printf("%03u,", src0[i] - src1[i]);
            if ( (i&0xf) == 0xf) {
                puts("");
            }
        }
        puts("");
        puts("src0");
        for (int i = 0; i < 32; i++) {
            printf("%03d,", src0[i]);
            if ( (i&0xf) == 0xf) {
                puts("");
            }
        }
        puts("");
        puts("src1");
        for (int i = 0; i < 32; i++) {
            printf("%03d,", src1[i]);
            if ( (i&0xf) == 0xf) {
                puts("");
            }
        }
        puts("");
#endif
        avxsatsub(src0, src1, dst1,size);
#if 0
        puts("dst1");
        for (int i = 0; i < 32; i++) {
            printf("%03d,", dst1[i]);
            if ( (i&0xf) == 0xf) {
                puts("");
            }
        }
        puts("");
#endif
        satsub(src0, src1, dst1,size);
#if 0
        puts("");
        puts("dst1");
        for (int i = 0; i < 32; i++) {
            printf("%03d,", dst1[i]);
            if ( (i&0xf) == 0xf) {
                puts("");
            }
        }
        puts("");
#endif
        uint8_t  mm = avx2max(memp, size);




        printf("mm1=%d\n", mm);
    uint8_t* maxp = std::max_element(std::execution::par_unseq, memp, memp+size);
        printf("mm2=%d\n", *maxp);
    uint8_t  u8 = maxvalarray((uint8_t*)memp, size);
        printf("mm3=%d\n", u8);
#if 0
        //memset(dst, 0, WIDTH * HEIGHT);
        memset(memp, 0, PITCH * HEIGHT);
#endif
#if 0
        uint64_t addr = uint64_t(dst);
        addr &= ~0xfff;
        VirtualLock((void*)addr, memsz);
#endif
        //cout << (void*)dst << endl;
        LARGE_INTEGER PerformanceCount0;
        QueryPerformanceCounter(&PerformanceCount0);
#if 0
        memcpy(dst, memp, WIDTH*HEIGHT);
#else
//#pragma omp parallel for num_threads(4)
        for (int h = 0; h < HEIGHT; h++) {
#if 0
            int srcidx = h * PITCH;
            int dstidx = h * WIDTH;
            memcpy(&dst[dstidx], &memp[srcidx], WIDTH);
#else
#pragma vector nontemporal

            for (int w = 0; w < WIDTH; w++) {
                int srcidx = h * PITCH + w;
                int dstidx = h * WIDTH + w;
                dst[dstidx] = memp[srcidx];
            }
#endif
        }
#endif
        LARGE_INTEGER PerformanceCount1;
        QueryPerformanceCounter(&PerformanceCount1);
        long long count0 = PerformanceCount0.QuadPart;
        long long count1 = PerformanceCount1.QuadPart;
        long long count01 = PerformanceCount01.QuadPart;
        elappsed[ii] = (count1 - count0)/10.0;
        elappsed1[ii] = (count1 - count01)/10.0;
        //std::cout << "count0," << count0 << ",count1," << count1 << ",diff," << count1 - count0 << endl;
        //VirtualFree(memp, memsz, MEM_DECOMMIT|MEM_RELEASE);
        //VirtualFree(dst, memsz, MEM_DECOMMIT|MEM_RELEASE);
        VirtualFree(memp, 0, MEM_RELEASE);
        VirtualFree(dst, 0, MEM_RELEASE);
    }

    double sum = elappsed.sum();
    std::sort(std::begin(elappsed), std::end(elappsed));
    auto ave = sum / elappsed.size();
    cout << "ave=" << ave << endl;
    auto err = elappsed - ave;
    auto err2 = err * err;
    cout << "err2size" << err2.size() << endl;
    auto var = err2.sum() / err2.size();
    auto sigma = std::sqrt(var);
    double medium = elappsed[NUMLOOP / 2];
    cout << "ave="<<ave << ",var=" << var << ",sigma=" << sigma << ",max=" << elappsed.max() << ",min=" << elappsed.min() <<",mid=" << medium << endl;


    double sum1 = elappsed1.sum();
    std::sort(std::begin(elappsed1), std::end(elappsed1));
    auto ave1 = sum1 / elappsed1.size();
    cout << "ave1=" << ave1 << endl;
    auto err1 = elappsed1 - ave1;
    auto err21 = err1 * err1;
    cout << "err21size" << err21.size() << endl;
    auto var1 = err21.sum() / err21.size();
    auto sigma1 = std::sqrt(var1);
    double medium1 = elappsed1[NUMLOOP / 2];
    cout << "ave1="<<ave1 << ",var1=" << var1 << ",sigma1=" << sigma1 << ",max1=" << elappsed1.max() << ",min1=" << elappsed1.min() <<",mid1=" << medium1 << endl;




}

// プログラムの実行: Ctrl + F5 または [デバッグ] > [デバッグなしで開始] メニュー
// プログラムのデバッグ: F5 または [デバッグ] > [デバッグの開始] メニュー

// 作業を開始するためのヒント: 
//    1. ソリューション エクスプローラー ウィンドウを使用してファイルを追加/管理します 
//   2. チーム エクスプローラー ウィンドウを使用してソース管理に接続します
//   3. 出力ウィンドウを使用して、ビルド出力とその他のメッセージを表示します
//   4. エラー一覧ウィンドウを使用してエラーを表示します
//   5. [プロジェクト] > [新しい項目の追加] と移動して新しいコード ファイルを作成するか、[プロジェクト] > [既存の項目の追加] と移動して既存のコード ファイルをプロジェクトに追加します
//   6. 後ほどこのプロジェクトを再び開く場合、[ファイル] > [開く] > [プロジェクト] と移動して .sln ファイルを選択します

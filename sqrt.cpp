// sqrt.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#ifdef __GNUC__
#include <cpuid.h>
#else
#include <windows.h>
#undef max
#undef min
#include <intrin.h>
#endif
#include <immintrin.h>
//#include <xmmintrin.h>
#include <mmintrin.h>

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <malloc.h>
using namespace std;
__declspec(noinline)  void convsqrtavx2(uint32_t* src0, int shiftval, uint32_t len, uint8_t* dst)
{
    uint32_t size = len * sizeof(uint32_t);
    uint8_t* src = (uint8_t*)src0;
    auto ts0 = chrono::system_clock::now();
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
    
    for (int i = 0; size >= 256; size -= 256, i++) {// per 256 bytes
        __m256i c0, c1, c2, c3, c4, c5, c6, c7;
        //load 32bytes x 8= 256bytes
        c0 = _mm256_loadu_si256(((const __m256i*)src) + 0); //load 256 BITS = 32bytes
        c1 = _mm256_loadu_si256(((const __m256i*)src) + 1);
        c2 = _mm256_loadu_si256(((const __m256i*)src) + 2);
        c3 = _mm256_loadu_si256(((const __m256i*)src) + 3);
        c4 = _mm256_loadu_si256(((const __m256i*)src) + 4);
        c5 = _mm256_loadu_si256(((const __m256i*)src) + 5);
        c6 = _mm256_loadu_si256(((const __m256i*)src) + 6);
        c7 = _mm256_loadu_si256(((const __m256i*)src) + 7);
        //shift
        __m256i scaled_0;
        __m256i scaled_1;
        __m256i scaled_2;
        __m256i scaled_3;
        __m256i scaled_4;
        __m256i scaled_5;
        __m256i scaled_6;
        __m256i scaled_7;
        scaled_0 = _mm256_slli_epi32(c0, shiftval);
        scaled_1 = _mm256_slli_epi32(c1, shiftval);
        scaled_2 = _mm256_slli_epi32(c2, shiftval);
        scaled_3 = _mm256_slli_epi32(c3, shiftval);
        scaled_4 = _mm256_slli_epi32(c4, shiftval);
        scaled_5 = _mm256_slli_epi32(c5, shiftval);
        scaled_6 = _mm256_slli_epi32(c6, shiftval);
        scaled_7 = _mm256_slli_epi32(c7, shiftval);
        // & 0xffff
        __m256i mask16 = _mm256_set1_epi32(0xffff);
        scaled_0 = _mm256_and_si256(scaled_0, mask16 );
        scaled_1 = _mm256_and_si256(scaled_1, mask16 );
        scaled_2 = _mm256_and_si256(scaled_2, mask16 );
        scaled_3 = _mm256_and_si256(scaled_3, mask16 );
        scaled_4 = _mm256_and_si256(scaled_4, mask16 );
        scaled_5 = _mm256_and_si256(scaled_5, mask16 );
        scaled_6 = _mm256_and_si256(scaled_6, mask16 );
        scaled_7 = _mm256_and_si256(scaled_7, mask16 );
        __m256 f_0;
        __m256 f_1;
        __m256 f_2;
        __m256 f_3;
        __m256 f_4;
        __m256 f_5;
        __m256 f_6;
        __m256 f_7;
        f_0 = _mm256_cvtepi32_ps(scaled_0);
        f_1 = _mm256_cvtepi32_ps(scaled_1);
        f_2 = _mm256_cvtepi32_ps(scaled_2);
        f_3 = _mm256_cvtepi32_ps(scaled_3);
        f_4 = _mm256_cvtepi32_ps(scaled_4);
        f_5 = _mm256_cvtepi32_ps(scaled_5);
        f_6 = _mm256_cvtepi32_ps(scaled_6);
        f_7 = _mm256_cvtepi32_ps(scaled_7);
        f_0 = _mm256_sqrt_ps(f_0);
        f_1 = _mm256_sqrt_ps(f_1);
        f_2 = _mm256_sqrt_ps(f_2);
        f_3 = _mm256_sqrt_ps(f_3);
        f_4 = _mm256_sqrt_ps(f_4);
        f_5 = _mm256_sqrt_ps(f_5);
        f_6 = _mm256_sqrt_ps(f_6);
        f_7 = _mm256_sqrt_ps(f_7);
        __m256i out32_0;
        __m256i out32_1;
        __m256i out32_2;
        __m256i out32_3;
        __m256i out32_4;
        __m256i out32_5;
        __m256i out32_6;
        __m256i out32_7;
        out32_0 = _mm256_cvttps_epi32(f_0);
        out32_1 = _mm256_cvttps_epi32(f_1);
        out32_2 = _mm256_cvttps_epi32(f_2);
        out32_3 = _mm256_cvttps_epi32(f_3);
        out32_4 = _mm256_cvttps_epi32(f_4);
        out32_5 = _mm256_cvttps_epi32(f_5);
        out32_6 = _mm256_cvttps_epi32(f_6);
        out32_7 = _mm256_cvttps_epi32(f_7);

        //make 64bit
        uint8_t dst0 = _mm256_extract_epi32(out32_0, 0) ;
        uint8_t dst1 = _mm256_extract_epi32(out32_1, 1) ;
        uint8_t dst2 = _mm256_extract_epi32(out32_2, 2) ;
        uint8_t dst3 = _mm256_extract_epi32(out32_3, 3) ;
        uint8_t dst4 = _mm256_extract_epi32(out32_4, 4) ;
        uint8_t dst5 = _mm256_extract_epi32(out32_5, 5) ;
        uint8_t dst6 = _mm256_extract_epi32(out32_6, 6) ;
        uint8_t dst7 = _mm256_extract_epi32(out32_7, 7) ;
        printf("%02x ", dst0);
        printf("%02x ", dst1);
        printf("%02x ", dst2);
        printf("%02x ", dst3);
        printf("%02x ", dst4);
        printf("%02x ", dst5);
        printf("%02x ", dst6);
        printf("%02x ", dst7);
        puts("");
        uint64_t dst64 = 
            ((uint64_t)dst0 )|
            ((uint64_t)dst1 << 8 )|
            ((uint64_t)dst2 << 16) |
            ((uint64_t)dst3 << 24) |
            ((uint64_t)dst4 << 32) |
            ((uint64_t)dst5 << 40) |
            ((uint64_t)dst6 << 48) |
            ((uint64_t)dst7 << 54)
            ;
        //__m64 d64;
        //d64 = _mm_set_pi8(dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7);
        //_mm_stream_pi((__m64*)dst, d64);
        printf("%16llx\n", dst64);
        *(uint64_t*)dst = dst64;
        dst += 8;

        _mm_prefetch((const char*)(src+512), _MM_HINT_NTA);
        src += 256;
    }
    for (; size > 128; size -= 128) {// per 128 bytes
        __m256i c0, c1, c2, c3;
        c0 = _mm256_loadu_si256(((const __m256i*)src) + 0); //load 256 BITS = 32bytes
        c1 = _mm256_loadu_si256(((const __m256i*)src) + 1);
        c2 = _mm256_loadu_si256(((const __m256i*)src) + 2);
        c3 = _mm256_loadu_si256(((const __m256i*)src) + 3);
        max_c0 = _mm256_max_epu32(max_c0, c0);
        max_c1 = _mm256_max_epu32(max_c1, c0);
        max_c2 = _mm256_max_epu32(max_c2, c0);
        max_c3 = _mm256_max_epu32(max_c3, c0);
        _mm_prefetch((const char*)(src+256), _MM_HINT_NTA);
        src += 128;
    }
    if(size >= 64){
        __m256i c0, c1;
        c0 = _mm256_loadu_si256(((const __m256i*)src) + 0); //load 256 BITS = 32bytes
        c1 = _mm256_loadu_si256(((const __m256i*)src) + 1);
        max_c0 = _mm256_max_epu32(max_c0, c0);
        max_c1 = _mm256_max_epu32(max_c1, c1);
        size -= 64;
        src += 64;
    }
    if(size >= 32){
        __m256i c0, c1, c2, c3, c4, c5, c6, c7;
        c0 = _mm256_loadu_si256(((const __m256i*)src) + 0); //load 256 BITS = 32bytes
        max_c0 = _mm256_max_epu32(max_c0, c0);
        src += 32;
        size -= 32;//  // per 32 bytes
    }
    max_c0 = _mm256_max_epu32(max_c0, max_c1);
    max_c2 = _mm256_max_epu32(max_c2, max_c3);
    max_c4 = _mm256_max_epu32(max_c4, max_c5);
    max_c6 = _mm256_max_epu32(max_c6, max_c7);

    max_c0 = _mm256_max_epu32(max_c0, max_c2);
    max_c4 = _mm256_max_epu32(max_c4, max_c6);

    max_c0 = _mm256_max_epu32(max_c0, max_c4);

    uint32_t ret = 0;
    for(int i = 0; i < 32/sizeof(uint32_t); i++){
        if (max_c0.m256i_u32[i] > ret) {
            ret = max_c0.m256i_u32[i];
        }
    }
    if (size > 0) {
        printf("size=%d\n", size);
#if 0
        for(int i = 0; i < size; i++){
            uint32_t v = *((uint32_t*)src)+i
            if (v > ret) {
                ret = v;
            }
        }
#endif
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
    //LARGE_INTEGER PerformanceCount01;
    //QueryPerformanceCounter(&PerformanceCount01);
    //printf("%s:", __func__);
    //::printf("max=%x,%d\n", ret, (uint32_t)(PerformanceCount01.QuadPart - PerformanceCount00.QuadPart));;

    return ;
}
__declspec(noinline)  void convsqrt(uint32_t* src, int shiftval, uint32_t len, uint8_t* dst)
{
    //cout << len << endl;

    for (uint32_t i = 0; i < len; i++) {
        uint32_t v0 = src[i];
        float v1 = (v0 >> shiftval)&0xffff;
        float vsqr = sqrt(v1);
        uint8_t vout = vsqr;
        dst[i] = vout;
        //cout << v0 << ","<<v1 << ","<<(int)vout << ","<<vsqr;  
        //cout << endl;
    }
}
#define WIDTH 2432
#define HEIGHT 140
#define NLINE 640
int main()
{
    std::cout << "Hello World!\n" << RAND_MAX;
    uint32_t* src = (uint32_t*)malloc((WIDTH * HEIGHT * NLINE) * sizeof(uint32_t));
    uint8_t* dst = (uint8_t*)malloc((WIDTH * HEIGHT * NLINE) * sizeof(uint32_t));
    if (dst == nullptr || src == nullptr) {
        cout << "error src or dst \n";
        return 1;
    }
    for (int i = 0; i < WIDTH * HEIGHT * NLINE; i++) {
        //src[i] = (uint32_t)rand();
        src[i] = (uint32_t)(4096 << 16);
    }
    uint32_t len = WIDTH * HEIGHT;
#define PITCH (WIDTH*HEIGHT)
    for (volatile int ii = 0; ii < 10; ii++) {
        int shift = rand()& 0xf;
        //int shift = 16;
        auto t0 = chrono::system_clock::now();
#pragma omp parallel for 
        for (int k = 0; k < NLINE; k++) {
            convsqrt(&src[k*PITCH], shift, len, &dst[k*PITCH]);
        }
        auto t1 = chrono::system_clock::now();
        std::chrono::duration<double> diff = t1 - t0;
        cout << diff.count() << endl;
        convsqrtavx2(src, shift, len, dst);
    }

    _exit(0);

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

#include <windows.h>
#undef max
#undef min
#include <intrin.h>
#include <immintrin.h>
//#include <ammintrin.h>
//#include <emmintrin.h>
//#include <xmmintrin.h>

#include <iostream>
#include <stdlib.h>   // For _MAX_PATH definition
#include <stdio.h>
#include <malloc.h>
#include <valarray>
#include <algorithm>
#include <execution>
#include <omp.h>

#if 1
#define PITCH 640
#define WIDTH 2432
#else
#define PITCH 2240
#define WIDTH 2240 
#endif

#define HEIGHT 220
//#define HEIGHT 110
//#define HEIGHT 221
//#define PITCH 64
//#define WIDTH 48
//#define HEIGHT 220
//#define NUMLOOP 10000
#define NUMLOOP 4

using namespace std;

void bin2_64_avx2(uint8_t* src0, uint8_t* src1, int16_t* dst)
{
    //load line0
    __m128i line00w = _mm_loadu_si128(((const __m128i*)src0) + 0); //read upper 16byte 8words 
    __m128i line01w = _mm_loadu_si128(((const __m128i*)src0) + 1); //read upper 16byte 8words 
        //convert bytes->words
    __m256i line00 = _mm256_cvtepu8_epi16(line00w); //line0 +0 byte -> short 16words
    __m256i line01 = _mm256_cvtepu8_epi16(line01w);//line0 +128 upper1 byte -> short 16words
    //load line1
    __m128i line10w = _mm_loadu_si128(((const __m128i*)src1) + 0); //read lower 16byte 8words 
    __m128i line11w = _mm_loadu_si128(((const __m128i*)src1) + 1); //read lower 16byte 8words 
        //convert bytes->words
    __m256i line10 = _mm256_cvtepu8_epi16(line10w); //line1 +0 byte -> short 16 words
    __m256i line11 = _mm256_cvtepu8_epi16(line11w);// line1 +128 -> short 16 words

    //sum of line0 and line1
    line00 = _mm256_add_epi16(line00, line10); //first 16 words
    line01 = _mm256_add_epi16(line01, line11); //first 16 words
#if 0
    puts("line sum ");
    for(int i = 0; i < 16; i++){
        printf("%04x ", line00.m256i_u16[i]);
    }
    for(int i = 0; i < 16; i++){
        printf("%04x ", line01.m256i_u16[i]);
    }
    puts("");
#endif

    //store
    _mm256_storeu_si256(((__m256i*)dst) + 0, line00);//first 16words
    _mm256_storeu_si256(((__m256i*)dst) + 1, line01);//second 16words
}
void bin2_avx2(uint8_t* src0, int16_t* dst)
{
#pragma omp parallel for num_threads(8)

    for (int line = 0; line < HEIGHT; line += 2 ) {
        int n64 = WIDTH / 64;
        for (int i64 = 0; i64 < n64; i64++) {
//void bin2_64_avx2(uint8_t* src0, uint8_t* src1, int32_t* dst)
            bin2_64_avx2 
            (
				&src0 [(line + 0) * WIDTH + i64*64],
				&src0 [(line + 1) * WIDTH + i64*64],
                &dst  [(line ) * WIDTH + i64*64] 
                );
        }
    }
}


void bin2satsub_line_avx2_64(uint8_t* src0, uint8_t* src1, int16_t *sub, int32_t* dst, float *gain)
{
    __m128i line00w = _mm_loadu_si128(((const __m128i*)src0) + 0); //read upper 16byte 8words 
    __m128i line01w = _mm_loadu_si128(((const __m128i*)src0) + 1); //read upper 16byte 8words 
    __m128i line10w = _mm_loadu_si128(((const __m128i*)src0) + 2); //read upper 16byte 8words 
    __m128i line11w = _mm_loadu_si128(((const __m128i*)src0) + 3); //read upper 16byte 8words 
    __m256i line00 = _mm256_cvtepu8_epi16(line00w); //line0 +0 byte -> short
    __m256i line01 = _mm256_cvtepu8_epi16(line01w);//line0 +128 upper1 byte -> short

    __m256i line10 = _mm256_cvtepu8_epi16(line10w); //line1 +0 byte -> short
    __m256i line11 = _mm256_cvtepu8_epi16(line11w);// line1 +128 -> short

    //__m256i line0, line1;// , c2, c3, c4, c5, c6, c7;
	//line0 = _mm256_loadu_si256(((const __m256i*)src0) + 0); //read upper 32byte 
	//line1 = _mm256_loadu_si256(((const __m256i*)src1) + 0); //read lower 32byte 
#if 0
    puts("load line0 32byte");
    for(int i = 0; i < 32; i++){
        printf("%02x ", line0.m256i_u8[i]);
    }
    puts("");
    puts("load line1 32byte");
    for(int i = 0; i < 32; i++){
        printf("%02x ", line1.m256i_u8[i]);
    }
    puts("");
#endif
    //byte -> short
#if 0
    __m128i line00w = _mm256_extracti128_si256(line0, 0); // upper split 2 128
    __m128i line01w = _mm256_extracti128_si256(line0, 1);
    __m256i line00 = _mm256_cvtepu8_epi16(line00w); //upper0 byte -> short
    __m256i line01 = _mm256_cvtepu8_epi16(line01w);//upper1 byte -> short

    __m128i line10w = _mm256_extracti128_si256(line1, 0); // lower0 split 2 128
    __m128i line11w = _mm256_extracti128_si256(line1, 1);
    __m256i line10 = _mm256_cvtepu8_epi16(line10w); //lower0 byte -> short
    __m256i line11 = _mm256_cvtepu8_epi16(line11w);//lower1 byte -> short
#endif
#if 0
    puts("line0 -> int16");
    for(int i = 0; i < 16; i++){
        printf("%04x ", line00.m256i_u16[i]);
    }
    for(int i = 0; i < 16; i++){
        printf("%04x ", line01.m256i_u16[i]);
    }
    puts("");
    puts("line1 -> int16");
    for(int i = 0; i < 16; i++){
        printf("%04x ", line10.m256i_u16[i]);
    }
    for(int i = 0; i < 16; i++){
        printf("%04x ", line11.m256i_u16[i]);
    }
    puts("");
#endif
//sum
    line00 = _mm256_add_epi16(line00, line10);
    line01 = _mm256_add_epi16(line01, line11);
#if 0
    puts("line sum ");
    for(int i = 0; i < 16; i++){
        printf("%04x ", line00.m256i_u16[i]);
    }
    for(int i = 0; i < 16; i++){
        printf("%04x ", line01.m256i_u16[i]);
    }
    puts("");
#endif

    //line 0 shift << 3
    line00 = _mm256_slli_epi16(line00, 3);
    line01 = _mm256_slli_epi16(line01, 3);
#if 0
    puts("line0 shift 3 input << 3");
    for(int i = 0; i < 16; i++){
        printf("%04x ", line00.m256i_u16[i]);
    }
    for(int i = 0; i < 16; i++){
        printf("%04x ", line01.m256i_u16[i]);
    }
    puts("");
#endif
    //load sub
	__m256i sub0 = _mm256_loadu_si256(((const __m256i*)sub) + 0); //read sub 16words
	__m256i sub1 = _mm256_loadu_si256(((const __m256i*)sub) + 1); //read sub 16words
#if 0
    puts("load sub int16 32");
    for(int i = 0; i < 16; i++){
        printf("%04x ", sub0.m256i_u16[i]);
    }
    for(int i = 0; i < 16; i++){
        printf("%04x ", sub1.m256i_u16[i]);
    }
    puts("");
#endif

    //subtract
    __m256i sub_result0 = _mm256_subs_epu16(line00, sub0);
    __m256i sub_result1 = _mm256_subs_epu16(line01, sub1);
#if 0
    puts("after subtract 32 -> sub_result");
    for(int i = 0; i < 16; i++){
        printf("%04x ", sub_result0.m256i_u16[i]);
    }
    for(int i = 0; i < 16; i++){
        printf("%04x ", sub_result1.m256i_u16[i]);
    }
    puts("");
#endif
    //short -> int32
        //split
    __m128i e00 = _mm256_extracti128_si256(sub_result0, 0); // split 2 128
    __m128i e01 = _mm256_extracti128_si256(sub_result0, 1); // split 2 128

    __m128i e10 = _mm256_extracti128_si256(sub_result1, 0);
    __m128i e11 = _mm256_extracti128_si256(sub_result1, 1);

    __m256i i32_00 = _mm256_cvtepu16_epi32(e00);
    __m256i i32_01 = _mm256_cvtepu16_epi32(e01);
    __m256i i32_10 = _mm256_cvtepu16_epi32(e10);
    __m256i i32_11 = _mm256_cvtepu16_epi32(e11);
#if 0
    puts("convert result -> i32");
    for(int i = 0; i < 8; i++){
        printf("%08x ", i32_00.m256i_u32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%08x ", i32_01.m256i_u32[i]);
    }
    //puts("");
    for(int i = 0; i < 8; i++){
        printf("%08x ", i32_10.m256i_u32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%08x ", i32_11.m256i_u32[i]);
    }
    puts("");
#endif
    //int32->float
    __m256 f32_00 = _mm256_cvtepi32_ps(i32_00);//32byte 8 float
    __m256 f32_01 = _mm256_cvtepi32_ps(i32_01);//32byte 8 float
    __m256 f32_10 = _mm256_cvtepi32_ps(i32_10);//32byte 8 float
    __m256 f32_11 = _mm256_cvtepi32_ps(i32_11);//32byte 8 float
#if 0
    puts("convert result -> float32");
    for(int i = 0; i < 8; i++){
        printf("%f ", f32_00.m256_f32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%f ", f32_01.m256_f32[i]);
    }
    //puts("");
    for(int i = 0; i < 8; i++){
        printf("%f ", f32_10.m256_f32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%f ", f32_11.m256_f32[i]);
    }
    puts("");
#endif

    //load gain
    __m256 g00 = _mm256_load_ps(gain + 0);//load 8 word
    __m256 g01 = _mm256_load_ps(gain + 1);//load 8
    __m256 g10 = _mm256_load_ps(gain + 2);//load 8
    __m256 g11 = _mm256_load_ps(gain + 3);//load 8
#if 0
    puts("load gain 32 ");
    for(int i = 0; i < 8; i++){
        printf("%.2f ", g00.m256_f32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%.2f ", g01.m256_f32[i]);
    }
    //puts("");
    for(int i = 0; i < 8; i++){
        printf("%.2f ", g10.m256_f32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%.2f ", g11.m256_f32[i]);
    }
    puts("");
#endif
    //mull
    f32_00 = _mm256_mul_ps(f32_00, g00);
    f32_01 = _mm256_mul_ps(f32_01, g01);
    f32_10 = _mm256_mul_ps(f32_10, g10);
    f32_11 = _mm256_mul_ps(f32_11, g11);
#if 0
    puts("mul result*gain 32 ");
    for(int i = 0; i < 8; i++){
        printf("%.3f ", f32_00.m256_f32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%.3f ", f32_01.m256_f32[i]);
    }
    //puts("");
    for(int i = 0; i < 8; i++){
        printf("%.3f ", f32_10.m256_f32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%.3f ", f32_11.m256_f32[i]);
    }
    puts("");

#endif
    //float -> int32
    __m256i o00 = _mm256_cvtps_epi32(f32_00);
    __m256i o01 = _mm256_cvtps_epi32(f32_01);
    __m256i o10 = _mm256_cvtps_epi32(f32_10);
    __m256i o11 = _mm256_cvtps_epi32(f32_11);
#if 0
    puts("mul float -> i32");
    for(int i = 0; i < 8; i++){
        printf("%08x ", o00.m256i_u32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%08x ", o01.m256i_u32[i]);
    }
    //puts("");
    for(int i = 0; i < 8; i++){
        printf("%08x ", o10.m256i_u32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%08x ", o11.m256i_u32[i]);
    }
    puts("");
    for(int i = 0; i < 8; i++){
        printf("%6d ", o00.m256i_u32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%6d ", o01.m256i_u32[i]);
    }
    //puts("");
    for(int i = 0; i < 8; i++){
        printf("%6d ", o10.m256i_u32[i]);
    }
    for(int i = 0; i < 8; i++){
        printf("%6d ", o11.m256i_u32[i]);
    }
    puts("");
#endif
    //store
    _mm256_storeu_si256(((__m256i*)dst) + 0, o00);
    _mm256_storeu_si256(((__m256i*)dst) + 1, o01);
    _mm256_storeu_si256(((__m256i*)dst) + 2, o10);
    _mm256_storeu_si256(((__m256i*)dst) + 3, o11);
#if 0
    puts("dst=");
    for (int i = 0; i < 32; i++) {
        printf("%08x ", dst[i]);
        if ((i & 0xf) == 0xf) {
            puts("");
        }
    }
    puts("");
#endif
}
void bin2satsub_avx2(uint8_t* src0, int16_t* sub, int32_t* dst, float *gain)
{
#pragma omp parallel for num_threads(8)

    for (int line = 0; line < HEIGHT/2; line++) {
        int n64 = WIDTH / 64;
        for (int i64 = 0; i64 < n64; i64++) {
#if 0
            bin2satsub_line64(line0p, line1p, 
                &sub[(line / 2) * WIDTH + i64], 
                &dst[(line / 2) * WIDTH + i64], 
                &gain[(line/2)*WIDTH+i64]);
#else
            //bin2satsub_line_avx2_64(line0p, line1p, 
                //int16_t *sub, int32_t* dst, float *gain)
            bin2satsub_line_avx2_64
            (
				&src0 [(line*2 + 0) * WIDTH + i64*64],
				&src0 [(line*2 + 1) * WIDTH + i64*64],
                &sub  [(line ) * WIDTH + i64*64], 
                &dst  [(line ) * WIDTH + i64*64] ,
                &gain [(line ) * WIDTH + i64*64]);
#endif
        }
    }
}
void bin2satsub_gain_line64(uint8_t* src0, uint8_t* src1, int16_t *sub, int16_t* dst, float *gain)
{
    for (int j = 0; j < 64; j++) {

		uint8_t c0 = *src0++;
		uint8_t c1 = *src1++;
		int16_t sum = (c0 + c1) << 3;
		int16_t s = *sub++;
		int16_t out = sum - s;
		if (out < 0) {
		    out = 0;
		}
        float g = *gain++;
        out *= g;
        *dst++ = out;
    }

}
void bin2satsub_line64(uint8_t* src0, uint8_t* src1, int16_t *sub, int16_t* dst)
{
    for (int j = 0; j < 64; j++) {

		uint8_t c0 = *src0++;
		uint8_t c1 = *src1++;
		int16_t sum = (c0 + c1) << 3;
		int16_t s = *sub++;
		int16_t out = sum - s;
		if (out < 0) {
		    out = 0;
		}
        //float g = *gain++;
        //out *= g;
        *dst++ = out;
    }

}
void bin2satsub(uint8_t* src0, int16_t* sub, int16_t* dst, float *gain)
{
//#pragma omp parallel for num_threads(4)

    for (int line = 0; line < HEIGHT/2; line++) {
        int n64 = WIDTH / 64;
        uint8_t* line0p = &src0 [ (line*2 +0) * WIDTH];
        uint8_t* line1p = &src0 [ (line*2 +1) * WIDTH];
        for (int i64 = 0; i64 < n64; i64++) {
#if 0
            bin2satsub_line64(line0p, line1p, 
                &sub[(line / 2) * WIDTH + i64], 
                &dst[(line / 2) * WIDTH + i64], 
                &gain[(line/2)*WIDTH+i64]);
#else
            bin2satsub_line64(line0p, line1p, 
                &sub[(line ) * WIDTH + i64], 
                &dst[(line ) * WIDTH + i64] 
                );
#endif
        }
    }
}
void satsub64(uint16_t* src0, uint16_t* src1, uint16_t* dst)
{
    for (int j = 0; j < 64; j++) {

       uint16_t c0 = *src0++;
       c0 <<= 3;
       uint16_t c1 = *src1++;
       int16_t diff = (int16_t)c0 - (int16_t)c1;
       if (diff < 0) {
        diff = 0;
       }
        *dst++ = diff;
     }

}
void satsub256(uint16_t* src0, uint16_t* src1, uint16_t* dst)
{
    for (int j = 0; j < 256; j++) {
		uint8_t c0 = src0[j];
        c0 <<= 3;
		uint16_t c1 = src1[j];
		int16_t diff = (int16_t)c0 - (int16_t)c1;
		if (diff < 0) {
			diff = 0;
		}
        *dst++ = diff;
     }

}
void satsub512(uint16_t* binout, uint16_t*sub , uint16_t* dst)
{
    for (int j = 0; j < 512; j++) {

       uint16_t c0 = *binout++;
       int16_t s0 = *sub++;
       int16_t diff = (int16_t)c0 - (int16_t)s0;
       if (diff < 0) {
        diff = 0;
       }
        *dst++ = diff;
     }

}
void satsub(uint16_t* binout, uint16_t* sub, uint16_t* dst)
{
    //LARGE_INTEGER PerformanceCount00;
    //QueryPerformanceCounter(&PerformanceCount00);
    uint8_t maxv = 0;
    int size = WIDTH * HEIGHT / 2;
    int n512 = size / 512;
#pragma omp parallel for //num_threads(4)
    for (int i = 0; i < n512; i++){
        satsub512(&binout[i * 512], &sub[i * 512], &dst[i * 512]);
    }
    int n256 = (size%512) / 256;
    uint16_t* src0_256 = &binout[n512 * 512];
    uint16_t* sub_256 = &sub[n512 * 512];
    uint16_t* dst_256 = &dst[n512 * 512];
    for (int i = 0; i < n256; i++) {
        satsub256(&src0_256[i * 256], &sub_256[i * 512], &dst_256[i * 512]);
    }
    int n64 = (size % 256) / 64;
    uint16_t* src0_64 = &src0_256[n256 * 256 + n512*512];
    uint16_t* src1_64 = &sub_256[n256 * 256 + n512*512];
    uint16_t* dst_64 = &dst_256[n256 * 256 + n512*512];
    for (int i = 0; i < n64; i++) {
        satsub64(&src0_64[i * 256], &src1_64[i * 64], &dst_64[i * 64]);
    }
    //LARGE_INTEGER PerformanceCount01;
    //QueryPerformanceCounter(&PerformanceCount01);
    //double diff = (PerformanceCount01.QuadPart - PerformanceCount00.QuadPart)/10.0;
    //::printf("satsub= %f,maxv=%d\n", diff, (int)maxv);
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
        //satsub(src0, src1, dst, size);
    }
    LARGE_INTEGER PerformanceCount01;
    QueryPerformanceCounter(&PerformanceCount01);
    float diff=(PerformanceCount01.QuadPart - PerformanceCount00.QuadPart)/10.0;;
    ::printf("avxsatsub= %.1f\n", diff);

}
uint32_t  avx2max32(uint32_t* src0, int size)
{
    uint8_t* src = (uint8_t*)src0;
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
        max_c0 = _mm256_max_epu32(max_c0, c0);
        max_c1 = _mm256_max_epu32(max_c1, c1);
        max_c2 = _mm256_max_epu32(max_c2, c2);
        max_c3 = _mm256_max_epu32(max_c3, c3);
        max_c4 = _mm256_max_epu32(max_c4, c4);
        max_c5 = _mm256_max_epu32(max_c5, c5);
        max_c6 = _mm256_max_epu32(max_c6, c6);
        max_c7 = _mm256_max_epu32(max_c7, c7);
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
    LARGE_INTEGER PerformanceCount01;
    QueryPerformanceCounter(&PerformanceCount01);
    printf("%s:", __func__);
    ::printf("max=%x,%d\n", ret, (uint32_t)(PerformanceCount01.QuadPart - PerformanceCount00.QuadPart));;

    return ret;

}
uint8_t  avx2max(uint8_t* src, int size)
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
    printf("%s:", __func__);
    ::printf("max=%d,%d\n", ret, (uint32_t)(PerformanceCount01.QuadPart - PerformanceCount00.QuadPart));;

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
void bin2_64(uint8_t* line0p, uint8_t * line1p, uint16_t* outp)
{
	for (int i = 0; i < 64; i++) {
		uint16_t c0 = (uint16_t)line0p[i];
		uint16_t c1 = (uint16_t)line1p[i];;
		uint16_t sum = c0 + c1;
		outp[i] = sum;
	}

}
void bin2_512(uint8_t* line0p, uint8_t * line1p, uint16_t* outp)
{
	for (int i = 0; i < 512; i++) {
		uint16_t c0 = (uint16_t)line0p[i];
		uint16_t c1 = (uint16_t)line1p[i];;
		uint16_t sum = c0 + c1;
		outp[i] = sum;
	}

}
void gain_64(uint16_t* subout, float* gain, uint32_t* out)
{
    for (int j = 0; j < 8; j++) {
        float* gainp = &gain[j * 8];
        uint16_t* inp = &subout[j * 8];
        uint32_t* outp = &out[j*8];
        //load gain
        __m256 g256 = _mm256_load_ps(gainp);//8 float 32bytes
        //load in
        __m128i in0 = _mm_loadu_si128(((const __m128i*)inp)+0);//read 8words 16bytes
        //converet to u32
		__m256i d00 = _mm256_cvtepu8_epi16(in0); // 8 dwords  32 bytes
        //to float
		__m256 if32_00 = _mm256_cvtepi32_ps(d00);//32byte 8 float 32bytes
        //mul
		__m256 of32_00 = _mm256_mul_ps(if32_00, g256);
        //to uint32
        __m256i o00 = _mm256_cvtps_epi32(of32_00);
		_mm256_storeu_si256 (((__m256i*)outp) + 0, o00);
    }
}
void gain_512(uint16_t* subout, float* gain, uint32_t* out)
{
	for (int i = 0; i < 512; i++) {
		uint16_t c0 = (uint16_t)subout[i];
        float c0f = (float)c0;
        float g = gain[i];
        float mul = c0f * g;
        uint32_t u32 = mul;
		out[i] = u32;
	}
}
void gain_mul(uint16_t* subout, float* gain, uint32_t *out)
{
#pragma omp parallel for num_threads(8)
    for (int h = 0; h < HEIGHT/2 ; h++ ) {
#if 0
        int n512 = WIDTH / 512;
        for (int i512 = 0; i512 < n512; i512++) {
            uint16_t* linep = &subout[h  * WIDTH + i512 * 512];
            uint32_t* outp = &out[h * WIDTH + i512 * 512];
            float* gainp = &gain[h * WIDTH + i512 * 512];
            gain_512(linep, gainp, outp);
        }
        int n64 = (WIDTH %512)/ 64;
#else
        int n64 = (WIDTH )/ 64;
        for (int i64 = 0; i64 < n64; i64++) {
            uint16_t* linep = &subout[h * WIDTH + i64*64];
            uint32_t* outp = &out[h * WIDTH + i64 * 64];
            float* gainp = &gain[h * WIDTH + i64 * 64];
            gain_64(linep, gainp, outp);
        }
#endif
    }
}
void bin2(uint8_t* src, uint16_t *out)
{
//#pragma omp parallel for num_threads(8)

    for (int h = 0; h < HEIGHT/2 ; h++ ) {
        int n512 = WIDTH / 512;
        for (int i512 = 0; i512 < n512; i512++) {
            uint8_t* line0p = &src[(h * 2 + 0) * WIDTH + i512 * 512];
            uint8_t* line1p = &src[(h * 2 + 1) * WIDTH + i512 * 512];
            uint16_t* outp = &out[h * WIDTH + i512 * 512];
            bin2_512(line0p, line1p, outp);
        }
        //uint8_t *src64_0 = &src[(h * 2 + 0) * WIDTH + n512*512];
        //uint8_t *src64_1 = &src[(h * 2 + 1) * WIDTH + n512*512];
        int n64 = (WIDTH %512)/ 64;
        for (int i64 = 0; i64 < n64; i64++) {
            uint8_t* line0p = &src[(h * 2 + 0) * WIDTH + i64*64+n512*512];
            uint8_t* line1p = &src[(h * 2 + 1) * WIDTH + i64*64+n512*512];
            uint16_t* outp = &out[h*WIDTH+i64*64 ];
            bin2_64(line0p, line1p, outp);
        }
    }

}
int main()
{
    //std::cout << "Hello World!\n";
    valarray<double> elappsed(NUMLOOP);
    valarray<double> elappsed1(NUMLOOP);
    size_t memsz = (((WIDTH * HEIGHT) + 0xfff) & ~0xfff);
    size_t memszp = (((PITCH * HEIGHT*WIDTH*sizeof(uint32_t)) + 0xfff) & ~0xfff);
    size_t memszp8 = (((PITCH * HEIGHT*WIDTH*sizeof(uint8_t)) + 0xfff) & ~0xfff);
    size_t memszsub = (((PITCH * HEIGHT*WIDTH)*sizeof(int16_t) + 0xfff) & ~0xfff);
    size_t memszout = (((WIDTH * HEIGHT)*sizeof(int32_t) + 0xfff) & ~0xfff);
    size_t memszgain = (((WIDTH * HEIGHT*PITCH)*sizeof(float) + 0xfff) & ~0xfff);
    size_t memszbinout = (((PITCH * HEIGHT*WIDTH)*sizeof(int16_t) + 0xfff) & ~0xfff);
    //int16_t *dst =   (int16_t*)VirtualAlloc(nullptr, memszout, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    //uint8_t *memp = (uint8_t*)VirtualAlloc(nullptr, memszp, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    //uint8_t *src0  = (uint8_t*)VirtualAlloc(nullptr, memszp,   MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    uint32_t *src0  = (uint32_t*)VirtualAlloc(nullptr, memszp,   MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    uint8_t *src08  = (uint8_t*)VirtualAlloc(nullptr, memszp,   MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    uint16_t *binout  = (uint16_t*)VirtualAlloc(nullptr, memszbinout ,   MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    if (binout) {
        memset(binout, 0, memszbinout);
    }
    int16_t *sub   = (int16_t*)VirtualAlloc(nullptr, memszsub, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    int32_t *dst1  = (int32_t*)VirtualAlloc(nullptr, memszout, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    float *gain  = (float*)VirtualAlloc(nullptr, memszgain , MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    uint32_t *dst32  = (uint32_t*)VirtualAlloc(nullptr, memszout, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    if (dst32) {
        memset(dst32, 0, memszout);
    }

    if (src0) {
	    VirtualLock(src0, memszp);
        memset(src0, 0, memszp);
    }
    if (src08) {
	    VirtualLock(src08, memszp8);
        memset(src08, 0, memszp8);
    }
	//VirtualLock(memp, memsz);
	//VirtualLock(src1, memsz);
    if (sub) {
        VirtualLock(sub, memszsub);
        memset(sub, 0, memszsub);
    }
    if (dst1) {
        VirtualLock(dst1, memszout);
        memset(dst1, 0, memszout);
    }
    if (gain) {
        VirtualLock(gain, memszgain);
        memset(gain, 0, memszgain);
    }
	//VirtualLock(dst, memsz);
    if (!src0 || !src08|| !sub || !dst1 || !gain || !binout ||!dst32) {
        return 1;
    }

	size_t size = PITCH * HEIGHT; 
	//for (int i = 0; i < PITCH * HEIGHT; i+=1)  {
	for (int i = 0; i < memszp8; i+=1)  {
		src08[i] = (i+1)&0x7f;
        gain[i] = 1.+ (i%100)/100.0;
	}
    for (int i = 0; i < memszp/sizeof(uint32_t); i += 1) {
        src0[i] = ((i + 1)<<24)& 0x7fffffff;
    }
	for (int i = 0; i < PITCH * HEIGHT; i+=1)  {
		sub[i] = (i+1)&0x7f;
	}

    printf("->avx2max size=%zd\n", memszp8);
    static uint8_t maxval[NUMLOOP];
    for (int ii = 0; ii < NUMLOOP; ii++) {
        maxval[ii] = avx2max(src08, memszp8);
    }
    printf("\n done ->avx2max %zd\n", memszp);
    printf("->avx2max32 size=%zd\n", memszp);
    static uint32_t maxval32[NUMLOOP];
    for (int ii = 0; ii < NUMLOOP; ii++) {
        maxval32[ii] = avx2max32(src0, memszp);
    }
    printf("\n done ->avx2max32 %zd\n", memszp);


    for (int ii = 0; ii < NUMLOOP; ii++) {
#pragma omp parallel for num_threads(8)
        for (int i = 0; i < memsz; i += 64) {
            _mm_clflush(src0 + i);
        }
        for (int i = 0; i < WIDTH*HEIGHT/2; i += 64 / 2) {
            _mm_clflush(sub + i );
        }
        for (int i = 0; i < WIDTH*HEIGHT/2; i += 64 / 4) {
            _mm_clflush(dst1 + i );
        }
        for (int i = 0; i < WIDTH*HEIGHT/2; i += 64 / 4) {
            _mm_clflush(dst32 + i );
        }
        for (int i = 0; i < WIDTH*HEIGHT/2; i += 64 / 4) {
            _mm_clflush(gain + i);
        }
        for (int i = 0; i < WIDTH*HEIGHT/2; i += 64 / 2) {
            _mm_clflush(binout + i);
        }

        if (1) {
            LARGE_INTEGER PerformanceCount00;
            QueryPerformanceCounter(&PerformanceCount00);
            //int size = WIDTH*HEIGHT;
            //bin2satsub(src0, sub, dst1, gain);
            uint8_t* src1 = src08 + WIDTH;
            //void bin2satsub_line_avx2_64(uint8_t* src0, uint8_t* src1, int16_t *sub, int32_t* dst, float *gain)
            //void bin2_avx2(uint8_t* src0, int16_t* dst)
            bin2_avx2(src08, (int16_t*)dst1);
            LARGE_INTEGER PerformanceCount01;
            QueryPerformanceCounter(&PerformanceCount01);
            double diff = (PerformanceCount01.QuadPart - PerformanceCount00.QuadPart) / 10.0;
            ::printf("bin2_avx2= %f\n", diff);
        }


        if (1) {
            LARGE_INTEGER PerformanceCount00;
            QueryPerformanceCounter(&PerformanceCount00);
            //int size = WIDTH*HEIGHT;
            //bin2satsub(src0, sub, dst1, gain);
            uint8_t* src1 = src08 + WIDTH;
            //void bin2satsub_line_avx2_64(uint8_t* src0, uint8_t* src1, int16_t *sub, int32_t* dst, float *gain)
#if 1
            bin2satsub_avx2(src08, sub, dst1, gain);//200 - 300 us
#endif
#if 0
            bin2(src0, (uint16_t*)dst1); //70us - 100 us
#endif
            LARGE_INTEGER PerformanceCount01;
            QueryPerformanceCounter(&PerformanceCount01);
            double diff = (PerformanceCount01.QuadPart - PerformanceCount00.QuadPart) / 10.0;
            ::printf("bin2satsub_avx2= %f\n", diff);
        }


		if(1){
            LARGE_INTEGER PerformanceCount00;
			QueryPerformanceCounter(&PerformanceCount00);
			satsub((uint16_t*)binout, (uint16_t*)sub, (uint16_t*)dst1); //150us
			LARGE_INTEGER PerformanceCount01;
			QueryPerformanceCounter(&PerformanceCount01);
			double diff = (PerformanceCount01.QuadPart - PerformanceCount00.QuadPart)/10.0;
			::printf("satsub= %f\n", diff);
		}

		if(1){
            LARGE_INTEGER PerformanceCount00;
			QueryPerformanceCounter(&PerformanceCount00);
			gain_mul((uint16_t*)dst1, gain, dst32); //200us
	//void gain_mul(uint16_t* subout, float* gain, uint32_t *out)
			LARGE_INTEGER PerformanceCount01;
			QueryPerformanceCounter(&PerformanceCount01);
			double diff = (PerformanceCount01.QuadPart - PerformanceCount00.QuadPart)/10.0;
			::printf("gainmul= %f\n", diff);
		}
        if (0) {
            LARGE_INTEGER PerformanceCount00;
			QueryPerformanceCounter(&PerformanceCount00);
//void avxsatsub(uint8_t* src0, uint8_t* src1, uint8_t* dst, int size)
            //avxsatsub(src0, src0, dst1,size);
			LARGE_INTEGER PerformanceCount01;
			QueryPerformanceCounter(&PerformanceCount01);
			double diff = (PerformanceCount01.QuadPart - PerformanceCount00.QuadPart)/10.0;
			::printf("gainmul= %f\n", diff);
        }

    }
        //VirtualFree(memp, 0, MEM_RELEASE);
        VirtualFree(src0, 0, MEM_RELEASE);
        VirtualFree(sub, 0, MEM_RELEASE);
        VirtualFree(dst1, 0, MEM_RELEASE);
        //VirtualFree(dst, 0, MEM_RELEASE);
#if 0
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

#endif


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

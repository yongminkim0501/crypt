/*
 * Copyright(c) 2020-2024 All rights reserved by Heekuck Oh.
 * 이 프로그램은 한양대학교 ERICA 컴퓨터학부 학생을 위한 교육용으로 제작되었다.
 * 한양대학교 ERICA 학생이 아닌 자는 이 프로그램을 수정하거나 배포할 수 없다.
 * 프로그램을 수정할 경우 날짜, 학과, 학번, 이름, 수정 내용을 기록한다.
 */
#include "miller_rabin.h"

/*
 * mod_add() - computes a+b mod m
 * a와 b가 m보다 작다는 가정하에서 a+b >= m이면 결과에서 m을 빼줘야 하므로
 * 오버플로가 발생하지 않도록 a-(m-b)를 계산하고, 그렇지 않으면 그냥 a+b를 계산하면 된다.
 * a+b >= m을 검사하는 과정에서 오버플로가 발생할 수 있으므로 a >= m-b를 검사한다.
 */
uint64_t mod_add(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t result = 0;
    if (a>=m-b){
        result = a -(m-b);
    }else{
        result = a+b;
    }
    return result % m;

}

/*
 * mod_sub() - computes a-b mod m
 * 만일 a < b이면 결과가 음수가 되므로 m을 더해서 양수로 만든다.
 */
uint64_t mod_sub(uint64_t a, uint64_t b, uint64_t m)
{
    if (a < b) {
        return (a + m - b) % m; // 결과가 m보다 작도록 조정
    }
    return (a - b) % m;
}

/*
 * mod_mul() - computes a*b mod m
 * a*b에서 오버플로가 발생할 수 있기 때문에 덧셈을 사용하여 빠르게 계산할 수 있는
 * "double addition" 알고리즘을 사용한다. 그 알고리즘은 다음과 같다.
 *     r = 0;
 *     while (b > 0) {
 *         if (b & 1)
 *             r = mod_add(r, a, m);
 *         b = b >> 1;
 *         a = mod_add(a, a, m);
 *     }
 */
uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t r = 0;
    a = a % m; // a를 m으로 나눈 나머지로 초기화
    while (b > 0) {
        if (b & 1) {
            r = mod_add(r, a, m); // b의 최하위 비트가 1인 경우
        }
        b = b >> 1;              // b를 2로 나눈다.
        a = mod_add(a, a, m);   // a를 2배로 만든다.
    }
    return r;
}

/*
 * mod_pow() - computes a^b mod m
 * a^b에서 오버플로가 발생할 수 있기 때문에 곱셈을 사용하여 빠르게 계산할 수 있는
 * "square multiplication" 알고리즘을 사용한다. 그 알고리즘은 다음과 같다.
 *     r = 1;
 *     while (b > 0) {
 *         if (b & 1)
 *             r = mod_mul(r, a, m);
 *         b = b >> 1;
 *         a = mod_mul(a, a, m);
 *     }
 */
uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t r = 1;
    a = a % m; // a를 m으로 나눈 나머지로 초기화
    while (b > 0) {
        if (b & 1) {
            r = mod_mul(r, a, m); // b의 최하위 비트가 1인 경우
        }
        b = b >> 1;              // b를 2로 나눈다.
        a = mod_mul(a, a, m);   // a의 제곱을 계산한다.
    }
    return r;
}

/*
 * Miller-Rabin Primality Testing against small sets of bases
 *
 * if n < 2^64,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
 *
 * if n < 3,317,044,064,679,887,385,961,981,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, and 41.
 */
const uint64_t a[BASELEN] = {2,3,5,7,11,13,17,19,23,29,31,37};

/*
 * miller_rabin() - Miller-Rabin Primality Test (deterministic version)
 *
 * n > 3, an odd integer to be tested for primality
 * It returns PRIME if n is prime, COMPOSITE otherwise.
 */
int miller_rabin(uint64_t n)
{
    if (n < 2) return 0; // 2보다 작은 수는 소수가 아님
    if (n == 2 || n == 3) return 1; // 2와 3은 소수
    if (n % 2 == 0) return 0; // 짝수는 소수가 아님

    // n-1을 2의 제곱수로 분해
    uint64_t r = 0; // 2의 몇 제곱?
    uint64_t d = n - 1;
    while (d % 2 == 0) {
        d /= 2;
        r++;
    }

    // 기본 집합을 사용하여 n의 소수 여부를 판별
    for (int i = 0; i < BASELEN; i++) {
        uint64_t a_i = a[i];
        if (a_i > n - 2) break; // a_i가 n보다 크면 중지

        uint64_t x = mod_pow(a_i, d, n); // a^d mod n
        if (x == 1 || x == n - 1) continue;

        // 반복 검사
        int composite = 1;
        for (uint64_t j = 0; j < r - 1; j++) {
            x = mod_mul(x, x, n); // x^2 mod n
            if (x == n - 1) {
                composite = 0;
                break;
            }
        }
        if (composite) return 0; // 합성수
    }
    return 1; // 소수
}

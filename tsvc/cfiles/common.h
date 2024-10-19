#ifndef TSVC_COMMON_HDR
#define TSVC_COMMON_HDR

#define iterations 100000
#define LEN_1D 32000
#define LEN_2D 256

#include <sys/time.h>

struct args_t {
    double c1;
    double c2;
    const char *name;
    int status;
    int seed;
    struct timeval t1;
    struct timeval t2;
    void * __restrict__ arg_info;
};


#define TSVC_DUMP_TARGET stderr

#define DATA_PRINTF_MODIFIER "%lf "

#if 1
typedef double real_t;
#define ABS fabs
#else
typedef float real_t;
#define ABS fabsf
#endif

int dummy(real_t[LEN_1D], real_t[LEN_1D], real_t[LEN_1D], real_t[LEN_1D], real_t[LEN_1D], real_t[LEN_2D][LEN_2D], real_t[LEN_2D][LEN_2D], real_t[LEN_2D][LEN_2D], real_t);

void init(int** ip, real_t* s1, real_t* s2);

int initialise_arrays(const char* name);
real_t calc_checksum(const char * name);
void save_array(const char* name);

#endif

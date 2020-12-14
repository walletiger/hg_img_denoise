#pragma once 

#define  int8_t	signed char
#define uint8_t unsigned char

#define int16_t signed short
#define uint16_t unsigned short

#define int32_t signed	int
#define uint32_t unsigned int

#define uint64_t unsigned long long
#define int64_t signed long long


/*time*/
#ifdef __linux
#include <time.h>
#include <sys/time.h>
#include <sys/select.h>
#include <stdlib.h>

inline static uint32_t	HG_GetTimeMs()
{
	struct timespec t;
	
	clock_gettime(CLOCK_MONOTONIC_RAW , &t); 
	//clock_gettime(CLOCK_MONOTONIC , &t);
	return (t.tv_sec * 1000 + t.tv_nsec / 1000000);
}

inline static void * hg_aligned_malloc(int align, int size)
{
	return memalign(align, size);
}

inline static void hg_aligned_free(void *ptr)
{
	free(ptr);
}

#else

#include <sys/timeb.h>
#include <windows.h>
inline static uint32_t	HG_GetTimeMs()
{
    struct timeb t;

    ftime(&t);

    return (t.time * 1000 + t.millitm);
}

inline static void * hg_aligned_malloc(int align, int size)
{

    return _aligned_malloc(size, align);
}

inline static void hg_aligned_free(void *ptr)
{
    _aligned_free(ptr);
}


#endif



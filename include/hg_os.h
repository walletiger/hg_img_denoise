#pragma once

#ifdef __cplusplus
extern "C"
{
#endif


enum {
    HG_SOK              = 0,
    HG_ERR_FAILED       = -1,
    HG_ERR_INVALID_ARG  = -2,
    HG_ERR_NO_MEMORY    = -3
};

#define HG_PI (3.141592653)

/*data types*/
#ifdef __linux 
#include <stdint.h>
#else
#if (_MSC_VER < 1300)
typedef signed char       int8_t;
typedef signed short      int16_t;
typedef signed int        int32_t;
typedef unsigned char     uint8_t;
typedef unsigned short    uint16_t;
typedef unsigned int      uint32_t;
#else
typedef signed __int8     int8_t;
typedef signed __int16    int16_t;
typedef signed __int32    int32_t;
typedef unsigned __int8   uint8_t;
typedef unsigned __int16  uint16_t;
typedef unsigned __int32  uint32_t;
#endif
typedef signed __int64       int64_t;
typedef unsigned __int64     uint64_t;

#endif //#ifdef __linux

/*log*/
#ifndef HG_LOG
#define HG_LOG	printf
#endif 

#ifdef HG_DEBUG_LOG
#define HG_DEBUG HG_LOG
#else
#define HG_DEBUG(...)
#endif

/*mutex*/

#ifdef __linux 
#include <pthread.h>

typedef pthread_mutex_t		hg_mutex_t;

inline static int hg_mutex_init(hg_mutex_t * mutex)
{
	return pthread_mutex_init(mutex , NULL);
}

inline static void hg_mutex_lock(hg_mutex_t * mutex)
{
	pthread_mutex_lock(mutex);
}

inline static void hg_mutex_unlock(hg_mutex_t * mutex)
{
	pthread_mutex_unlock(mutex);
}

inline static void hg_mutex_destroy(hg_mutex_t * mutex)
{
	pthread_mutex_destroy(mutex);
}

#else 

#include <windows.h>
#include <process.h>

typedef HANDLE 	hg_mutex_t;

inline static int hg_mutex_init(hg_mutex_t * mutex)
{
	*mutex = CreateMutex(NULL , FALSE , NULL);
	return 0;		
}

inline static void hg_mutex_lock(hg_mutex_t * mutex)
{
	WaitForSingleObject(*mutex , INFINITE);
}

inline static void hg_mutex_unlock(hg_mutex_t *mutex)
{
	ReleaseMutex(*mutex);
}

inline static void hg_mutex_destroy(hg_mutex_t *mutex)
{
	CloseHandle(*mutex);
}

#endif //end #ifdef __linux


/*time*/
#ifdef __linux 
#include <time.h>
#include <sys/time.h>
#include <sys/select.h>

inline static uint32_t	hg_gettime_ms()
{
	struct timespec t;
	
	clock_gettime(CLOCK_MONOTONIC_RAW , &t); 
	//clock_gettime(CLOCK_MONOTONIC , &t);
	return (t.tv_sec * 1000 + t.tv_nsec / 1000000);
}

inline static void hg_wait_msec(int ms)
{
	struct timeval t;

	t.tv_sec  = ms / 1000;
	t.tv_usec = (ms % 1000) * 1000 ; 

	while(1){
		select(0 , NULL , NULL , NULL , &t);
		if (0 == t.tv_sec && 0 == t.tv_usec)
			break;
	}	
}

#else

#include <sys/time.h>

inline static uint32_t	hg_gettime_ms()
{
	struct timeb t;

	ftime(&t);

	return (t.time * 1000 + t.millitm);	
}

inline static void hg_wait_msec(int ms)
{
	SleepEx(ms , FALSE);	
}

#endif  //#ifdef __linux

#define HG_MIN(a , b)   ((a) < (b) ? (a) : (b))
#define HG_MAX(a , b)   ((a) > (b) ? (a) : (b))

inline static  uint8_t CLIP_UINT8(int32_t a)
{
    if( a&(~255) ) return (-a)>>31;
    else           return a;
}


/*all color spaces support*/
typedef enum HG_CS_FMT
{
	HG_CS_UNKNOWN = -1,
	HG_CS_RAW    = 0,
	HG_CS_GRAY8  = 0,
	HG_CS_RGB24,
	HG_CS_ARGB,
	HG_CS_YUYV,
	HG_CS_UYVY,
	HG_CS_I420,
	HG_CS_NV12 ,
	HG_CS_MAX
}HG_CS_FMT;

typedef struct HG_Point
{
	int x;
	int y;
}HG_Point;

typedef struct HG_Rect
{
	int start_x;
	int start_y;
	int cols;
	int rows;
}HG_Rect;

/*time*/
#ifdef __linux
#include <time.h>
#include <sys/time.h>
#include <sys/select.h>
#include <malloc.h>
#include <stdlib.h>

inline static void * hg_aligned_malloc(int align, int size)
{
    return memalign(align, size);
}

inline static void hg_aligned_free(void *ptr)
{
    free(ptr);
}

#else

#include <windows.h>

inline static void * hg_aligned_malloc(int align, int size)
{

    return _aligned_malloc(size, align);
}

inline static void hg_aligned_free(void *ptr)
{
    _aligned_free(ptr);
}


#endif

#ifdef __cplusplus
}
#endif


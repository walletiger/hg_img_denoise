#pragma once 


/*
    局部均方差滤波算法

 */

#define HG_LOCAL_VARIANCE_DENOISE_FILTER_R_MAX	(127)
#define HG_LOCAL_VARIANCE_DENOISE_FILTER_IMG_W_MAX	(4096)

typedef void * HG_FILTER_HANDLE;


#ifdef __cplusplus
extern "C"
{
#endif 

/**
 * @brief hg_local_variance_denoise_filter_create
 * @param w
 * @param h
 * @param r
 * @param lamda
 * @return
 */
HG_FILTER_HANDLE    hg_local_variance_denoise_filter_create(int w , int h , int r , float lamda /*0 ~ 1*/);/*r 最大15*/

/**
 * @brief hg_local_variance_denoise_filter_process
 * @param h
 * @param pIP
 * @param pitch_ip
 * @param pQ
 * @param pitch_q
 * @return
 */
int		hg_local_variance_denoise_filter_process(HG_FILTER_HANDLE h , unsigned char * pIP, int pitch_ip, unsigned char *pQ, int pitch_q);

/**
 * @brief hg_local_variance_denoise_filter_destroy
 * @param h
 */
void	hg_local_variance_denoise_filter_destroy(HG_FILTER_HANDLE h);


#ifdef __cplusplus
}
#endif 






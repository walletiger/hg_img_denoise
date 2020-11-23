#pragma once 

/*
	引导滤波降噪算法，原理来源于何凯明的 <Guided Image FIlter> , 	http://kaiminghe.com/eccv10/

	引导滤波是一种局部线程滤波模型，用于降噪时功能与双边滤波类似(保边降噪)，效果要比双波滤波略好。
	理论速度不受滤边窗口大小影响，大窗口时速度远远快于双边滤波。

	本算法实现针对降噪场景做了优化，小滤波半径(r <= 15) 在 Intel(R) Core(TM) i7-3632QM CPU @ 2.20GHz 主机系统ubuntu16.04 64bit，
	处理 720p yuv 一帧数据 单线程平均 6~7ms， cv::ximgproc::guidedFilter() 需要90ms.
	实时视频yuv处理时，仅处理y通道。
	
	
	假定输出图像Q与与引导图像I满足局部线性模型，
		Q(i) = I(i)   * a + b
		
	基于Q与输入P图像最小误差求得得到原作者公式:
		A = Cov(I, P) / (Delta(I, r) + lamda)
		B = Mean(P, r) - A .* Mean(I, r)
		Q = Mean(A, r) .* I + Mean(B, r)  

	
	用于图像平滑的引导图像I == 输入图像P,得到如下简化公式:
		Q = Mean(A, r) .* I + Mean(B, r)  
		A = Delta(I, r) / (Delta(I, r) + lamda)
		B = Mean(I, r) - A .* Mean(I, r)

	其中:
		(1)		Q 为输出图像(二维矩阵)
		(2)		Mean(A, r) 表示二维矩阵 A 以r 为半径的均值矩阵
		(3)   	.* 每个对应元素依次相乘
		(4)		lamda 	平滑因子，越大越平滑（lamda越大平坦区域 A 越小, B 越接近均值, Q 越接近B即均值）


	算法包含两种实现：
		算法1 先计算 Mean(I), Delta(I), Mean(A), Mean(B)，再算出Q
		算法2    Mean(I), Delta(I), Mean(A), Mean(B) 计算过程做了合并，省去了四个矩阵的转存。仅开辟 (2*r+1) * width 
		的 A、B ring buffer, 内存空间开销小，速度比较稳定

	算法速度:
	64bit windows , 64bit linux 要快于 32bit 平台 ，前者720p 一帧 i7约 6ms, 后者约 8ms
	大半径(r>15)用了64位乘法，某些32位平台速度会有减慢

 */

#define HG_GUIDED_DENOISE_FILTER_R_MAX	(127)
#define HG_GUIDED_DENOISE_FILTER_IMG_W_MAX	(4096)

typedef void * HG_FILTER_HANDLE;


#ifdef __cplusplus
extern "C"
{
#endif 

/**
 * @brief hg_guided_denoise_filter_create
 * @param w
 * @param h
 * @param r
 * @param lamda
 * @return
 */
HG_FILTER_HANDLE    hg_guided_denoise_filter_create(int w , int h , int r , float lamda /*0 ~ 1*/);/*r 最大15*/

/**
 * @brief hg_guided_denoise_filter_process
 * @param h
 * @param pIP
 * @param pitch_ip
 * @param pQ
 * @param pitch_q
 * @return
 */
int		hg_guided_denoise_filter_process(HG_FILTER_HANDLE h , unsigned char * pIP, int pitch_ip, unsigned char *pQ, int pitch_q);

/**
 * @brief hg_guided_denoise_filter_destroy
 * @param h
 */
void	hg_guided_denoise_filter_destroy(HG_FILTER_HANDLE h);

/**
 * @brief hg_guided_denoise_filter2_create
 * @param w
 * @param h
 * @param r
 * @param lamda
 * @return
 */
HG_FILTER_HANDLE    hg_guided_denoise_filter2_create(int w , int h , int r , float lamda);/*r 最大HG_GUIDED_DENOISE_FILTER_R_MAX=127*/


/**
 * @brief hg_guided_denoise_filter2_process
 * @param h
 * @param pIP input buf
 * @param pitch_ip
 * @param pQ output buf, pQ can't be the same buffer with pIP
 * @param pitch_q
 * @return
 */
int		hg_guided_denoise_filter2_process(HG_FILTER_HANDLE h , unsigned char * pIP, int pitch_ip, unsigned char *pQ, int pitch_q);

/**
 * @brief hg_guided_denoise_filter2_destroy
 * @param h
 */
void	hg_guided_denoise_filter2_destroy(HG_FILTER_HANDLE h);


#ifdef __cplusplus
}
#endif 






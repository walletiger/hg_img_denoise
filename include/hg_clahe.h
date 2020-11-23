#pragma once

#ifdef __cplusplus
extern "C"
{
#endif 

#define CLAHE_NBLOCK_W        (8)
#define CLAHE_NBLOCK_H        (8)
#define HIST_MAX_BIN         (256)


int hg_gray8_clahe(uint8_t *plane, int w, int h, int pitch, float hist_peak_crop);

#if 0
long HG_I420_CLAHE_Create(int w , int h , float hist_peak_crop , int i_use_autolevel , int i_fast_level /*0 ~ 9*/);

int HG_I420_CLAHE_Process(long h , struct HG_Mat_U8 *img);

int HG_I420_CLAHE_Destroy(long h);
#endif

#ifdef __cplusplus
}
#endif 




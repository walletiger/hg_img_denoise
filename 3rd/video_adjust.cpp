
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "video_adjust.h"

#ifndef NULL
#define NULL 0
#endif

#define EXTRA_WIDTH     (32)

#define MAX_WIDTH       (1920 + EXTRA_WIDTH)
#define MAX_HEIGHT      (1080)
#if 0
static AdjustParams g_AdjustParam = {1.0, 1.0, 2.0, 1.0, 0};
static unsigned char *g_Ydata = NULL;
static unsigned char *g_Udata = NULL;
static unsigned char *g_Vdata = NULL;
int VideoAdjustInit()
{
    g_Ydata = (unsigned char *)malloc(1920 * 1080);
    if (NULL == g_Ydata)
    {
        printf("VideoAdjustInit:malloc g_Ydata failed!\n");
        assert(0);
        return -1;
    }
    g_Udata = (unsigned char *)malloc(1920 * 1080 / 4);
    if (NULL == g_Ydata)
    {
        printf("VideoAdjustInit:malloc g_Ydata failed!\n");
        assert(0);
        return -1;
    }
    g_Vdata = (unsigned char *)malloc(1920 * 1080 / 4);
    if (NULL == g_Ydata)
    {
        printf("VideoAdjustInit:malloc g_Ydata failed!\n");
        assert(0);
        return -1;
    }
    return 0;
}
#endif
#if 0
int SetVideoAdjust(AdjustParams *param)
{
    /* legal to do */
    if (param)
    {
        g_AdjustParam = *param;
    }
    return 0;
}
#endif
void GetAdjustParamsDefault(AdjustParams *param)
{
    if (NULL == param)
    {
        printf("GetAdjustParamsDefault:NULL == param\n");
        assert(0);
        return;
    }
    param->bEnable      = false;
    param->f_contrast   = 1.0;
    param->f_brightness = 1.0;
    param->f_saturation = 1.0;
    param->f_gamma      = 1.0;
    param->i_hue        = 0;
}

int planar_sat_hue_clip_C(YUVFrame *pstFrameIn, YUVFrame *pstFrameOut, int i_sin, int i_cos,
                         int i_sat, int i_x, int i_y );
int planar_sat_hue_C(YUVFrame *pstFrameIn, YUVFrame *pstFrameOut, int i_sin, int i_cos,
                         int i_sat, int i_x, int i_y );

unsigned char clip_uint8_vlc(int a)
{
    if (a & (~255))
        return (-a)>>31;
    else
        return a;
}
#if 0
void *VideoAdjust(YUVFrame *pstFrameIn, YUVFrame *pstFrameOut)
{
    int pi_luma[256];
	int pi_gamma[256];

//    picture_t *p_outpic;
    unsigned char *p_in, *p_in_end, *p_line_end;
    unsigned char *p_out;

//    bool b_thres;
    double  f_hue;
    double  f_gamma;
    int i_cont, i_lum;
    int i_sat, i_sin, i_cos, i_x, i_y;
    int i;

    /* Get variables */
//    vlc_mutex_lock( &p_sys->lock );
    i_cont = (int)(g_AdjustParam.f_contrast * 255);
    i_lum = (int)((g_AdjustParam.f_brightness - 1.0) * 255);
    f_hue = (float)(g_AdjustParam.i_hue * PI / 180);
    i_sat = (int)(g_AdjustParam.f_saturation * 256);
	f_gamma = (int)(1.0 / g_AdjustParam.f_gamma);
//    vlc_mutex_unlock( &p_sys->lock );
#if 1
    pstFrameOut->apucData[0] = g_Ydata;
    pstFrameOut->apucData[1] = g_Udata;
    pstFrameOut->apucData[2] = g_Vdata;
    pstFrameOut->apucData[3] = NULL;

    pstFrameOut->iLineSize[0] = pstFrameIn->iLineSize[0];
    pstFrameOut->iLineSize[1] = pstFrameIn->iLineSize[1];
    pstFrameOut->iLineSize[2] = pstFrameIn->iLineSize[2];
    pstFrameOut->iLineSize[3] = pstFrameIn->iLineSize[3];

    pstFrameOut->iHeight = pstFrameIn->iHeight;
    pstFrameOut->iWidth = pstFrameIn->iWidth;
#else
    memcpy(pstFrameOut, pstFrameIn, sizeof(YUVFrame));
#endif
    /* Contrast is a fast but kludged function, so I put this gap to be cleaner :) */
    i_lum += 128 - i_cont / 2;

    /* Fill the gamma lookup table */
    for (i = 0 ; i < 256; i++)
    {
        pi_gamma[ i ] = clip_uint8_vlc(pow(i / 255.0, f_gamma) * 255.0);
    }

    /* Fill the luma lookup table */
    for (i = 0; i < 256; i++)
    {
        pi_luma[ i ] = pi_gamma[clip_uint8_vlc(i_lum + i_cont * i / 256)];
    }
    

    /*
     * Do the Y plane
     */

    p_in = pstFrameIn->apucData[0];
    p_in_end = p_in + pstFrameIn->iHeight * pstFrameIn->iLineSize[Y_PLANE] - 8;

    p_out = pstFrameOut->apucData[0];

    for( ; p_in < p_in_end; )
    {
        p_line_end = p_in + pstFrameIn->iWidth - 8;

        for( ; p_in < p_line_end; )
        {
            /* Do 8 pixels at a time */
            *p_out++ = pi_luma[ *p_in++ ]; *p_out++ = pi_luma[ *p_in++ ];
            *p_out++ = pi_luma[ *p_in++ ]; *p_out++ = pi_luma[ *p_in++ ];
            *p_out++ = pi_luma[ *p_in++ ]; *p_out++ = pi_luma[ *p_in++ ];
            *p_out++ = pi_luma[ *p_in++ ]; *p_out++ = pi_luma[ *p_in++ ];
        }

        p_line_end += 8;

        for( ; p_in < p_line_end; )
        {
            *p_out++ = pi_luma[ *p_in++ ];
        }
#if 1
        p_in += pstFrameIn->iLineSize[Y_PLANE]
              - pstFrameIn->iWidth;
        p_out += pstFrameOut->iLineSize[Y_PLANE]
              - pstFrameOut->iWidth;
#endif
    }

    /*
     * Do the U and V planes
     */

    i_sin = sin(f_hue) * 256;
    i_cos = cos(f_hue) * 256;

    i_x = (cos(f_hue) + sin(f_hue)) * 32768;
    i_y = (cos(f_hue) - sin(f_hue)) * 32768;

    if (i_sat > 256)
    {
        /* Currently no errors are implemented in the function, if any are added
         * check them here */
        planar_sat_hue_clip_C(pstFrameIn, pstFrameOut, i_sin, i_cos, i_sat,
                                        i_x, i_y );
    }
    else
    {
        /* Currently no errors are implemented in the function, if any are added
         * check them here */
        planar_sat_hue_C(pstFrameIn, pstFrameOut, i_sin, i_cos, i_sat,
                                        i_x, i_y );
    }

    return 0;
}
#else
void *VideoAdjust(YUVFrame *pstFrameIn, YUVFrame *pstFrameOut, AdjustParams *param)
{
    int pi_luma[256];
    int pi_gamma[256];

//    picture_t *p_outpic;
    unsigned char *p_in, *p_in_end, *p_line_end;
    unsigned char *p_out;

//    bool b_thres;
    double  f_hue;
    double  f_gamma;
    int i_cont, i_lum;
    int i_sat, i_sin, i_cos, i_x, i_y;
    int i;

    assert(NULL != pstFrameIn);
    assert(NULL != pstFrameOut);
    assert(NULL != param);

    /* Get variables */
    i_cont = (int)(param->f_contrast * 255);
    i_lum = (int)((param->f_brightness - 1.0) * 255);
    f_hue = (float)(param->i_hue * PI / 180);
    i_sat = (int)(param->f_saturation * 256);
	f_gamma = (float)(1.0 / param->f_gamma);
#if 0
    pstFrameOut->apucData[0] = g_Ydata;
    pstFrameOut->apucData[1] = g_Udata;
    pstFrameOut->apucData[2] = g_Vdata;
    pstFrameOut->apucData[3] = NULL;

    pstFrameOut->iLineSize[0] = pstFrameIn->iLineSize[0];
    pstFrameOut->iLineSize[1] = pstFrameIn->iLineSize[1];
    pstFrameOut->iLineSize[2] = pstFrameIn->iLineSize[2];
    pstFrameOut->iLineSize[3] = pstFrameIn->iLineSize[3];

    pstFrameOut->iHeight = pstFrameIn->iHeight;
    pstFrameOut->iWidth = pstFrameIn->iWidth;
#endif
    /* Contrast is a fast but kludged function, so I put this gap to be cleaner :) */
    i_lum += 128 - i_cont / 2;

    /* Fill the gamma lookup table */
    for (i = 0 ; i < 256; i++)
    {
        pi_gamma[ i ] = clip_uint8_vlc(pow(i / 255.0, f_gamma) * 255.0);
    }

    /* Fill the luma lookup table */
    for (i = 0; i < 256; i++)
    {
        pi_luma[ i ] = pi_gamma[clip_uint8_vlc(i_lum + i_cont * i / 256)];
    }
    

    /*
     * Do the Y plane
     */

    p_in = pstFrameIn->apucData[0];
    p_in_end = p_in + pstFrameIn->iHeight * pstFrameIn->iLineSize[Y_PLANE] - 8;

    p_out = pstFrameOut->apucData[0];

    for( ; p_in < p_in_end; )
    {
        p_line_end = p_in + pstFrameIn->iWidth - 8;

        for( ; p_in < p_line_end; )
        {
            /* Do 8 pixels at a time */
            *p_out++ = pi_luma[ *p_in++ ]; *p_out++ = pi_luma[ *p_in++ ];
            *p_out++ = pi_luma[ *p_in++ ]; *p_out++ = pi_luma[ *p_in++ ];
            *p_out++ = pi_luma[ *p_in++ ]; *p_out++ = pi_luma[ *p_in++ ];
            *p_out++ = pi_luma[ *p_in++ ]; *p_out++ = pi_luma[ *p_in++ ];
        }

        p_line_end += 8;

        for( ; p_in < p_line_end; )
        {
            *p_out++ = pi_luma[ *p_in++ ];
        }
#if 1
        p_in += pstFrameIn->iLineSize[Y_PLANE]
              - pstFrameIn->iWidth;
        p_out += pstFrameOut->iLineSize[Y_PLANE]
              - pstFrameOut->iWidth;
#endif
    }

    /*
     * Do the U and V planes
     */

    i_sin = sin(f_hue) * 256;
    i_cos = cos(f_hue) * 256;

    i_x = (cos(f_hue) + sin(f_hue)) * 32768;
    i_y = (cos(f_hue) - sin(f_hue)) * 32768;

    if (i_sat > 256)
    {
        /* Currently no errors are implemented in the function, if any are added
         * check them here */
        planar_sat_hue_clip_C(pstFrameIn, pstFrameOut, i_sin, i_cos, i_sat,
                                        i_x, i_y );
    }
    else
    {
        /* Currently no errors are implemented in the function, if any are added
         * check them here */
        planar_sat_hue_C(pstFrameIn, pstFrameOut, i_sin, i_cos, i_sat,
                                        i_x, i_y );
    }

    return 0;
}

#endif


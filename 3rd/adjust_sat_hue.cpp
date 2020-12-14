
#include "video_adjust.h"


#define PLANAR_WRITE_UV_CLIP() \
    i_u = *p_in++ ; i_v = *p_in_v++ ; \
    *p_out++ = clip_uint8_vlc( (( ((i_u * i_cos + i_v * i_sin - i_x) >> 8) \
                           * i_sat) >> 8) + 128); \
    *p_out_v++ = clip_uint8_vlc( (( ((i_v * i_cos - i_u * i_sin - i_y) >> 8) \
                           * i_sat) >> 8) + 128)

#define PLANAR_WRITE_UV() \
    i_u = *p_in++ ; i_v = *p_in_v++ ; \
    *p_out++ = (( ((i_u * i_cos + i_v * i_sin - i_x) >> 8) \
                       * i_sat) >> 8) + 128; \
    *p_out_v++ = (( ((i_v * i_cos - i_u * i_sin - i_y) >> 8) \
                       * i_sat) >> 8) + 128

#define PACKED_WRITE_UV_CLIP() \
    i_u = *p_in; p_in += 4; i_v = *p_in_v; p_in_v += 4; \
    *p_out = clip_uint8_vlc( (( ((i_u * i_cos + i_v * i_sin - i_x) >> 8) \
                           * i_sat) >> 8) + 128); \
    p_out += 4; \
    *p_out_v = clip_uint8_vlc( (( ((i_v * i_cos - i_u * i_sin - i_y) >> 8) \
                           * i_sat) >> 8) + 128); \
    p_out_v += 4

#define PACKED_WRITE_UV() \
    i_u = *p_in; p_in += 4; i_v = *p_in_v; p_in_v += 4; \
    *p_out = (( ((i_u * i_cos + i_v * i_sin - i_x) >> 8) \
                       * i_sat) >> 8) + 128; \
    p_out += 4; \
    *p_out_v = (( ((i_v * i_cos - i_u * i_sin - i_y) >> 8) \
                       * i_sat) >> 8) + 128; \
    p_out_v += 4

#define ADJUST_2_TIMES(x) x; x
#define ADJUST_4_TIMES(x) x; x; x; x
#define ADJUST_8_TIMES(x) x; x; x; x; x; x; x; x

/*****************************************************************************
 * Hue and saturation adjusting routines
 *****************************************************************************/

int planar_sat_hue_clip_C(YUVFrame *pstFrameIn, YUVFrame *pstFrameOut, int i_sin, int i_cos,
                         int i_sat, int i_x, int i_y )
{
    unsigned char *p_in, *p_in_v, *p_in_end, *p_line_end;
    unsigned char *p_out, *p_out_v;
    unsigned char i_u, i_v;
#if 0

    p_in = p_pic->p[U_PLANE].p_pixels;
    p_in_v = p_pic->p[V_PLANE].p_pixels;
    p_in_end = p_in + p_pic->p[U_PLANE].i_visible_lines
                      * p_pic->p[U_PLANE].i_pitch - 8;

    p_out = p_outpic->p[U_PLANE].p_pixels;
    p_out_v = p_outpic->p[V_PLANE].p_pixels;
#endif
    p_in = pstFrameIn->apucData[1];
    p_in_v = pstFrameIn->apucData[2];
    p_in_end = pstFrameIn->apucData[1] + pstFrameIn->iHeight / 2 * pstFrameIn->iLineSize[U_PLANE] - 8;

    p_out = pstFrameOut->apucData[1];
    p_out_v = pstFrameOut->apucData[2];


    for( ; p_in < p_in_end; )
    {
        p_line_end = p_in + pstFrameIn->iWidth / 2 - 8;

        for( ; p_in < p_line_end; )
        {
            /* Do 8 pixels at a time */
            ADJUST_8_TIMES(PLANAR_WRITE_UV_CLIP());
        }

        p_line_end += 8;

        for( ; p_in < p_line_end; )
        {
            PLANAR_WRITE_UV_CLIP();
        }
#if 1
        p_in += pstFrameIn->iLineSize[U_PLANE]
                - pstFrameIn->iWidth / 2;
        p_in_v += pstFrameIn->iLineSize[V_PLANE]
                - pstFrameIn->iWidth / 2;
        p_out += pstFrameOut->iLineSize[U_PLANE]
                - pstFrameOut->iWidth / 2;
        p_out_v += pstFrameOut->iLineSize[V_PLANE]
                - pstFrameOut->iWidth / 2;
#endif
    }

    return 0;
}

int planar_sat_hue_C(YUVFrame *pstFrameIn, YUVFrame *pstFrameOut, int i_sin, int i_cos,
                         int i_sat, int i_x, int i_y )
{
    unsigned char *p_in, *p_in_v, *p_in_end, *p_line_end;
    unsigned char *p_out, *p_out_v;
    unsigned char i_u, i_v;
#if 0
    p_in = p_pic->p[U_PLANE].p_pixels;
    p_in_v = p_pic->p[V_PLANE].p_pixels;
    p_in_end = p_in + p_pic->p[U_PLANE].i_visible_lines
                      * p_pic->p[U_PLANE].i_pitch - 8;

    p_out = p_outpic->p[U_PLANE].p_pixels;
    p_out_v = p_outpic->p[V_PLANE].p_pixels;
#endif
    p_in = pstFrameIn->apucData[1];
    p_in_v = pstFrameIn->apucData[2];
    p_in_end = pstFrameIn->apucData[1] + pstFrameIn->iHeight / 2 * pstFrameIn->iLineSize[U_PLANE] - 8;

    p_out = pstFrameOut->apucData[1];
    p_out_v = pstFrameOut->apucData[2];


    for( ;p_in < p_in_end; )
    {
        p_line_end = p_in + pstFrameIn->iWidth / 2 - 8;

        for( ; p_in < p_line_end; )
        {
            /* Do 8 pixels at a time */
            ADJUST_8_TIMES(PLANAR_WRITE_UV());
        }

        p_line_end += 8;

        for( ; p_in < p_line_end; )
        {
            PLANAR_WRITE_UV();
        }
#if 1
        p_in += pstFrameIn->iLineSize[U_PLANE]
                - pstFrameIn->iWidth / 2;
        p_in_v += pstFrameIn->iLineSize[V_PLANE]
                - pstFrameIn->iWidth / 2;
        p_out += pstFrameOut->iLineSize[U_PLANE]
                - pstFrameOut->iWidth / 2;
        p_out_v += pstFrameOut->iLineSize[V_PLANE]
                - pstFrameOut->iWidth / 2;
#endif
    }

    return 0;
}



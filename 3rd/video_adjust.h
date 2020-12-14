#ifndef __VIDEO_ADJUST__
#define __VIDEO_ADJUST__

#define PI 3.14159265358979323846

typedef struct
{
   unsigned char * apucData[4];
   int iLineSize[4];
   int iHeight;
   int iWidth;
}YUVFrame;

typedef struct
{
   float f_contrast;        /* range from (0.0, 2.0) defalut 1.0 */
   float f_brightness;      /* range from (0.0, 2.0) defalut 1.0 */
   float f_saturation;      /* range from (0.0, 3.0) defalut 1.0 */
   float f_gamma;           /* range from (0.1, 10.0) defalut 1.0 */
   int i_hue;               /* range from (0, 360) defalut 0 */
   bool bEnable;            /* enable flag */
}AdjustParams;


enum
{
    Y_PLANE = 0,
    U_PLANE,
    V_PLANE
};
#if 0
int VideoAdjustInit();

void *VideoAdjust(YUVFrame *pstFrameIn, YUVFrame *pstFrameOut);
#endif
void GetAdjustParamsDefault(AdjustParams *param);

void *VideoAdjust(YUVFrame *pstFrameIn, YUVFrame *pstFrameOut, AdjustParams *param);

unsigned char clip_uint8_vlc(int a);

#endif

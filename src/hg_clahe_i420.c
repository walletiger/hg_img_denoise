#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <malloc.h>
#include <hg_os.h>
#include <hg_mat.h>
#include <hg_clahe.h>
#include "hg_clahe_priv.h"

typedef struct HG_I420_CLAHE_Context
{
    int         w;
    int         h;
    int         block_w;
    int         block_h;
    int         block_size;
    float       hist_peak_crop;

    uint32_t    bin_learn[CLAHE_NBLOCK_W][CLAHE_NBLOCK_H][HIST_MAX_BIN];
    uint32_t    bin_cur[CLAHE_NBLOCK_W][CLAHE_NBLOCK_H][HIST_MAX_BIN];
    uint32_t    bin_crop[CLAHE_NBLOCK_W][CLAHE_NBLOCK_H][HIST_MAX_BIN];
    uint8_t     pix_map[CLAHE_NBLOCK_W][CLAHE_NBLOCK_H][HIST_MAX_BIN];

    uint8_t     *buf_rgb;
    HG_Mat_U8   img_rgb;
    uint32_t    n_frames;
    int         i_use_autolevel;
    uint32_t    i_fast_level;
}HG_I420_CLAHE_Context;

long HG_I420_CLAHE_Create(int w , int h , float hist_peak_crop , int i_use_autolevel , int i_fast_level)
{
    HG_I420_CLAHE_Context *ctx;

    if (w % CLAHE_NBLOCK_W != 0 || h % CLAHE_NBLOCK_H != 0)
        return 0;

    ctx = (HG_I420_CLAHE_Context *)calloc(1 , sizeof(HG_I420_CLAHE_Context));
    assert(ctx != 0);

    ctx->w = w;
    ctx->h = h;
    ctx->block_w = w / CLAHE_NBLOCK_W;
    ctx->block_h = h / CLAHE_NBLOCK_H;

    ctx->i_fast_level = i_fast_level;
    ctx->hist_peak_crop = hist_peak_crop;
    ctx->i_use_autolevel = i_use_autolevel;

    ctx->buf_rgb = (uint8_t *)memalign(16 , w * h * 3);
    assert(ctx->buf_rgb != 0);

    HG_Mat_U8_CreateSimplely(&ctx->img_rgb , HG_CS_RGB24 , w , h , ctx->buf_rgb);

    return (long)ctx;
}

int HG_I420_CLAHE_Destroy(long h)
{
    HG_I420_CLAHE_Context *ctx;

    ctx = (HG_I420_CLAHE_Context *)h;

    if (ctx->buf_rgb){
        free(ctx->buf_rgb);
        ctx->buf_rgb = 0;
    }

    free(ctx);
    return 0;
}

static void CalcSubHistsRGB(HG_I420_CLAHE_Context *ctx , struct HG_Mat_U8 *img)
{
    uint8_t *p_pix;
    uint32_t *p_bin_cur;
    uint32_t *p_bin_learn;
    int idx_block_h , idx_block_w;
    int row , col;
    int skip;
    int i;

    skip = img->pitch[0] - 3 * ctx->block_w;

    memset(ctx->bin_cur , 0 , sizeof(ctx->bin_cur));

    for(idx_block_h = 0 ; idx_block_h < CLAHE_NBLOCK_H ; ++idx_block_h){
        for(idx_block_w = 0 ; idx_block_w < CLAHE_NBLOCK_W ; ++idx_block_w){
            p_pix = img->planes[0] + idx_block_h * ctx->block_h * img->pitch[0] + idx_block_w * ctx->block_w * 3/*RGB*/;
            p_bin_cur = ctx->bin_cur[idx_block_h][idx_block_w];

            for(row = 0 ; row < ctx->block_h ; ++row){
                for(col = 0 ; col < ctx->block_w ; ++col) {
                    ++p_bin_cur[*p_pix++];
                    ++p_bin_cur[*p_pix++];
                    ++p_bin_cur[*p_pix++];
                }
                p_pix += skip;
            }
        }
    }

    if (ctx->n_frames > 0) {
        for(idx_block_h = 0 ; idx_block_h < CLAHE_NBLOCK_H ; ++idx_block_h){
            for(idx_block_w = 0 ; idx_block_w < CLAHE_NBLOCK_W ; ++idx_block_w){
                p_bin_cur = ctx->bin_cur[idx_block_h][idx_block_w];
                p_bin_learn = ctx->bin_learn[idx_block_h][idx_block_w];
                for(i = 0 ; i < HIST_MAX_BIN ; ++i){
                    p_bin_learn[i] = (1 - LEARN_RATIO) * p_bin_learn[i] + LEARN_RATIO * p_bin_cur[i];
                }
            }
        }
    }else {
        for(idx_block_h = 0 ; idx_block_h < CLAHE_NBLOCK_H ; ++idx_block_h){
            for(idx_block_w = 0 ; idx_block_w < CLAHE_NBLOCK_W ; ++idx_block_w){
                p_bin_cur = ctx->bin_cur[idx_block_h][idx_block_w];
                p_bin_learn = ctx->bin_learn[idx_block_h][idx_block_w];
                for(i = 0 ; i < HIST_MAX_BIN ; ++i){
                    p_bin_learn[i] = p_bin_cur[i];
                }
            }
        }
    }

    memcpy(ctx->bin_crop , ctx->bin_learn , sizeof(ctx->bin_learn));

}

static void BilinearInterPoSubRGB(struct HG_Mat_U8 *img , struct HG_Rect *rect ,
                              int U , int D , int L , int R ,
                              uint8_t *p_mapUL , uint8_t *p_mapUR , uint8_t *p_mapDL , uint8_t * p_mapDR)
{
    uint8_t v;
    uint8_t *p_pix;
    int skip;
    int col , row;
    int end_row , start_col , end_col;
    int32_t vy0 , vy1;
    int step_w , step_h;
    int32_t factor_y , factor_x , factor_x_step;

    p_pix       = img->planes[0] + rect->start_y * img->pitch[0] + rect->start_x * 3 /*RGB*/;

    skip        = img->pitch[0]  - rect->cols * 3 /*RGB*/;
    end_row     = rect->start_y + rect->rows - 1;
    start_col   = rect->start_x;
    end_col     = rect->start_x + rect->cols - 1;

    step_w  = (R - L + 1);
    step_h  = (D - U + 1);

    factor_x_step = (1 << 16) / step_w;

    for(row = rect->start_y ; row <= end_row  ; ++row){
        factor_y            = ((row       - U) << 16) / step_h;
        factor_x            = ((start_col - L) << 16) / step_w;
        for(col = start_col ; col  <= end_col ; ++col){
            v               = *p_pix;
            vy0             = p_mapUL[v] + ((factor_x * (p_mapUR[v] - p_mapUL[v])) >> 16);
            vy1             = p_mapDL[v] + ((factor_x * (p_mapDR[v] - p_mapDL[v])) >> 16);
            *p_pix++        = vy0        + ((factor_y * (vy1        - vy0       )) >> 16);

            v               = *p_pix;
            vy0             = p_mapUL[v] + ((factor_x * (p_mapUR[v] - p_mapUL[v])) >> 16);
            vy1             = p_mapDL[v] + ((factor_x * (p_mapDR[v] - p_mapDL[v])) >> 16);
            *p_pix++        = vy0        + ((factor_y * (vy1        - vy0       )) >> 16);

            v               = *p_pix;
            vy0             = p_mapUL[v] + ((factor_x * (p_mapUR[v] - p_mapUL[v])) >> 16);
            vy1             = p_mapDL[v] + ((factor_x * (p_mapDR[v] - p_mapDL[v])) >> 16);
            *p_pix++        = vy0        + ((factor_y * (vy1        - vy0       )) >> 16);

            factor_x       += factor_x_step;
        }
        p_pix += skip;
    }
}

static void BilinearInterPoRGB(struct HG_Mat_U8 *img , int block_w , int block_h ,
                               uint8_t  (*pix_map)[CLAHE_NBLOCK_H][HIST_MAX_BIN] )
{

    int  U , D , L , R ; /*center point*/
    int U_Block , D_Block , L_Block , R_Block; /*pixel map index*/
    int col_block , row_block;
    int half_block_w , half_block_h;
    struct HG_Rect rect_proc;

    half_block_h    = block_h / 2;
    half_block_w    = block_w / 2;

    //HG_DEBUG("block_w = %d , block_h = %d\n" , block_w , block_h);

    for(row_block = 0 ; row_block < CLAHE_NBLOCK_H + 1 ; ++row_block){
        if (row_block == 0){
            U_Block = 0;
            D_Block = 0;

            U = half_block_h;
            D = half_block_h;

            rect_proc.start_y = 0;
            rect_proc.rows    = half_block_h;

        }else if (row_block == CLAHE_NBLOCK_H){
            U_Block = CLAHE_NBLOCK_H - 1;
            D_Block = CLAHE_NBLOCK_H - 1;

            U = U_Block * block_h + half_block_h;
            D = D_Block * block_h + half_block_h;

            rect_proc.start_y = U;
            rect_proc.rows    = img->rows - rect_proc.start_y;
        }else {
            U_Block = row_block - 1;
            D_Block = row_block;

            U = U_Block * block_h + half_block_h;
            D = D_Block * block_h + half_block_h;

            rect_proc.start_y = U;
            rect_proc.rows    = D - U;
        }

        for(col_block = 0 ; col_block < CLAHE_NBLOCK_W + 1 ; ++col_block){
            if (col_block == 0){
                L_Block = 0;
                R_Block = 0;

                L = half_block_w;
                R = half_block_w;

                rect_proc.start_x = 0;
                rect_proc.cols    = half_block_w;

            }else if (col_block == CLAHE_NBLOCK_W){
                L_Block = CLAHE_NBLOCK_W - 1;
                R_Block = CLAHE_NBLOCK_W - 1;

                L = L_Block * block_w + half_block_w;
                R = L_Block * block_w + half_block_w;

                rect_proc.start_x = L;
                rect_proc.cols    = img->cols - rect_proc.start_x;
            }else {
                L_Block = col_block - 1;
                R_Block = col_block;

                L = L_Block * block_w + half_block_w;
                R = R_Block * block_w + half_block_w;

                rect_proc.start_x = L;
                rect_proc.cols    = R - L;
            }

            BilinearInterPoSubRGB(img ,  &rect_proc , U , D , L , R ,
                    pix_map[U_Block][L_Block] , pix_map[U_Block][R_Block] ,
                    pix_map[D_Block][L_Block] , pix_map[D_Block][R_Block]
            );
        }
    }
}

int HG_I420_CLAHE_Process(long h , struct HG_Mat_U8 *img)
{
    uint32_t t0 , t1 , t2 , t3 , t4 , t5;
    HG_I420_CLAHE_Context *ctx;
    int row_block , col_block;

    ctx = (HG_I420_CLAHE_Context *)h;

    if (img->cols != ctx->w || img->rows != ctx->h)
        return HG_ERR_INVALID_ARG;

    t0 = HG_GetTimeMs();
    HG_Mat_U8_Convert(img , &ctx->img_rgb);
    t1 = HG_GetTimeMs();

    if (ctx->n_frames % (ctx->i_fast_level + 1) == 0)
        CalcSubHistsRGB(ctx , &ctx->img_rgb);
    t2 = HG_GetTimeMs();

    /*crop hist bin*/
    if (ctx->n_frames % (ctx->i_fast_level + 1) == 0) { 
        for(row_block = 0 ; row_block < CLAHE_NBLOCK_H ; ++row_block){
            for(col_block = 0 ; col_block < CLAHE_NBLOCK_W ; ++col_block){
                if (ctx->i_use_autolevel){
                    PixMap_AutoLevel(ctx->bin_crop[row_block][col_block] , ctx->pix_map[row_block][col_block]);
                }else {
                    CropHistogram(ctx->bin_crop[row_block][col_block] , ctx->hist_peak_crop);
                    MapPixelHistEq(ctx->bin_crop[row_block][col_block] , ctx->pix_map[row_block][col_block]);
                }
            }
        }
    }
    t3 = HG_GetTimeMs();

    /*output image*/
    BilinearInterPoRGB(&ctx->img_rgb  , ctx->block_w , ctx->block_h , ctx->pix_map);
    t4 = HG_GetTimeMs();

    /*recover to yuv420*/
    HG_Mat_U8_Convert(&ctx->img_rgb , img);

    t5 = HG_GetTimeMs();
   
    ++ctx->n_frames;

    if (ctx->n_frames >= (1 << 30)){
        ctx->n_frames = LEARN_FRAMES;
    }

    if (0)
    printf("delay = %d,%d,%d,%d,%d , all = %d ms \n" , (int)(t1 - t0) , (int)(t2 - t1) , (int)(t3 - t2) , (int)(t4 - t3) , (int)(t5 - t4) , 
        (int)(t5 - t0)
    );

    return 0;
}




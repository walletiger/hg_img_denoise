//
// Created by hugo on 19-1-10.
//

#include <malloc.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <hg_os.h>

#include "hg_local_variance_denoise.h"


typedef struct _HG_LocalVarianceDenoise_Context
{
    int w;
    int h;
    int r;

    uint16_t	*sum_line_I;
    uint32_t	*sum_sq_line_I;
    uint16_t	*u_line_I;
    uint32_t	*delta_line_I;
    float		lamda;
    int 		lamda_div;
    uint8_t     delta_a_map[65536];
}HG_LocalVarianceDenoise_Context;

HG_FILTER_HANDLE    hg_local_variance_denoise_filter_create(int w , int h , int r, float lamda)
{
    HG_LocalVarianceDenoise_Context *ctx;

    if (r > HG_LOCAL_VARIANCE_DENOISE_FILTER_R_MAX)
        return 0;

    ctx = (HG_LocalVarianceDenoise_Context *)calloc(1 , sizeof(HG_LocalVarianceDenoise_Context));
    assert(ctx != NULL);


    ctx->w = w;
    ctx->h = h;
    ctx->r = r;
    ctx->lamda = lamda;
    ctx->lamda_div = lamda * 255 * 255;

    ctx->sum_line_I	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);
    ctx->sum_sq_line_I	= (uint32_t *)hg_aligned_malloc(16, sizeof(uint32_t) * w);

    ctx->u_line_I	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);
    ctx->delta_line_I	= (uint32_t *)hg_aligned_malloc(16, sizeof(uint32_t) * w);

    assert(ctx->sum_line_I != 0);
    assert(ctx->sum_sq_line_I != 0);
    {
        int delta = 0;

        for(delta = 0; delta < 65536; delta++){
            ctx->delta_a_map[delta] = 0.5 + (1.0f * (delta << 16)) / ((delta << 8) + ctx->lamda_div);
        }

    }
    return (HG_FILTER_HANDLE)ctx;
}


inline static void hg_local_variance_denoise_output_q_line_s(HG_LocalVarianceDenoise_Context *ctx, int wnd_h, uint8_t *I, uint8_t *Q)
{
    uint32_t sum_box;
    uint32_t sum_box_sq;
    int col , wnd_w , wnd_size ;
    //uint16_t u; // u * u -> signed32
    uint32_t u;
    uint32_t delta;
    int r;
    uint16_t *sum_col;
    uint32_t *sum_col_sq;

    uint32_t coeff_mul;
    register uint32_t coeff_mul_h16, coeff_mul_l8;
    uint16_t *u_line;
    uint32_t *delta_line;
    uint8_t k;

    u_line = ctx->u_line_I;
    delta_line = ctx->delta_line_I;

    r = ctx->r;

    sum_col = ctx->sum_line_I;
    sum_col_sq = ctx->sum_sq_line_I;

    /*calc the first col*/
    wnd_w	 = r + 1;
    sum_box = 0;
    sum_box_sq = 0;

    for(col = 0 ; col <= r ; ++col){
        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];
    }

    wnd_size = wnd_w * wnd_h;

    coeff_mul = (1u << 24) / wnd_size;
    coeff_mul_h16 = (coeff_mul >> 8);
    coeff_mul_l8 = (coeff_mul & 0xffu);

    u_line[0] = u  = (sum_box * coeff_mul) >> 16;

    /*w * 2^16 * (2^16 / w) <  2^32*/
    delta_line[0]  = (((sum_box_sq * coeff_mul_h16) + ((sum_box_sq >> 8) * coeff_mul_l8) - (u * u)) >> 16);

    /*calc the left*/
    for(; col <= 2 * r ; ++col){
        ++wnd_w;

        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];

        wnd_size 	= wnd_w * wnd_h;

        coeff_mul = (1u << 24) / wnd_size;
        coeff_mul_h16 = (coeff_mul >> 8);
        coeff_mul_l8 = (coeff_mul & 0xffu);

        u_line[col - r]	= u = (sum_box * coeff_mul) >> 16;
        delta_line[col - r]	  = (((sum_box_sq * coeff_mul_h16) + ((sum_box_sq >> 8) * coeff_mul_l8) - (u * u)) >> 16);

    }

    /*calc the middle*/
    wnd_size = wnd_w * wnd_h;
    coeff_mul = (1u << 24) / wnd_size;
    coeff_mul_h16 = (coeff_mul >> 8);
    coeff_mul_l8 = (coeff_mul & 0xffu);

    for(; col < ctx->w	; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        sum_box_sq += sum_col_sq[col] - sum_col_sq[col - wnd_w];

        u_line[col - r]	 = u =  (sum_box * coeff_mul) >> 16;
        delta_line[col - r]	  = (((sum_box_sq * coeff_mul_h16) + ((sum_box_sq >> 8) * coeff_mul_l8) - (u * u)) >> 16);
    }

    /*calc the right*/
    for( ; col <  ctx->w + r ; ++col){
        --wnd_w;
        sum_box	 -= sum_col[ctx->w - wnd_w - 1];
        sum_box_sq  -= sum_col_sq[ctx->w - wnd_w - 1];

        wnd_size -= wnd_h;

        coeff_mul = (1u << 24) / wnd_size;
        coeff_mul_h16 = (coeff_mul >> 8);
        coeff_mul_l8 = (coeff_mul & 0xffu);

        u_line[col - r]	 = u =  (sum_box * coeff_mul) >> 16;
        delta_line[col - r]	  = (((sum_box_sq * coeff_mul_h16) + ((sum_box_sq >> 8) * coeff_mul_l8) - (u * u)) >> 16);
    }

    for(col = 0; col < ctx->w; ++col){
        u = u_line[col] >> 8;
        delta = delta_line[col];
        k = ctx->delta_a_map[delta];
        Q[col] = ((255 - k) * u +  k * I[col]) >> 8;
    }
}

inline static void hg_local_variance_denoise_output_q_line_l(HG_LocalVarianceDenoise_Context *ctx, int wnd_h, uint8_t *I, uint8_t *Q)
{
    uint32_t sum_box;
    uint32_t sum_box_sq;
    int col , wnd_w , wnd_size ;
    uint32_t u;
    uint32_t delta;
    int r;
    uint16_t *sum_col;
    uint32_t *sum_col_sq;

    uint64_t coeff_mul;
    uint16_t *u_line;
    uint32_t *delta_line;
    uint8_t k;

    u_line = ctx->u_line_I;
    delta_line = ctx->delta_line_I;

    r = ctx->r;

    sum_col = ctx->sum_line_I;
    sum_col_sq = ctx->sum_sq_line_I;

    /*calc the first col*/
    wnd_w	 = r + 1;
    sum_box = 0;
    sum_box_sq = 0;

    for(col = 0 ; col <= r ; ++col){
        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];
    }

    wnd_size = wnd_w * wnd_h;

    coeff_mul = (1ull << 32) / wnd_size;
    u_line[0] = u  = (sum_box * coeff_mul) >> 24;

    delta_line[0]  = (((sum_box_sq * coeff_mul) >> 16) - u * u) >> 16;

    /*calc the left*/
    for(; col <= 2 * r ; ++col){
        ++wnd_w;

        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];

        wnd_size 	= wnd_w * wnd_h;

        coeff_mul = (1ull << 32) / wnd_size;

        u_line[col - r]	= u = (sum_box * coeff_mul) >> 24;
        delta_line[col - r]	  =  (((sum_box_sq * coeff_mul) >> 16) - u * u) >> 16;

    }

    /*calc the middle*/
    wnd_size = wnd_w * wnd_h;
    coeff_mul = (1ull << 32) / wnd_size;

    for(; col < ctx->w	; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        sum_box_sq += sum_col_sq[col] - sum_col_sq[col - wnd_w];

        u_line[col - r]	 = u =  (sum_box * coeff_mul) >> 24;
        delta_line[col - r]	  =  (((sum_box_sq * coeff_mul) >> 16) - u * u) >> 16;
    }

    /*calc the right*/
    for( ; col <  ctx->w + r ; ++col){
        --wnd_w;
        sum_box	 -= sum_col[ctx->w - wnd_w - 1];
        sum_box_sq  -= sum_col_sq[ctx->w - wnd_w - 1];

        wnd_size -= wnd_h;

        coeff_mul = (1ull << 32) / wnd_size;

        u_line[col - r]	= u = (sum_box * coeff_mul) >> 24;
        delta_line[col - r]	  =  (((sum_box_sq * coeff_mul) >> 16) - u * u) >> 16;
    }

    for(col = 0; col < ctx->w; ++col){
        u = u_line[col] >> 8;
        delta = delta_line[col];

        k = ctx->delta_a_map[delta];
        Q[col] = ((255 - k) * u +  k * I[col]) >> 8;
    }
}

inline static void hg_local_variance_denoise_output_q_line(HG_LocalVarianceDenoise_Context *ctx, int wnd_h, uint8_t *I, uint8_t *Q)
{
    if (ctx->r <= 15)
        hg_local_variance_denoise_output_q_line_s(ctx, wnd_h, I, Q);
    else
        hg_local_variance_denoise_output_q_line_l(ctx, wnd_h, I, Q);
}


int hg_local_variance_denoise_filter_process(HG_FILTER_HANDLE h , unsigned char * pIP, int pitch_ip, unsigned char *pQ, int pitch_q)
{
    HG_LocalVarianceDenoise_Context *ctx;

    uint8_t *p_src , *p_last , v, v1;

    int col , row;
    int  wnd_h;
    int offset;
    int r, w;


    ctx = (HG_LocalVarianceDenoise_Context *)h;

    r = ctx->r;
    w = ctx->w;

    memset(ctx->sum_line_I, 0, sizeof(uint16_t) * w);
    memset(ctx->sum_sq_line_I, 0, sizeof(uint32_t) * w);

    p_src  = pIP;

    /*calc the first row sum and sum square*/
    for(row = 0 ; row <= r ; ++row , p_src += pitch_ip){
        for(col = 0; col < w ; ++col){
            v = p_src[col];

            ctx->sum_line_I[col] += v;
            ctx->sum_sq_line_I[col] += v * v;
        }
    }
    wnd_h = r + 1;


    hg_local_variance_denoise_output_q_line(ctx, wnd_h, pIP,  pQ);

    /*calc the top r */
    for( ; row <= 2 * r ; ++row , p_src += pitch_ip){
        ++wnd_h;

        for(col = 0; col < w ; ++col){
            v = p_src[col];
            ctx->sum_line_I[col] += v;
            ctx->sum_sq_line_I[col] += v * v;
        }

        hg_local_variance_denoise_output_q_line(ctx, wnd_h, pIP + (row - r) * pitch_ip, pQ + (row - r) * pitch_q);
    }


    /*calc vertical middle*/
    offset = pitch_ip * wnd_h;

    for(; row < ctx->h ; ++row , p_src += pitch_ip){
        for(col = 0; col < w ; ++col){
            v =  p_src[col];
            v1 = (p_src - offset)[col];

            ctx->sum_line_I[col] +=  v - v1;
            ctx->sum_sq_line_I[col] += v * v - v1 * v1;
            //ctx->auto_level_bin[v]++;
        }

        hg_local_variance_denoise_output_q_line(ctx, wnd_h, pIP + (row - r) * pitch_ip, pQ + (row - r) * pitch_q);
    }

    /*calc the bottom sum and sum_sq*/
    p_last = p_src - pitch_ip; // last line of pix buf
    for(; row < ctx->h + r ; ++row){
        --wnd_h;
        offset = pitch_ip * wnd_h;

        for(col = 0 ; col < w ; ++col){
            v = (p_last - offset)[col];
            ctx->sum_line_I[col] -= v;
            ctx->sum_sq_line_I[col] -= v * v;
        }

        hg_local_variance_denoise_output_q_line(ctx, wnd_h, pIP + (row - r) * pitch_ip, pQ + (row - r) * pitch_q);

    }
    return 0;
}


void    hg_local_variance_denoise_filter_destroy(HG_FILTER_HANDLE h)
{
    HG_LocalVarianceDenoise_Context *ctx;

    ctx = (HG_LocalVarianceDenoise_Context *)h;
    hg_aligned_free(ctx->sum_line_I);
    hg_aligned_free(ctx->sum_sq_line_I);
    hg_aligned_free(ctx->u_line_I);
    hg_aligned_free(ctx->delta_line_I);

    free(ctx);
}



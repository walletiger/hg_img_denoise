//
// Created by hugo on 19-1-10.
//

#include <malloc.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <hg_os.h>

#ifdef __arm__
#include <arm_neon.h>


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
    uint8_t    delta_a_map[65536];
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
    uint32_t u;
    uint32_t delta;
    int r;
    uint16_t *sum_col;
    uint32_t *sum_col_sq;
    uint32_t coeff_mul;
    uint16_t k;
    uint32_t *p_sum_box;
    uint32_t *p_sum_box_sq;
    
    
    uint32_t sum_box_line[HG_LOCAL_VARIANCE_DENOISE_FILTER_IMG_W_MAX];
    uint32_t sum_box_sq_line[HG_LOCAL_VARIANCE_DENOISE_FILTER_IMG_W_MAX];   

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

    sum_box_line[0] = sum_box;
    sum_box_sq_line[0] = sum_box_sq;

    /*calc the left*/
    for(; col <= 2 * r ; ++col){
        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];

        sum_box_line[col - r] = sum_box;
        sum_box_sq_line[col - r] = sum_box_sq;
    }

    /*calc the middle*/
    wnd_w = 2 * ctx->r + 1;
        
    for(; (col < ctx->w) && col % 8 != 0 ; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        sum_box_sq += sum_col_sq[col] - sum_col_sq[col - wnd_w];
        sum_box_line[col - r] = sum_box;
        sum_box_sq_line[col - r] = sum_box_sq;
    }
    if (1)
    {
        uint16x8_t sum_col_16x8, sum_col_p16x8, sum_col_diff;
        uint32x4_t sum_col_sq_32x4, sum_col_sq_p32x4, sum_col_sq_diff32x4;
        
        uint16_t *p_col, *p_col_p;
        uint32_t *p_col_sq, *p_col_p_sq;
        uint32_t  *p_sum_box_line, *p_sum_box_sq_line;
        
        int16_t *p_diff;
        int32_t *p_diff32;
        
        p_col             = sum_col + col;
        p_col_p           = sum_col + col - wnd_w;
        p_col_sq          = sum_col_sq + col;
        p_col_p_sq        = sum_col_sq + col - wnd_w;
        p_sum_box_line    = sum_box_line + col - r;
        p_sum_box_sq_line = sum_box_sq_line + col - r;
        
        for(; col + 8 < ctx->w ; col += 8){
            sum_col_16x8 = vld1q_u16(p_col);
            sum_col_p16x8 = vld1q_u16(p_col_p);  

            sum_col_diff = vsubq_u16(sum_col_16x8, sum_col_p16x8);
            p_diff = (int16_t *)&sum_col_diff;
            sum_box += p_diff[0];  p_sum_box_line[0] = sum_box;
            sum_box += p_diff[1];  p_sum_box_line[1] = sum_box;
            sum_box += p_diff[2];  p_sum_box_line[2] = sum_box;
            sum_box += p_diff[3];  p_sum_box_line[3] = sum_box;
            sum_box += p_diff[4];  p_sum_box_line[4] = sum_box;
            sum_box += p_diff[5];  p_sum_box_line[5] = sum_box;
            sum_box += p_diff[6];  p_sum_box_line[6] = sum_box;
            sum_box += p_diff[7];  p_sum_box_line[7] = sum_box;

            p_col += 8;
            p_col_p += 8;
            p_sum_box_line += 8;
            

            sum_col_sq_32x4     = vld1q_u32(p_col_sq);
            sum_col_sq_p32x4    = vld1q_u32(p_col_p_sq);

            sum_col_sq_diff32x4 = vsubq_u32(sum_col_sq_32x4, sum_col_sq_p32x4);

            p_diff32 = (int32_t *)&sum_col_sq_diff32x4;

            sum_box_sq += p_diff32[0]; p_sum_box_sq_line[0] = sum_box_sq;
            sum_box_sq += p_diff32[1]; p_sum_box_sq_line[1] = sum_box_sq;
            sum_box_sq += p_diff32[2]; p_sum_box_sq_line[2] = sum_box_sq;
            sum_box_sq += p_diff32[3]; p_sum_box_sq_line[3] = sum_box_sq;

            p_col_sq += 4;
            p_col_p_sq += 4;
            p_sum_box_sq_line += 4;
            
            sum_col_sq_32x4     = vld1q_u32(p_col_sq);
            sum_col_sq_p32x4    = vld1q_u32(p_col_p_sq);

            sum_col_sq_diff32x4 = vsubq_u32(sum_col_sq_32x4, sum_col_sq_p32x4);

            p_diff32 = (int32_t *)&sum_col_sq_diff32x4;

            sum_box_sq += p_diff32[0]; p_sum_box_sq_line[0] = sum_box_sq;
            sum_box_sq += p_diff32[1]; p_sum_box_sq_line[1] = sum_box_sq;
            sum_box_sq += p_diff32[2]; p_sum_box_sq_line[2] = sum_box_sq;
            sum_box_sq += p_diff32[3]; p_sum_box_sq_line[3] = sum_box_sq;

            p_col_sq += 4;
            p_col_p_sq += 4;
            p_sum_box_sq_line += 4;
        }
    }

    for(; col < ctx->w	; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        sum_box_sq += sum_col_sq[col] - sum_col_sq[col - wnd_w];
        sum_box_line[col - r] = sum_box;
        sum_box_sq_line[col - r] = sum_box_sq;
    }    

    /*calc the right*/
    for( ; col < ctx->w + r ; ++col){
        --wnd_w;
        sum_box	 -= sum_col[ctx->w  - 1 - wnd_w];
        sum_box_sq  -= sum_col_sq[ctx->w - 1 - wnd_w];
        sum_box_line[col - r] = sum_box;
        sum_box_sq_line[col - r] = sum_box_sq;
    }

    /**********output q*/
    
    wnd_w = r + 1;
    p_sum_box = sum_box_line;
    p_sum_box_sq = sum_box_sq_line;

    for(col = 0; col < r; col++){
        wnd_size = wnd_w * wnd_h;
        coeff_mul = (1u<<16) / wnd_size;
    
        u = ((*p_sum_box) * coeff_mul) >> 8;
        delta = (*p_sum_box_sq *coeff_mul - u * u) >> 16;
        k = ctx->delta_a_map[delta];

        Q[col] = ((((255 - k) * u) >> 8)  +  k * I[col]) >> 8;

        ++p_sum_box;
        ++p_sum_box_sq;
        ++wnd_w;
    }
    
    wnd_w = 2 * r + 1;
    wnd_size = wnd_w * wnd_h;
    coeff_mul = (1u<<16) / wnd_size;
    
    for(; col < ctx->w && (col % 8 != 0); ++col){
        u = ((*p_sum_box) * coeff_mul) >> 8;
        delta = (*p_sum_box_sq *coeff_mul - u * u) >> 16;
        k = ctx->delta_a_map[delta];

        Q[col] = ((((255 - k) * u) >> 8)  +  k * I[col]) >> 8;

        ++p_sum_box;
        ++p_sum_box_sq;

    }
    if (1)
    {
        uint32x4_t  delta_arg0_32x4, delta_arg1_32x4, *p_sum_box_sq_32x4, sum_box0_32x4, sum_box1_32x4, *p_sum_box_u32x4;
        uint16x8_t  usq_u16x8, u_16x8, delta_16x8, delta_arg_16x8, kxor_16x8, kxorXu, kxi_16x8;
        uint8x8_t k_8x8, mask_8x8, *p_i_8x8, *p_q_8x8;
        uint16_t *p_delta;
        uint8_t *p_k;

        
        p_sum_box_u32x4 = (uint32x4_t *) p_sum_box;
        p_sum_box_sq_32x4 = (uint32x4_t *)p_sum_box_sq;
        p_i_8x8 = (uint8x8_t *)(I + col);
        p_q_8x8 = (uint8x8_t *)(Q + col);
        
        mask_8x8 = vcreate_u8(~0ull);
        
        for(; col + 8 < ctx->w - r; col += 8){
            sum_box0_32x4 = vmulq_n_u32 (*p_sum_box_u32x4, coeff_mul); 
            sum_box1_32x4 = vmulq_n_u32 (*(p_sum_box_u32x4 + 1), coeff_mul);
            ((uint16x4_t *)&u_16x8)[0] = vrshrn_n_u32(sum_box0_32x4, 16);
            ((uint16x4_t *)&u_16x8)[1] = vrshrn_n_u32(sum_box1_32x4, 16);

            usq_u16x8 = vmulq_u16(u_16x8, u_16x8); // u * u 

            delta_arg0_32x4 = vmulq_n_u32 (*p_sum_box_sq_32x4, coeff_mul);
            delta_arg1_32x4 = vmulq_n_u32 (*(p_sum_box_sq_32x4 + 1), coeff_mul);
            ((uint16x4_t *)&delta_arg_16x8)[0] = vrshrn_n_u32(delta_arg0_32x4, 16); // >>16 and rounded value
            ((uint16x4_t *)&delta_arg_16x8)[1] = vrshrn_n_u32(delta_arg1_32x4, 16); // >>16 and rounded value

            delta_16x8 = vqsubq_u16(delta_arg_16x8, usq_u16x8); // Mean(sq) - u * u,  saturating subtract
            p_delta = (uint16_t *)&delta_16x8;

            p_k = (uint8_t *)&k_8x8;
            p_k[0] = ctx->delta_a_map[p_delta[0]]; p_k[1] = ctx->delta_a_map[p_delta[1]]; 
            p_k[2] = ctx->delta_a_map[p_delta[2]]; p_k[3] = ctx->delta_a_map[p_delta[3]];
            p_k[4] = ctx->delta_a_map[p_delta[4]]; p_k[5] = ctx->delta_a_map[p_delta[5]];
            p_k[6] = ctx->delta_a_map[p_delta[6]]; p_k[7] = ctx->delta_a_map[p_delta[7]]; 


            kxor_16x8 = vsubl_u8 (mask_8x8, k_8x8);// 255 - k
            kxorXu = vmulq_u16 (kxor_16x8, u_16x8); // (255 - k) * u            
            kxi_16x8 = vmull_u8(k_8x8, *p_i_8x8);
            
            *p_q_8x8 = vraddhn_u16(kxorXu, kxi_16x8); // ((255 - k) * u + k * i) >> 8

            p_sum_box_u32x4 += 2;
            p_sum_box_sq_32x4 += 2;
            ++p_i_8x8;
            ++p_q_8x8;
        }

        p_sum_box = (uint32_t *)p_sum_box_u32x4;
        p_sum_box_sq = (uint32_t *)p_sum_box_sq_32x4;
    }

    for(; col < ctx->w - r; ++col){
        u = ((*p_sum_box) * coeff_mul) >> 8;
        delta = (*p_sum_box_sq *coeff_mul - u * u) >> 16;
        k = ctx->delta_a_map[delta];

        Q[col] = ((((255 - k) * u) >> 8)  +  k * I[col]) >> 8;


        ++p_sum_box;
        ++p_sum_box_sq;
    }

    for(; col < ctx->w; col++){
        --wnd_w;
        wnd_size = wnd_w * wnd_h;
        coeff_mul = (1u<<16) / wnd_size;
        
        u = ((*p_sum_box) * coeff_mul) >> 8;
        delta = (*p_sum_box_sq *coeff_mul - u * u) >> 16;
        k = ctx->delta_a_map[delta];

        Q[col] = ((((255 - k) * u) >> 8)  +  k * I[col]) >> 8;

        ++p_sum_box;
        ++p_sum_box_sq;
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
        uint16x8_t *p_sum;
        uint8x8_t *p_pix, *p_pix_prev, pix, pix_prev;
        uint16x8_t sum_sq, sum_sq_prev;
        uint16x4_t *p_sum_sq_cur, *p_sum_sq_prev;
        uint32x4_t *p_sum_sq_dst;

        p_pix = (uint8x8_t *)p_src;
        p_pix_prev = (uint8x8_t *)(p_src - offset);
        p_sum = (uint16x8_t *)ctx->sum_line_I;
        p_sum_sq_dst = (uint32x4_t *)ctx->sum_sq_line_I;

        for(col = 0; col + 8 < w ; col += 8){
        	pix = *p_pix;
        	pix_prev = *p_pix_prev;

        	*p_sum = vaddw_u8(*p_sum, pix);
        	*p_sum = vsubw_u8(*p_sum, pix_prev);

        	sum_sq = vmull_u8(pix, pix);
        	sum_sq_prev = vmull_u8(pix_prev, pix_prev);

        	p_sum_sq_cur = (uint16x4_t *)&sum_sq;
        	p_sum_sq_prev = (uint16x4_t *)&sum_sq_prev; 
        	
        	*p_sum_sq_dst = vaddw_u16(*p_sum_sq_dst, *p_sum_sq_cur);
        	*(p_sum_sq_dst + 1) = vaddw_u16(*(p_sum_sq_dst + 1), *(p_sum_sq_cur + 1));

        	*p_sum_sq_dst = vsubw_u16(*p_sum_sq_dst, *p_sum_sq_prev);
        	*(p_sum_sq_dst + 1) = vsubw_u16(*(p_sum_sq_dst + 1), *(p_sum_sq_prev + 1));

        	++p_sum;
        	++p_pix;
        	++p_pix_prev;

        	p_sum_sq_dst += 2;
        }
        
        for(; col < w ; ++col){
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

#endif

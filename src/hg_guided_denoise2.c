#include <malloc.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <hg_os.h>
#include "hg_guided_denoise_filter.h"

typedef struct _HG_GuidedDenoise_Context2
{
    int w;
    int h;
    int r;
    float       min_pcnt;
    float       max_pcnt;

    uint16_t	*ringbuf_ab;
    uint16_t	*sum_line_I; 
    uint32_t	*sum_sq_line_I;
    uint16_t	*sum_line_A;	
    uint16_t	*sum_line_B;	
	uint16_t	*u_line_I;
	uint32_t	*delta_line_I;	
    float		lamda;
    int 		idx_ab_wnd;
    int			ab_wnd_size;
    int 		lamda_div;
    uint8_t    delta_a_map[65536]; 
}HG_GuidedDenoise_Context2;

HG_FILTER_HANDLE    hg_guided_denoise_filter2_create(int w , int h , int r, float lamda)
{
    HG_GuidedDenoise_Context2 *ctx;

    if (r > HG_GUIDED_DENOISE_FILTER_R_MAX)
        return 0;

    ctx = (HG_GuidedDenoise_Context2 *)calloc(1 , sizeof(HG_GuidedDenoise_Context2));
    assert(ctx != NULL);


    ctx->w = w;
    ctx->h = h;
    ctx->r = r;
    ctx->min_pcnt = 4.0/255;
    ctx->max_pcnt = 250.0/255;
    ctx->lamda = lamda;
    ctx->lamda_div = lamda * 255 * 255;
    ctx->idx_ab_wnd = 0;
    ctx->ab_wnd_size = 2 * r + 1;

    ctx->ringbuf_ab = (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * (2 * r + 1) * w);
    ctx->sum_line_I	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);
    ctx->sum_sq_line_I	= (uint32_t *)hg_aligned_malloc(16, sizeof(uint32_t) * w);
    ctx->sum_line_A	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);
    ctx->sum_line_B	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);

    ctx->u_line_I	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);
    ctx->delta_line_I	= (uint32_t *)hg_aligned_malloc(16, sizeof(uint32_t) * w);

    assert(ctx->ringbuf_ab != 0);
    assert(ctx->sum_line_I != 0);
    assert(ctx->sum_sq_line_I != 0);
    assert(ctx->sum_line_A != 0);
    assert(ctx->sum_line_B != 0);

    {
        int delta = 0;
        
        for(delta = 0; delta < 65536; delta++){
            ctx->delta_a_map[delta] = 0.5 + (1.0f * (delta << 16)) / ((delta << 8) + ctx->lamda_div);
        }

    }
    return (HG_FILTER_HANDLE)ctx;
}



inline static void hg_guided_denoise_output_q_line_s(HG_GuidedDenoise_Context2 *ctx,  int wnd_h, uint8_t *I, uint8_t *Q)
{
    int col;
    int wnd_w;
    uint32_t wnd_size;
    uint32_t sum_a, sum_b;
    uint32_t coeff_mul;
    int r;
    uint16_t q;

    r = ctx->r;

    /*1 calc the first pix*/
    wnd_w = r + 1;
    sum_a = 0;
    sum_b = 0;

    for(col = 0 ; col <= ctx->r; ++col){
        sum_a += ctx->sum_line_A[col];
        sum_b += ctx->sum_line_B[col];
    }

    wnd_size = wnd_w * wnd_h;
    coeff_mul = (1u << 16) / wnd_size;

    // (*I++) * sum_a + (sum_b << 8) < wnd_size * 2^16 (max pixel value 255 < 256)
    // coeff_mul  <= 2 ^ 16 / wnd_size
    // so , must ,  (((*I++) * sum_a + (sum_b << 8)) * coeff_mul) < 2^32
    *Q++ = q =  ((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 24);

    /*2 calc the left cols*/
    for(; col <= 2 * ctx->r ; ++col){
        ++wnd_w;
        wnd_size += wnd_h;

        sum_a += ctx->sum_line_A[col];
        sum_b += ctx->sum_line_B[col];

        coeff_mul = (1u << 16) / wnd_size;
		*Q++ = q =	((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 24);


    }

    /*calc the middle*/
    coeff_mul = (1u << 16) / wnd_size;

    for( ; col < ctx->w  ; ++col){
        sum_a += ctx->sum_line_A[col] - ctx->sum_line_A[col - wnd_w];
        sum_b += ctx->sum_line_B[col] - ctx->sum_line_B[col - wnd_w];
        *Q++ = q =  ((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 24);

    }
    
    /*calc the right*/
    for( ; col < ctx->w + r ; ++col){
        --wnd_w;
        wnd_size -= wnd_h;
        coeff_mul = (1u << 16) / wnd_size;

        sum_a -= ctx->sum_line_A[ctx->w - wnd_w -1];
        sum_b -= ctx->sum_line_B[ctx->w - wnd_w -1];
		
        *Q++ = q = ((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 24);

    }

}

inline static void hg_guided_denoise_output_q_line_l(HG_GuidedDenoise_Context2 *ctx,  int wnd_h, uint8_t *I, uint8_t *Q)
{
    int col;
    int wnd_w;
    uint32_t wnd_size;
    uint32_t sum_a, sum_b;
    uint64_t coeff_mul;
    int r;
	uint8_t q;

    r = ctx->r;

    /*1 calc the first pix*/
    wnd_w = r + 1;
    sum_a = 0;
    sum_b = 0;

    for(col = 0 ; col <= ctx->r; ++col){
        sum_a += ctx->sum_line_A[col];
        sum_b += ctx->sum_line_B[col];
    }

    wnd_size = wnd_w * wnd_h;
    coeff_mul = (1ull << 24)/ wnd_size;

    *Q++ = q =  ((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 32);

    /*2 calc the left cols*/
    for(; col <= 2 * ctx->r ; ++col){
        ++wnd_w;
        wnd_size += wnd_h;

        sum_a += ctx->sum_line_A[col];
        sum_b += ctx->sum_line_B[col];

        coeff_mul = (1ull << 24)/ wnd_size;
        *Q++ = q = ((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 32);
	
    }

    /*calc the middle*/
    coeff_mul = (1ull << 24)/ wnd_size;

    for( ; col < ctx->w  ; ++col){
        sum_a += ctx->sum_line_A[col] - ctx->sum_line_A[col - wnd_w];
        sum_b += ctx->sum_line_B[col] - ctx->sum_line_B[col - wnd_w];

        *Q++ = q =  ((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 32);
		
    }

    /*calc the right*/
    for( ; col < ctx->w + r ; ++col){
        --wnd_w;
        wnd_size -= wnd_h;
        coeff_mul = (1ull << 24)/ wnd_size;

        sum_a -= ctx->sum_line_A[ctx->w - wnd_w -1];
        sum_b -= ctx->sum_line_B[ctx->w - wnd_w -1];

        *Q++ = q = ((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 32);
	
    }

}

inline static void hg_guided_denoise_output_q_line(HG_FILTER_HANDLE h,  int wnd_h, uint8_t *I, uint8_t *Q)
{
    HG_GuidedDenoise_Context2 *ctx;

	ctx = (HG_GuidedDenoise_Context2 *)h;

	if (ctx->r <= 15)
		hg_guided_denoise_output_q_line_s(ctx, wnd_h, I, Q);
	else
		hg_guided_denoise_output_q_line_l(ctx, wnd_h, I, Q);

}


inline static void hg_guided_denoise_learn_ab_line_l(HG_GuidedDenoise_Context2 *ctx, int wnd_h)
{
    uint32_t sum_box , sum_box_sq;
    int col , wnd_w , wnd_size ;
    uint32_t u, delta;
    int r;
    uint16_t *sum_col;
    uint32_t *sum_col_sq; 
    uint16_t prev_ab;

    uint16_t *p_ab;
    uint8_t a, b;

    uint16_t *p_sum_a, *p_sum_b;

    uint64_t coeff_mul;
    uint16_t *u_line;
    uint32_t *delta_line;

	u_line = ctx->u_line_I;
	delta_line = ctx->delta_line_I;
	
    r = ctx->r;

    sum_col = ctx->sum_line_I;
    sum_col_sq = ctx->sum_sq_line_I;


    p_ab 		= ctx->ringbuf_ab + ctx->idx_ab_wnd * ctx->w;
    p_sum_a 	= ctx->sum_line_A;
    p_sum_b 	= ctx->sum_line_B;

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
    delta_line[0]  = ((((sum_box_sq * coeff_mul) >> 16)  - (u * u)) >> 16);

    /*calc the left*/
    for(; col <= 2 * r ; ++col){
        ++wnd_w;

        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];

        wnd_size 	= wnd_w * wnd_h;

        coeff_mul = (1ull << 32) / wnd_size;

        u_line[col - r]	= u = (sum_box * coeff_mul) >> 24;
        delta_line[col - r]	= ((((sum_box_sq * coeff_mul) >> 16)  - (u * u)) >> 16);

    }

    /*calc the middle*/
    wnd_size = wnd_w * wnd_h;
    coeff_mul = (1ull << 32) / wnd_size;
    for(; col < ctx->w	; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        sum_box_sq += sum_col_sq[col] - sum_col_sq[col - wnd_w];

        u_line[col - r]	= u = (sum_box * coeff_mul) >> 24;
        delta_line[col - r]	= ((((sum_box_sq * coeff_mul) >> 16)  - (u * u)) >> 16);
    }

    /*calc the right*/
    for( ; col <  ctx->w + r ; ++col){
        --wnd_w;
        sum_box	 -= sum_col[ctx->w - wnd_w - 1];
        sum_box_sq  -= sum_col_sq[ctx->w - wnd_w - 1];

        wnd_size -= wnd_h;


        coeff_mul = (1ull << 32) / wnd_size;

        u_line[col - r]	= u = (sum_box * coeff_mul) >> 24;
        delta_line[col - r]	= ((((sum_box_sq * coeff_mul) >> 16)  - (u * u)) >> 16);
    }


    for(col = 0; col < ctx->w; ++col){
        u = u_line[col];
        delta = delta_line[col];
        a = ctx->delta_a_map[delta];
        b = ((255 - a) * u) >> 16;	  /*0 ~ 255*/
        prev_ab = p_ab[col];
        p_sum_a[col] += a - (prev_ab & 0x00ffu);
        p_sum_b[col] += b - (prev_ab >> 8);
        p_ab[col] = (a | (b << 8));
    }

    ctx->idx_ab_wnd = (ctx->idx_ab_wnd + 1) % ctx->ab_wnd_size;
}



inline static void hg_guided_denoise_learn_ab_line_s(HG_GuidedDenoise_Context2 *ctx, int wnd_h)
{
    uint32_t sum_box;
    uint32_t sum_box_sq;
    int col , wnd_w , wnd_size ;
    uint32_t u;
    uint32_t delta;
    int r;
    uint16_t *sum_col;
    uint32_t *sum_col_sq; 
    uint16_t prev_ab;

    uint16_t *p_ab;
    uint16_t a, b;

    uint16_t *p_sum_a, *p_sum_b;
    uint32_t coeff_mul;
    uint16_t *u_line;
    uint32_t *delta_line;
    register uint32_t coeff_mul_h16, coeff_mul_l8;

    u_line = ctx->u_line_I;
    delta_line = ctx->delta_line_I;
	
    r = ctx->r;

    sum_col = ctx->sum_line_I;
    sum_col_sq = ctx->sum_sq_line_I;

    p_ab 		= ctx->ringbuf_ab + ctx->idx_ab_wnd * ctx->w;
    p_sum_a 	= ctx->sum_line_A;
    p_sum_b 	= ctx->sum_line_B;

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

        u_line[col - r]	= u = (sum_box * coeff_mul) >> 16;
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

        u_line[col - r]	= u = (sum_box * coeff_mul) >> 16;
        delta_line[col - r]	  = (((sum_box_sq * coeff_mul_h16) + ((sum_box_sq >> 8) * coeff_mul_l8) - (u * u)) >> 16);
    }

    for(col = 0; col < ctx->w; ++col){
        u = u_line[col];
        delta = delta_line[col];
        a = ctx->delta_a_map[delta];
        b = ((255 - a) * u) >> 16;	  /*0 ~ 255*/
        prev_ab = p_ab[col];
        p_sum_a[col] += a - (prev_ab & 0x00ffu);
        p_sum_b[col] += b - (prev_ab >> 8);
        p_ab[col] = (a | (b << 8));
    }

    ctx->idx_ab_wnd = (ctx->idx_ab_wnd + 1) % ctx->ab_wnd_size;
}

inline static void hg_guided_denoise_learn_ab_line(HG_FILTER_HANDLE h, int wnd_h)
{
    HG_GuidedDenoise_Context2 *ctx;

	ctx = (HG_GuidedDenoise_Context2 *)h;

	if (ctx->r <= 15)
		hg_guided_denoise_learn_ab_line_s(ctx, wnd_h);
	else
		hg_guided_denoise_learn_ab_line_l(ctx, wnd_h);

}


int hg_guided_denoise_filter2_process(HG_FILTER_HANDLE h , unsigned char * pIP, int pitch_ip, unsigned char *pQ, int pitch_q)
{
    HG_GuidedDenoise_Context2 *ctx;

    uint8_t *p_src , *p_last , v, v1;

    int col , row;
    int  wnd_h;
    int offset;
    int r, w;


    ctx = (HG_GuidedDenoise_Context2 *)h;

    r = ctx->r;
    w = ctx->w;

    memset(ctx->ringbuf_ab, 0, sizeof(uint16_t) * (2 * r + 1) * w);
    memset(ctx->sum_line_I, 0, sizeof(uint16_t) * w);
    memset(ctx->sum_sq_line_I, 0, sizeof(uint32_t) * w);

    memset(ctx->sum_line_A, 0, sizeof(uint16_t) * w);
    memset(ctx->sum_line_B, 0, sizeof(uint16_t) * w);

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


    hg_guided_denoise_learn_ab_line(h, wnd_h);

    /*calc the top r */
    for( ; row <= 2 * r ; ++row , p_src += pitch_ip){
        ++wnd_h;

        for(col = 0; col < w ; ++col){
            v = p_src[col];
            ctx->sum_line_I[col] += v;
            ctx->sum_sq_line_I[col] += v * v;
        }

        hg_guided_denoise_learn_ab_line(h, wnd_h);
    }

    /*calc first line of Q*/
    hg_guided_denoise_output_q_line(h, wnd_h - r,  pIP, pQ);


    /*calc vertical middle*/
    offset = pitch_ip * wnd_h;

    for(; row < ctx->h ; ++row , p_src += pitch_ip){
        for(col = 0; col < w ; ++col){
            v =  p_src[col];
            v1 = (p_src - offset)[col];

            ctx->sum_line_I[col] +=  v - v1;
            ctx->sum_sq_line_I[col] += v * v - v1 * v1;
        }

        {
            int ab_wnd_h;
            
            ab_wnd_h =  row < ctx->ab_wnd_size + r ? row - r + 1 : ctx->ab_wnd_size;

            hg_guided_denoise_learn_ab_line(h, wnd_h);
            hg_guided_denoise_output_q_line(h, ab_wnd_h,  pIP + (row - r * 2) * pitch_ip , pQ + (row - r * 2) * pitch_q);
        }
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

        hg_guided_denoise_learn_ab_line(h, wnd_h);
        hg_guided_denoise_output_q_line(h, ctx->ab_wnd_size,  pIP + (row - r * 2) * pitch_ip , pQ + (row - r * 2) * pitch_q);
        
        
    }
    /*calc the bottom r line of  Q*/
	ctx->idx_ab_wnd = (ctx->idx_ab_wnd + ctx->ab_wnd_size - 1) % ctx->ab_wnd_size;
    
    uint16_t *p_prev_ab;
    uint16_t prev_ab;
    
    wnd_h = ctx->r * 2;
    for(; row < ctx->h + 2 * r ; ++row, --wnd_h){
        p_prev_ab = ctx->ringbuf_ab + w * ((ctx->idx_ab_wnd + ctx->ab_wnd_size - wnd_h) % ctx->ab_wnd_size);
        for(col = 0; col < w; col++){
            prev_ab = p_prev_ab[col];
            ctx->sum_line_A[col] -= prev_ab &(0xffu);
            ctx->sum_line_B[col] -= (prev_ab >> 8);
        }
        hg_guided_denoise_output_q_line(h, wnd_h,  pIP + (row - r * 2) * pitch_ip , pQ + (row - r * 2) * pitch_q);
    }
	
    return 0;

}


void    hg_guided_denoise_filter2_destroy(HG_FILTER_HANDLE h)
{
    HG_GuidedDenoise_Context2 *ctx;

    ctx = (HG_GuidedDenoise_Context2 *)h;

    hg_aligned_free(ctx->ringbuf_ab);
    hg_aligned_free(ctx->sum_line_I);
    hg_aligned_free(ctx->sum_sq_line_I);
    hg_aligned_free(ctx->sum_line_A);
    hg_aligned_free(ctx->sum_line_B);
    hg_aligned_free(ctx->u_line_I);
    hg_aligned_free(ctx->delta_line_I);

    free(ctx);
}



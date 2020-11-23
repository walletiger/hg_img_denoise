#include <malloc.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <hg_os.h>
#include "hg_guided_denoise_filter.h"

#ifdef __arm__

#include <arm_neon.h>

typedef struct _HG_GuidedDenoise_Context2
{
    int w;
    int h;
    int r;
    float       min_pcnt;
    float       max_pcnt;

    uint8_t	    *ringbuf_a;
    uint8_t	    *ringbuf_b;
    
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
    uint8_t     delta_a_map[65536]; 
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

    ctx->ringbuf_a = (uint8_t *)hg_aligned_malloc(16, sizeof(uint8_t) * (2 * r + 1) * w);
    ctx->ringbuf_b = (uint8_t *)hg_aligned_malloc(16, sizeof(uint8_t) * (2 * r + 1) * w);    
    ctx->sum_line_I	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);
    ctx->sum_sq_line_I	= (uint32_t *)hg_aligned_malloc(16, sizeof(uint32_t) * w);
    ctx->sum_line_A	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);
    ctx->sum_line_B	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);

    ctx->u_line_I	= (uint16_t *)hg_aligned_malloc(16, sizeof(uint16_t) * w);
    ctx->delta_line_I	= (uint32_t *)hg_aligned_malloc(16, sizeof(uint32_t) * w);

    assert(ctx->ringbuf_a != 0);
    assert(ctx->ringbuf_b != 0);    
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
    uint32_t sum_a_box[HG_GUIDED_DENOISE_FILTER_IMG_W_MAX];
    uint32_t sum_b_box[HG_GUIDED_DENOISE_FILTER_IMG_W_MAX];
    uint8_t mean_a, mean_b;
    
    r = ctx->r;

    sum_a = 0;
    sum_b = 0;

    for(col = 0 ; col <= ctx->r; ++col){
        sum_a += ctx->sum_line_A[col];
        sum_b += ctx->sum_line_B[col];
    }

    sum_a_box[0] = sum_a;
    sum_b_box[0] = sum_b;
    
    /*2 calc the left cols*/
    for(; col <= 2 * ctx->r ; ++col){
        sum_a += ctx->sum_line_A[col];
        sum_b += ctx->sum_line_B[col];

        sum_a_box[col - r] = sum_a;
        sum_b_box[col - r] = sum_b;
    }
    
    /*3 calc the middle*/
    wnd_w = 2 * r + 1;
    for( ; col < ctx->w  ; ++col){
        sum_a += ctx->sum_line_A[col] - ctx->sum_line_A[col - wnd_w];
        sum_b += ctx->sum_line_B[col] - ctx->sum_line_B[col - wnd_w];
    
        sum_a_box[col - r] = sum_a;
        sum_b_box[col - r] = sum_b;

    }    

    /*4 calc the right*/
    for( ; col < ctx->w + r ; ++col){
        --wnd_w;
        wnd_size -= wnd_h;

        sum_a -= ctx->sum_line_A[ctx->w - wnd_w - 1];
        sum_b -= ctx->sum_line_B[ctx->w - wnd_w - 1];

        sum_a_box[col - r] = sum_a;
        sum_b_box[col - r] = sum_b;

   }

    //wnd_size = wnd_w * wnd_h;
    //coeff_mul = (1u << 16) / wnd_size;

    // (*I++) * sum_a + (sum_b << 8) < wnd_size * 2^16 (max pixel value 255 < 256)
    // coeff_mul  <= 2 ^ 16 / wnd_size
    // so , must ,  (((*I++) * sum_a + (sum_b << 8)) * coeff_mul) < 2^32

    
    wnd_w = r + 1;
    for(col = 0; col < r ; col++){
        wnd_size = wnd_w * wnd_h;
        coeff_mul = (1u << 16) / wnd_size;
    
        mean_a = (sum_a_box[col] * coeff_mul + 32768) >> 16;
        mean_b = (sum_b_box[col] * coeff_mul + 32768) >> 16;

        *Q++ = ((*I++) * mean_a + (mean_b << 8) + 127) >> 8;

        ++wnd_w ;
    }
    
    wnd_w = 2 * r + 1;
    wnd_size = wnd_w * wnd_h;
    coeff_mul = (1u << 16) / wnd_size;    
    
    for(; col < ctx->w && (col % 8 != 0) ; col++){
        mean_a = (sum_a_box[col] * coeff_mul + 32768) >> 16;
        mean_b = (sum_b_box[col] * coeff_mul + 32768) >> 16;
        *Q++ = ((*I++) * mean_a + (mean_b << 8) + 127) >> 8;
    }

    if (1)
    {
        uint8x8_t *p_i, *p_q, u_a_8x8;
        uint32x4_t *p_sum_a_32x4, *p_sum_b_32x4;
        uint16x8_t u_a_16x8, u_b_16x8;
        uint16x8_t axi_plus_b;
        uint32x4_t sum_0_32x4, sum_1_32x4;
        
        p_i = (uint8x8_t *)I;
        p_q = (uint8x8_t *)Q;
        
        p_sum_a_32x4 = (uint32x4_t *)(&sum_a_box[col]);
        p_sum_b_32x4 = (uint32x4_t *)(&sum_b_box[col]);

        for(; col + 8 < ctx->w - r ; col += 8){
            sum_0_32x4 = vmulq_n_u32 (p_sum_a_32x4[0], coeff_mul);
            sum_1_32x4 = vmulq_n_u32 (p_sum_a_32x4[1], coeff_mul); 

            ((uint16x4_t *)&u_a_16x8)[0] = vrshrn_n_u32(sum_0_32x4, 16);
            ((uint16x4_t *)&u_a_16x8)[1] = vrshrn_n_u32(sum_1_32x4, 16);
            
            sum_0_32x4 = vmulq_n_u32 (p_sum_b_32x4[0], coeff_mul);
            sum_1_32x4 = vmulq_n_u32 (p_sum_b_32x4[1], coeff_mul); 

            ((uint16x4_t *)&u_b_16x8)[0] = vrshrn_n_u32(sum_0_32x4, 8);
            ((uint16x4_t *)&u_b_16x8)[1] = vrshrn_n_u32(sum_1_32x4, 8);

            /*printf("1b = %d,  %d, %d, %d, %d, %d, %d, %d\n", 
                vgetq_lane_u16(u_b_16x8, 0), vgetq_lane_u16(u_b_16x8, 1), vgetq_lane_u16(u_b_16x8, 2), vgetq_lane_u16(u_b_16x8, 3),
                sum_b_box[col],  sum_b_box[col + 1], sum_b_box[col + 2], sum_b_box[col + 3]
            );
            exit(1);*/
            
            u_a_8x8  = vmovn_u16 (u_a_16x8);

            axi_plus_b = vmlal_u8(u_b_16x8, u_a_8x8, *p_i);
            *p_q = vqrshrn_n_u16(axi_plus_b, 8); // >> 8

            p_i++;
            p_q++;
            p_sum_a_32x4 += 2;
            p_sum_b_32x4 += 2;
        }
        
        Q = (uint8_t *)p_q;
        I = (uint8_t *)p_i;
    }

    
    for(; col < ctx->w - r ; col++){
        mean_a = (sum_a_box[col] * coeff_mul + 32768) >> 16;
        mean_b = (sum_b_box[col] * coeff_mul + 32768) >> 16;

        *Q++ = ((*I++) * mean_a + (mean_b << 8) + 127)>> 8;
    }

    for(; col < ctx->w; ++col){
        --wnd_w;
        wnd_size = wnd_w * wnd_h;
        coeff_mul = (1u << 16) / wnd_size;   
        
        mean_a = (sum_a_box[col] * coeff_mul + 32768) >> 16;
        mean_b = (sum_b_box[col] * coeff_mul + 32768) >> 16;
        *Q++ = ((*I++) * mean_a + (mean_b << 8) + 127) >> 8;
    }
    
    //TODO

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

    *Q++ = ((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 32);

    /*2 calc the left cols*/
    for(; col <= 2 * ctx->r ; ++col){
        ++wnd_w;
        wnd_size += wnd_h;

        sum_a += ctx->sum_line_A[col];
        sum_b += ctx->sum_line_B[col];

        coeff_mul = (1ull << 24)/ wnd_size;
        *Q++ = ((((*I++) * sum_a + (sum_b << 8)) * coeff_mul) >> 32);
		
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

	if (ctx->r <= 127)
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
    uint8_t prev_a;
    uint8_t prev_b;

    uint8_t *p_a, *p_b;
    uint8_t a, b;

    uint16_t *p_sum_col_a, *p_sum_col_b;
    uint64_t coeff_mul;
    uint16_t *u_line;
    uint32_t *delta_line;

	u_line = ctx->u_line_I;
	delta_line = ctx->delta_line_I;
	
    r = ctx->r;

    sum_col = ctx->sum_line_I;
    sum_col_sq = ctx->sum_sq_line_I;


    p_a 		= ctx->ringbuf_a + ctx->idx_ab_wnd * ctx->w;
    p_b 		= ctx->ringbuf_b + ctx->idx_ab_wnd * ctx->w;
    p_sum_col_a 	= ctx->sum_line_A;
    p_sum_col_b 	= ctx->sum_line_B;

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
        prev_a = p_a[col];
        prev_b = p_b[col];

        p_sum_col_a[col] += a - prev_a;
        p_sum_col_b[col] += b - prev_b;
        p_a[col] = a;
        p_b[col] = b;
        
    }



    ctx->idx_ab_wnd = (ctx->idx_ab_wnd + 1) % ctx->ab_wnd_size;
}


inline static void hg_guided_denoise_learn_ab_line_s(HG_GuidedDenoise_Context2 *ctx, int wnd_h)
{
    uint32_t sum_box;
    uint32_t sum_box_sq;
    int col , wnd_w , wnd_size ;
    uint16_t u;
    uint16_t delta;
    int r;
    uint16_t *sum_col;
    uint32_t *sum_col_sq; 
    uint8_t prev_a;
    uint8_t prev_b;

    uint8_t *p_a, *p_b;
    uint8_t a, b;

    uint16_t *p_sum_col_a, *p_sum_col_b;

    uint32_t coeff_mul;
    uint32_t sum_box_line[HG_GUIDED_DENOISE_FILTER_IMG_W_MAX];
    uint32_t sum_box_sq_line[HG_GUIDED_DENOISE_FILTER_IMG_W_MAX];    
    uint32_t *p_sum_box_sq;
    uint32_t *p_sum_box;
	
    r = ctx->r;

    sum_col = ctx->sum_line_I;
    sum_col_sq = ctx->sum_sq_line_I;

    /*calc the first col*/

    sum_box = 0;
    sum_box_sq = 0;

    p_sum_box = sum_box_line;
    p_sum_box_sq = sum_box_sq_line;
    

    for(col = 0 ; col <= r ; ++col){
        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];
    }

    *p_sum_box++ = sum_box;
    *p_sum_box_sq++ = sum_box_sq;
        
    /*calc the left*/
    for(; col <= 2 * r ; ++col){
        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];
        
        *p_sum_box++ = sum_box;
        *p_sum_box_sq++ = sum_box_sq;
    }

    /*calc the middle*/
    wnd_w	 = 2 * r + 1;    
    for(; col < ctx->w ; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        sum_box_sq += sum_col_sq[col] - sum_col_sq[col - wnd_w];

        *p_sum_box++ = sum_box;
        *p_sum_box_sq++ = sum_box_sq;

    }

    /*calc the right*/
    for( ; col <  ctx->w + r ; ++col){
        --wnd_w;
        sum_box	 -= sum_col[ctx->w - wnd_w - 1];
        sum_box_sq  -= sum_col_sq[ctx->w - wnd_w - 1];

        *p_sum_box++ = sum_box;
        *p_sum_box_sq++ = sum_box_sq;
    }

    p_sum_box = sum_box_line;
    p_sum_box_sq = sum_box_sq_line;
    
    p_a 	= ctx->ringbuf_a + ctx->idx_ab_wnd * ctx->w;
    p_b 	= ctx->ringbuf_b + ctx->idx_ab_wnd * ctx->w;    
    p_sum_col_a = ctx->sum_line_A;
    p_sum_col_b = ctx->sum_line_B;

    wnd_w = r + 1;
    p_sum_box = sum_box_line;
    p_sum_box_sq = sum_box_sq_line;

    for(col = 0; col < r; col++){
        wnd_size = wnd_w * wnd_h;
        coeff_mul = (1u<<16) / wnd_size;
    
        u = ((*p_sum_box) * coeff_mul) >> 8;
        delta = (*p_sum_box_sq *coeff_mul - u * u)>>16;

        a = ctx->delta_a_map[delta];
        b = ((255 - a) * u) >> 16;

        prev_a = *p_a;
        prev_b = *p_b;        

        *p_sum_col_a += a - prev_a;
        *p_sum_col_b += b - prev_b;

        *p_a++ = a;
        *p_b++ = b;
        ++p_sum_box;
        ++p_sum_box_sq;
        ++p_sum_col_a;
        ++p_sum_col_b;
        ++wnd_w;
    }
    
    wnd_w = 2 * r + 1;
    wnd_size = wnd_w * wnd_h;
    coeff_mul = (1u<<16) / wnd_size;
    
    for(; col < ctx->w && (col % 8 != 0); ++col){
        u = ((*p_sum_box) * coeff_mul) >> 8;
        delta = (*p_sum_box_sq *coeff_mul - u * u)>>16;

        a = ctx->delta_a_map[delta];
        b = ((255 - a) * u) >> 16;


        prev_a = *p_a;
        prev_b = *p_b;        
        
        *p_sum_col_a += a - prev_a;
        *p_sum_col_b += b - prev_b;

        *p_a++ = a;
        *p_b++ = b;
        ++p_sum_box;
        ++p_sum_box_sq;
        ++p_sum_col_a;
        ++p_sum_col_b;
    }

    if(1)
    {

        uint16x8_t  usq_u16x8, u_16x8, delta_16x8, delta_arg_16x8, axor_16x8, axorXu;
        uint32x4_t  delta_arg0_32x4, delta_arg1_32x4, *p_sum_box_sq_32x4, sum_box0_32x4, sum_box1_32x4, *p_sum_box_u32x4;
        uint16x4_t  delta_arg0_16x4, delta_arg1_16x4;
        uint8_t *p_a1;
        uint16_t *p_delta;
        uint8x8_t   mask_8x8, a_8x8, b_8x8;
        uint16x8_t  *p_sum_a_16x8, *p_sum_b_16x8, sum_a_16x8, sum_b_16x8;
        uint8x8_t   *p_a_8x8, *p_b_8x8;
        mask_8x8 = vcreate_u8(~0ull);
        
        p_sum_box_u32x4 = (uint32x4_t *) p_sum_box;
        p_sum_box_sq_32x4 = (uint32x4_t *) p_sum_box_sq;

        p_sum_a_16x8 = (uint16x8_t *)p_sum_col_a;
        p_sum_b_16x8 = (uint16x8_t *)p_sum_col_b;
        p_a_8x8      = (uint8x8_t *)p_a;
        p_b_8x8      = (uint8x8_t *)p_b;

        for(; col + 8 < ctx->w - r; col += 8){ 
            sum_box0_32x4 = vmulq_n_u32 (*p_sum_box_u32x4, coeff_mul); 
            sum_box1_32x4 = vmulq_n_u32 (*(p_sum_box_u32x4 + 1), coeff_mul);
            ((uint16x4_t *)&u_16x8)[0] = vrshrn_n_u32(sum_box0_32x4, 16);
            ((uint16x4_t *)&u_16x8)[1] = vrshrn_n_u32(sum_box1_32x4, 16);

            usq_u16x8 = vmulq_u16(u_16x8, u_16x8); // u * u 

            delta_arg0_32x4 = vmulq_n_u32 (*p_sum_box_sq_32x4, coeff_mul);
            delta_arg1_32x4 = vmulq_n_u32 (*(p_sum_box_sq_32x4 + 1), coeff_mul);

            delta_arg0_16x4 = vrshrn_n_u32(delta_arg0_32x4, 16); // >>16 and rounded value
            delta_arg1_16x4 = vrshrn_n_u32(delta_arg1_32x4, 16); // >>16 and rounded value
          
            ((uint16x4_t *)&delta_arg_16x8)[0] = delta_arg0_16x4; 
            ((uint16x4_t *)&delta_arg_16x8)[1] = delta_arg1_16x4;


            delta_16x8 = vqsubq_u16(delta_arg_16x8, usq_u16x8); // Mean(sq) - u * u,  saturating subtract
            p_delta = (uint16_t *)&delta_16x8;

            p_a1 = (uint8_t *)&a_8x8;
            p_a1[0] = ctx->delta_a_map[p_delta[0]]; p_a1[1] = ctx->delta_a_map[p_delta[1]]; p_a1[2] = ctx->delta_a_map[p_delta[2]]; p_a1[3] = ctx->delta_a_map[p_delta[3]];
            p_a1[4] = ctx->delta_a_map[p_delta[4]]; p_a1[5] = ctx->delta_a_map[p_delta[5]]; p_a1[6] = ctx->delta_a_map[p_delta[6]]; p_a1[7] = ctx->delta_a_map[p_delta[7]]; 
            
            axor_16x8 = vsubl_u8 (mask_8x8, a_8x8);// 255 - a
            axorXu = vmulq_u16 (axor_16x8, u_16x8); // (255 - a) * u
            b_8x8 = vrshrn_n_u16(axorXu, 8); // ((255 - a) * u )>> 8  , round 
            
            sum_a_16x8 = vaddw_u8(*p_sum_a_16x8, a_8x8);
            *p_sum_a_16x8 = vsubw_u8(sum_a_16x8, *p_a_8x8);
            
            sum_b_16x8 = vaddw_u8(*p_sum_b_16x8, b_8x8);
            *p_sum_b_16x8 = vsubw_u8(sum_b_16x8, *p_b_8x8);

            *p_a_8x8 = a_8x8;
            *p_b_8x8 = b_8x8;
            
            p_sum_box_u32x4 += 2;
            p_sum_box_sq_32x4 += 2;
            ++p_a_8x8;
            ++p_b_8x8;
            ++p_sum_a_16x8;
            ++p_sum_b_16x8;
        }

        
        p_sum_box = (uint32_t *)p_sum_box_u32x4;
        p_sum_box_sq = (uint32_t *)p_sum_box_sq_32x4;
        
        p_sum_col_a = (uint16_t *)p_sum_a_16x8;
        p_sum_col_b = (uint16_t *)p_sum_b_16x8;
        p_a     = (uint8_t *)p_a_8x8;
        p_b     = (uint8_t *)p_b_8x8;
    }

    for(; col < ctx->w - r; ++col){

        u = ((*p_sum_box) * coeff_mul + 32768) >> 16;
        delta = ((*p_sum_box_sq *coeff_mul + 32768 )>>16)  -  u * u;
       
        a = ctx->delta_a_map[delta];
        //printf("delta = %d, a = %d\n", delta, a);        
        b = ((255 - a) * u + 127) >> 8; // u -  ((a * u) >> 8)

        prev_a = *(p_a);
        prev_b = *(p_b);        
        
        *p_sum_col_a += a - prev_a;
        *p_sum_col_b += b - prev_b;
        
        *p_a++ = a;
        *p_b++ = b;
        ++p_sum_box;
        ++p_sum_box_sq;
        ++p_sum_col_a;
        ++p_sum_col_b;
    }

    for(; col < ctx->w; col++){
        --wnd_w;
        wnd_size = wnd_w * wnd_h;
        coeff_mul = (1u<<16) / wnd_size;
        
        u = ((*p_sum_box) * coeff_mul) >> 8;
        delta = (*p_sum_box_sq *coeff_mul - u * u)>>16;
        
        a = ctx->delta_a_map[delta];
        b = ((255 - a) * u + 127) >> 16;

        
        prev_a = *p_a;
        prev_b = *p_b;        
        
        *p_sum_col_a += a - prev_a;
        *p_sum_col_b += b - prev_b;
        
        *p_a++ = a;
        *p_b++ = b;
        ++p_sum_box;
        ++p_sum_box_sq;
        ++p_sum_col_a;
        ++p_sum_col_b;
    }    

    p_a 	= ctx->ringbuf_a + ctx->idx_ab_wnd * ctx->w;
    p_b 	= ctx->ringbuf_b + ctx->idx_ab_wnd * ctx->w;   

    ctx->idx_ab_wnd = (ctx->idx_ab_wnd + 1) % ctx->ab_wnd_size;
}

inline static void hg_guided_denoise_learn_ab_line(HG_FILTER_HANDLE h, int wnd_h)
{
    HG_GuidedDenoise_Context2 *ctx;

	ctx = (HG_GuidedDenoise_Context2 *)h;

	if (ctx->r <= 127)
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

    memset(ctx->ringbuf_a, 0, sizeof(uint8_t) * (2 * r + 1) * w);
    memset(ctx->ringbuf_b, 0, sizeof(uint8_t) * (2 * r + 1) * w);    
    
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
    int ab_wnd_size = ctx->ab_wnd_size - 1;

	ctx->idx_ab_wnd = (ctx->idx_ab_wnd + ctx->ab_wnd_size - 1) % ctx->ab_wnd_size;
    for(; row < ctx->h + 2 * r ; row++, --ab_wnd_size){
        uint8_t *p_prev_a = ctx->ringbuf_a + w * ((ctx->idx_ab_wnd + ctx->ab_wnd_size - ab_wnd_size) % ctx->ab_wnd_size);
        uint8_t *p_prev_b = ctx->ringbuf_b + w * ((ctx->idx_ab_wnd + ctx->ab_wnd_size - ab_wnd_size) % ctx->ab_wnd_size);        
        uint8_t prev_a, prev_b;
        
        for(col = 0; col < w; col++){
            prev_a = p_prev_a[col];
            prev_b = p_prev_b[col];            
            ctx->sum_line_A[col] -= prev_a;
            ctx->sum_line_B[col] -= prev_b;
        }
        hg_guided_denoise_output_q_line(h, ab_wnd_size,  pIP + (row - r * 2) * pitch_ip , pQ + (row - r * 2) * pitch_q);
    }

    return 0;

}


void    hg_guided_denoise_filter2_destroy(HG_FILTER_HANDLE h)
{
    HG_GuidedDenoise_Context2 *ctx;

    ctx = (HG_GuidedDenoise_Context2 *)h;

    hg_aligned_free(ctx->ringbuf_a);
    hg_aligned_free(ctx->ringbuf_b);    
    hg_aligned_free(ctx->sum_line_I);
    hg_aligned_free(ctx->sum_sq_line_I);
    hg_aligned_free(ctx->sum_line_A);
    hg_aligned_free(ctx->sum_line_B);
    hg_aligned_free(ctx->u_line_I);
    hg_aligned_free(ctx->delta_line_I);

    free(ctx);
}



#endif // __arm__

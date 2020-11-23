#include <malloc.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "hg_os.h"
#include "hg_integral_like_filter.h"
#include "hg_guided_denoise_filter.h"


typedef struct HG_GuidedFilter_Context
{
    int w;
    int h;
    int r;
    int         pitch;
    float       min_pcnt;
    float       max_pcnt;

    uint16_t    *buf_mean_IP;
    uint32_t    *buf_var_IP;
    uint8_t     *buf_a;
    uint8_t     *buf_b;
    uint16_t    *buf_mean_a;
    uint16_t    *buf_mean_b;
    uint8_t     delta_a_map[65536];
	float		lamda;
}HG_GuidedFilter_Context;


static int hg_guided_denoise_filter_learnAB(HG_GuidedFilter_Context *ctx)
{
    //const float lamda = ctx->lamda * 256;
    int col , row , max_row , max_col;
    int16_t a , b;
    uint8_t *p_a , *p_b;
    uint16_t  *p_mean_IP;
    uint32_t *p_var_I;
    int skip_a , skip_b , skip_var , skip_mean_IP;
    uint16_t v_uIP;
    uint32_t v_delta;
    //uint32_t lamda1;
    
    static double min_a = 1e10 , max_a = -1e10 , min_b = 1e10 , max_b = -1e10;

    p_a         = ctx->buf_a;
    p_b         = ctx->buf_b;
    p_mean_IP   = ctx->buf_mean_IP;
    p_var_I     = ctx->buf_var_IP;

    skip_a      = ctx->pitch - ctx->w;
    skip_b      = ctx->pitch - ctx->w;
    skip_var    = ctx->pitch - ctx->w;
    skip_mean_IP = ctx->pitch - ctx->w;

    max_row = ctx->h;
    max_col = ctx->w;

    //lamda1   = lamda * 255;
    for(row = 0 ; row < max_row ; ++row){
        for(col = 0 ; col < max_col ; ++col){
                v_uIP       = *p_mean_IP++;
                v_delta     = *p_var_I++;

                a           = ctx->delta_a_map[v_delta >> 8]; // (v_delta << 8) / (v_delta  + lamda1);
                b           = ((256 - a) * v_uIP) >> 16;

            if (0)
            {
                if (min_a > v_delta) min_a = v_delta;
                if (min_b > b) min_b = b;
                if (max_a < v_delta) max_a = v_delta;
                if (max_b < b) max_b = b;
            }

            *p_a++    = a; 
            *p_b++    = b; 
        }
        p_a      += skip_a;
        p_b      += skip_b;
        p_mean_IP += skip_mean_IP;
        p_var_I  += skip_var;
    }

    if (0)
    printf("a = [%.6f , %.6f] , b = [%.6f , %.6f]\n" , 
            min_a , max_a , min_b , max_b
          );

    hg_gray8_meanfilter(ctx->buf_a, ctx->w, ctx->h, ctx->pitch, ctx->r, ctx->buf_mean_a, ctx->pitch);
    hg_gray8_meanfilter(ctx->buf_b, ctx->w, ctx->h, ctx->pitch, ctx->r, ctx->buf_mean_b, ctx->pitch);

    if (0)
    {
        int ret;

        ret = hg_gray8_mean_filter_verify(ctx->buf_a, ctx->w, ctx->h, ctx->pitch, ctx->r, ctx->buf_mean_a, ctx->pitch);

        if (ret == 0)
            printf("vertify mean filter to a matrix ok !\n");

        ret = hg_gray8_mean_filter_verify(ctx->buf_b, ctx->w, ctx->h, ctx->pitch, ctx->r, ctx->buf_mean_b, ctx->pitch);

        if (ret == 0)
            printf("vertify mean filter to b matrix ok !\n");
    }

    return 0;
}

static int hg_guided_denoise_filter_output(HG_GuidedFilter_Context *ctx , uint8_t *p_I , int pitch_i, uint8_t *p_Q, int pitch_q)
{
    uint16_t *p_a , *p_b;
    int col , row;
    int skip_I , skip_Q , skip_a , skip_b;
    int32_t q_min , q_max;
    int max_row , max_col;
	int32_t a_U , b_U;

    max_row = ctx->h;
    max_col = ctx->w;

    p_a = ctx->buf_mean_a;
    p_b = ctx->buf_mean_b;
    
    skip_I = pitch_i - ctx->w;
    skip_Q = pitch_q - ctx->w;

    skip_a = ctx->pitch - ctx->w;
    skip_b = ctx->pitch - ctx->w;

    q_min = ctx->min_pcnt * 255;
    q_max = ctx->max_pcnt * 255;

    (void)q_min;
    (void)q_max;

    for(row = 0 ; row < max_row ; ++row){
        for(col = 0 ; col < max_col ; ++col){

            a_U =  (*p_a);
            b_U = (*p_b);
            *p_Q   =  ((((*p_I) * a_U) + (b_U << 8)) >> 16) + 8;

            /*if (q < q_min)
				q = q_min;
            else if (q > q_max) 
				q = q_max;*/
			
            //*p_Q = q;

            ++p_I ; ++p_Q ; ++p_a; ++p_b;
        }

        p_I += skip_I;
        p_Q += skip_Q;
        p_a += skip_a;
        p_b += skip_b;
    }

    return 0;
}


HG_FILTER_HANDLE    hg_guided_denoise_filter_create(int w , int h , int r, float lamda)
{
    HG_GuidedFilter_Context *ctx;

    ctx = (HG_GuidedFilter_Context *)calloc(1 , sizeof(HG_GuidedFilter_Context));
    assert(ctx != NULL);

    ctx->w = w;
    ctx->h = h;
    ctx->r = r;
    ctx->min_pcnt = 16.0/255;
    ctx->max_pcnt = 235.0/255;
	ctx->lamda = lamda;
    ctx->pitch = (ctx->w + 15) & (~15);

    ctx->buf_mean_IP = (uint16_t *)hg_aligned_malloc(16, ctx->pitch * ctx->h * sizeof(uint16_t));
    ctx->buf_var_IP = (uint32_t *)hg_aligned_malloc(16, ctx->pitch * ctx->h * sizeof(uint32_t));
    ctx->buf_a = (uint8_t *)hg_aligned_malloc(16, ctx->pitch * ctx->h * sizeof(uint8_t));
    ctx->buf_b = (uint8_t *)hg_aligned_malloc(16, ctx->pitch * ctx->h * sizeof(uint8_t));
    ctx->buf_mean_a = (uint16_t *)hg_aligned_malloc(16, ctx->pitch * ctx->h * sizeof(uint16_t));
    ctx->buf_mean_b = (uint16_t *)hg_aligned_malloc(16, ctx->pitch * ctx->h * sizeof(uint16_t));

    // use lookup table to speed up div
    {
        int delta = 0;
        int lamda_div = ctx->lamda * 65536 ;

        for(delta = 0; delta < 65536; delta++){
            ctx->delta_a_map[delta] = (delta << 16) / ((delta << 8) + lamda_div);
        }

    }
    return (HG_FILTER_HANDLE)ctx;
}


int     hg_guided_denoise_filter_process(HG_FILTER_HANDLE h , unsigned char * pIP, int pitch_ip, unsigned char *pQ, int pitch_q)
{
    HG_GuidedFilter_Context *ctx;
    unsigned int t0 , t1 , t2 , t3;
    //int ret;

    ctx = (HG_GuidedFilter_Context *)h;

    t0 = hg_gettime_ms();

    hg_gray8_mean_variance_filter(pIP, ctx->w, ctx->h, pitch_ip, ctx->r, ctx->buf_mean_IP, ctx->pitch, ctx->buf_var_IP, ctx->pitch);

#if 0
    ret = hg_gray8_mean_filter_verify(pIP, ctx->w, ctx->h, pitch_ip, ctx->r, ctx->buf_mean_IP, ctx->pitch);

    if (ret == 0)
        printf("vertify mean filter ok !\n");

    ret = hg_gray8_variance_verify(pIP, ctx->w, ctx->h, pitch_ip, ctx->r, ctx->buf_var_IP, ctx->pitch);
    if (ret == 0)
        printf("vertify variance  filter ok !\n");
#endif

    t1 = hg_gettime_ms();

    hg_guided_denoise_filter_learnAB(ctx);
    t2 = hg_gettime_ms();

    hg_guided_denoise_filter_output(ctx, pIP, pitch_ip, pQ, pitch_q);

    t3 = hg_gettime_ms();

    if (1)
    printf("hg_guided_denoise_filter_process : ==========time cost = %d %d %d , all = %d ms\n" , 
        (t1 - t0) , (t2 - t1) , (t3 - t2) , 
        (t3 - t0)
    );

    return 0;
}

#define CHECK_AND_FREE(a) do{ if ((a)) {free((a)) ; (a) = 0 ; } ; }while(0)

void    hg_guided_denoise_filter_destroy(HG_FILTER_HANDLE h)
{
    HG_GuidedFilter_Context *ctx;

    ctx = (HG_GuidedFilter_Context *)h;

    CHECK_AND_FREE(ctx->buf_mean_IP);
    CHECK_AND_FREE(ctx->buf_a);
    CHECK_AND_FREE(ctx->buf_b);
    CHECK_AND_FREE(ctx->buf_mean_a);
    CHECK_AND_FREE(ctx->buf_mean_b);
    CHECK_AND_FREE(ctx->buf_var_IP);
    CHECK_AND_FREE(ctx->buf_mean_IP);

    free(ctx);
}




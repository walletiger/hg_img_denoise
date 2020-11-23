#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>

#include "hg_os.h"
#include "hg_integral_like_filter.h"

/**
 * @brief hg_meanfilter_line_s
 * @param max_cols
 * @param p_dst
 * @param r <= 15
 * @param wnd_h
 * @param sum_col
 */
inline static void hg_meanfilter_line_s(int max_cols, uint16_t *p_dst,
                                     int r, int wnd_h, uint16_t *sum_col)
{
    uint32_t sum_box;
    int col , wnd_w ;
    uint32_t wnd_size;
    uint32_t coeff_mul;

    /*calc the first col*/
    wnd_w	= r + 1;
    sum_box = 0;

    /*calc the first pixel*/
    for(col = 0 ; col <= r ; ++col){
        sum_box += sum_col[col];
    }

    wnd_size = wnd_w * wnd_h;

    *p_dst++ = (sum_box << 8) / wnd_size;

    /*calc the left*/
    for(; col <= 2 * r ; ++col){
        ++wnd_w;
        wnd_size += wnd_h;
        sum_box += sum_col[col];
        *p_dst++ = (sum_box << 8) / wnd_size;
    }

    /*calc the middle*/
    coeff_mul = (1u << 24) / wnd_size;

    for( ; col < max_cols  ; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        *p_dst++ = ((sum_box * coeff_mul) >> 16);
    }

    /*calc the right*/
    for( ; col < max_cols + r ; ++col){
        --wnd_w;
        wnd_size -= wnd_h;
        sum_box -= sum_col[max_cols - 1 - wnd_w];
        *p_dst++ = (sum_box << 8) / wnd_size;
    }
}

/**
 * @brief hg_meanfilter_line_l
 * @param max_cols
 * @param p_dst
 * @param r > 15
 * @param wnd_h
 * @param sum_col
 */
inline static void hg_meanfilter_line_l(int max_cols, uint16_t *p_dst,
                                     int r, int wnd_h, uint16_t *sum_col)
{
    uint32_t sum_box;
    int col , wnd_w ;
    uint32_t wnd_size;
    uint64_t coeff_mul;

    /*calc the first col*/
    wnd_w	= r + 1;
    sum_box = 0;

    /*calc the first pixel*/
    for(col = 0 ; col <= r ; ++col){
        sum_box += sum_col[col];
    }

    wnd_size = wnd_w * wnd_h;

    *p_dst++ = (sum_box << 8) / wnd_size;

    /*calc the left*/
    for(; col <= 2 * r ; ++col){
        ++wnd_w;
        wnd_size += wnd_h;
        sum_box += sum_col[col];
        *p_dst++ = (sum_box << 8) / wnd_size;
    }

    /*calc the middle*/
    coeff_mul = (1ull << 32) / wnd_size;

    for( ; col < max_cols  ; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        *p_dst++ = ((sum_box * coeff_mul) >> 24);
    }

    /*calc the right*/
    for( ; col < max_cols + r ; ++col){
        --wnd_w;
        wnd_size -= wnd_h;

        sum_box -= sum_col[max_cols - wnd_w - 1];
        *p_dst++ = (sum_box << 8) / wnd_size;
    }
}

inline static void hg_meanfilter_line(int max_cols, uint16_t *p_dst,
                                         int r, int wnd_h, uint16_t *sum_col)
{
    if (r > 15)
        hg_meanfilter_line_l(max_cols, p_dst, r, wnd_h, sum_col);
    else
        hg_meanfilter_line_s(max_cols, p_dst, r, wnd_h, sum_col);
}

/**
 * @brief hg_gray8_meanfilter
 * @param p_src
 * @param width
 * @param height
 * @param src_pitch
 * @param r
 * @param p_dst
 * @param dst_pitch
 * @return
 */
int hg_gray8_meanfilter(const uint8_t *p_src, int width, int height, int src_pitch, int r, uint16_t *p_dst, int dst_pitch)
{
    const uint8_t *p_last;
    int col , row;
    int  wnd_h;
    uint16_t sum_col[HG_MEAN_FILTER_MAX_WIDTH];

    int offset;

    if (r > HG_MEAN_FILTER_MAX_R)
        return -1;

    if (width < 2 * r + 1 || height < 2 * r + 1)
        return -1;

    memset(sum_col , 0 , sizeof(uint16_t) * width);

    /*calc the first row*/
    for(row = 0 ; row <= r ; ++row , p_src += src_pitch){
        for(col = 0 ; col < width ; ++col){
            sum_col[col] += p_src[col];
        }
    }
    wnd_h = r + 1;

    hg_meanfilter_line(width , p_dst , r , wnd_h , sum_col);
    p_dst += dst_pitch;

    /*calc the top*/
    for( ; row <= 2 * r ; ++row , p_src += src_pitch , p_dst += dst_pitch){
        ++wnd_h;

        for(col = 0 ; col < width ; ++col){
            sum_col[col] += p_src[col];
        }
        hg_meanfilter_line(width , p_dst , r , wnd_h , sum_col);
    }

    /*calc vertical middle*/
    offset = src_pitch * wnd_h;

    for(; row < height ; ++row , p_src += src_pitch , p_dst += dst_pitch){
        for(col = 0 ; col < width ; ++col){
            sum_col[col] += p_src[col] - (p_src - offset)[col];
        }
        hg_meanfilter_line(width , p_dst , r , wnd_h , sum_col);
    }

    /*calc the bottom*/
    p_last = p_src - src_pitch;
    for(; row < height + r ; ++row , p_src += src_pitch , p_dst += dst_pitch){
        const uint8_t *p_prev;
        --wnd_h;
        offset = src_pitch * wnd_h;

        p_prev = p_last - offset;
        for(col = 0 ; col < width ; ++col){
            sum_col[col] -= p_prev[col];
        }
        hg_meanfilter_line(width , p_dst , r , wnd_h , sum_col);
    }
    return 0;
}

/**
 * @brief hg_mean_and_variance_line_s
 * @param max_cols
 * @param p_dst_mean
 * @param p_dst_var
 * @param r
 * @param wnd_h
 * @param sum_col
 * @param sum_col_sq
 */
inline static void hg_mean_and_variance_line_s(int max_cols , uint16_t *p_dst_mean , uint32_t *p_dst_var ,
                                     int r , int wnd_h , uint16_t *sum_col , uint32_t *sum_col_sq)
{
    uint32_t sum_box , sum_box_sq;
    int col , wnd_w , wnd_size ;
    uint32_t u;
    uint32_t coeff_mul;
    register uint32_t coeff_mul_h16, coeff_mul_l8;

    /*calc the first col*/
    wnd_w   = r + 1;
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

    u			 = ((sum_box * coeff_mul) >> 16);
    *p_dst_mean++ = u;
    *p_dst_var++ = (((sum_box_sq * coeff_mul_h16) + ((sum_box_sq >> 8) * coeff_mul_l8) - (u * u)) >> 8);

    /*calc the left*/
    for(; col <= 2 * r ; ++col){
        ++wnd_w;
        wnd_size += wnd_h;

        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];

        coeff_mul = (1u << 24) / wnd_size;
        coeff_mul_h16 = (coeff_mul >> 8);
        coeff_mul_l8 = (coeff_mul & 0xffu);

        u			 = ((sum_box * coeff_mul) >> 16);
        *p_dst_mean++ = u;
        *p_dst_var++ = (((sum_box_sq * coeff_mul_h16) + ((sum_box_sq >> 8) * coeff_mul_l8) - (u * u)) >> 8);

    }

    /*calc the middle*/
    wnd_size        = wnd_w * wnd_h;
    coeff_mul       = (1u << 24) / wnd_size;
    coeff_mul_h16   = (coeff_mul >> 8);
    coeff_mul_l8    = (coeff_mul & 0xffu);

    for( ; col < max_cols  ; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        sum_box_sq += sum_col_sq[col] - sum_col_sq[col - wnd_w];

        u			 = ((sum_box * coeff_mul) >> 16);
        *p_dst_mean++ = u;
        *p_dst_var++ = (((sum_box_sq * coeff_mul_h16) + ((sum_box_sq >> 8) * coeff_mul_l8) - (u * u)) >> 8);
    }

    /*calc the right*/
    for( ; col < max_cols + r ; ++col){
        --wnd_w;
        wnd_size    -= wnd_h;
        coeff_mul       = (1u << 24) / wnd_size;
        coeff_mul_h16   = (coeff_mul >> 8);
        coeff_mul_l8    = (coeff_mul & 0xffu);

        sum_box     -= sum_col[max_cols - wnd_w - 1];
        sum_box_sq  -= sum_col_sq[max_cols - wnd_w - 1];

        u			 = ((sum_box * coeff_mul) >> 16);
        *p_dst_mean++ = u;
        *p_dst_var++ = (((sum_box_sq * coeff_mul_h16) + ((sum_box_sq >> 8) * coeff_mul_l8) - (u * u)) >> 8);
     }
}

/**
 * @brief hg_mean_and_variance_line_l
 * @param max_cols
 * @param p_dst_mean
 * @param p_dst_var
 * @param r
 * @param wnd_h
 * @param sum_col
 * @param sum_col_sq
 */
inline static void hg_mean_and_variance_line_l(int max_cols , uint16_t *p_dst_mean , uint32_t *p_dst_var ,
                                     int r , int wnd_h , uint16_t *sum_col , uint32_t *sum_col_sq)
{
    uint32_t sum_box , sum_box_sq;
    int col , wnd_w , wnd_size ;
    uint32_t u; // must u32!!!
    uint64_t coeff_mul;

    /*calc the first col*/
    wnd_w   = r + 1;
    sum_box = 0;
    sum_box_sq = 0;

    for(col = 0 ; col <= r ; ++col){
        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];
    }

    wnd_size = wnd_w * wnd_h;
    coeff_mul = (1ull << 32) / wnd_size;

    u			  = ((sum_box * coeff_mul) >> 24);
    *p_dst_mean++ = u;
    *p_dst_var++  = ((((sum_box_sq * coeff_mul) >> 16) - (u * u)) >> 8);

    /*calc the left*/
    for(; col <= 2 * r ; ++col){
        ++wnd_w;
        wnd_size += wnd_h;

        sum_box += sum_col[col];
        sum_box_sq += sum_col_sq[col];

        coeff_mul = (1ull << 32) / wnd_size;

        u            = ((sum_box * coeff_mul) >> 24);
        *p_dst_mean++ = u;
        *p_dst_var++ = ((((sum_box_sq * coeff_mul) >> 16) - (u * u)) >> 8);

    }

    /*calc the middle*/
    wnd_size = wnd_w * wnd_h;
    coeff_mul    = (1ull << 32) / wnd_size;

    for( ; col < max_cols  ; ++col){
        sum_box += sum_col[col] - sum_col[col - wnd_w];
        sum_box_sq += sum_col_sq[col] - sum_col_sq[col - wnd_w];

        u            = ((sum_box * coeff_mul) >> 24);
        *p_dst_mean++ = u;
        *p_dst_var++ = ((((sum_box_sq * coeff_mul) >> 16) - (u * u)) >> 8);
    }

    /*calc the right*/
    for( ; col < max_cols + r ; ++col){
        --wnd_w;
        wnd_size    -= wnd_h;
        coeff_mul    = (1ull << 32) / wnd_size;

        sum_box     -= sum_col[max_cols - wnd_w - 1];
        sum_box_sq  -= sum_col_sq[max_cols - wnd_w - 1];

        u   = ((sum_box * coeff_mul) >> 24);
        *p_dst_mean++ = u;
        *p_dst_var++ =  ((((sum_box_sq * coeff_mul) >> 16) - (u * u)) >> 8);
     }
}

inline static void hg_mean_and_variance_line(int max_cols , uint16_t *p_dst_mean , uint32_t *p_dst_var ,
                                     int r , int wnd_h , uint16_t *sum_col , uint32_t *sum_col_sq)
{

    if (r <= 15)
        hg_mean_and_variance_line_s(max_cols, p_dst_mean, p_dst_var, r, wnd_h, sum_col, sum_col_sq);
    else
        hg_mean_and_variance_line_l(max_cols, p_dst_mean, p_dst_var, r, wnd_h, sum_col, sum_col_sq);
}

/**
 * @brief hg_gray8_mean_variance_filter
 * @param p_src
 * @param width
 * @param height
 * @param src_pitch
 * @param r
 * @param p_dst_mean output mean matrix , (mean << 8)
 * @param pitch_mean
 * @param p_dst_variance output variance matrix , (variance << 8)
 * @param pitch_variance
 * @return
 */
int hg_gray8_mean_variance_filter(const uint8_t *p_src, int width, int height, int src_pitch, int r, uint16_t *p_dst_mean, int pitch_mean, uint32_t *p_dst_variance, int pitch_variance)
{
    uint8_t v , v1;
    const uint8_t * p_last;
    int col , row;
    int  wnd_h;
    uint16_t sum_col[HG_MEAN_FILTER_MAX_WIDTH];
    uint32_t sum_col_sq[HG_MEAN_FILTER_MAX_WIDTH];
    int offset;


    if (r > HG_MEAN_FILTER_MAX_R)
        return -1;

    if (width < 2 * r + 1 || height < 2 * r + 1)
        return -1;

    memset(sum_col, 0 , sizeof(uint16_t) * width);
    memset(sum_col_sq , 0 , sizeof(uint32_t) * width);

    /*calc the first row*/
    for(row = 0 ; row <= r ; ++row , p_src += src_pitch){
        for(col = 0 ; col < width; ++col){
            v = p_src[col];
            sum_col[col] += v;
            sum_col_sq[col] += v * v;
        }
    }
    wnd_h = r + 1;

    hg_mean_and_variance_line(width , p_dst_mean , p_dst_variance ,  r , wnd_h , sum_col , sum_col_sq);

    p_dst_variance += pitch_variance;
    p_dst_mean += pitch_mean;

    /*calc the top*/
    for( ; row <= 2 * r ; ++row , p_src += src_pitch , p_dst_variance += pitch_variance , p_dst_mean += pitch_mean){
        ++wnd_h;

        for(col = 0 ; col < width ; ++col){
            v = p_src[col];
            sum_col[col] += v;
            sum_col_sq[col] += v * v;

        }
        hg_mean_and_variance_line(width , p_dst_mean , p_dst_variance ,  r , wnd_h , sum_col , sum_col_sq);
    }

    /*calc vertical middle*/
    offset = src_pitch * wnd_h;

    for(; row < height ; ++row , p_src += src_pitch ,  p_dst_variance += pitch_variance , p_dst_mean += pitch_mean){
        for(col = 0 ; col < width ; ++col){
            v =  p_src[col];
            v1 = (p_src - offset)[col];
            sum_col[col] +=  v - v1;
            sum_col_sq[col] += v * v - v1 * v1;
        }
        hg_mean_and_variance_line(width , p_dst_mean , p_dst_variance ,  r , wnd_h , sum_col , sum_col_sq);
    }

    /*calc the bottom*/

    p_last = p_src - src_pitch;

    for(; row < height + r ; ++row,  p_dst_variance += pitch_variance , p_dst_mean += pitch_mean){
        const uint8_t *p_prev;

        --wnd_h;
        offset = src_pitch * wnd_h;
        p_prev = p_last - offset;

        for(col = 0 ; col < width ; ++col){
            v = p_prev[col];
            sum_col[col] -= v;
            sum_col_sq[col] -= v * v;
        }

        hg_mean_and_variance_line(width, p_dst_mean, p_dst_variance,  r, wnd_h, sum_col, sum_col_sq);
    }


    return 0;
}

/**
 * @brief hg_gray8_mean_filter_verify
 * @param p_src
 * @param width
 * @param height
 * @param src_pitch
 * @param r
 * @param p_dst
 * @param dst_pitch
 * @return
 */
int hg_gray8_mean_filter_verify(const uint8_t *p_src, int width, int height, int src_pitch, int r , uint16_t *p_dst, int dst_pitch)
{
    int col , row , i , j;
    uint32_t sum , npix;
    int diff;
    double u1, u2;

    for(row = 0 ; row < height ; ++row){
        for(col = 0 ; col < width ; ++col){
            sum = 0;
            npix = 0;
            for(i = -r ; i <= r ; ++i){
                if (i + row < 0 || i + row >= height)
                    continue;

                for(j = -r ; j <= r ; ++j){
                    if (j + col < 0 || j + col >= width)
                        continue;

                    sum += p_src[ (i + row) * src_pitch + j + col];
                    npix++;
                }
            }

            u1 = (sum * 1.0) / npix;
            u2 = (p_dst[row * dst_pitch + col] / 256.0);
            diff =  u1 - u2 ;
            if (diff > 2 || diff < -2){
                printf("Mean Not identfy at <%d,%d> , val = %.3f (correct), %.3f, wnd_size = %d\n" ,
                       row , col, u1, u2, npix
                       );
                return -1;
            }
        }
    }
    return 0;
}

/**
 * @brief hg_gray8_variance_verify
 * @param p_src
 * @param width
 * @param height
 * @param src_pitch
 * @param r
 * @param p_dst
 * @param dst_pitch
 * @return
 */
int hg_gray8_variance_verify(const uint8_t *p_src, int width, int height, int src_pitch, int r, uint32_t *p_dst, int dst_pitch)
{
    int col , row , i , j;
    uint8_t v;
    double u;
    int max_diff_abs = 0, diff_abs;
    uint64_t sum_sq;
    uint32_t sum, npix;

    double variance;

    int diff;

    (void)sum_sq;

    for(row = 0 ; row < height ; ++row){
        for(col = 0 ; col < width ; ++col){
            sum = 0;
            npix = 0;
            sum_sq = 0;
#if 1
            for(i = -r ; i <= r ; ++i){
                if (i + row < 0 || i + row >= height)
                    continue;

                for(j = -r ; j <= r ; ++j){
                    if (j + col < 0 || j + col >= width)
                        continue;

                    v = p_src[ (i + row) * src_pitch + j + col];
                    sum += v;
                    sum_sq += v * v;
                    npix++;
                }
            }

            u = ((double)sum) / npix;
            variance = (double)sum_sq / npix - u * u;
#else
            for(i = -r ; i <= r ; ++i){
                if (i + row < 0 || i + row >= height)
                    continue;

                for(j = -r ; j <= r ; ++j){
                    if (j + col < 0 || j + col >= width)
                        continue;

                    v = p_src[ (i + row) * src_pitch + j + col];
                    sum += v;
                    npix++;
                }
            }

            u = ((double)sum) / npix;
            variance = 0;
            for(i = -r ; i <= r ; ++i){
                if (i + row < 0 || i + row >= height)
                    continue;

                for(j = -r ; j <= r ; ++j){
                    if (j + col < 0 || j + col >= width)
                        continue;

                    v = p_src[ (i + row) * src_pitch + j + col];

                   variance += (u - v) * (u - v);
                }
            }

            variance /= npix;

#endif
            diff = (int)(variance - p_dst[row * dst_pitch + col] / 256.0);

            if (diff > 4 || diff < -4){

                diff_abs = fabs(diff);

                if (diff_abs > max_diff_abs){
                    max_diff_abs = diff_abs;

                    printf("Variance Not identfy at <%d,%d> , wnd_size = %d, diff = %d , u = %.3f , sum = %d, sum_sq = %ld , variance = %.3f (correct), %.3f\n" ,
                           row, col, npix, diff ,  u , sum , sum_sq,  (variance) , p_dst[row * dst_pitch + col] / 256.0

                           );
                }

                if (diff > 36 ||diff < -36)
                    return -1;
            }
        }
    }
    return 0;
}

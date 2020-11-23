#ifndef __HG_CLAHE_PRIV_H
#define __HG_CLAHE_PRIV_H

#define COLOR_AUTOLEVEL_CLIP (0.01f)
#define PIX_LOW_THRESH  (80)
#define PIX_MAX_THRESH  (130)

#define LEARN_RATIO     (0.85)
#define LEARN_FRAMES    (100)

static int PixMap_AutoLevel(uint32_t *bin , uint8_t *pix_map)
{
    int pix_max , pix_min ;
    uint32_t clip_thresh , sum_all , sum;
    int i;

    //clip the lowest and the highest 0.005 pixels vals
    sum_all = 0;
    for(i = 0 ; i < HIST_MAX_BIN ; ++i){
        sum_all += bin[i];
    }

    clip_thresh = sum_all * COLOR_AUTOLEVEL_CLIP;
    sum = 0;
    pix_min = 0 ; pix_max = HIST_MAX_BIN - 1;
    for(i = HIST_MAX_BIN - 1 ; i >= 0 ; --i){
        sum += bin[i];
        if (sum >= clip_thresh){
            pix_max = i < (HIST_MAX_BIN - 1) ? i + 1 : i;
            break;
        }
    }
    sum = 0;
    for(i = 0 ; i < HIST_MAX_BIN ; ++i){
        sum += bin[i];
        if (sum >= clip_thresh){
            pix_min = i > 0 ? i - 1 : 0;
            break;
        }
    }

    if (pix_min  > PIX_LOW_THRESH) pix_min = PIX_LOW_THRESH;
    if (pix_max < PIX_MAX_THRESH)  pix_max = PIX_MAX_THRESH;

    //printf("pix_min = %d , pix_max = %d , clip_thresh = %d , sum = %d\n" , pix_min , pix_max , clip_thresh , sum);

    for(i = 0 ; i <= pix_min ; ++i){
        bin[i] = 0;
        pix_map[i] = 0;
    }

    for( ; i < pix_max ; ++i){
        pix_map[i] = (i - pix_min) * (HIST_MAX_BIN - 1) / (pix_max - pix_min);
    }

    for( ; i < HIST_MAX_BIN ; ++i){
        pix_map[i] = HIST_MAX_BIN -1;
    }

    return 0;
}

static void MapPixelHistEq(uint32_t *bin , uint8_t *pix_map)
{
    int i;
    uint32_t sum_all , cum_sum;

    sum_all = 0;
    for(i = 0 ; i < HIST_MAX_BIN ; ++i){
        sum_all += bin[i];
    }

    cum_sum = 0;
    for(i = 0 ;  i < HIST_MAX_BIN ; ++i){
        cum_sum += bin[i];
        pix_map[i] = (((float)cum_sum) / sum_all) * (HIST_MAX_BIN - 1);
    }
}

static void CropHistogram(uint32_t *bin , float crop_peak_limit)
{
    int i , j , n_excess;
    int n_ave_dis;
    int i_step;
    uint32_t peak_limit_min , peak_limit;
    uint32_t sum_all;

    sum_all = 0;
    for(i = 0 ; i < HIST_MAX_BIN ; ++i){
        sum_all += bin[i];
    }

    peak_limit_min = (int)(0.5f + ((float)sum_all) / HIST_MAX_BIN);

    /*crop all pixel bin values that excess peak_limit*/
    peak_limit = peak_limit_min + crop_peak_limit * (sum_all - peak_limit_min);

    n_excess = 0;
    /**calc excess pixels*/
    for(i = 0 ; i < HIST_MAX_BIN ; ++i){
        if (bin[i] > peak_limit)    n_excess += (bin[i] - peak_limit);
    }

    /** clip and averagely redistribute excess pixels int each bin*/
    n_ave_dis = (int)(0.5f + ((float)(n_excess)) / HIST_MAX_BIN);

    for(i = 0 ; i < HIST_MAX_BIN ; ++i){
        if (bin[i] > peak_limit)    bin[i] = peak_limit;
        else if (bin[i] + n_ave_dis > peak_limit){
            n_excess -= peak_limit - bin[i];
            bin[i] = peak_limit;
        }else {
            n_excess -= n_ave_dis;
            bin[i] += n_ave_dis;
        }
    }
    /*averagely redistribute the left*/

    i = 0;

    while (n_excess > 0) {
        i_step = HIST_MAX_BIN / n_excess;
        if (i_step < 1)
            i_step = 1;

        for(j = i ; n_excess > 0 && j < HIST_MAX_BIN ; j += i_step){
            if (bin[j] < peak_limit){
                ++bin[j];
                --n_excess;
            }
        }

        i = i + 1;
        if (i >= HIST_MAX_BIN)  i = 0;
    }
}

#endif

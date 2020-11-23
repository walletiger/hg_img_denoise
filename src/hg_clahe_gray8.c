#include <string.h>
#include <stdlib.h>

#include <hg_os.h>
#include <hg_clahe.h>
#include "hg_clahe_priv.h"

static void BilinearInterPoSub(uint8_t *plane, int pitch, struct HG_Rect *rect ,
                              int U , int D , int L , int R , uint8_t * p_mapUL , 
                              uint8_t *p_mapUR , uint8_t *p_mapDL , uint8_t * p_mapDR)
{
    uint8_t v;
    uint8_t *p_pix;
    int skip;
    int col , row;
    int end_row , end_col;
    int32_t vx0 , vx1 , vy0 , vy1;
    int step_w , step_h;
    uint64_t coeff_step_w, coeff_step_h;

    p_pix   = plane + rect->start_y * pitch + rect->start_x;
    skip    = pitch - rect->cols;
    end_row = rect->start_y + rect->rows - 1;
    end_col = rect->start_x + rect->cols - 1;

    step_w  = (R - L + 1);
    step_h  = (D - U + 1);

    coeff_step_w =(1ull << 24) / step_w;
    coeff_step_h =(1ull << 24) / step_h;

   // printf("step_w = %d , step_h = %d\n" , step_w , step_h);

    for(row = rect->start_y ; row <= end_row  ; ++row){
        for(col = rect->start_x ; col <= end_col ; ++col){
            v       = *p_pix;
                
            vx0     = p_mapUL[v] ; vx1 = p_mapUR[v];
            //vy0     = vx0 + (vx1 - vx0) * (col - L) / step_w;
            vy0     = vx0 + (((vx1 - vx0) * (col - L) * coeff_step_w) >> 24);
            vx0     = p_mapDL[v] ; vx1 = p_mapDR[v];
            vy1     = vx0 + (((vx1 - vx0) * (col - L) * coeff_step_w) >> 24);
            *p_pix  = vy0 + (((vy1 - vy0) * (row - U) * coeff_step_h) >> 24);

            p_pix++;
        }
        p_pix += skip; 
    }
}


static void BilinearInterPo(uint8_t *plane, int w, int h, int pitch , int block_w , int block_h , uint8_t  (*pix_map)[CLAHE_NBLOCK_H][HIST_MAX_BIN])
{

    uint8_t  *p_mapUL , *p_mapUR , *p_mapDL , *p_mapDR;
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
            rect_proc.rows    = h - rect_proc.start_y;
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
                rect_proc.cols    = w - rect_proc.start_x;
            }else {
                L_Block = col_block - 1;
                R_Block = col_block;
            
                L = L_Block * block_w + half_block_w;
                R = R_Block * block_w + half_block_w;

                rect_proc.start_x = L; 
                rect_proc.cols    = R - L;
            }

            p_mapUL = pix_map[U_Block][L_Block];
            p_mapUR = pix_map[U_Block][R_Block];
            p_mapDL = pix_map[D_Block][L_Block];
            p_mapDR = pix_map[D_Block][R_Block];

            BilinearInterPoSub(plane, pitch,  &rect_proc , U , D , L , R ,
                p_mapUL , p_mapUR , p_mapDL , p_mapDR);

        }
    } 
}

int hg_gray8_clahe(uint8_t *plane, int w, int h, int pitch, float hist_peak_crop)
{
    int32_t block_w , block_h;
    uint8_t *p_pix , gray;
    int skip;
    int idx_block_h;

    uint32_t bin[CLAHE_NBLOCK_W][CLAHE_NBLOCK_H][HIST_MAX_BIN];
    uint8_t  pix_map[CLAHE_NBLOCK_W][CLAHE_NBLOCK_H][HIST_MAX_BIN];

    int col, row;

    if (w % CLAHE_NBLOCK_W != 0)
        return HG_ERR_INVALID_ARG;

    if (h % CLAHE_NBLOCK_H != 0)
        return HG_ERR_INVALID_ARG;

    memset(bin , 0 , sizeof(bin));

    HG_DEBUG("hg_gray8_clahe , size = %d , %d\n" ,
        w , h
    );

    block_w     = w / CLAHE_NBLOCK_W;
    block_h     = h / CLAHE_NBLOCK_H;

    /*calc hist bin*/
    p_pix   = plane;
    skip    = pitch - w;

    for(row = 0 ; row < h ; ++row){
        idx_block_h = row / block_h;
        for(col = 0 ; col < w ; ++col){
            gray = *p_pix++;
            ++bin[idx_block_h][col / block_w][gray];
        }
        p_pix += skip;
    }

    /*crop hist bin*/
    for(row = 0 ; row < CLAHE_NBLOCK_H ; ++row){
        for(col = 0 ; col < CLAHE_NBLOCK_W ; ++col){
            if (0){
                PixMap_AutoLevel(bin[row][col] , pix_map[row][col]);
            }else {
                CropHistogram(bin[row][col] , hist_peak_crop);
                MapPixelHistEq(bin[row][col] , pix_map[row][col]);
            }
        }
    }
    /*output image*/
    BilinearInterPo(plane, w, h, pitch, block_w , block_h , pix_map);

    return 0;
}


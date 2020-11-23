#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <hg_os.h>
#include "hg_clahe.h"
#include "hg_guided_denoise_filter.h"
#include "hg_local_variance_denoise.h"

int main(int ac, char *av[])
{
	FILE *fp_r, *fp_w;
	unsigned char  *buf_in, *buf_out;

	HG_FILTER_HANDLE h_filter;
	int ret;
	int size;
	int filter_type;
	int r;
	float lamda;

    unsigned int t0, t1, t2;
	
	if (ac < 6){
		printf("usage: %s filter_type(0-1) radius(1~127)  lamda (0 ~ 10.0) in_yuv out_yuv (in 720)\n", av[0]);
		return 0;
	}

	filter_type = atoi(av[1]);
	r = atoi(av[2]);
	lamda = atof(av[3]);
	
	fp_r = fopen(av[4], "rb");
	if (!fp_r)
		return -1;

	fp_w = fopen(av[5], "wb");

	if (!fp_w)
		return -1;

	if (r < 1)
		r = 1;
	else if (r > HG_GUIDED_DENOISE_FILTER_R_MAX)
		r = HG_GUIDED_DENOISE_FILTER_R_MAX;

	if (lamda < 1e-4)
		lamda = 1e-4;
	else if (lamda > 1e4)
		lamda = 1e4;
	

	buf_in = (unsigned char *)hg_aligned_malloc(16, 1280 * 720 * 3 / 2);
	buf_out = (unsigned char *)hg_aligned_malloc(16, 1280 * 720 * 3 / 2);

	if (r > 15 &&0 == filter_type){
        printf("denoise filter 1 only support max radius of 15, switch to filter 1\n");
        filter_type = 0;
	}

	if (0 == filter_type)
		h_filter = hg_guided_denoise_filter_create(1280, 720, r, lamda);
	else if (1 == filter_type)
		h_filter = hg_guided_denoise_filter2_create(1280, 720, r, lamda);
    else if (2 == filter_type)
		h_filter = hg_local_variance_denoise_filter_create(1280, 720, r, lamda);
    else {
        return -1;
    }

	size = 1280 * 720 * 3 / 2;

	while(1){
        ret = fread(buf_in, 1, size, fp_r);
        if (ret < size)
            break;
        
	    t0 = hg_gettime_ms();
        //hg_gray8_clahe(buf_in, 1280, 720, 1280, 0.001);
        t1 = hg_gettime_ms();

		if (0 == filter_type)
			hg_guided_denoise_filter_process(h_filter, buf_in, 1280, buf_out, 1280);
		else if (1 == filter_type)
			hg_guided_denoise_filter2_process(h_filter, buf_in, 1280, buf_out, 1280);
        else if (2 == filter_type)
            hg_local_variance_denoise_filter_process(h_filter, buf_in, 1280, buf_out, 1280);

        t2 = hg_gettime_ms();

        printf("delay = %d ms, %d\n", (t1 - t0), (t2 - t1));

		memcpy(buf_out + 1280 * 720 , buf_in + 1280 * 720, 1280 * 720 / 2);

		if (1){
			fwrite(buf_out, 1, 1280 * 720 * 3 / 2 , fp_w);
		}else {
			/*Y*/
			fwrite(buf_out, 1, 1280 * 360, fp_w);
			fwrite(buf_in + 1280 * 360 , 1, 1280 * 360, fp_w);
			
			/*U*/
			fwrite(buf_out + 1280 * 720, 1, 640 * 180, fp_w);
			fwrite(buf_in + 1280 * 720 + 640 * 180, 1, 640 * 180, fp_w);
			
			/*V*/
			fwrite(buf_out + 1280 * 720 * 5 / 4, 1, 640 * 180, fp_w);
			fwrite(buf_in + 1280 * 720 * 5 / 4 + 640 * 180, 1, 640 * 180, fp_w);
		}
	}

	fclose(fp_r);
	fclose(fp_w);

	hg_aligned_free(buf_in);
	hg_aligned_free(buf_out);

	if (0 == filter_type)
		hg_guided_denoise_filter_destroy(h_filter);
	else if (1 == filter_type)
		hg_guided_denoise_filter2_destroy(h_filter);
    else 
	    hg_local_variance_denoise_filter_destroy(h_filter);
	return 0;
	
}




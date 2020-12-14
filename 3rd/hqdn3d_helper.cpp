#include "hqdn3d.h"
#include "hqdn3d_helper.h"
#include "hg_denoise_priv.h"

typedef struct  
{
	int w;
	int h;
	int i_coeff_changed;

	float luma_spat;
	float luma_temp;
	float chroma_spat;
	float chroma_temp;
	vf_priv_s cfg;

}HQDN3D_Context;

typedef void * HQDN3D_HANDLE;

HQDN3D_HANDLE hqdn3d_create(int w, int h, float luma_spat, float luma_temp, float chroma_spat, float chroma_temp)
{
	HQDN3D_Context *ctx;

	ctx = (HQDN3D_Context *)calloc(1, sizeof(HQDN3D_Context));
	ctx->w = w;
	ctx->h = h;
	ctx->cfg.Line = (unsigned int  *)hg_aligned_malloc(16, sizeof(unsigned int) * w);

	if (!ctx->cfg.Line){
		free(ctx);
		return NULL;
	}
	ctx->luma_spat = luma_spat;
	ctx->luma_temp = luma_temp;
	ctx->chroma_spat = chroma_spat;
	ctx->chroma_temp = chroma_temp;
	ctx->i_coeff_changed = 1;

	return (HQDN3D_Context *)ctx;
}


void hqdn3d_set_coeff(HQDN3D_HANDLE h_filter, float luma_spat, float luma_temp, float chroma_spat, float chroma_temp)
{
	HQDN3D_Context *ctx;

	ctx = (HQDN3D_Context *)h_filter;

	ctx->luma_spat = luma_spat;
	ctx->luma_temp = luma_temp;
	ctx->chroma_spat = chroma_spat;
	ctx->chroma_temp = chroma_temp;

	ctx->i_coeff_changed = 1;
}


void hqdn3d_process_Y(HQDN3D_HANDLE h_filter, unsigned char *src_Y, int pitch_src, unsigned char *dst_Y, int pitch_dst)
{
	HQDN3D_Context *ctx;

	ctx = (HQDN3D_Context *)h_filter;

	if (ctx->i_coeff_changed){
		PrecalcCoefs(ctx->cfg.Coefs[0],  ctx->luma_spat);
		PrecalcCoefs(ctx->cfg.Coefs[1],  ctx->luma_temp);
		PrecalcCoefs(ctx->cfg.Coefs[2],  ctx->chroma_spat);
		PrecalcCoefs(ctx->cfg.Coefs[3],  ctx->chroma_temp);
		ctx->i_coeff_changed = 0;
	}

	deNoise(src_Y, dst_Y, ctx->cfg.Line,  &ctx->cfg.Frame[0],  ctx->w, ctx->h,  pitch_src, pitch_dst, ctx->cfg.Coefs[0], ctx->cfg.Coefs[0], ctx->cfg.Coefs[1]);
}

void hqdn3d_process_U(HQDN3D_HANDLE h_filter, unsigned char *src_U, int pitch_src, unsigned char *dst_U, int pitch_dst)
{
	HQDN3D_Context *ctx;

	ctx = (HQDN3D_Context *)h_filter;

	deNoise(src_U, dst_U, ctx->cfg.Line,  &ctx->cfg.Frame[1],  ctx->w / 2 , ctx->h / 2 ,  pitch_src, pitch_dst, ctx->cfg.Coefs[2], ctx->cfg.Coefs[2], ctx->cfg.Coefs[3]);

}

void hqdn3d_process_V(HQDN3D_HANDLE h_filter, unsigned char *src_V, int pitch_src, unsigned char *dst_V, int pitch_dst)
{
	HQDN3D_Context *ctx;

	ctx = (HQDN3D_Context *)h_filter;

	deNoise(src_V,  dst_V, ctx->cfg.Line,  &ctx->cfg.Frame[2],  ctx->w / 2 , ctx->h / 2 ,  pitch_src, pitch_dst, ctx->cfg.Coefs[2], ctx->cfg.Coefs[2], ctx->cfg.Coefs[3]);
}

void hqdn3d_destroy(HQDN3D_HANDLE h_filter)
{
	HQDN3D_Context *ctx;
	int i;
	ctx = (HQDN3D_Context *)h_filter;

	if (ctx){
		hg_aligned_free(ctx->cfg.Line);
	}

	for(i = 0 ; i < 3; ++i){
		if (ctx->cfg.Frame[i]){
			free(ctx->cfg.Frame[i]);
		}
	}
}


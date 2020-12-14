#pragma once 

#ifdef __cplusplus
extern "C"
{
#endif 


typedef void * HQDN3D_HANDLE;

HQDN3D_HANDLE hqdn3d_create(int w, int h, float luma_spat, float luma_temp, float chroma_spat, float chroma_temp);

void hqdn3d_set_coeff(HQDN3D_HANDLE h_filter, float luma_spat, float luma_temp, float chroma_spat, float chroma_temp);

void hqdn3d_process_Y(HQDN3D_HANDLE h_filter, unsigned char *src_Y, int pitch_src, unsigned char *dst_Y, int pitch_dst);

void hqdn3d_process_U(HQDN3D_HANDLE h_filter, unsigned char *src_U, int pitch_src, unsigned char *dst_U, int pitch_dst);

void hqdn3d_process_V(HQDN3D_HANDLE h_filter, unsigned char *src_V, int pitch_src, unsigned char *dst_V, int pitch_dst);

void hqdn3d_destroy(HQDN3D_HANDLE h);


#ifdef __cplusplus
};
#endif 


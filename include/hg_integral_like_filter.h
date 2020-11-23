#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

#define HG_MEAN_FILTER_MAX_R    (127)
#define HG_MEAN_FILTER_MAX_WIDTH    (4096)

/**
 * @brief hg_gray8_meanfilter
 * @param p_src
 * @param width
 * @param height
 * @param src_pitch
 * @param r
 * @param p_dst_mean output matrix, (mean << 8)
 * @param dst_pitch
 * @return
 */
int hg_gray8_meanfilter(const uint8_t *p_src, int width, int height, int src_pitch, int r, uint16_t *p_dst_mean, int dst_pitch);

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
int hg_gray8_mean_variance_filter(const uint8_t *p_src, int width, int height, int src_pitch, int r, uint16_t *p_dst_mean, int pitch_mean, uint32_t *p_dst_variance, int pitch_variance);

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
int hg_gray8_mean_filter_verify(const uint8_t *p_src, int width, int height, int src_pitch, int r , uint16_t *p_dst, int dst_pitch);

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
int hg_gray8_variance_verify(const uint8_t *p_src, int width, int height, int src_pitch, int r, uint32_t *p_dst, int dst_pitch);


#ifdef __cplusplus
}
#endif


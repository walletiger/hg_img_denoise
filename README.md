# hg_img_denoise
some implements of image denoising alg

local  variance denoise and guided denoise filter alg . 

1 局部均方差滤波和 引导滤波的算法实现， 算法做了一定优化。
局部均方差滤波 在 i7-3632QM CPU @ 2.20GHz 上可以 到 720p@200fps . 在 hi3798 的 arm neno 优化版本可以到 720p@60fps 
引导滤波在 i7-3632QM CPU @ 2.20GHz 可以到 720p@150fps 

这两种算法处理静态图片都是比双边滤波效率和效果上比较好的。 可以很好的起到 降噪磨皮效果。
但用在实时视频 最好的应该是3d的算法。 实时的 3d 降噪算法有很多。 




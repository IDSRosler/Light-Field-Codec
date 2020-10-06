import cv2
import os
import glob
import sys
import numpy as np
import pandas as pd
from scipy import signal
from skimage.metrics import structural_similarity as ssim


def immse(x, y):
    return np.mean((x - y).astype('double')**2)


def _ssim(img1, img2, n_bits):
    K = [0.01, 0.03]
    L = 2 ** n_bits - 1
    kernelX = cv2.getGaussianKernel(11, 1.5)
    window = kernelX * kernelX.T

    M, N = np.shape(img1)

    C1 = (K[0]*L)**2
    C2 = (K[1]*L)**2

    mu1 = signal.convolve2d(img1, window, 'valid')
    mu2 = signal.convolve2d(img2, window, 'valid')

    mu1_sq = mu1*mu1
    mu2_sq = mu2*mu2
    mu1_mu2 = mu1*mu2

    sigma1_sq = signal.convolve2d(img1*img1, window, 'valid') - mu1_sq
    sigma2_sq = signal.convolve2d(img2*img2, window, 'valid') - mu2_sq
    sigma12 = signal.convolve2d(img1*img2, window, 'valid') - mu1_mu2

    ssim_map = ((2*mu1_mu2 + C1)*(2*sigma12 + C2)) / \
        ((mu1_sq + mu2_sq + C1)*(sigma1_sq + sigma2_sq + C2))
    mssim = np.mean(ssim_map)
    return mssim, ssim_map


def rgb2ycbcrn(rgb, n_yuv):
    M = np.r_[
        [[0.212600,   0.715200,  0.072200],
         [-0.114572,  -0.385428,  0.500000],
         [0.500000,  -0.454153, -0.045847]]
    ]
    ycbcr = rgb.reshape(-1, 3) @ M.T
    ycbcr = ycbcr.reshape(rgb.shape)
    ycbcr[:, :, 0] = (219 * ycbcr[:, :, 0] + 16)*2**(n_yuv-8)
    ycbcr[:, :, 1:3] = (224 * ycbcr[:, :, 1:3] + 128)*2**(n_yuv-8)
    if n_yuv == 8:
        return ycbcr.astype('int8')
    elif n_yuv in [10, 16]:
        return ycbcr.astype('int16')
    raise ValueError('Unkown/Invalid n_yuv')


def quantitative_metrics(ref_img, rec_img, n_rgb, n_yuv, *results):
    ref = ref_img.astype('double') / (2**n_rgb - 1)
    rec = rec_img.astype('double') / (2**n_rgb - 1)
    ref_ycbcr = rgb2ycbcrn(ref, n_yuv)
    rec_ycbcr = rgb2ycbcrn(rec, n_yuv)

    y1 = ref_ycbcr[:, :, 0]
    u1 = ref_ycbcr[:, :, 1]
    v1 = ref_ycbcr[:, :, 2]

    y2 = rec_ycbcr[:, :, 0]
    u2 = rec_ycbcr[:, :, 1]
    v2 = rec_ycbcr[:, :, 2]

    y_mse = immse(y1, y2)
    u_mse = immse(u1, u2)
    v_mse = immse(v1, v2)
    y_psnr = 10 * np.log10((2**n_rgb-1)*(2**n_yuv-1)/y_mse)
    u_psnr = 10 * np.log10((2**n_rgb-1)*(2**n_yuv-1)/u_mse)
    v_psnr = 10 * np.log10((2**n_rgb-1)*(2**n_yuv-1)/v_mse)

    yuv_psnr = (6*y_psnr + u_psnr + v_psnr)/8
    y_ssim = ssim(y1, y2)
    lc = locals()
    dict_results = {r: lc[r.lower()] for r in results}

    return dict_results


def calculate(ref_folder, rec_folder, *metrics):
    headers = ['VIEW', *metrics]
    lst = []
    ppms = glob.glob(os.path.join(rec_folder, '*.ppm'))

    for view in (ppm.split('/')[-1] for ppm in ppms):
        # print('[%s' % view, end=' ')
        # sys.stdout.flush()
        ref_path = os.path.join(ref_folder, view)
        rec_path = os.path.join(rec_folder, view)
        ref = cv2.cvtColor(cv2.imread(ref_path, -1), cv2.COLOR_BGR2RGB)
        rec = cv2.cvtColor(cv2.imread(rec_path, -1), cv2.COLOR_BGR2RGB)
        result = quantitative_metrics(ref, rec, 10, 10, *metrics)
        entry = [result[m] for m in metrics]
        lst.append([view, *entry])
        # print('Done]', end=' ')
        # sys.stdout.flush()
    # print('\nGenerating data frame...', end=' ')
    # sys.stdout.flush()
    df = pd.DataFrame(lst, columns=headers)
    # print('exporting csv...', end=' ')
    # sys.stdout.flush()
    df.to_csv(os.path.join(rec_folder, 'metrics.csv'))
    # print('Done.')
    # sys.stdout.flush()

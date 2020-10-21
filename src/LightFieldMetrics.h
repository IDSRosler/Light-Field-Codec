

#include "Point4D.h"
#include "utils.h"
#include <opencv2/core/core.hpp> // Basic OpenCV structures (cv::Mat, cv::Scalar)
#include <opencv2/imgproc/imgproc.hpp> // Gaussian Blur



  cv::Mat convertViewToMat(const float *lf, const Point4D &shape, int view_u, int view_v) {
    auto stride = make_stride(shape);

    cv::Mat mat = cv::Mat::zeros(cv::Size(shape.x, shape.y), CV_32F);
    for (Point4D::value_type y = 0; y < shape.y; y++)
      for (Point4D::value_type x = 0; x < shape.x; x++) {
        int index = calc_offset(x, y, view_u, view_v, stride);
        float value = lf[index] + 512;
        float adjusted_value = value;
        mat.at<float>(y, x) = adjusted_value;
      }
    return mat;
  }

  double getPSNR(const cv::Mat &I1, const cv::Mat &I2, const unsigned bit_depth = 8) {
    cv::Mat s1;
    cv::absdiff(I1, I2, s1);  // |I1 - I2|
    s1.convertTo(s1, CV_32F); // cannot make a square on 8 bits
    s1 = s1.mul(s1);          // |I1 - I2|^2

    cv::Scalar s = sum(s1); // sum elements per channel

    double sse = s.val[0]; // sum channels

    if (sse <= 1e-10) // for small values return zero
      return 0;
    else {
      double mse = sse / (double)(I1.channels() * I1.total());
      double dynamic_range = (double)((1 << bit_depth) - 1);
      double psnr = 10.0 * log10((dynamic_range * dynamic_range) / mse);
      return psnr;
    }
  }

  auto getMSSIM(const cv::Mat &i1, const cv::Mat &i2, const unsigned bit_depth = 8) {
    const double dynamic_range = (double)((1 << bit_depth) - 1); // 2 ^ bit_depth - 1
    const double C1 = cv::pow(0.01, 2);
    const double C2 = cv::pow(0.03, 2);

    /***************************** INITS **********************************/
    int d = CV_32F;

    cv::Mat I1, I2;
    i1.convertTo(I1, d); // cannot calculate on one byte large values
    i2.convertTo(I2, d);

    I1 /= dynamic_range;
    I2 /= dynamic_range;

    cv::Mat I2_2 = I2.mul(I2);  // I2^2
    cv::Mat I1_2 = I1.mul(I1);  // I1^2
    cv::Mat I1_I2 = I1.mul(I2); // I1 * I2

    /***********************PRELIMINARY COMPUTING ******************************/

    cv::Mat mu1, mu2; //
    cv::GaussianBlur(I1, mu1, cv::Size(11, 11), 1.5);
    cv::GaussianBlur(I2, mu2, cv::Size(11, 11), 1.5);

    cv::Mat mu1_2 = mu1.mul(mu1);
    cv::Mat mu2_2 = mu2.mul(mu2);
    cv::Mat mu1_mu2 = mu1.mul(mu2);

    cv::Mat sigma1_2, sigma2_2, sigma12;

    cv::GaussianBlur(I1_2, sigma1_2, cv::Size(11, 11), 1.5);
    sigma1_2 -= mu1_2;

    cv::GaussianBlur(I2_2, sigma2_2, cv::Size(11, 11), 1.5);
    sigma2_2 -= mu2_2;

    cv::GaussianBlur(I1_I2, sigma12, cv::Size(11, 11), 1.5);
    sigma12 = cv::abs(sigma12 - mu1_mu2);

    ///////////////////////////////// FORMULA ////////////////////////////////
    cv::Mat t1, t2, t3;

    t1 = 2 * mu1_mu2 + C1;
    t2 = 2 * sigma12 + C2;
    t3 = t1.mul(t2); // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))

    t1 = mu1_2 + mu2_2 + C1;
    t2 = sigma1_2 + sigma2_2 + C2;
    t1 = t1.mul(t2); // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))

    cv::Mat ssim_map;
    cv::divide(t3, t1, ssim_map); // ssim_map =  t3./t1;

    cv::Scalar mssim = cv::mean(ssim_map); // mssim = average of ssim map
    return mssim[0];
  }


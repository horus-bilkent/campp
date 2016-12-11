#ifndef HORUS_TRANSFORM_HPP
#define HORUS_TRANSFORM_HPP
#include <opencv2/imgproc.hpp>
using namespace cv;

void polarTranform(Mat& src, Mat &dst, float cx, float cy, float cr);
void skeletonTransform(Mat src, Mat &dst);

#endif

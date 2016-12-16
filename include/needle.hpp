#ifndef HORUS_NEEDLE_HPP
#define HORUS_NEEDLE_HPP
#include <opencv2/imgproc.hpp>
using namespace cv;

double locateNeedle(Mat img, Mat &img_bgr, Point &center, double cr);

#endif

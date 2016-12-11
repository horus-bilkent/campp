#ifndef HORUS_CENTER_HPP
#define HORUS_CENTER_HPP
#include <opencv2/imgproc.hpp>
using namespace cv;

void locateCenter(Mat img, Mat &img_bgr, Point &center, int &center_range);

#endif

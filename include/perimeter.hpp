#ifndef HORUS_PERIMETER_HPP
#define HORUS_PERIMETER_HPP
#include <opencv2/imgproc.hpp>
#include <vector>
using namespace std;
using namespace cv;

double locateScalemarks(Mat img, vector<double> &theta, Point center, int center_range);

#endif

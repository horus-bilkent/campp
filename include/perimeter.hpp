#ifndef HORUS_PERIMETER_HPP
#define HORUS_PERIMETER_HPP
#include <opencv2/imgproc.hpp>
#include <vector>
using namespace std;
using namespace cv;

double locateScalemarks(Mat img, vector<double> &theta, Point center, int center_range);
double scoreScalemarkModel(const Mat& img_bin, const vector<vector<Point> > &contours, vector<double>& tempt, double small_range, double min_lim);
bool isScaleMark(const vector<Point>& contour, float imgY, float limY);

#endif

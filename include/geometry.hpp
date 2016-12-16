#ifndef CAMPP_GEOMETRY_HPP
#define CAMPP_GEOMETRY_HPP
#include <opencv2/imgproc.hpp>
using namespace cv;

bool lineCircleIntersect(Point2f a, Point2f b, Point2f center, double range);
double calcRange2(double big_range, double small_range, double thickness, double cols);
bool lineLineIntersect(Vec4i a, Vec4i b, Point &res);
void minCovCircleIt(Point &center, double &rad, Point a, Point b);
double angleDiff(double a, double b);

#endif

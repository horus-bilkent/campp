#include <vector>
#include "needle.hpp"
#include "cv_algo.hpp"
#include "geometry.hpp"
#include "transform.hpp"
#include <cmath>
using namespace std;

double locateNeedle(Mat img, Mat &img_bgr, Point center, double cr) {
	double cx = center.x;
	double cy = center.y;
	Mat img_bin;
    skeletonTransform(img, img_bin);
	//adaptiveThreshold(img, img_bin, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 51, 15);
	//thinning(img_bin);

	vector<Vec4i> lines;
	HoughLinesP(img_bin, lines, 1, CV_PI/180, 50, 50, 10);
	Vec4i res;
	double resSize;
	bool isRes = false;

	for(size_t i = 0; i < lines.size(); i++) {
		Vec4i l = lines[i];

		bool flag = lineCircleIntersect(Point2f(l[0], l[1]), Point2f(l[2], l[3]), Point2f(cx, cy), cr);
		double size = (l[0] - l[2]) * (l[0] - l[2]) + (l[1] - l[3]) * (l[1] - l[3]);

		if(flag && (!isRes || resSize < size)) {
			isRes = true;
			resSize = size;
			res = l;
			//cout << "POINT 1: " << l[0] << " " << l[1] << " - POINT 2: " << l[2] << " " << l[3] << endl;
			//cout << "SIZE: " << size << endl;
		}
	}

	line(img_bgr, Point(cx, cy), Point(cx - res[2] + res[0], cy - res[3] + res[1]), Scalar(0, 255, 255), 3, CV_AA);
	circle(img_bgr, Point(cx, cy), cr, Scalar(0, 0, 255), 3, CV_AA);
    return atan2(res[2] - res[0], res[3] - res[1]);
}

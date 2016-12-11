#include "center.hpp"
#include "geometry.hpp"
#include "transform.hpp"
#include <vector>
#include <iostream>
using namespace std;

void locateCenter(Mat img, Mat &img_bgr, Point &center, int &center_range) {
	Mat img_bin;
	skeletonTransform(img, img_bin);
	cvtColor(img_bin, img_bgr, COLOR_GRAY2BGR);
	vector<Vec4i> lines;
	HoughLinesP(img_bin, lines, 1, CV_PI/180, 30, 5, 3);
	Vec4i res;
	double resSize;
	bool isRes = false;

	for(size_t i = 0; i < lines.size(); i++) {
		Vec4i l = lines[i];

		double size = (l[0] - l[2]) * (l[0] - l[2]) + (l[1] - l[3]) * (l[1] - l[3]);
		//line(img_bgr, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(255, 0, 255), 2, CV_AA);

		if(!isRes || resSize < size) {
			isRes = true;
			resSize = size;
			res = l;
			//cout << "POINT 1: " << l[0] << " " << l[1] << " - POINT 2: " << l[2] << " " << l[3] << endl;
			//cout << "SIZE: " << size << endl;
		}
	}

	vector<Point> ints;

	double lim = ((img.rows + img.cols) / 2) * 0.1;
	for(size_t i = 0; i < lines.size(); i++)
		for(size_t j = i + 1; j < lines.size(); j++)
			for(size_t k = j + 1; k < lines.size(); k++) {
				Point p1, p2, p3;
				bool isInter = lineLineIntersect(lines[i], lines[j], p1);
				if(!isInter) continue;
				isInter = lineLineIntersect(lines[j], lines[k], p2);
				if(!isInter) continue;
				isInter = lineLineIntersect(lines[i], lines[k], p3);
				if(!isInter) continue;
				Point center;
				double rad = -1;
				minCovCircleIt(center, rad, p1, p2);
				minCovCircleIt(center, rad, p1, p3);
				minCovCircleIt(center, rad, p3, p2);
				if(rad <= lim) {
					//cout << "CENTER: " << center << " p1: " << p1 << " p2: " << p2 << " p3: " << p3 << endl;
					ints.push_back(center);
				}
			}

	//for(size_t i = 0; i < ints.size(); i++)
	//	circle(img_bgr, ints[i], 1, Scalar(0, 255, 0), 1, CV_AA);

	printf("NUMBER OF INTERSECTIONS FOUND: %lu\n", ints.size());

	double scale = 3;
	double max_range = min(img.rows,img.cols) / 2.0;
	double range_scale = max_range / 10;
	const int center_threshold = ints.size() * 0.1;
	printf("MAX_RANGE: %.3f\n CENTER_THRESH: %d\n", max_range, center_threshold);
	for(int range = range_scale; range < max_range; range += range_scale) {
		int max_count = -1;
		Point centerMax;

		for(int i = 0; i < img.cols / scale; i++)
			for(int j = 0; j < img.rows / scale; j++) {
				double x = i * scale;
				double y = j * scale;
				int count = 0;
				for(size_t t = 0; t < ints.size(); t++) {
					double dist = (x - ints[t].x) * (x - ints[t].x) + (y - ints[t].y) * (y - ints[t].y);
					count += (dist <= range * range);
					//cout << "RANGE*RANGE: " << range * range << " DIST: " << dist << " cent: " << Point(x, y) << " CHECK: " << ints[t] << " = COUNT( " << count << ")" << endl;
				}

				if(count > max_count) {
					cout << "Range: " << range << "\tCENT: " << Point(x, y) << "\t= count( " << count << " )" << endl;
					max_count = count;
					centerMax = Point(x, y);
				}
			}

		if(max_count > center_threshold) {
			circle(img_bgr, centerMax, range, Scalar(255, 255, 0), 3, CV_AA);
			center = centerMax;
			center_range = range;
			return;
		}
	}
}

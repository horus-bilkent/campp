#include "center.hpp"
#include "perimeter.hpp"
#include "needle.hpp"
#include "transform.hpp"
#include <unistd.h>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace cv;
using namespace std;

char winName[] = "DEBUG";
const float thr = 0.25;

void showImageAndExit(const Mat &img) {
	namedWindow(winName, WINDOW_AUTOSIZE);
	imshow(winName, img);

	waitKey(0);
	exit(0);
}

void drawScaleMark(Mat &im, double theta, double cx, double cy, double cr) {
	double rho = cr;
	double a = cos(theta), b = sin(theta);
	double x0 = a * rho, y0 = b * rho;
	//cout << "THETA: " << (theta * 360) / (2 * CV_PI) << endl;
	Point pt1, pt2;
   pt1.x = cvRound(x0 + 10*(a) + cx);
   pt1.y = cvRound(y0 + 10*(b) + cy);
   pt2.x = cvRound(x0 - 10*(a) + cx);
   pt2.y = cvRound(y0 - 10*(b) + cy);
   line(im, pt1, pt2, Scalar(0,0,255), 2, CV_AA);
	line(im, Point(cx, cy), Point(x0 + cx, y0 + cy), Scalar(0, 255, 0), 1);
}

int main(int argc, char** argv) {
	if(argc != 2) {
		cout << "Usage: ./a.out image_path" << endl;
		return -1;
	}

	Mat img, img_tmp, img_bin, img_bgr;
	img_tmp = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	if(!img_tmp.data) {
		cout << "ACAMADIM YA LA" << endl;
	}
	Size size(800, img_tmp.rows * (800.0/img_tmp.cols));
	resize(img_tmp, img, size, 0, 0, INTER_LANCZOS4);
	//img_tmp.copyTo(img);
	cvtColor(img, img_bgr, COLOR_GRAY2BGR);

	Point center;
	int center_range;
	locateCenter(img, img_bgr, center, center_range);
	vector<double> thetas;
	/*circle(img_bgr, center, perimeter, Scalar(0, 255, 255), 1, CV_AA);
	double perimeter = locateScalemarks(img, thetas, center, center_range);
	locateNeedle(img, img_bgr, center, center_range);
	for(size_t i = 0; i < thetas.size(); i++) {
		drawScaleMark(img_bgr, thetas[i], center.x, center.y, perimeter);
	}*/

	namedWindow(winName, WINDOW_AUTOSIZE);
	imshow(winName, img_bgr);

	while(1)
		waitKey(0);
	destroyWindow(winName);

	return 0;
}

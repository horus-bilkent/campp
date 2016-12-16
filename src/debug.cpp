#include "center.hpp"
#include "perimeter.hpp"
#include "needle.hpp"
#include "transform.hpp"
#include "geometry.hpp"
#include <unistd.h>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace cv;
using namespace std;

char winName[] = "DEBUG";
char winName2[] = "DEBUG2";
const float thr = 0.25;
int center_range;
Mat img, img_tmp, img_bin, img_bgr, img_bgr2;
const int start_slider_max = 1;
int start_slider, old_start_slider = 0;
const int perimeter_slider_max = 4000;
int perimeter_slider;
const int small_range_slider_max = 4000;
int small_range_slider;
double r;
Point mcenter;

void showImageAndExit(const Mat &img) {
	namedWindow(winName, WINDOW_AUTOSIZE);
	imshow(winName, img);

	waitKey(0);
	exit(0);
}

void drawScaleMark(Mat &im, double theta, double cx, double cy, double cr, bool isNeedle) {
	double rho = cr;
	double a = cos(theta), b = sin(theta);
	double x0 = a * rho, y0 = b * rho;
	cout << "THETA: " << (theta * 360) / (2 * CV_PI) << endl;
	Point pt1, pt2;
    pt1.x = cvRound(x0 + 10*(a) + cx);
    pt1.y = cvRound(y0 + 10*(b) + cy);
    pt2.x = cvRound(x0 - 10*(a) + cx);
    pt2.y = cvRound(y0 - 10*(b) + cy);
    if(isNeedle)
	   line(im, Point(cx, cy), Point(x0 + cx, y0 + cy), Scalar(255, 255, 0), 3);
    else {
        line(im, pt1, pt2, Scalar(0,0,255), 2, CV_AA);
    	line(im, Point(cx, cy), Point(x0 + cx, y0 + cy), Scalar(0, 255, 0), 1);
    }
}

void drawShit() {
    img_bgr2.copyTo(img_bgr);
    circle(img_bgr, mcenter, center_range, Scalar(0, 255, 255), 3, CV_AA);
    circle(img_bgr, mcenter, perimeter_slider, Scalar(255, 0, 0), 3, CV_AA);;
    if(small_range_slider != 0) {
        r = calcRange2(perimeter_slider, small_range_slider, img.cols * 0.05, img.cols);
        if(r > 0)
            circle(img_bgr, mcenter, r, Scalar(255, 0, 255), 3, CV_AA);
    }
    imshow(winName, img_bgr);
}

void start_debug() {
    cout << "DEBUG STARTED: CENTER(" << mcenter.x << ", " << mcenter.y << ")" << " - PERIMETER: " << perimeter_slider << endl;
    Mat img_pro;
    polarTranform(img, img_pro, mcenter.x, mcenter.y, perimeter_slider);
	vector<vector<Point> > contours;
	vector<Vec4i> hie;
    vector<double> tempt;
	cvtColor(img_pro, img_tmp, COLOR_GRAY2BGR);
    double slid_max = small_range_slider;
    line(img_tmp, Point(small_range_slider - img.cols * 0.03, 0), Point(small_range_slider - img.cols * 0.03, img.rows), Scalar(0, 255, 0), 2);
    line(img_tmp, Point(small_range_slider, 0), Point(small_range_slider, img.rows), Scalar(0, 255, 0), 2);
	findContours(img_pro.clone(), contours, hie, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, Point(0, 0));
    double p = scoreScalemarkModel(img_pro, contours, tempt, small_range_slider, img.cols * 0.05);
    cout << "SCORE: " << p << endl;
    cout << "TEMPT: ";
    for(size_t i = 0; i < tempt.size(); i++)
        cout << tempt[i] << " ";
    cout << endl;
    for(size_t t = 0; t < contours.size(); t++)
	   if(isScaleMark(contours[t], small_range_slider, img.cols)) {
           drawContours(img_tmp, contours, t, Scalar(255, 0, 255), 1, 8, hie, 0, Point());
	   }
    for(size_t i = 0; i < tempt.size(); i++) {
		drawScaleMark(img_bgr, tempt[i], mcenter.x, mcenter.y, r, false);
	}
    imshow(winName, img_bgr);
    imshow(winName2, img_tmp);
}

void mouseCB(int event, int x, int y, int flags, void* userdata) {
    if(flags & EVENT_FLAG_LBUTTON) {
        mcenter = Point(x, y);
        cout << "CENTER:" << mcenter << endl;
        drawShit();
    }
}

void perimeter_on(int, void *) {
    drawShit();
    if(start_slider == 1)
        start_debug();
}

void start_on(int, void *) {
    if(start_slider == 1 && old_start_slider == 0)
        start_debug();
    old_start_slider = start_slider;
}

void createDebugWin() {
    namedWindow(winName2, WINDOW_NORMAL);
    setMouseCallback(winName, mouseCB, NULL);
    createTrackbar("BASLAT", winName, &start_slider, start_slider_max, start_on);
    start_on(start_slider, 0);
    createTrackbar("PERIMETER", winName, &perimeter_slider, perimeter_slider_max, perimeter_on);
    createTrackbar("SMALL_RANGE", winName, &small_range_slider, small_range_slider_max, perimeter_on);
    perimeter_on(perimeter_slider, 0);
}

void workflow() {
    Point center;
	cvtColor(img, img_bgr, COLOR_GRAY2BGR);
	locateCenter(img, img_bgr, center, center_range);
	vector<double> thetas;
	double perimeter = locateScalemarks(img, thetas, center, center_range);
	double needleAngle = locateNeedle(img, img_bgr, center, center_range);
    drawScaleMark(img_bgr, needleAngle, center.x, center.y, perimeter, true);
    drawScaleMark(img_bgr, CV_PI / 2, center.x, center.y, perimeter, true);
	circle(img_bgr, center, perimeter, Scalar(0, 255, 255), 1, CV_AA);
	/*for(size_t i = 0; i < thetas.size(); i++) {
		drawScaleMark(img_bgr, thetas[i], center.x, center.y, perimeter, false);
	}*/
}

int main(int argc, char** argv) {
	/*if(argc != 2) {
		cout << "Usage: ./a.out image_path" << endl;
		return -1;
	}*/

	img_tmp = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	if(!img_tmp.data) {
		cout << "ACAMADIM YA LA" << endl;
	}

    if(argc == 6) {
        mcenter = Point(stoi(argv[2]), stoi(argv[3]));
        perimeter_slider = stoi(argv[4]);
        small_range_slider = stoi(argv[5]);
    }

	Size size(800, img_tmp.rows * (800.0/img_tmp.cols));
	resize(img_tmp, img, size, 0, 0, INTER_LANCZOS4);
	//img_tmp.copyTo(img);
	cvtColor(img, img_bgr, COLOR_GRAY2BGR);
    img_bgr.copyTo(img_bgr2);

    //skeletonTransform(img, img_bgr);
    workflow();
	namedWindow(winName, WINDOW_AUTOSIZE);
    //createDebugWin();
	imshow(winName, img_bgr);

	while(1)
		waitKey(0);
	destroyWindow(winName);

	return 0;
}

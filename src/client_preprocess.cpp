#include "center.hpp"
#include "perimeter.hpp"
#include "needle.hpp"
#include "transform.hpp"
#include "geometry.hpp"
#include <unistd.h>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <raspicam/raspicam_cv.h>

using namespace cv;
using namespace std;

const string winName = "PREPROCESS";
Mat img, img_bgr, img_tmp;
double CENTER_RANGE = 40;
int request_slider = 0;
int preprocess_slider = 0;
int rec_width = 0;
int rec_height = 0;
bool request_new = false;
Point center, rec_center;
double center_range;

vector<int> roi(4);
vector<double> thetas;

void writeOutput(const string& output_path) {
	
		ofstream outfile(output_path);
		// first read center
		outfile << center.x;
		outfile << "\n";
		outfile << center.y;
		outfile << "\n" ;
		outfile << center_range;
		outfile << "\n";
		
		for (auto theta : thetas) {
			outfile << theta;
			outfile << "\n";
		}
		
		outfile << "X";
		outfile << "\n";
		for (auto point : roi) {
			outfile << point;
			outfile << "\n";
		}	
}


void drawScaleMark(Mat &im, double theta, double cx, double cy, double cr, bool isNeedle) {
	double rho = cr;
	double a = cos(theta), b = sin(theta);
	double x0 = a * rho, y0 = b * rho;
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

void preprocess(int, void *) {
	if (preprocess_slider == 1) {
		request_new = true;
		cout << "PREPROCESS STARTED: CENTER(" << center.x << ", " << center.y << ")"  << endl;
		cout << "Locate Needle" << endl;
		double needleAngle = locateNeedle(img, center, CENTER_RANGE);
		cout << "Locate Scale Marks" << endl;
		double perimeter = locateScalemarks(img, thetas, center, center_range);
		drawScaleMark(img_bgr, needleAngle, center.x, center.y, perimeter, true);
		for(size_t i = 0; i < thetas.size(); i++) {
			drawScaleMark(img_bgr, thetas[i], center.x, center.y, perimeter, false);
		}
		sort(thetas.begin(), thetas.end());
		imshow(winName, img_bgr);
	}
}

void drawImage(int, void *) {
	img.copyTo(img_bgr);
    circle(img_bgr, center, 5, Scalar(212, 35, 78), 5, CV_AA);
    rectangle(img_bgr, Point(rec_center.x - rec_width / 2, rec_center.y - rec_height/2), Point(rec_center.x + rec_width / 2, rec_center.y + rec_height/2), Scalar(0,0,255), 5, CV_AA,0);
    circle(img_bgr, rec_center, 5, Scalar(67,13, 214), 5, CV_AA);
    imshow(winName, img_bgr);
    roi.at(0) = rec_center.x - rec_width / 2; 
    roi.at(1) = rec_center.y - rec_height / 2;
    roi.at(2) = rec_center.x + rec_width / 2;
    roi.at(3) = rec_center.y + rec_height / 2;
}

void mouseCB(int event, int x, int y, int flags, void* userdata) {
    if(flags & EVENT_FLAG_LBUTTON) {
        center = Point(x, y);
        drawImage(0,0);
    } else if (flags & EVENT_FLAG_RBUTTON)  {
        rec_center = Point(x, y);
        drawImage(0,0);
	}
}

int main(int argc, char** argv) {
	
	if (argc != 3) {
		cout << "GIVE THE IMAGE FILE AND CONFIGURATION FILE OUTPUT!" << endl;
		return -1;
	}
	
	string image_file = argv[1];
	string config_file = argv[2];
	
	
	img_tmp = imread(image_file, CV_LOAD_IMAGE_GRAYSCALE);
	Size size(750, img_tmp.rows * (750.0/img_tmp.cols));
	resize(img_tmp, img, size, 0, 0, INTER_LANCZOS4);
	// img_tmp.copyTo(img);
	cvtColor(img, img_bgr, COLOR_GRAY2BGR);
	namedWindow(winName, WINDOW_AUTOSIZE);
	createTrackbar("PREPROCESS", winName, &preprocess_slider, 2, preprocess);
	createTrackbar("Rec Width", winName, &rec_width, 2000, drawImage);
	createTrackbar("Rec Height", winName, &rec_height, 2000, drawImage);

	imshow(winName, img_bgr);
	setMouseCallback(winName, mouseCB, NULL);


	waitKey(0);
	
	if (!request_new) {
		writeOutput(config_file);
		return 0;
	} else {
		return 1;
	}
}

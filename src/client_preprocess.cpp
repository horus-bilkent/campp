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
Mat img, img_bgr, img_tmp, img_pro;
double CENTER_RANGE = 50;
int request_slider = 0;
int preprocess_slider = 0;
int rec_width = 0;
int rec_height = 0;
bool request_new = true;
bool cropped = false;
bool center_assigned = false;
Point center, rec_center;

vector<int> roi(4);
vector<double> thetas;

void writeOutput(const string& output_path) {
	
		ofstream outfile(output_path);
		// first read center
		outfile << center.x;
		outfile << "\n";
		outfile << center.y;
		outfile << "\n" ;
		outfile << CENTER_RANGE;
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
	if (preprocess_slider == 1 && !cropped) {
			cout << "CROPPING" << endl;
			cout << roi.at(0) << " " << roi.at(1) <<  " " << roi.at(2) << " " << roi.at(3) << endl;
			cropped = true;	
			cv::Rect myROI(roi.at(0), roi.at(1), roi.at(2), roi.at(3));
			cv::Mat cropped(img, myROI);
			cropped.copyTo(img);
			imshow(winName, img);
	}
	if (preprocess_slider == 2) {
		center_assigned = true;
		img.copyTo(img_bgr);
		cout << "PREPROCESS STARTED: CENTER(" << center.x << ", " << center.y << ")"  << endl;
		cout << "Locate Needle" << endl;
		double needleAngle = locateNeedle(img, center, CENTER_RANGE);
		cout << "Locate Scale Marks" << endl;
		double perimeter = locateScalemarks(img, thetas, center, CENTER_RANGE);
		drawScaleMark(img_bgr, needleAngle, center.x, center.y, perimeter, true);
		circle(img_bgr, center, perimeter, Scalar(0, 255, 255), 3, CV_AA);
		for(size_t i = 0; i < thetas.size(); i++) {
			drawScaleMark(img_bgr, thetas[i], center.x, center.y, perimeter, false);
		}
		sort(thetas.begin(), thetas.end());
	}
	
	if (preprocess_slider == 3) {
		cout << "Accepting the image." << endl;
		request_new = false;
	}
	
	imshow(winName, img_bgr);
}

void drawImage(int event, void *) {
	img.copyTo(img_bgr);
	if (cropped == false) {
		circle(img_bgr, rec_center, 5, Scalar(67,13, 214), 5, CV_AA);
		rectangle(img_bgr, Point(rec_center.x - rec_width / 2, rec_center.y - rec_height/2), Point(rec_center.x + rec_width / 2, rec_center.y + rec_height/2), Scalar(0,0,255), 5, CV_AA,0);
		roi.at(0) = rec_center.x - rec_width / 2; 
		roi.at(1) = rec_center.y - rec_height / 2;
		roi.at(2) = rec_width;
		roi.at(3) = rec_height;	
	}
	
	if (center_assigned == false) {
		circle(img_bgr, center, 5, Scalar(212, 35, 78), 5, CV_AA);
	}
	imshow(winName, img_bgr);    
}

void mouseCB(int event, int x, int y, int flags, void* userdata) {
    if(flags & EVENT_FLAG_LBUTTON) {
		if (center_assigned == false) {
        center = Point(x, y);
        drawImage(0,0);
		}
    } else if (flags & EVENT_FLAG_RBUTTON)  {
		if (cropped == false) {
			rec_center = Point(x, y);
			drawImage(0,0);
		}
	}
}

int main(int argc, char** argv) {
	
	if (argc != 3) {
		cout << "GIVE THE IMAGE FILE AND CONFIGURATION FILE OUTPUT!" << endl;
		return -1;
	}
	
	string image_file = argv[1];
	string config_file = argv[2];
	
	
	img = imread(image_file, CV_LOAD_IMAGE_GRAYSCALE);
	cvtColor(img, img_bgr, COLOR_GRAY2BGR);
	namedWindow(winName, WINDOW_NORMAL);
	createTrackbar("PREPROCESS", winName, &preprocess_slider, 3, preprocess);
	createTrackbar("Rec Width", winName, &rec_width, 2000, drawImage);
	createTrackbar("Rec Height", winName, &rec_height, 2000, drawImage);

	imshow(winName, img_bgr);
	setMouseCallback(winName, mouseCB, NULL);


	waitKey(0);
	
	if (!request_new) {
		writeOutput(config_file);
		exit(0);
	} else {
		exit(1);
	}
}

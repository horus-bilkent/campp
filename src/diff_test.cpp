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
#include <raspicam/raspicam_cv.h>

using namespace cv;
using namespace std;

char winName[] = "DEBUG";
char winName2[] = "DEBUG2";

void initCamera(raspicam::RaspiCam_Cv& Camera) {
	    Camera.set( CV_CAP_PROP_FORMAT, CV_8UC1 );
		if (!Camera.open()) {
			cerr<<"Error opening the camera"<<endl;return;
		}
}


void destructCamera(raspicam::RaspiCam_Cv& Camera) {
	    Camera.release();
}

Mat captureImage(raspicam::RaspiCam_Cv& Camera) {
	  cv::Mat image;
	  Camera.grab();
      Camera.retrieve(image);
	  return image;
}


void showImageAndExit(const Mat &img) {
	namedWindow(winName, WINDOW_AUTOSIZE);
	imshow(winName, img);

    while(1)
	   waitKey(0);
	exit(0);
}

Point locateControl(const Mat &img, const Mat &templ) {
    Mat result;
    int result_cols =  img.cols - templ.cols + 1;
    int result_rows = img.rows - templ.rows + 1;
    result.create( result_rows, result_cols, CV_32FC1 );

    matchTemplate( img, templ, result, CV_TM_SQDIFF);
    normalize( result, result, 0, 1, NORM_MINMAX, -1, Mat() );

    double minVal; double maxVal; Point minLoc; Point maxLoc;
    Point matchLoc;

    minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, Mat() );
    matchLoc = minLoc;
    cout << "MATCH:" << matchLoc << endl;
    return matchLoc;
	/*namedWindow(winName2, WINDOW_AUTOSIZE);
	imshow(winName2, templ);
    rectangle( img, matchLoc, Point( matchLoc.x + templ.cols , matchLoc.y + templ.rows ), Scalar::all(0), 2, 8, 0 );
    showImageAndExit(img);*/
}

void diffTransform(Mat img1, Mat img2, Mat &img_bgr) {
    /*Point newC = locateControl(img2, templ);
    int roiL = min(img1.cols - cPoint.x, img2.cols - newC.x);
    int roiH = min(img1.rows - cPoint.y, img2.rows - newC.y);
    cout << img1.cols << " " << img1.rows << endl;
    cout << img2.cols << " " << img2.rows << endl;
    cout << "ROI " << roiL << " " << roiH << endl;
    cout << "cPoint: " << cPoint << endl;
    cout << "newC: " << newC << endl;
    img1 = img1(Rect(cPoint.x, cPoint.y, roiL, roiH));
    img2 = img2(Rect(newC.x, newC.y, roiL, roiH));*/

    Mat img_diff;
    absdiff(img1, img2, img_diff);
    Mat fg_mask = Mat::zeros(img_diff.rows, img_diff.cols, CV_8UC1);
    float threshold = 30.0f;
    int threshX = img1.cols / 10;
    int threshY = img1.rows / 10;
    float dist;
    for(int j = 30; j < img_diff.rows - 30; j++)
        for(int i = 30; i < img_diff.cols - 30; i++) {
            Vec3b pix = img_diff.at<Vec3b>(j, i);
            dist = (pix[0]*pix[0] + pix[1] * pix[1] + pix[2] * pix[2]);
            dist = sqrt(dist);
            if(dist > threshold) {
                fg_mask.at<unsigned char>(j, i) = 255;
            }
        }

    fg_mask.copyTo(img_bgr);
}

Mat getTemplate(const Mat &img) {
    return img(Rect(30, 30, 70, 70));
}

int main(int argc, char** argv) {
    VideoCapture cap(0);
    if (!cap.isOpened()) {
		std::cout << "cam not open" << std::endl;
		return 0;
	}
    Mat frame, img1, img2, img_bgr;
    bool flag = false;
	namedWindow(winName, WINDOW_AUTOSIZE);
    while(1) {
        cap >> frame;
            imshow(winName, frame);
        waitKey(20);
    }

    /*Mat img1, img2, img_bgr;
    img1 = imread(argv[1], CV_LOAD_IMAGE_COLOR);
	if(!img1.data) {
		cout << "ACAMADIM YA LA" << endl;
        return 0;
	}

    img2 = imread(argv[2], CV_LOAD_IMAGE_COLOR);
	if(!img2.data) {
		cout << "ACAMADIM YA LA" << endl;
        return 0;
	}

    Mat templ = getTemplate(img1);
    diffTransform(img1, img2, img_bgr, templ, Point(30, 30));

	namedWindow(winName, WINDOW_AUTOSIZE);
	//imshow(winName, templ);
	imshow(winName, img_bgr);

	while(1)
		waitKey(0);
	destroyWindow(winName);*/

	return 0;
}

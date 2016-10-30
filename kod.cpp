#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include <cmath>
#include <vector>
using namespace cv;
using namespace std;

char winName1[] = "YARRAK";
char winName2[] = "DASSAK";
char winName3[] = "AM";
const double thr = 0.05;

void thinningIteration(cv::Mat& im, int iter)
{
    cv::Mat marker = cv::Mat::zeros(im.size(), CV_8UC1);

    for (int i = 1; i < im.rows-1; i++)
    {
        for (int j = 1; j < im.cols-1; j++)
        {
            uchar p2 = im.at<uchar>(i-1, j);
            uchar p3 = im.at<uchar>(i-1, j+1);
            uchar p4 = im.at<uchar>(i, j+1);
            uchar p5 = im.at<uchar>(i+1, j+1);
            uchar p6 = im.at<uchar>(i+1, j);
            uchar p7 = im.at<uchar>(i+1, j-1);
            uchar p8 = im.at<uchar>(i, j-1);
            uchar p9 = im.at<uchar>(i-1, j-1);

            int A  = (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) + 
                     (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) + 
                     (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
                     (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1);
            int B  = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
            int m1 = iter == 0 ? (p2 * p4 * p6) : (p2 * p4 * p8);
            int m2 = iter == 0 ? (p4 * p6 * p8) : (p2 * p6 * p8);

            if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
                marker.at<uchar>(i,j) = 1;
        }
    }

    im &= ~marker;
}

void thinning(cv::Mat& im)
{
    im /= 255;

    cv::Mat prev = cv::Mat::zeros(im.size(), CV_8UC1);
    cv::Mat diff;

    do {
        thinningIteration(im, 0);
        thinningIteration(im, 1);
        cv::absdiff(im, prev, diff);
        im.copyTo(prev);
    } 
    while (cv::countNonZero(diff) > 0);

    im *= 255;
}

bool circleLineIntersect(float x1, float y1, float x2, float y2, float 
	cx, float cy, float cr ) {
	float dx = x2 - x1;
	float dy = y2 - y1;
	float a = dx * dx + dy * dy;
	float b = 2 * (dx * (x1 - cx) + dy * (y1 - cy));
	float c = cx * cx + cy * cy;
	c += x1 * x1 + y1 * y1;
	c -= 2 * (cx * x1 + cy * y1);
	c -= cr * cr;
	float bb4ac = b * b - 4 * a * c;
	
	if(bb4ac<0) {
		return false;    // No collision
	}
	else {
		return true;      //Collision
	}
}

void myhough(cv::Mat& im, vector<Vec4i>& lines, float cx, float cy, float cr) {
	const float theta = CV_PI / 180;
	const float rho = 5;
	int acc[500][500];

	for(int i = 0; i < im.rows; i++)
		for(int j = 0; j < im.cols; j++) {
			uchar p = im.at<uchar>(i, j);
			if(p != 0) {
				int x = i - cx;
				int y = j - cy;
				
				for(int a = 0; a * theta <= CV_PI; a++) {
					float b = x * cos(a * theta) + y * sin(a * theta);
					if(b <= cr && b >= -cr) {
						int bint = (b + cr) / rho;
						acc[a][bint]++;
					}
				}
			}
		}
	
		
	for(int i = 0; i * theta <= CV_PI; i++)
		for(int j = 0; j * rho <= cr + cr; j++) {
			if(acc[i][j] > 200) {
				printf("ACC: %d, %d - %g, %g = %d\n", i, j, i * theta, (j * rho) - cr, acc[i][j]);
				int x1 = -cx;
				int y1 = ((j * rho) - cr - x1 * cos(i * theta)) / sin(i * theta);
				int x2 = im.rows - 1 -cx;
				int y2 = ((j * rho) - cr - x2 * cos(i * theta)) / sin(i * theta);
				lines.push_back(Vec4i(x1 + cx, y1 + cy, x2 + cx, y2 + cy));
			}
		}
}

int main(int argc, char** argv) {
	if(argc != 2) {
		cout << "Usage: ./a.out image_path" << endl;
		return -1;
	}

	Mat img, img_bin, img_bgr;
	img = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	if(!img.data) {
		cout << "ACAMADIM YA LA" << endl;
	}


	//GaussianBlur(img, img, Size(11, 11), 6.0);
	Mat img_equ;
	equalizeHist(img, img_equ);
	GaussianBlur(img_equ, img_equ, Size(11, 11), 6.0);
	adaptiveThreshold(img, img_bin, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 51, 15);
	adaptiveThreshold(img_equ, img_equ, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 51, 15);
	thinning(img_equ);
	thinning(img_bin);

	cvtColor(img_equ, img_bgr, CV_GRAY2BGR);
	float cx = img_bgr.rows / 2, cy = img_bgr.cols / 2;
	float cr = img_bgr.rows / 25;
	float minLength = sqrt(img_bgr.rows * img_bgr.rows + img_bgr.cols + img_bgr.cols) * thr;

	vector<Vec4i> lines;
	//myhough(img_bin, lines, cx, cy, cr);
	HoughLinesP(img_equ, lines, 1, 1 * (CV_PI/180), 80, minLength, 10);

	for(size_t i = 0; i < lines.size(); i++) {
		Vec4i l = lines[i];
		cout << "POINT 1: " << l[0] << " " << l[1] << " - POINT 2: " << l[2] << " " << l[3] << endl;

		bool flag = circleLineIntersect(l[0], l[1], l[2], l[3], cx, cy, cr);
		
		//if(flag)
			line(img_bgr, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(0, 255, 255), 1, CV_AA);
	}

	circle(img_bgr, Point(img_bgr.rows / 2, img_bgr.cols / 2), cr, Scalar(0, 0, 255), -cr, CV_AA);

	namedWindow(winName2, WINDOW_AUTOSIZE);
	imshow(winName2, img_bgr);
	//imshow(winName2, img_equ);

	waitKey(0);
	//destroyWindow(winName1);
	destroyWindow(winName2);

	return 0;
}

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <cmath>
using namespace cv;
using namespace std;

char winName[] = "YARRAK";
const float thr = 0.25;

/**
 * Perform one thinning iteration.
 * Normally you wouldn't call this function directly from your code.
 *
 * @param  im    Binary image with range = 0-1
 * @param  iter  0=even, 1=odd
 */
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

/**
 * Function for thinning the given binary image
 *
 * @param  im  Binary image with range = 0-255
 */
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

void precalc(Mat& src, Mat &dst, float cx, float cy, float cr) {
	Mat temp;
	equalizeHist(src, dst);
	medianBlur(dst, dst, 5);

	linearPolar(dst, dst, Point2f(cx, cy), cr, INTER_LINEAR + WARP_FILL_OUTLIERS); 
	adaptiveThreshold(dst, dst, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 51, 15);
	Sobel(dst, dst, dst.type(), 0, 1, 9, 1, 0, BORDER_DEFAULT);
	//Canny(dst, dst, 50, 150, 3);
	thinning(dst);

	//dilate(dst, dst, element);
	//erode(dst, dst, element);
}

bool isScaleMark(const vector<Point>& contour, float imgY) {
	float minX, maxX;
	bool flagX = true, flagY = true;

	float thresY = imgY * (9.5 / 10); 

	cout << "########### START CONTOUR #############" << endl;
	for(size_t i = 0; i < contour.size(); i++) {
		//cout << contour[i] << endl;
		if(contour[i].x > thresY) {
			if(flagX || contour[i].y < minX) {
				minX = contour[i].y;
				flagX = false;
			}
			if(flagY || contour[i].y > maxX) {
				maxX = contour[i].y;
				flagY = false;
			}
		}
	}	

	cout << "########### END CONTOUR[" << (!(flagX || flagY) && (maxX - minX < 10)) << "] #############" << endl;
	if(flagX || flagY)
		return false;

	return maxX - minX < 10;
}

double getPolarScaleMark(Mat& im, const vector<Point> & contour, int cx, int cy, int cr) {
	float avgY = 0;
	int numY = 0;

	for(size_t i = 0; i < contour.size(); i++)
		if(contour[i].x > im.cols * 0.9) {
			numY++;
			avgY += contour[i].y;
		}
	avgY /= numY;

	double theta = (avgY * 2 * CV_PI) / (im.rows - 1);
	return theta;
}

void drawScaleMark(Mat &im, double theta, double cx, double cy, double cr) {
	double rho = cr;
	double a = cos(theta), b = sin(theta);
	double x0 = a * rho, y0 = b * rho;
	cout << "THETA: " << (theta * 360) / (2 * CV_PI) << endl;
	Point pt1, pt2;
   pt1.x = cvRound(x0 + 10*(a) + cx);
   pt1.y = cvRound(y0 + 10*(b) + cy);
   pt2.x = cvRound(x0 - 10*(a) + cx);
   pt2.y = cvRound(y0 - 10*(b) + cy);
   line(im, pt1, pt2, Scalar(0,0,255), 2, CV_AA);
	line(im, Point(cx, cy), Point(x0 + cx, y0 + cy), Scalar(0, 255, 0), 1);
}

void findNeedle(Mat img, Mat &img_bgr, double cx, double cy, double cr) {
	Mat img_bin;
	adaptiveThreshold(img, img_bin, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 51, 15);
	thinning(img_bin);

	vector<Vec4i> lines;
	HoughLinesP(img_bin, lines, 1, CV_PI/180, 50, 50, 10);
	Vec4i res;
	double resSize;
	bool isRes = false;

	for(size_t i = 0; i < lines.size(); i++) {
		Vec4i l = lines[i];

		bool flag = circleLineIntersect(l[0], l[1], l[2], l[3], cx, cy, cr);
		double size = (l[0] - l[2]) * (l[0] - l[2]) + (l[1] - l[3]) * (l[1] - l[3]);
		
		if(flag && (!isRes || resSize < size)) {
			isRes = true;
			resSize = size;
			res = l;
			cout << "POINT 1: " << l[0] << " " << l[1] << " - POINT 2: " << l[2] << " " << l[3] << endl;
			cout << "SIZE: " << size << endl;
		}
	}

	line(img_bgr, Point(res[0], res[1]), Point(res[2], res[3]), Scalar(0, 255, 255), 3, CV_AA);
	circle(img_bgr, Point(img_bgr.rows / 2, img_bgr.cols / 2), 50, Scalar(0, 0, 255), 3, CV_AA);	 
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

	adaptiveThreshold(img, img_bin, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 51, 15);
	thinning(img_bin);
	Mat img_pro;
	float cx = img.rows / 2, cy = img.cols / 2;
	float cr = img.rows / 30;
	//float minLength = sqrt(img_bgr.rows * img_bgr.rows + img_bgr.cols * img_bgr.cols) * thr;
	precalc(img, img_pro, cx, cy, 190);

	cvtColor(img, img_bgr, COLOR_GRAY2BGR);
	findNeedle(img, img_bgr, cx, cy, cr);
	img_pro.copyTo(img_bin);

	vector<vector<Point> > contours;
	vector<Vec4i> hie;
	findContours(img_bin, contours, hie, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, Point(0, 0)); 
	for(size_t i = 0; i < contours.size(); i++)
		if(isScaleMark(contours[i], img_bin.cols)) {
			double theta = getPolarScaleMark(img_bgr, contours[i], cx, cy, 190);
			drawScaleMark(img_bgr, theta, cx, cy, 190);
		}
	
	
	cout << "POLAR IMAGE SIZE: " << img_bgr.rows << ", " << img_bgr.cols << ";" << endl;

	namedWindow(winName, WINDOW_AUTOSIZE);
	imshow(winName, img_bgr);

	waitKey(0);
	destroyWindow(winName);

	return 0;
}

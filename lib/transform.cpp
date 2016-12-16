#include "transform.hpp"
#include "cv_algo.hpp"
#include <iostream>
using namespace std;

void polarTranform(Mat& src, Mat &dst, float cx, float cy, float cr) {
	Mat temp;
	equalizeHist(src, dst);
	medianBlur(dst, dst, 5);

	linearPolar(dst, dst, Point2f(cx, cy), cr, INTER_LINEAR + WARP_FILL_OUTLIERS);
	adaptiveThreshold(dst, dst, 255, ADAPTIVE_THRESH_GAUSSIAN_C, THRESH_BINARY_INV, 51, 15);
	//Canny(dst, dst, 50, 150, 3);
	//thinning(dst);
	//Sobel(dst, dst, dst.type(), 0, 1, 9, 1, 0, BORDER_DEFAULT);
	Mat element = getStructuringElement(MORPH_RECT, Size(17, 1));
	erode(dst, dst, element, Point(-1, -1));
	dilate(dst, dst, element, Point(-1, -1));
	//thinning(dst);
	//showImageAndExit(dst);

	//erode(dst, dst, element);
}

void skeletonTransform(Mat src, Mat &dst) {
	equalizeHist(src, dst);
	medianBlur(dst, dst, 7);

    //cout << mean(src)[0] << endl;
    //int bin_thresh = mean(src)[0] < 70 ? THRESH_BINARY: THRESH_BINARY_INV;
	adaptiveThreshold(src, dst, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 41, 15);
	//Sobel(dst, dst, dst.type(), 0, 1, 9, 1, 0, BORDER_DEFAULT);
	//Canny(dst, dst, 50, 150, 3);

	/*GaussianBlur(src, dst, Size(3, 3), 0, 0);
	adaptiveThreshold(dst, dst, 255, ADAPTIVE_THRESH_GAUSSIAN_C, THRESH_BINARY_INV, 51, 15);
	//threshold(dst, dst, 0, 255, CV_THRESH_BINARY_INV | CV_THRESH_OTSU);*/
	Mat element = getStructuringElement(MORPH_RECT, Size(7, 7), Point(3, 3));
	morphologyEx(dst, dst, MORPH_CLOSE, element);
	//Mat element2 = getStructuringElement(MORPH_RECT, Size(3, 3), Point(1, 1));
	//morphologyEx(dst, dst, MORPH_OPEN, element2);
	//Canny(dst, dst, 50, 150, 3);

	thinning(dst);

}

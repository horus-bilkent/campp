#include "perimeter.hpp"
#include "transform.hpp"
#include "geometry.hpp"
#include <iostream>

bool isScaleMark(const vector<Point>& contour, float imgY, float limY) {
	float minX, maxX;
	double minY = contour[0].x, maxY = contour[0].x;
	bool flagX = true, flagY = true;

	float thresY = imgY  - limY * 0.05;

	//cout << "imgY: " << imgY << " thresY: " << thresY << " maxY: " << limY << endl;
	//cout << "########### START CONTOUR #############" << endl;
	for(size_t i = 0; i < contour.size(); i++) {
		//cout << contour[i] << endl;
		minY = min(minY, (double)contour[i].x);
		maxY = max(maxY, (double)contour[i].x);
		if(contour[i].x >= thresY && contour[i].x <= imgY) {
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

	if(minY > thresY || maxY < imgY)
		return false;
	//cout << "########### END CONTOUR[" << (!(flagX || flagY) && (maxX - minX < 10)) << "] #############" << endl;
	if(flagX || flagY)
		return false;

	return maxX - minX < 10;
}

double getScalemarkAngle(Mat& im, const vector<Point> & contour, int xMax) {
	float avgY = 0;
	int numY = 0;

	for(size_t i = 0; i < contour.size(); i++)
		if(contour[i].x > xMax - im.cols * 0.1 && contour[i].x <= xMax) {
			numY++;
			avgY += contour[i].y;
		}
	avgY /= numY;

	double theta = (avgY * 2 * CV_PI) / (im.rows - 1);
	return theta;
}

void getScoreIt(double diff, int &lSub, int &cSub, double &prevDiff) {
	const double epsilon = (2.0 / 360) * 2 * CV_PI;

	if(abs(diff - prevDiff) < epsilon) {
		cSub++;
	} else {
		cSub = 1;
		prevDiff = diff;
	}
	if(cSub > lSub) {
		lSub = cSub;
	}
}

double getScore(vector<double> &thetas) {
	if(thetas.size() < 1)
		return -1;
	sort(thetas.begin(), thetas.end());
	vector<double> diff;

	diff.push_back(thetas[0] + 2 * CV_PI - thetas[thetas.size() - 1]);
	for(size_t i = 1; i < thetas.size(); i++)
		diff.push_back(thetas[i] - thetas[i - 1]);


	/*cout << "############## DIFFS START [SIZE: " << diff.size() << "] #################### " << endl;
	for(size_t i = 0; i < diff.size(); i++)
		cout << "#" << i << " Theta: " << thetas[i] << " Diff: " << (diff[i] * 360) / (2 * CV_PI)  << endl;*/

	int lSub = 0, cSub = 0;
	double prevDiff;
	prevDiff = diff[0];
	lSub = 1;
	cSub = 1;

	for(size_t i = 1; i < diff.size(); i++)
		getScoreIt(diff[i], lSub, cSub, prevDiff);

	for(size_t i = 0; i < diff.size(); i++)
		getScoreIt(diff[i], lSub, cSub, prevDiff);
	//cout << "############## DIFFS END [MAX: " << lSub << "] #################### " << endl;

	return lSub;
}

double locateScalemarks(Mat img, vector<double> &theta, Point center, int center_range) {
	Mat img_pro, img_bin, img_bgr;
	double maxPoint;
	bool isMax = false;
	vector<double> tempt;
	cvtColor(img, img_bgr, COLOR_GRAY2BGR);
	double res = -1;
	double max_range = sqrt(img.cols * img.cols + img.rows * img.rows) * 2;
	double min_range = sqrt(img.cols * img.cols + img.rows * img.rows) * 0.1;
	double big_rate = (max_range - min_range) / 10;
	double min_lim = 0.05 * img.cols;
	Mat tempIm;
	double bigTemp;
	double smallTemp;

	for(int i = 0; i < center_range; i++)
		for(double j = 0; j < 2 * CV_PI; j += CV_PI / 180) {
			double x = i * cos(j) + center.x;
			double y = j * sin(j) + center.y;

			for(double big_range = min_range; big_range <= max_range; big_range += big_rate) {
				polarTranform(img, img_pro, x, y, big_range);
				img_pro.copyTo(img_bin);
				vector<vector<Point> > contours;
				vector<Vec4i> hie;
				findContours(img_bin.clone(), contours, hie, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, Point(0, 0));

				for(int small_range = min_lim; small_range <= img_bin.cols; small_range++)  {
					tempt.clear();
					double r = calcRange2(big_range, small_range, min_lim, img.cols);

					for(size_t t = 0; t < contours.size(); t++)
						if(isScaleMark(contours[t], small_range, img.cols)) {
							tempt.push_back(getScalemarkAngle(img_bin, contours[t], small_range));
						}

					double p = getScore(tempt);
					if(!isMax || p > maxPoint) {
						cout << "#### Trying (" << x << ", " << y << ", " << r << ", " << big_range << ", "  << small_range << ") #########" << endl;
						img_bin.copyTo(tempIm);
						isMax = true;
						maxPoint = p;
						theta = tempt;
						bigTemp = big_range;
						smallTemp = small_range;
						res = r;
						cout << "SCORE: " << p << " MAX: " << maxPoint << endl;
					}
				}
			}

			/*vector<vector<Point> > contours;
			vector<Vec4i> hie;
			cout << "DEBUG#########################################################" << endl;
			findContours(tempIm.clone(), contours, hie, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, Point(0, 0));
			for(size_t t = 0; t < contours.size(); t++)
				if(isScaleMark(contours[t], smallTemp, img.cols)) {
					drawContours(tempIm, contours, t, Scalar(255, 0, 0), 2, 8, hie, 0, Point());
				}
			showImageAndExit(tempIm);*/

			return res;
		}

	return res;
}

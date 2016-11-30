#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <cmath>
#include <algorithm>
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

	//cout << "########### START CONTOUR #############" << endl;
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

	//cout << "########### END CONTOUR[" << (!(flagX || flagY) && (maxX - minX < 10)) << "] #############" << endl;
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
	//cout << "THETA: " << (theta * 360) / (2 * CV_PI) << endl;
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
			//cout << "POINT 1: " << l[0] << " " << l[1] << " - POINT 2: " << l[2] << " " << l[3] << endl;
			//cout << "SIZE: " << size << endl;
		}
	}

	line(img_bgr, Point(res[0], res[1]), Point(res[2], res[3]), Scalar(0, 255, 255), 3, CV_AA);
	circle(img_bgr, Point(cx, cy), cr, Scalar(0, 0, 255), 3, CV_AA);	 
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

Vec3i findCenter(Mat img, vector<double> &theta) {
	Mat img_pro, img_bin, img_bgr;
	double maxPoint;
	bool isMax = false;
	Vec3i res;
	vector<double> tempt;
	cvtColor(img, img_bgr, COLOR_GRAY2BGR);

	for(int x = 0; x < img.cols - 1; x++)
		for(int y = 0; y < img.rows - 1; y++) {
			for(int r = img.cols / 5; r < img.cols * 0.95; r++) {
				/*x = img.cols / 2 + 1;
				y = img.rows / 2 + 1;
				r = 190;*/
				cout << "#### Trying (" << x << ", " << y << ", " << r << ") #########" << endl;
				precalc(img, img_pro, x, y, r);
				img_pro.copyTo(img_bin);

				vector<vector<Point> > contours;
				tempt.clear();
				vector<Vec4i> hie;
				findContours(img_bin, contours, hie, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, Point(0, 0));
				for(size_t i = 0; i < contours.size(); i++)
					if(isScaleMark(contours[i], img_bin.cols)) {
						tempt.push_back(getPolarScaleMark(img_bin, contours[i], x, y, r));
					}

				double p = getScore(tempt);
				if(!isMax || p > maxPoint) {
					maxPoint = p;
					theta = tempt;
				}

				cout << "SCORE: " << p << endl;
	/*			for(size_t i = 0; i < theta.size(); i++) {
					drawScaleMark(img_bgr, theta[i], x, y, r);
				}

				namedWindow(winName, WINDOW_AUTOSIZE);
				imshow(winName, img_bgr);

				waitKey(0);
				destroyWindow(winName);

				return maxPoint;*/
			}
		}
	
	return maxPoint;
}

bool getIntersection(Vec4i a, Vec4i b, Point &res) {
	double A1 = a[3] - a[1], B1 = a[0] - a[2], C1 = A1*a[0] + B1 * a[1];
	double A2 = b[3] - b[1], B2 = b[0] - b[2], C2 = A2*b[0] + B2 * b[1];
	double det = A1 * B2 - A2 * B1;
	if(det == 0)
		return false;
	res.x = (B2 * C1 - B1 * C2) / det;
	res.y = (A1 * C2 - A2 * C1) / det;
	return true;
}

void calcCircle(Point &center, double &rad, Point a, Point b) {
	double dist = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
	if(rad == -1 || dist >= rad) {
		rad = dist;
		center.x = (a.x + b.x) / 2.0;
		center.y = (a.y + b.y) / 2.0;
	}
}

void test_img(Mat img, Mat &img_bgr) {
	Mat img_bin;
	adaptiveThreshold(img, img_bin, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, 51, 15);
	thinning(img_bin);
	cvtColor(img_bin, img_bgr, COLOR_GRAY2BGR);

	vector<Vec4i> lines;
	HoughLinesP(img_bin, lines, 1, CV_PI/180, 30, 5, 10);
	Vec4i res;
	double resSize;
	bool isRes = false;

	for(size_t i = 0; i < lines.size(); i++) {
		Vec4i l = lines[i];

		double size = (l[0] - l[2]) * (l[0] - l[2]) + (l[1] - l[3]) * (l[1] - l[3]);
		line(img_bgr, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(255, 0, 255), 2, CV_AA);
		
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
				bool isInter = getIntersection(lines[i], lines[j], p1);
				if(!isInter) continue;
				isInter = getIntersection(lines[j], lines[k], p2);
				if(!isInter) continue;
				isInter = getIntersection(lines[i], lines[k], p3);
				if(!isInter) continue;
				Point center;
				double rad = -1;
				calcCircle(center, rad, p1, p2);
				calcCircle(center, rad, p1, p3);
				calcCircle(center, rad, p3, p2);
				if(rad <= lim) {
					//cout << "CENTER: " << center << " p1: " << p1 << " p2: " << p2 << " p3: " << p3 << endl;
					ints.push_back(center);
				}
			}
	
	for(size_t i = 0; i < ints.size(); i++)
		circle(img_bgr, ints[i], 1, Scalar(0, 255, 0), 1, CV_AA);	 
	
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
			return;	
		}
	}

	//line(img_bgr, Point(res[0], res[1]), Point(res[2], res[3]), Scalar(0, 255, 255), 3, CV_AA);
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

	vector<double> thetas;
	//Vec3i center = findCenter(img, thetas);
	//circle(img_bgr, Point(center[0], center[1]), center[2], Scalar(0, 0, 255), 3, CV_AA);	 

	test_img(img, img_bgr);
	//findNeedle(img, img_bgr, cx, cy, cr); 

	namedWindow(winName, WINDOW_AUTOSIZE);
	imshow(winName, img_bgr);

	waitKey(0);
	destroyWindow(winName);

	return 0;
}

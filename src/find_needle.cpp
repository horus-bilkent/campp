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

using namespace cv;
using namespace std;


int readInput(const string& input_path, vector<double>& thetas, vector<int>& roi, Point& center, double& center_range) {
		std::ifstream infile(input_path);
		std::string line;
		// first read center
		
		double center_x, center_y;
		infile >> center_x >> center_y >> center_range;
		center = Point(center_x, center_y);
		
		while  (getline(infile, line)) {

			if (line.empty()) {
				continue;
			} 
			if (line == "X") {
				break;
			}
			double theta_val = std::stod(line);
			thetas.push_back(theta_val);
		}
				
		while (getline(infile, line)) {
			if (line.empty()) {
				continue;
			} 
			roi.push_back(std::stoi(line));
		}
		return 0;	
	}

int main(int argc, char** argv) {
	const string winName = "APPLIANCE";
	bool debug = false;
	if (argc < 3) {
		cout << "GIVE THE CONFIGURATION FILE AND THE IMAGE!" << endl;
		return -1;
	}
	
	if (argc == 4) {
		debug = true;
	}
	
	vector<double> thetas;
	Point center;
	double center_range;
	vector<int> roi;
	
	if (readInput(argv[2], thetas, roi, center, center_range) != 0) {
		cout << "Error!" << endl;
		return -1;
	}


	string line;
	Mat img;
	Mat img_tmp = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	Size size(750, img_tmp.rows * (750.0/img_tmp.cols));
	resize(img_tmp, img, size, 0, 0, INTER_LANCZOS4);
	cv::Rect myROI(roi.at(0), roi.at(1), roi.at(2), roi.at(3));
	cv::Mat cropped(img, myROI);
	double needleAngle = locateNeedle(cropped, center, center_range);
	int needle_value = 100 * (needleAngle - thetas.front()) / (thetas.back() - thetas.front());
	if (needle_value > 100) {
		needle_value = 100;
	}
	
	if (needle_value < 0) {
		needle_value = 0;
	}
	
	if (debug) {
		const string winName = "DEBUG";
		namedWindow(winName, WINDOW_AUTOSIZE);
		imshow(winName, cropped);
		waitKey(3000);
	}

	std::cout << "Value read: " << needle_value << std::endl;
	return needle_value;
}

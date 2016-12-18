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

void writeOutput(const string& output_path, vector<double> thetas, Point center, double center_range) {
	
		ofstream outfile(output_path);
		// first read center
		outfile << center.x;
		outfile << "\n";
		outfile << center.y;
		outfile << "\n" ;
		outfile << center_range;
		outfile << "\n";
		
		for (auto theta : thetes) {
			outfile << theta;
			outfile << "\n";
		}	
}

int main(int argc, char** argv) {
	
	if (argc != 3) {
		cout << "GIVE THE IMAGE FILE AND CONFIGURATION FILE OUTPUT!" << endl;
		return -1;
	}
	
	string image_file = argv[1];
	string config_file = argv[2];
	
	Mat img = imread(image_file, CV_LOAD_IMAGE_GRAYSCALE);
	namedWindow(winName, WINDOW_AUTOSIZE);
	imshow(winName, img);
	
	Point center(1,2);
	double center_range = 5;
	vector<double> thetas = {1.1, 2.2, 3.3};
	waitKey(0);
	writeOutput(config_file, thetas, center, center_range);
}

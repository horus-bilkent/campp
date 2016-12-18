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

using namespace std;

int readInput(const string& input_path, vector<double>& thetas, Point& center, double& center_range) {
	
	try {
		std::ifstream infile(input_path);
		std::string line;
		// first read center
		double center_x, center_y;
		std::getline(infile, line);
		center_x = std::stod(line);
		std::getline(infile, line);
		center_y = std::stod(line);
		center.x = center_x;
		center.y = center_y;
		
		// read the center range
		std::getline(infile, line);
		center_range = std::stod(line);
		
		// read the thetas
		while  (getline(infile, line) ) {
			if (line.empty()) {
				continue;
			}
			double theta_val = std::stod(line);
			thetas.push_back(theta_val);
		}
		return 0;	
	} catch (std::exception& e) {
		cout << string(e.what()) << endl;
		return -1;
	}

}

int main(int argc, char** argv) {
	
	if (argc != 2) {
		cout << "GIVE THE CONFIGURATION FILE!" << endl;
		return -1;
	}
	
	vector<double> thetas;
	Point center;
	double center_range;
	
	if (readInput(argv[1], thetas, center, center_range) != 0) {
		cout << "Error!" << endl;
		return -1;
	}
	
	// omer add your code
	
	// hersey bitince, python -> stdin -> readNeedlevalue -> stdout -> python -> kafka (loop)

	string line;
	while (std::cin >> line) {
		// double needle_value = readNeedleValue();  // ya da int hangisi olursa
		
		// to try
		double needle_value = center_range;
		std::cout << needle_value << std::endl;
		std::cout.flush();
	}


}

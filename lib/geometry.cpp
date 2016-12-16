#include "geometry.hpp"

double angleDiff(double a, double b) {
    double res = a - b;
    if(a - b < 0)
        return res + CV_PI * 2;
    return res;
}

bool lineCircleIntersect(Point2f a, Point2f b, Point2f center, double range) {
	float dx = b.x - a.x;
	float dy = b.y - a.y;
	float A = dx * dx + dy * dy;
	float B = 2 * (dx * (a.x - center.x) + dy * (a.y - center.y));
	float C = center.x * center.y + center.y * center.y;
	C += a.x * a.x + a.y * a.y;
	C -= 2 * (center.x * a.x + center.y * a.y);
	C -= range * range;
	float bb4ac = B * B - 4 * A * C;

	if(bb4ac<0) {
		return false;    // No collision
	}
	else {
		return true;      //Collision
	}
}

double calcRange2(double big_range, double small_range, double thickness,
                       double cols) {
	double temp = small_range - thickness / 2;
	return big_range * (temp / cols);
}

bool lineLineIntersect(Vec4i a, Vec4i b, Point &res) {
	double A1 = a[3] - a[1], B1 = a[0] - a[2], C1 = A1*a[0] + B1 * a[1];
	double A2 = b[3] - b[1], B2 = b[0] - b[2], C2 = A2*b[0] + B2 * b[1];
	double det = A1 * B2 - A2 * B1;
	if(det == 0)
		return false;
	res.x = (B2 * C1 - B1 * C2) / det;
	res.y = (A1 * C2 - A2 * C1) / det;
	return true;
}

void minCovCircleIt(Point &center, double &rad, Point a, Point b) {
	double dist = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
	if(rad == -1 || dist >= rad) {
		rad = dist;
		center.x = (a.x + b.x) / 2.0;
		center.y = (a.y + b.y) / 2.0;
	}
}

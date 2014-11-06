#include "CubicInterpolation.h"
#include <vector>

using namespace std;

int main()
{

	CubicInterpolation *cubicInterpolation = new CubicInterpolation();

	vector<float> xi;
	vector<float> yi;

	//横坐标
	xi.push_back(1);
	xi.push_back(2);
	xi.push_back(3);

	//纵坐标
	yi.push_back(3);
	yi.push_back(6);
	yi.push_back(3);

	xiSize = xi.size();

	cubicInterpolation->initVector(xi, yi, 0, 0, xiSize);
	cubicInterpolation->calcCoefs();
	//计算横坐标为2.5的值
	float result = cubicInterpolation->evaluate(2.5);
	
	delete cubicInterpolation;
	return 0;
}
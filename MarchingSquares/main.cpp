#include "MarchingSquares.h"
#include <vector>

using namespace std;

int main()
{
	//用于存储等值线生成结果，相同值的等值线会放在一起。
	vector<vector<isotools::Isoline>> pathLinesV;

	//pData是网格数据
	//isoValue是需要生成的等值线的等值
	//maxValue和minValue是格点数据中的最大值和最小值
	//处理完成后，结果会存放在pathLinesV
	marchingsquares::doMarchingSquaresAccelerateOMP(pData, isoValue, pathLinesV, startLongitude,
		longitudeGridSpace, startLatitude, latitudeGridSpace, maxValue, minValue);

	return 0;
}
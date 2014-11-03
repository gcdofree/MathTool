#pragma once
#include <vector>

/*
* Author: gcdofree
* Date: 2014.11.3
* Description: 三次样条插值的实现
*/

class CubicInterpolation
{
public:
	CubicInterpolation();
	~CubicInterpolation();

	float Ml, Mr;
	int originSize;
	std::vector<float> xi, yi;
	std::vector<std::vector<float>> coefs;

	/* 初始化插值参数
	* xi代表所有节点的x坐标，必须依次递增，不得递减或相等
	* yi代表所有节点的y坐标
	* ml为左侧的二阶导数
	* mr为右侧的二阶导数
	* oSize是初始节点的个数，此处至少为3个
	*/
	void initVector(vector<float> &xi, vector<float> &yi, float ml, float mr, int oSize);

	/* 计算各种矩阵系数 */
	void calcCoefs();

	/* 计算某点在三次样条曲线上的坐标
	* x为该点的x坐标
	*/
	float evaluate(float x);

	/* 计算2阶导数 */
	void derivative2(vector<float> &dx, vector<float> &d1, vector<float> &d2);
};


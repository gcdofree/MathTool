#include "CubicInterpolation.h"
#include <iostream>

/*
* Author: gcdofree
* Date: 2014.11.3
* Description: 三次样条插值的实现
*/

CubicInterpolation::CubicInterpolation()
{
	Ml = 0;
	Mr = 0;
	originSize = 0;
}


CubicInterpolation::~CubicInterpolation()
{
}

/* 初始化插值参数
* xi代表所有节点的x坐标，必须依次递增，不得递减或相等
* yi代表所有节点的y坐标
* ml为左侧的二阶导数
* mr为右侧的二阶导数
* oSize是初始节点的个数，此处至少为3个
*/
void CubicInterpolation::initVector(vector<float> &xxi, vector<float> &yyi, float ml, float mr, int oSize)
{
	xi.clear();
	yi.clear();
	xi.resize(oSize);
	yi.resize(oSize);
	for (int i = 0; i < oSize; ++i)
	{
		xi[i] = xxi[i];
		yi[i] = yyi[i];
	}
	originSize = oSize;
	this->Ml = ml;
	this->Mr = mr;
}

/* 计算各种矩阵系数 */
void CubicInterpolation::calcCoefs()
{
	int N = originSize,
		M = N - 1;

	vector<float> m(N),
		h(M),
		d(M);

	m[0] = Ml;
	m[M] = Mr;
	for (int i = 0; i<M; ++i)
	{
		h[i] = xi[i + 1] - xi[i];
		d[i] = (yi[i + 1] - yi[i]) / h[i];
	}

	derivative2(h, d, m);

	coefs.resize(M);
	for (int i = 0; i < M; ++i)
	{
		coefs[i].resize(4);
 		coefs[i][0] = yi[i];
		coefs[i][1] = d[i] - h[i] * (2 * m[i] + m[i + 1]) / 6;
		coefs[i][2] = m[i] / 2;
		coefs[i][3] = (m[i + 1] - m[i]) / (6 * h[i]);
	}
}

/* 计算某点在三次样条曲线上的坐标
* x为该点的x坐标
*/
float CubicInterpolation::evaluate(float x)
{
	int k = -1,
		N = originSize,
		M = N - 1;

	float dx,
		y;

	for (int i = 0; i < M; ++i)
	{
		if ((xi[i] <= x) && (xi[i + 1] >= x))
		{
			k = i;
			dx = x - xi[i];
			break;
		}
	}
	if (k != -1)
	{
		y = ((coefs[k][3] * dx + coefs[k][2]) * dx + coefs[k][1]) * dx
			+ coefs[k][0];
		return y;
	}
	else
	{
		cout << "The x value is out of range!" << endl;
		return float(0);
	}
}

/* 计算2阶导数 */
void CubicInterpolation::derivative2(vector<float> &dx, vector<float> &d1, vector<float> &d2)
{
	int N = originSize,
		M = N - 1;
	vector<float> b(M),
		v(M),
		y(M),
		alpha(M),
		beta(M - 1);

	for (int i = 1; i < M; ++i)
		b[i] = 2 * (dx[i] + dx[i - 1]);

	v[1] = 6 * (d1[1] - d1[0]) - dx[0] * d2[0];
	for (int i = 1; i < M - 1; ++i)
		v[i] = 6 * (d1[i] - d1[i - 1]);
		
	v[M - 1] = 6 * (d1[M - 1] - d1[M - 2]) - dx[M - 1] * d2[M];

	alpha[1] = b[1];
	for (int i = 2; i < M; ++i)
		alpha[i] = b[i] - dx[i] * dx[i - 1] / alpha[i - 1];
		

	for (int i = 1; i < M - 1; ++i)
		beta[i] = dx[i] / alpha[i];
		

	y[1] = v[1] / alpha[1];
	for (int i = 2; i<M; ++i)
		y[i] = (v[i] - dx[i] * y[i - 1]) / alpha[i];
		

	d2[M - 1] = y[M - 1];
	for (int i = M - 2; i>0; --i)
		d2[i] = y[i] - beta[i] * d2[i + 1];
}
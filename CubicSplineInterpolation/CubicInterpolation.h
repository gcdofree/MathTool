#pragma once
#include <vector>

/*
* Author: gcdofree
* Date: 2014.11.3
* Description: ����������ֵ��ʵ��
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

	/* ��ʼ����ֵ����
	* xi�������нڵ��x���꣬�������ε��������õݼ������
	* yi�������нڵ��y����
	* mlΪ���Ķ��׵���
	* mrΪ�Ҳ�Ķ��׵���
	* oSize�ǳ�ʼ�ڵ�ĸ������˴�����Ϊ3��
	*/
	void initVector(vector<float> &xi, vector<float> &yi, float ml, float mr, int oSize);

	/* ������־���ϵ�� */
	void calcCoefs();

	/* ����ĳ�����������������ϵ�����
	* xΪ�õ��x����
	*/
	float evaluate(float x);

	/* ����2�׵��� */
	void derivative2(vector<float> &dx, vector<float> &d1, vector<float> &d2);
};


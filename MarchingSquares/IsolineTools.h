#pragma once
#include <list>

using namespace std;

namespace isotools
{
	class Point2D 
	{
	public:
		bool operator==( const Point2D &right ) const
		{
			return ( *this ).x == right.x && ( *this ).y == right.y;
		}
		float x;//横轴坐标值
		float y;//纵轴坐标值
	};
	
	/**  等值线数据结构，保存一条等值线 **/
	struct Isoline
	{
		float isovalue;//等值线的值
		bool isCircle;//等值线是否构成环
		bool isBorder;//等值线是否在边界
		list< Point2D > points;//该等值线上所有的点
		Point2D startPoint;//该等值线 头结点 所在边的中点数组索引下标
		Point2D endPoint;//该等值线 尾结点 所在边的中点数组索引下标
	};

	/**  等值线数据结构，保存一条初始的边 **/
	struct Edge
	{
		int edgeIndex1;//第一个交点
		int edgeIndex2;//第二个交点
		int i;//所在网格位置的横坐标
		int j;//所在网格位置的纵坐标
		int m;//对应的等值线值的index
		float isovalue;//对应的等值线值
	};
}
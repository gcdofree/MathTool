#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include "IsolineTools.h"

using namespace std;
/************************************************************************/ 
/* Author: gcdofree
 * Date: 2014.11.3
 * Description: Marching Squares方法的实现
/************************************************************************/
namespace marchingsquares
{
/*   1     切点所在边的标号
  0     2
	 3
*/   
const int EdgeTable[ 16 ] = {
	0x0,     //0000
	0x9,     //1001
	0x3,     //0011
	0xa,     //1010

	0x6,     //0110,
	0xf,     //1111
	0x5,     //0101
	0xc,     //1100

	0xc,     //1100
	0x5,     //0101
	0xf,     //1111
	0x6,     //0110

	0xa,     //1010
	0x3,     //0011
	0x9,     //1001
	0x0,     //0000
};

 // 与上面的EdgeTable值相对应
const int SegmentTable[ 16 ][ 5 ] = {
	{ -1, -1, -1, -1, -1 },
	{ 0, 3, -1, -1, -1 },
	{ 2, 3, -1, -1, -1 },
	{ 0, 2, -1, -1, -1 },

	{ 1, 2, -1, -1, -1 },
	{ 0, 1, 2, 3, -1 }, //index =5 ,  0-1  2-3 
	{ 1, 3, -1, -1, -1 },
	{ 0, 1, -1, -1, -1 },

	{ 0, 1, -1, -1, -1 },
	{ 1, 3, -1, -1, -1 },
	{ 0, 3, 1, 2, -1 },//index =10 ,  0-3  1-2  
	{ 1, 2, -1, -1, -1 },

	{ 0, 2, -1, -1, -1 },
	{ 2, 3, -1, -1, -1 },
	{ 3, 0, -1, -1, -1 },
	{ -1, -1, -1, -1, -1 }
};

const float V_EPSILON = 1.0f / 256.0f;//阈值，用来判断切点和格点的距离

/************************************************************************/
/* Funciton:   VertexInterp     
 * Description: 插值函数
 * Input:
	isovalue: 等值线值
	x1: 第一个顶点的数组 x 下标
	y1: 第一个顶点的数组 y 下标
	v1:  第一个顶点值
	x2: 第二个顶点的数组 x 下标
	y2: 第二个顶点的数组 y 下标
	v2:  第二个顶点值
 * Output: 包含插值点的坐标及等值线值的Point2D对象
 * Author: gcdofree
 * Date: 2014.11.3
************************************************************************/
static isotools::Point2D VertexInterp( float isovalue, int x1, int y1, float v1, int x2, int y2, float v2, float startLongitude, float longitudeGridSpace, float startLatitude, float latitudeGridSpace )
{
	float mu;
	isotools::Point2D p;

	if ( fabs( v1 - isovalue ) < V_EPSILON )
	{
		p.x = y1 * longitudeGridSpace + startLongitude;
		p.y = x1 * latitudeGridSpace + startLatitude;
		return p;
	}
	else if ( fabs( v2 - isovalue ) < V_EPSILON )
	{
		p.x = y2 * longitudeGridSpace + startLongitude;
		p.y = x2 * latitudeGridSpace + startLatitude;
		return p;
	}

	//下面是边顶点值一个比isovalue大，一个比isovalue小的情况
	mu = ( isovalue - v1 ) / ( v2 - v1) ;
	p.x = ( y1 + mu * ( y2 - y1 ) ) * longitudeGridSpace + startLongitude;//由数组下标换成经纬度坐标系
	p.y = ( x1 + mu * ( x2 - x1 ) ) * latitudeGridSpace + startLatitude;//由数组下标换成经纬度坐标系

	return p;
}
 
 /************************************************************************/
/* Funciton:   getMiddlePoint     
 * Description: 获得交点所在边的中点（其坐标相当于该边的全局编号）
 * Input:
	edgeIndex: 交点所在边在cell中的局部编号（0,1,2,3）
	i: 交点所在cell中左上角点的数组x值下标
	j:交点所在cell中左上角点的数组y值下标
 * Output: isotools::Point2D  返回该边的中点
 * Author: gcdofree
 * Date: 2014.11.3
/************************************************************************/
static isotools::Point2D getMiddlePoint( int edgeIndex, int i,int j )
{
	isotools::Point2D p;
	switch ( edgeIndex )
	{
	case 0:
		p.x = i + 0.5;
		p.y = j;
		break;
	case 1:
		p.x = i;
		p.y = j + 0.5;
		break;
	case 2:
		p.x = i + 0.5;
		p.y = j + 1;
		break;
	case 3:
		p.x = i + 1;
		p.y = j + 0.5;
		break;
	default:
		break;
	}
	return p;
}

 /************************************************************************/
/* Funciton: getCutPoint     
 * Description: 通过插值获得交点的坐标
 * Input:
	edgeIndex: 交点所在边在cell中的局部编号（0,1,2,3）
	i:交点所在cell中左上角点的数组x值下标
	j:交点所在cell中左上角点的数组y值下标
	data:天气数据值二维数组
	isovalue:等值
	startLongitude: 起始经度（起始x坐标）
	longitudeGridSpace: 经度间隔（x坐标间隔）
	startLatitude: 起始纬度（起始y坐标）
	latitudeGridSpace: 纬度间隔（y坐标间隔）
 * Output: isotools::Point2D  返回该交点的坐标
 * Author: gcdofree
 * Date: 2014.11.3
/************************************************************************/
static isotools::Point2D getCutPoint( int edgeIndex, int i, int j, vector< vector< float > > &data, float isovalue,
	float startLongitude, float longitudeGridSpace, float startLatitude, float latitudeGridSpace )
{
	switch( edgeIndex )
	{
	case 0:
		return VertexInterp( isovalue, i + 1, j, data[ i + 1 ][ j ], i, j, data[ i ][ j ], startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace );
	case 1:
		return VertexInterp( isovalue, i, j, data[ i ][ j ], i, j + 1, data[ i ][ j + 1 ], startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace );
	case 2:
		return VertexInterp( isovalue, i, j + 1, data[ i ][ j + 1 ], i + 1, j + 1, data[ i + 1 ][ j + 1 ], startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace );
	case 3:
		return VertexInterp( isovalue, i + 1, j + 1, data[ i + 1 ][ j + 1 ], i + 1, j, data[ i + 1 ][ j ], startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace );
	}
}

/************************************************************************/
/* Funciton: isMergeIsoLineAccelerate
 * Description: 进行线段合并操作，
 * Input:
	mid: 当前需要处理的节点
	type: 0表示是头结点 ,1表示尾节点
	pathLines: 等值线集合
	m: 当前节点所在的等值线编号
	pathlineHit: 当删除线段后，pathlineHit需要进行更新
 * Output: 若有线段合并，则返回true，否则返回false
 * Author: gcdofree
 * Date: 2014.11.3
*/
static bool isMergeIsoLineAccelerate(isotools::Point2D mid, int type, vector<isotools::Isoline> &pathLines, int &m, vector<int> &pathlineHit)
{
	int pathlineSize = pathLines.size();

	for (int i = 0; i<pathlineSize; ++i)
	{
		if (i == m)continue;

		if (mid == pathLines[i].startPoint && type == 0)//两条线的 头结点 重合，则将短的线段加到长的线段，并删除短的线段
		{
			if (pathLines[m].points.size() > pathLines[i].points.size())//线段m较长，将i中的点移动到m中
			{
				// list不支持随机读取，不能通过下标访问
				list<isotools::Point2D>::iterator pIt = pathLines[i].points.begin();//指向头结点
				list<isotools::Point2D>::iterator pIt_end = pathLines[i].points.end();//指向头结点
				++pIt;//跳过头结点
				isotools::Point2D point;
				while (pIt != pIt_end)
				{
					point = *pIt;//取值
					pathLines[m].points.push_front(point);//在头部依次插入
					++pIt;
				}
				//修改m的头结点
				pathLines[m].startPoint = pathLines[i].endPoint;
				pathLines.erase(pathLines.begin() + i);//删除线段pathLines[i]
				if (i < m)
					--m;
				int pathlineHitSize = pathlineHit.size();
				for (int j = 0; j < pathlineHitSize; ++j)
				{
					if (pathlineHit[j] > i)
					{
						pathlineHit[j]--;
					}
					else if (pathlineHit[j] == i)
					{
						pathlineHit[j] = m;
					}

				}
				
			}
			else //线段i较长，将m中的点移动到i中
			{
				list<isotools::Point2D>::iterator pIt = pathLines[m].points.begin();//指向头结点
				list<isotools::Point2D>::iterator pIt_end = pathLines[m].points.end();//指向头结点
				++pIt;//跳过头结点
				isotools::Point2D point;
				while (pIt != pIt_end)
				{
					point = *pIt;//取值
					pathLines[i].points.push_front(point);//在头部依次插入
					++pIt;
				}
				//修改i的头结点
				pathLines[i].startPoint = pathLines[m].endPoint;
				pathLines.erase(pathLines.begin() + m);//删除线段pathLines[m]

				int pathlineHitSize = pathlineHit.size();
				for (int j = 0; j < pathlineHitSize; ++j)
				{
					if (pathlineHit[j] > m)
					{
						pathlineHit[j]--;
					}
					else if (pathlineHit[j] == m)
					{
						if (i < m)
							pathlineHit[j] = i;
						else
							pathlineHit[j] = i - 1;
					}

				}

				if (i < m)
					m = i;
				else
					m = i - 1;
			}
			return true;
		}
		else if (mid == pathLines[i].endPoint && type == 1)//两条线的 尾结点 重合，则将短的线段加到长的线段，并删除短的线段
		{
			if (pathLines[m].points.size() > pathLines[i].points.size())//线段m较长，将i中的点移动到m中
			{
				list<isotools::Point2D>::iterator pIt = pathLines[i].points.end();//指向尾结点
				list<isotools::Point2D>::iterator pIt_beg = pathLines[i].points.begin();
				--pIt;//跳过end
				--pIt;//跳过尾结点
				isotools::Point2D point;
				while (pIt != pIt_beg)
				{
					point = *pIt;//取值
					pathLines[m].points.push_back(point);//在尾部依次插入
					--pIt;
				}
				//再push头结点
				point = *pIt;//取值
				pathLines[m].points.push_back(point);//在尾部插入
				//修改m的尾结点
				pathLines[m].endPoint = pathLines[i].startPoint;
				pathLines.erase(pathLines.begin() + i);//删除线段pathLines[i]
				if (i < m)
					--m;
				int pathlineHitSize = pathlineHit.size();
				for (int j = 0; j < pathlineHitSize; ++j)
				{
					if (pathlineHit[j] > i)
					{
						pathlineHit[j]--;
					}
					else if (pathlineHit[j] == i)
					{
						pathlineHit[j] = m;
					}

				}
			}
			else //线段i较长，将m中的点移动到i中
			{
				list<isotools::Point2D>::iterator pIt = pathLines[m].points.end();//指向尾结点
				list<isotools::Point2D>::iterator pIt_beg = pathLines[m].points.begin();
				--pIt;//跳过end
				--pIt;//跳过尾结点
				isotools::Point2D point;
				while (pIt != pIt_beg)
				{
					point = *pIt;//取值
					pathLines[i].points.push_back(point);//在尾部依次插入
					--pIt;
				}
				//再push头结点
				point = *pIt;//取值
				pathLines[i].points.push_back(point);//在尾部插入
				//修改i的尾结点
				pathLines[i].endPoint = pathLines[m].startPoint;
				pathLines.erase(pathLines.begin() + m);//删除线段pathLines[m]

				int pathlineHitSize = pathlineHit.size();
				for (int j = 0; j < pathlineHitSize; ++j)
				{
					if (pathlineHit[j] > m)
					{
						pathlineHit[j]--;
					}
					else if (pathlineHit[j] == m)
					{
						if (i < m)
							pathlineHit[j] = i;
						else
							pathlineHit[j] = i - 1;
					}

				}

				if (i < m)
					m = i;
				else
					m = i - 1;
			}
			return true;
		}
		else if (mid == pathLines[i].endPoint && type == 0)//m线的头结点与i的尾结点 重合，则将m加到i后
		{
			//将m中的点移动到i中
			list<isotools::Point2D>::iterator pIt = pathLines[m].points.begin();//指向头结点
			list<isotools::Point2D>::iterator pIt_end = pathLines[m].points.end();
			++pIt;//跳过头结点
			isotools::Point2D point;
			while (pIt != pIt_end)
			{
				point = *pIt;//取值
				pathLines[i].points.push_back(point);//在尾部依次插入
				++pIt;
			}
			//修改i的尾结点
			pathLines[i].endPoint = pathLines[m].endPoint;
			pathLines.erase(pathLines.begin() + m);//删除线段pathLines[m]

			int pathlineHitSize = pathlineHit.size();
			for (int j = 0; j < pathlineHitSize; ++j)
			{
				if (pathlineHit[j] > m)
				{
					pathlineHit[j]--;
				}
				else if (pathlineHit[j] == m)
				{
					if (i < m)
						pathlineHit[j] = i;
					else
						pathlineHit[j] = i - 1;
				}

			}

			if (i < m)
				m = i;
			else
				m = i - 1;
			return true;
		}
		else if (mid == pathLines[i].startPoint && type == 1)//m线的尾结点与i的头结点 重合，则将i加到m后
		{
			//将m中的点移动到i中
			list<isotools::Point2D>::iterator pIt = pathLines[i].points.begin();//指向头结点
			list<isotools::Point2D>::iterator pIt_end = pathLines[i].points.end();
			++pIt;//跳过头结点
			isotools::Point2D point;
			while (pIt != pIt_end)
			{
				point = *pIt;//取值
				pathLines[m].points.push_back(point);//在尾部依次插入
				++pIt;
			}
			//修改m的头结点
			pathLines[m].endPoint = pathLines[i].endPoint;
			pathLines.erase(pathLines.begin() + i);//删除线段pathLines[i]
			if (i < m)
				--m;
			int pathlineHitSize = pathlineHit.size();
			for (int j = 0; j < pathlineHitSize; ++j)
			{
				if (pathlineHit[j] > i)
				{
					pathlineHit[j]--;
				}
				else if (pathlineHit[j] == i)
				{
					pathlineHit[j] = m;
				}

			}
			return true;
		}
	}
	return false;
}

/************************************************************************/
/* Funciton: setHitTableAccelerate
 * Description:  更新查找表
 * Input:
	hitnum: 命中的等值线编号
	hitsize: 查找表大小
	pathLinesHit:该等值线对应的查找表
 * Output: void
 * Author: gcdofree
 * Date: 2014.11.3
/************************************************************************/
static void setHitTableAccelerate(int hitnum, int hitsize, vector< int > &pathLinesHit)
{
	if (pathLinesHit.size() == 0)
	{
		pathLinesHit.push_back(hitnum);
		return;
	}
	if (pathLinesHit.size() < hitsize)
	{
		if (hitnum == pathLinesHit[0])
		{
			return;
		}
		pathLinesHit.push_back(hitnum);
		for (int i = pathLinesHit.size() - 1; i > 0; i--)
		{
			pathLinesHit[i] = pathLinesHit[i - 1];
		}
		pathLinesHit[0] = hitnum;
	}
	else
	{
		if (hitnum == pathLinesHit[0])
		{
			return;
		}
		for (int i = hitsize - 1; i > 0; i--)
		{
			pathLinesHit[i] = pathLinesHit[i - 1];
		}
		pathLinesHit[0] = hitnum;
	}
}

/************************************************************************/
/* Funciton: addPointToLineAccelerate
 * Description:  判断当前两边上的点是否与已有等值线的首尾点重合，若是，则忽略此重合点，再将另一个点插入到首或尾并判断加入后是否与已有其他线段端点重合，以进行合并操作；否则，新增一条等值线（增加最近查找表）
 * Input:
	edgeIndex1: 交点1所在边在cell中的局部编号（0,1,2,3）
	edgeIndex2: 交点2所在边在cell中的局部编号（0,1,2,3） （这两个交点构成了cell中的一条等值线）
	i: 交点所在cell中左上角点的数组x值下标
	j: 交点所在cell中左上角点的数组y值下标
	data: 天气数据值二维数组
	pathLinesHit: 等值对应的查找表
	isovalue: 等值
	pathLines: 返回该等值对应的各条等值线vector （每条保存在一个Isoline）
	startLongitude: 起始经度（起始x坐标）
	longitudeGridSpace: 经度间隔（x坐标间隔）
	startLatitude: 起始纬度（起始y坐标）
	latitudeGridSpace: 纬度间隔（y坐标间隔）
 * Output: void
 * Author: gcdofree
 * Date: 2014.11.3
/************************************************************************/
static void addPointToLineAccelerate(int edgeIndex1, int edgeIndex2, int i, int j, vector< vector< float > > &data, vector< int > &pathLinesHit,
	float isovalue, vector< isotools::Isoline > &pathLines, float startLongitude, float longitudeGridSpace, float startLatitude, float latitudeGridSpace)
{
	//判断当前边上的点是否与已有等值线的首尾点重合，若是，则忽略此重合点，再将另一个点插入到首或尾，否则，新增一条等值线
	isotools::Point2D mid1 = getMiddlePoint(edgeIndex1, i, j);//存放数组索引
	isotools::Point2D mid2 = getMiddlePoint(edgeIndex2, i, j);//存放数组索引

	int hitSize = 5;
	int pathlineHitSize = pathLinesHit.size();
	for (int k = 0; k < pathlineHitSize; ++k)
	{
		int m = pathLinesHit[k];

		if (pathLines[m].isCircle)
		{
			continue;//对于构成环的不予判断
		}
		//若边上的点与首尾结点都相同
		if ((mid1 == pathLines[m].startPoint && mid2 == pathLines[m].endPoint) || (mid2 == pathLines[m].startPoint && mid1 == pathLines[m].endPoint))
		{
			pathLines[m].isCircle = true;//形成回路
			return;
		}
		else if (mid1 == pathLines[m].startPoint)//若m1与头结点重合
		{
			//则在头结点前插入另一个点，并更新头结点
			isotools::Point2D point = getCutPoint(edgeIndex2, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
			pathLines[m].points.push_front(point);
			pathLines[m].startPoint = mid2;
			isMergeIsoLineAccelerate(mid2, 0, pathLines, m, pathLinesHit);//进行合并操作
			setHitTableAccelerate(m, hitSize, pathLinesHit);
			return;
		}
		else if (mid1 == pathLines[m].endPoint)//若m1与尾结点重合
		{
			isotools::Point2D point = getCutPoint(edgeIndex2, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
			pathLines[m].points.push_back(point);
			pathLines[m].endPoint = mid2;
			isMergeIsoLineAccelerate(mid2, 1, pathLines, m, pathLinesHit);//进行合并操作
			setHitTableAccelerate(m, hitSize, pathLinesHit);
			return;
		}
		else if (mid2 == pathLines[m].startPoint)//若m2与头结点重合
		{
			//则在头结点前插入另一个点，并更新头结点
			isotools::Point2D point = getCutPoint(edgeIndex1, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
			pathLines[m].points.push_front(point);
			pathLines[m].startPoint = mid1;
			isMergeIsoLineAccelerate(mid1, 0, pathLines, m, pathLinesHit);//进行合并操作
			setHitTableAccelerate(m, hitSize, pathLinesHit);
			return;
		}
		else if (mid2 == pathLines[m].endPoint)//若m2与尾结点重合
		{
			isotools::Point2D point = getCutPoint(edgeIndex1, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
			pathLines[m].points.push_back(point);
			pathLines[m].endPoint = mid1;
			isMergeIsoLineAccelerate(mid1, 1, pathLines, m, pathLinesHit);//进行合并操作
			setHitTableAccelerate(m, hitSize, pathLinesHit);
			return;
		}
	}
	
	int m = 0;
	int pathLinesSize = pathLines.size();
	for (m = 0; m < pathLinesSize; ++m)
	{
		if (pathLines[m].isCircle)
		{
			continue;//对于构成环的不予判断
		}
		//若边上的点与首尾结点都相同
		if ((mid1 == pathLines[m].startPoint && mid2 == pathLines[m].endPoint) || (mid2 == pathLines[m].startPoint && mid1 == pathLines[m].endPoint))
		{
			//isotools::Point2D point = pathLines[ m ].points.front();//返回第一个元素  【不再多插一个元素,而是只改isCircle标记
			//pathLines[ m ].points.push_back( point );//多插入第一个元素至末尾，形成回路  【如点1,2,3构成回路，则只保存点1,2,3，不再保存1,2,3,1】
			pathLines[m].isCircle = true;//形成回路
			return;
		}
		else if (mid1 == pathLines[m].startPoint)//若m1与头结点重合
		{
			//则在头结点前插入另一个点，并更新头结点
			isotools::Point2D point = getCutPoint(edgeIndex2, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
			pathLines[m].points.push_front(point);
			pathLines[m].startPoint = mid2;
			isMergeIsoLineAccelerate(mid2, 0, pathLines, m, pathLinesHit);//进行合并操作
			setHitTableAccelerate(m, hitSize, pathLinesHit);
			return;
		}
		else if (mid1 == pathLines[m].endPoint)//若m1与尾结点重合
		{
			isotools::Point2D point = getCutPoint(edgeIndex2, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
			pathLines[m].points.push_back(point);
			pathLines[m].endPoint = mid2;
			isMergeIsoLineAccelerate(mid2, 1, pathLines, m, pathLinesHit);//进行合并操作
			setHitTableAccelerate(m, hitSize, pathLinesHit);
			return;
		}
		else if (mid2 == pathLines[m].startPoint)//若m2与头结点重合
		{
			//则在头结点前插入另一个点，并更新头结点
			isotools::Point2D point = getCutPoint(edgeIndex1, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
			pathLines[m].points.push_front(point);
			pathLines[m].startPoint = mid1;
			isMergeIsoLineAccelerate(mid1, 0, pathLines, m, pathLinesHit);//进行合并操作
			setHitTableAccelerate(m, hitSize, pathLinesHit);
			return;
		}
		else if (mid2 == pathLines[m].endPoint)//若m2与尾结点重合
		{
			isotools::Point2D point = getCutPoint(edgeIndex1, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
			pathLines[m].points.push_back(point);
			pathLines[m].endPoint = mid1;
			isMergeIsoLineAccelerate(mid1, 1, pathLines, m, pathLinesHit);//进行合并操作
			setHitTableAccelerate(m, hitSize, pathLinesHit);
			return;
		}
	}//end for m
	if (m == pathLinesSize)//说明这条边没有与任何一条等值线相接，则新创建一条等值线
	{
		isotools::Point2D point1 = getCutPoint(edgeIndex1, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
		isotools::Point2D point2 = getCutPoint(edgeIndex2, i, j, data, isovalue, startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);

		isotools::Isoline isoList;
		isoList.points.push_front(point1);
		isoList.startPoint = mid1;
		isoList.points.push_back(point2);
		isoList.endPoint = mid2;
		isoList.isovalue = isovalue;
		isoList.isCircle = false;
		isoList.isBorder = false;
		pathLines.push_back(isoList);//新增一条等值线
	}
}

/************************************************************************/
/* Funciton: doMarchingSquaresAccelerate 【普通CPU串行算法】
 * Description: Marching Squares 算法的实现
 * Input:
	data: 天气数据值二维数组
	isovalues: 等值线值数组
	pathLinesV: 返回该等值数组下的所有各条等值线数组的集合
	startLongitude: 起始经度（起始x坐标）
	longitudeGridSpace: 经度间隔（x坐标间隔）
	startLatitude: 起始纬度（起始y坐标）
	latitudeGridSpace: 纬度间隔（y坐标间隔）
	maxGridValue: 网格点中的最大值
	minGridValue: 网格点中的最小值
* Output: void
* Author: gcdofree
* Date: 2014.11.3
/************************************************************************/
static void doMarchingSquaresAccelerate(vector<vector<float> > &data, vector<float> &isovalues, vector< vector<isotools::Isoline>> &pathLinesV,
	float startLongitude, float longitudeGridSpace, float startLatitude, float latitudeGridSpace, float maxGridValue, float minGridValue)
{

	pathLinesV.clear();

	vector<vector<int>> pathLinesHit;

	for (int m = 0; m < isovalues.size(); ++m)
	{
		if (isovalues[m] < minGridValue || isovalues[m] > maxGridValue)
		{
			isovalues.erase(isovalues.begin() + m);//对vector进行增删元素后，begin，end，size会失效，必须即时调用
			m--;
		}
	}
	int isovaluesNum = isovalues.size();
	pathLinesV.resize(isovaluesNum);
	pathLinesHit.resize(isovaluesNum);

	int dataSize_i = data.size() - 1;
	for (int i = 0; i<dataSize_i; ++i)//逐行扫
	{
		int dataSize_j = data[i].size() - 1;
		for (int j = 0; j<dataSize_j; ++j)//逐列扫
		{

			//预先算好邻域四个格点的最大最小值，避免重复计算squareIndexs[m]
			float maxValue = 0, minValue = FLT_MAX;
			if (data[i][j] >maxValue) maxValue = data[i][j];
			if (data[i][j+1] >maxValue) maxValue = data[i][j+1];
			if (data[i+1][j+1] >maxValue) maxValue = data[i+1][j+1];
			if (data[i+1][j] >maxValue) maxValue = data[i+1][j];
			if (data[i][j] < minValue) minValue = data[i][j];
			if (data[i][j + 1] < minValue) minValue = data[i][j + 1];
			if (data[i + 1][j + 1] < minValue) minValue = data[i + 1][j + 1];
			if (data[i + 1][j] < minValue) minValue = data[i + 1][j];

			for (int m = 0; m<isovaluesNum; ++m)
			{
				if (isovalues[m] < minValue || isovalues[m] > maxValue)
				{
					continue;
				}

				int squareIndex = 0;
				if (data[i][j] >= isovalues[m]) squareIndex |= 8;
				if (data[i][j + 1] >= isovalues[m]) squareIndex |= 4;
				if (data[i + 1][j + 1] >= isovalues[m]) squareIndex |= 2;
				if (data[i + 1][j] >= isovalues[m]) squareIndex |= 1;

				if (EdgeTable[squareIndex] != 0)
				{
					//消除二义性
					if (squareIndex == 5)// index = 0x5 
					{
						float centerValue = (1 / 4.0) *(data[i][j] + data[i][j + 1] + data[i + 1][j + 1] + data[i + 1][j]);
						if (centerValue < isovalues[m])
						{
							squareIndex = 10;
						}
					}
					else if (squareIndex == 10)//0x10
					{
						float centerValue = (1 / 4.0) *(data[i][j] + data[i][j + 1] + data[i + 1][j + 1] + data[i + 1][j]);
						if (centerValue < isovalues[m])
						{
							squareIndex = 5;
						}
					}
					int edgeIndex1 = 0, edgeIndex2 = 0;
					//将生成的该线段添加到某一等值线上或自己建立某条等值线
					for (int k = 0; SegmentTable[squareIndex][k] != -1; k = k + 2)
					{
						edgeIndex1 = SegmentTable[squareIndex][k];
						edgeIndex2 = SegmentTable[squareIndex][k + 1];
						addPointToLineAccelerate(edgeIndex1, edgeIndex2, i, j, data, pathLinesHit[m], isovalues[m], pathLinesV[m], startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
					}
				}
			}
		}
	}

	int pathLinesV_i = pathLinesV.size();
	for (int i = 0; i<pathLinesV_i; ++i)
	{
		int pathLinesV_j = pathLinesV[i].size();
		for (int j = 0; j<pathLinesV_j; ++j)
		{
			if (abs(pathLinesV[i][j].startPoint.x - pathLinesV[i][j].endPoint.x) <= 0.5 && abs(pathLinesV[i][j].startPoint.y - pathLinesV[i][j].endPoint.y) <= 0.5)
			{
				//近似首尾相连
				pathLinesV[i][j].isCircle = true;
				pathLinesV[i][j].endPoint = pathLinesV[i][j].startPoint;

			}
		}
	}
}

/************************************************************************/
/* Funciton: doGridCalcOMP 【多核CPU并行加速算法】
 * Description: Marching Squares 算法的实现，首先计算单个格点上短等值线
 * Input:
	data: 天气数据值二维数组
	isovalues: 等值线值数组
	i: 格点所在cell中左上角点的数组x值下标
	j: 格点所在cell中左上角点的数组y值下标
	isovaluesNum: 等值数量
	dataSize_j: 纬度方向有多少个格点
	edgeArray: 存放所有生成的等值线临时边
 * Output: void
 * Author: gcdofree
 * Date: 2014.11.3
/************************************************************************/
static void doGridCalcOMP(vector<vector<float> > &data, vector<float> &isovalues, int i, int j, int isovaluesNum, int dataSize_j, isotools::Edge * &edgeArray)
{
	//预先算好邻域四个格点的最大最小值，避免重复计算squareIndex
	float maxValue = 0, minValue = FLT_MAX;
	if (data[i][j] > maxValue) maxValue = data[i][j];
	if (data[i][j + 1] > maxValue) maxValue = data[i][j + 1];
	if (data[i + 1][j + 1] > maxValue) maxValue = data[i + 1][j + 1];
	if (data[i + 1][j] > maxValue) maxValue = data[i + 1][j];
	if (data[i][j] < minValue) minValue = data[i][j];
	if (data[i][j + 1] < minValue) minValue = data[i][j + 1];
	if (data[i + 1][j + 1] < minValue) minValue = data[i + 1][j + 1];
	if (data[i + 1][j] < minValue) minValue = data[i + 1][j];

	for (int m = 0; m < isovaluesNum; ++m)
	{
		if (isovalues[m] < minValue || isovalues[m] > maxValue)
		{
			continue;
		}

		int squareIndex = 0;
		if (data[i][j] >= isovalues[m]) squareIndex |= 8;
		if (data[i][j + 1] >= isovalues[m]) squareIndex |= 4;
		if (data[i + 1][j + 1] >= isovalues[m]) squareIndex |= 2;
		if (data[i + 1][j] >= isovalues[m]) squareIndex |= 1;

		if (EdgeTable[squareIndex] != 0)
		{
			//消除二义性
			if (squareIndex == 5)// index = 0x5 
			{
				float centerValue = (1 / 4.0) *(data[i][j] + data[i][j + 1] + data[i + 1][j + 1] + data[i + 1][j]);
				if (centerValue < isovalues[m])
				{
					squareIndex = 10;
				}
			}
			else if (squareIndex == 10)//0x10
			{
				float centerValue = (1 / 4.0) *(data[i][j] + data[i][j + 1] + data[i + 1][j + 1] + data[i + 1][j]);
				if (centerValue < isovalues[m])
				{
					squareIndex = 5;
				}
			}
			//将生成的该线段添加到某一等值线上或自己建立某条等值线
			for (int k = 0; SegmentTable[squareIndex][k] != -1; k = k + 2)
			{
				int currentIndex = ((i*dataSize_j + j) * isovaluesNum + m) * 2 + k / 2;
				//isotools::Edge tempEdge;
				edgeArray[currentIndex].edgeIndex1 = SegmentTable[squareIndex][k];
				edgeArray[currentIndex].edgeIndex2 = SegmentTable[squareIndex][k + 1];
				edgeArray[currentIndex].i = i;
				edgeArray[currentIndex].j = j;
				edgeArray[currentIndex].m = m;//第m个等值线
				edgeArray[currentIndex].isovalue = isovalues[m];//等值线的值
			}
		}
	}
}

/************************************************************************/
/* Funciton: isMergeTwoIsoLine
 * Description: 进行线段合并操作
 * Input:
	i: 第i个等值
	j1: 表示第一条等值线的位置
	j2: 表示第二条等值线的位置
	size1: 表示第一条等值线上点的数量
	size2: 表示第二条等值线上点的数量
	pathLines: 等值线集合
 * Output: 若有线段合并，则返回true，否则返回false
 * Author: gcdofree
 * Date: 2014.11.3
/************************************************************************/
static bool isMergeTwoIsoLine(int i, int j1, int j2, int size1, int size2, vector< vector<isotools::Isoline> > &pathLines)
{
	if (size1 == 0 || size2 == 0)
	{
		return false;
	}

	if (pathLines[i][j1].isCircle || pathLines[i][j2].isCircle)
	{
		return false;
	}

	if (pathLines[i][j1].startPoint == pathLines[i][j2].startPoint && pathLines[i][j1].endPoint == pathLines[i][j2].endPoint)//头头，尾尾
	{
		if (size1 < size2)//把j1的点添加到j2
		{
			pathLines[i][j2].isCircle = true;
			pathLines[i][j2].endPoint = pathLines[i][j2].startPoint;
			for (int kk = 0; kk < size1; ++kk)
			{
				pathLines[i][j2].points.push_back(pathLines[i][j1].points.back());
				pathLines[i][j1].points.pop_back();
			}
			pathLines[i].erase(pathLines[i].begin() + j1);
		}
		else//把j2的点添加到j1
		{
			pathLines[i][j1].isCircle = true;
			pathLines[i][j1].endPoint = pathLines[i][j1].startPoint;
			for (int kk = 0; kk < size2; ++kk)
			{
				pathLines[i][j1].points.push_back(pathLines[i][j2].points.back());
				pathLines[i][j2].points.pop_back();
			}
			pathLines[i].erase(pathLines[i].begin() + j2);
		}
		return true;
	}
	else if (pathLines[i][j1].startPoint == pathLines[i][j2].endPoint && pathLines[i][j1].endPoint == pathLines[i][j2].startPoint)//头尾，尾头
	{
		if (size1 < size2)//把j1的点添加到j2
		{
			pathLines[i][j2].isCircle = true;
			pathLines[i][j2].endPoint = pathLines[i][j2].startPoint;
			for (int kk = 0; kk < size1; ++kk)
			{
				pathLines[i][j2].points.push_back(pathLines[i][j1].points.front());
				pathLines[i][j1].points.pop_front();
			}
			pathLines[i].erase(pathLines[i].begin() + j1);
		}
		else//把j2的点添加到j1
		{
			pathLines[i][j1].isCircle = true;
			pathLines[i][j1].endPoint = pathLines[i][j1].startPoint;
			for (int kk = 0; kk < size2; ++kk)
			{
				pathLines[i][j1].points.push_back(pathLines[i][j2].points.front());
				pathLines[i][j2].points.pop_front();
			}
			pathLines[i].erase(pathLines[i].begin() + j2);
		}
		return true;
	}
	else if (pathLines[i][j1].startPoint == pathLines[i][j2].startPoint)//头头
	{
		if (size1 < size2)//把j1的点添加到j2
		{
			pathLines[i][j2].startPoint = pathLines[i][j1].endPoint;
			for (int kk = 0; kk < size1; ++kk)
			{
				pathLines[i][j2].points.push_front(pathLines[i][j1].points.front());
				pathLines[i][j1].points.pop_front();
			}
			pathLines[i].erase(pathLines[i].begin() + j1);
		}
		else//把j2的点添加到j1
		{
			pathLines[i][j1].startPoint = pathLines[i][j2].endPoint;
			for (int kk = 0; kk < size2; ++kk)
			{
				pathLines[i][j1].points.push_front(pathLines[i][j2].points.front());
				pathLines[i][j2].points.pop_front();
			}
			pathLines[i].erase(pathLines[i].begin() + j2);
		}
		return true;
	}
	else if (pathLines[i][j1].startPoint == pathLines[i][j2].endPoint)//头尾
	{
		if (size1 < size2)//把j1的点添加到j2
		{
			pathLines[i][j2].endPoint = pathLines[i][j1].endPoint;
			for (int kk = 0; kk < size1; ++kk)
			{
				pathLines[i][j2].points.push_back(pathLines[i][j1].points.front());
				pathLines[i][j1].points.pop_front();
			}
			pathLines[i].erase(pathLines[i].begin() + j1);
		}
		else//把j2的点添加到j1
		{
			pathLines[i][j1].startPoint = pathLines[i][j1].startPoint;
			for (int kk = 0; kk < size2; ++kk)
			{
				pathLines[i][j1].points.push_front(pathLines[i][j2].points.back());
				pathLines[i][j2].points.pop_back();
			}
			pathLines[i].erase(pathLines[i].begin() + j2);
		}
		return true;
	}
	else if (pathLines[i][j1].endPoint == pathLines[i][j2].endPoint)//尾尾
	{
		if (size1 < size2)//把j1的点添加到j2
		{
			pathLines[i][j2].endPoint = pathLines[i][j1].startPoint;
			for (int kk = 0; kk < size1; ++kk)
			{
				pathLines[i][j2].points.push_back(pathLines[i][j1].points.back());
				pathLines[i][j1].points.pop_back();
			}
			pathLines[i].erase(pathLines[i].begin() + j1);
		}
		else//把j2的点添加到j1
		{
			pathLines[i][j1].endPoint = pathLines[i][j2].startPoint;
			for (int kk = 0; kk < size2; ++kk)
			{
				pathLines[i][j1].points.push_back(pathLines[i][j2].points.back());
				pathLines[i][j2].points.pop_back();
			}
			pathLines[i].erase(pathLines[i].begin() + j2);
		}
		return true;
	}
	else if (pathLines[i][j1].endPoint == pathLines[i][j2].startPoint)//尾头
	{
		if (size1 < size2)//把j1的点添加到j2
		{
			pathLines[i][j2].startPoint = pathLines[i][j1].startPoint;
			for (int kk = 0; kk < size1; ++kk)
			{
				pathLines[i][j2].points.push_front(pathLines[i][j1].points.back());
				pathLines[i][j1].points.pop_back();
			}
			pathLines[i].erase(pathLines[i].begin() + j1);
		}
		else//把j2的点添加到j1
		{
			pathLines[i][j1].endPoint = pathLines[i][j2].endPoint;
			for (int kk = 0; kk < size2; ++kk)
			{
				pathLines[i][j1].points.push_back(pathLines[i][j2].points.front());
				pathLines[i][j2].points.pop_front();
			}
			pathLines[i].erase(pathLines[i].begin() + j2);
		}
		return true;
	}
	return false;
}

/************************************************************************/
/* Funciton: isMergeTwoArea
 * Description: 进行区域合并
 * Input:
	mergePos: 表示两个区域相邻的边界区域
	pathLinesV: 第一个区域的等值线集合
	pathLinesV1: 第二个区域的等值线集合
	pathLinesTemp: 临时存放等值线的集合
 * Output: void
 * Author: gcdofree
 * Date: 2014.11.3
/************************************************************************/
static void isMergeTwoArea(int mergePos, vector< vector<isotools::Isoline>> &pathLinesV,
	vector< vector<isotools::Isoline>> &pathLinesV1, vector< vector<isotools::Isoline>> &pathLinesTemp)
{
	pathLinesTemp.clear();
	//对局部拼接结果进行合并
	//pathLinesV + pathLinesV1
	//pathLinesV查找边界片段，并加入到pathLinesTemp中
	int pathLinesSize1_i = pathLinesV.size();
	pathLinesTemp.resize(pathLinesSize1_i);
	for (int i = 0; i < pathLinesSize1_i; ++i)
	{
		int pathLinesSize1_j = pathLinesV[i].size();
		for (int j = 0; j < pathLinesSize1_j; j++)
		{
			if ((pathLinesV[i][j].startPoint.x > mergePos - 1 && pathLinesV[i][j].startPoint.x <= mergePos) ||
				pathLinesV[i][j].endPoint.x > mergePos - 1 && pathLinesV[i][j].endPoint.x <= mergePos)
			{
				if (!pathLinesV[i][j].isCircle)
				{
					pathLinesV[i][j].isBorder = true;
					pathLinesTemp[i].push_back(pathLinesV[i][j]);
					pathLinesV[i].erase(pathLinesV[i].begin() + j);
					--j;
					--pathLinesSize1_j;
				}
			}
		}
	}
	//pathLinesV1查找边界片段，并加入到pathLinesTemp中
	for (int i = 0; i < pathLinesSize1_i; ++i)
	{
		int pathLinesSize2_j = pathLinesV1[i].size();
		for (int j = 0; j < pathLinesSize2_j; j++)
		{
			if ((pathLinesV1[i][j].startPoint.x > mergePos - 1 && pathLinesV1[i][j].startPoint.x <= mergePos) ||
				pathLinesV1[i][j].endPoint.x > mergePos - 1 && pathLinesV1[i][j].endPoint.x <= mergePos)
			{
				if (!pathLinesV1[i][j].isCircle)
				{
					pathLinesV1[i][j].isBorder = true;
					pathLinesTemp[i].push_back(pathLinesV1[i][j]);
					pathLinesV1[i].erase(pathLinesV1[i].begin() + j);
					--j;
					--pathLinesSize2_j;
				}
			}
		}
	}

	//对pathLinesTemp中的片段进行拼接
	for (int i = 0; i < pathLinesSize1_i; ++i)
	{
		int pathLinesSizej = pathLinesTemp[i].size();
		for (int j = 0; j < pathLinesSizej - 1; ++j)
		{
			for (int k = j + 1; k < pathLinesSizej; ++k)
			{
				bool isMerge = isMergeTwoIsoLine(i, j, k, pathLinesTemp[i][j].points.size(), pathLinesTemp[i][k].points.size(), pathLinesTemp);
				if (isMerge)
				{
					--j;
					--pathLinesSizej;
					break;
				}
			}
		}
		//将拼接结果加入到pathLinesV
		for (int j = 0; j < pathLinesSizej; ++j)
		{
			pathLinesV[i].push_back(pathLinesTemp[i][j]);
		}
		//将pathLinesV1中的等值线加入到pathLinesV
		pathLinesSizej = pathLinesV1[i].size();
		for (int j = 0; j < pathLinesSizej; ++j)
		{
			pathLinesV[i].push_back(pathLinesV1[i][j]);
		}
	}
}

/************************************************************************/
/* Funciton: doMarchingSquaresAccelerateOMP 【此方法适用于 多核CPU】
 * Description: Marching Squares 算法的实现，首先计算生成所有的格点上短等值线，之后再进行合并拼接
 * Input:
	data: 天气数据值二维数组
	isovalues: 等值线值数组
	pathLinesV: 返回该等值数组下的所有各条等值线数组的集合
	startLongitude: 起始经度（起始x坐标）
	longitudeGridSpace: 经度间隔（x坐标间隔）
	startLatitude: 起始纬度（起始y坐标）
	latitudeGridSpace: 纬度间隔（y坐标间隔）
	maxGridValue: 网格点中的最大值
	minGridValue: 网格点中的最小值
* Output: void
* Author: gcdofree
* Date: 2014.11.3
/************************************************************************/
static void doMarchingSquaresAccelerateOMP(vector<vector<float> > &data, vector<float> &isovalues, vector< vector<isotools::Isoline>> &pathLinesV,
	float startLongitude, float longitudeGridSpace, float startLatitude, float latitudeGridSpace, float maxGridValue, float minGridValue)
{
	int edgeNum = 0;

	//利用OpenMP加速。在拼接阶段，分四个section同步进行
	vector< vector<isotools::Isoline>> pathLinesV1, pathLinesV2, pathLinesV3, pathLinesTemp;
	pathLinesV.clear();
	pathLinesV1.clear();
	pathLinesV2.clear();
	pathLinesV3.clear();
	pathLinesTemp.clear();

	isotools::Edge *edgeArray;//用于存放生成的等值线片段
	vector<vector<int>> pathLinesHit, pathLinesHit1, pathLinesHit2, pathLinesHit3;

	int isovaluesNum = isovalues.size();
	for (int m = 0; m < isovaluesNum; ++m)
	{
		if (isovalues[m] < minGridValue || isovalues[m] > maxGridValue)
		{
			isovalues.erase(isovalues.begin() + m);//对vector进行增删元素后，begin，end，size会失效，必须即时调用
			m--;
			isovaluesNum--;
		}
	}
	
	pathLinesV.resize(isovaluesNum);
	pathLinesV1.resize(isovaluesNum);
	pathLinesV2.resize(isovaluesNum);
	pathLinesV3.resize(isovaluesNum);
	pathLinesTemp.resize(isovaluesNum);
	pathLinesHit.resize(isovaluesNum);
	pathLinesHit1.resize(isovaluesNum);
	pathLinesHit2.resize(isovaluesNum);
	pathLinesHit3.resize(isovaluesNum);

	//首先生成所有格点上短的等值线

	int dataSize_i = data.size() - 1;
	int dataSize_j = data[0].size() - 1;

	int edgeSize = dataSize_i*dataSize_j*isovaluesNum * 2;

	edgeArray = new isotools::Edge[edgeSize];

#pragma omp parallel for	
	for (int i = 0; i < edgeSize; ++i)
	{
		edgeArray[i].isovalue = 999;//初始化为999
	}

#pragma omp parallel for	
	for (int i = 0; i<dataSize_i; ++i)//逐行扫
	{
		for (int j = 0; j<dataSize_j; ++j)//逐列扫
		{
			doGridCalcOMP(data, isovalues, i, j, isovaluesNum, dataSize_j, edgeArray);
		}
	}

	int edgeSize1 = edgeSize / 4;
	int edgeSize2 = edgeSize / 2;
	int edgeSize3 = edgeSize / 4 * 3;
	//开始进行局部拼接，分为4个section，使用OpenMP同步进行
#pragma omp parallel sections
	{
#pragma omp section
		{
			for (int i = 0; i < edgeSize1; ++i)
			{
				if (edgeArray[i].isovalue != 999)
				{
					addPointToLineAccelerate(edgeArray[i].edgeIndex1, edgeArray[i].edgeIndex2, edgeArray[i].i, edgeArray[i].j,
						data, pathLinesHit[edgeArray[i].m], edgeArray[i].isovalue, pathLinesV[edgeArray[i].m],
						startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
				}
			}
		}

#pragma omp section
		{
		for (int i = edgeSize1; i < edgeSize2; ++i)
			{
				if (edgeArray[i].isovalue != 999)
				{
					addPointToLineAccelerate(edgeArray[i].edgeIndex1, edgeArray[i].edgeIndex2, edgeArray[i].i, edgeArray[i].j,
						data, pathLinesHit1[edgeArray[i].m], edgeArray[i].isovalue, pathLinesV1[edgeArray[i].m],
						startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
				}
			}
		}
#pragma omp section
		{
			for (int i = edgeSize2; i < edgeSize3; ++i)
			{
				if (edgeArray[i].isovalue != 999)
				{
					addPointToLineAccelerate(edgeArray[i].edgeIndex1, edgeArray[i].edgeIndex2, edgeArray[i].i, edgeArray[i].j,
						data, pathLinesHit2[edgeArray[i].m], edgeArray[i].isovalue, pathLinesV2[edgeArray[i].m],
						startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
				}
			}
		}
#pragma omp section
		{
			for (int i = edgeSize3; i < edgeSize; ++i)
			{
				if (edgeArray[i].isovalue != 999)
				{
					addPointToLineAccelerate(edgeArray[i].edgeIndex1, edgeArray[i].edgeIndex2, edgeArray[i].i, edgeArray[i].j,
						data, pathLinesHit3[edgeArray[i].m], edgeArray[i].isovalue, pathLinesV3[edgeArray[i].m],
						startLongitude, longitudeGridSpace, startLatitude, latitudeGridSpace);
				}
			}
		}
	}

	//区域拼接
	isMergeTwoArea(dataSize_i / 4, pathLinesV, pathLinesV1, pathLinesTemp);
	isMergeTwoArea(dataSize_i / 4 * 3, pathLinesV2, pathLinesV3, pathLinesTemp);
	isMergeTwoArea(dataSize_i / 2, pathLinesV, pathLinesV2, pathLinesTemp);

	int pathLinesV_i = pathLinesV.size();
	for (int i = 0; i < pathLinesV_i; ++i)
	{
		int pathLinesV_j = pathLinesV[i].size();
		for (int j = 0; j < pathLinesV_j; ++j)
		{
			if (abs(pathLinesV[i][j].startPoint.x - pathLinesV[i][j].endPoint.x) <= 0.5 && abs(pathLinesV[i][j].startPoint.y - pathLinesV[i][j].endPoint.y) <= 0.5)
			{
				//近似首尾相连
				pathLinesV[i][j].isCircle = true;
				pathLinesV[i][j].endPoint = pathLinesV[i][j].startPoint;
			}
		}
	}

	delete[] edgeArray;
}

}
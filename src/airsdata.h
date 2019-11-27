/*!
 * @file airsdata.h 声明AIRS工作路程中使用的数据结构
 * @version 0.1
 * @date 2019-10-16
 * @note
 */

#ifndef _AIRS_DATA_H_
#define _AIRS_DATA_H_

#include <string>
#include <deque>
#include <vector>
#include <boost/smart_ptr.hpp>
#include <string.h>
#include <longnam.h>
#include <fitsio.h>
#include "ADefine.h"

using std::string;
using namespace AstroUtil;

#define FAIL_IMGREDUCT		-1		//< 图像处理失败
#define FAIL_ASTROMETRY		-2		//< 天文定位失败
#define FAIL_PHOTOMETRY		-3		//< 流量定标失败
#define SUCCESS_INIT		0x00	//< 初始化
#define SUCCESS_COMPLETE	0x01	//< 完成处理流程
#define SUCCESS_IMGREDUCT	0x02	//< 完成图像处理
#define SUCCESS_ASTROMETRY	0x04	//< 完成天文定位
#define SUCCESS_PHOTOMETRY	0x08	//< 完成流量定标
#define PROCESS_IMGREDUCT	0x11	//< 执行图像处理
#define PROCESS_ASTROMETRY	0x12	//< 执行天文定位
#define PROCESS_PHOTOMETRY	0x13	//< 执行流量定标

/*!
 * @struct ObjectInfo 单个天体的特征信息
 */
struct ObjectInfo {
	/* 图像处理结果 */
	int id;		//< 在图像中的编号
	double x;	//< 图像坐标
	double y;
	double flux;	//< 积分流量
	double mag_img;	//< 仪器星等
	double magerr_img;	//< 仪器星等误差
	double fwhm;	//< FWHM
	double ellip;	//< 椭率. 0: 圆; 1: 线
	/* 天文定位 */
	double ra_fit;	//< 赤道坐标, J2000, 量纲: 角度
	double dec_fit;
	/* 流量定标 */
	/*!
	 * @var matched 恒星匹配结果
	 * 0 : 未匹配
	 * 1 : 与星表匹配成功
	 * 2 : 坏像素/热点
	 * 3 : 前后关联匹配成功
	 */
	int matched;
	double ra_cat;	//< 星表坐标+自行改正, J2000, 量纲: 角度
	double dec_cat;
	double mag_cat;	//< V星等: 星表. 无效值: 20
	double mag_fit;	//< V星等: 拟合

public:
	ObjectInfo() {
		memset(this, 0, sizeof(ObjectInfo));
		mag_img = 20;
		mag_cat = 20;
		mag_cat = 20;
	}
};
typedef boost::shared_ptr<ObjectInfo> NFObjPtr;
typedef std::vector<NFObjPtr> NFObjVec;

struct wcsinfo {
	int x1, y1, x2, y2;	//< 在全图中的XY左上角与右下角坐标
	int sizew, sizeh;	//< 参与生成WCS信息的区域大小...未来Lengdre多项式模型需要
	double x0, y0;	//< XY参考点
	double r0, d0;	//< RA/DEC参考点, 量纲: 弧度
	double cd[2][2];	//< 转换矩阵
	int orderA, orderB;	//< SIP改正阶数
	int ncoefA, ncoefB;	//< SIP系数数量
	double *A, *B;	//< 线性改正系数

public:
	wcsinfo() {
		x1 = y1 = 0;
		x2 = y2 = 0;
		sizew = sizeh = 0;
		x0 = y0 = 0.0;
		r0 = d0 = 0.0;
		orderA = orderB = 0;
		ncoefA = ncoefB = 0;
		A = B = NULL;
	}

	virtual ~wcsinfo() {
		if (A) {
			delete[] A;
			A = NULL;
		}
		if (B) {
			delete[] B;
			B = NULL;
		}
	}

	wcsinfo &operator=(const wcsinfo &other) {
		if (this != &other) {
			x1   = other.x1;
			y1   = other.y1;
			x2   = other.x2;
			y2   = other.y2;
			x0   = other.x0;
			y0   = other.y0;
			r0   = other.r0;
			d0   = other.d0;
			orderA = other.orderA;
			orderB = other.orderB;
			ncoefA = other.ncoefA;
			ncoefB = other.ncoefB;
			memcpy(&cd[0][0], &other.cd[0][0], sizeof(double) * 4);
			if (ncoefA) {
				A = new double[ncoefA];
				memcpy(A, other.A, sizeof(double) * ncoefA);
			}
			if (ncoefB) {
				B = new double[ncoefB];
				memcpy(B, other.B, sizeof(double) * ncoefB);
			}
		}

		return *this;
	}

protected:
	/*!
	 * @brief 计算SIP改正模型中与阶数对应的系数数量
	 * @return
	 * 系数数量
	 */
	int term_count(int order) {
		return (order + 1) * (order + 2) / 2;
	}

	/*!
	 * @brief 为SIP改正系数分配存储空间
	 * @param n   系数数量
	 * @param ptr 系数存储地址
	 */
	void alloc_coef(int n, double **ptr) {
		if ((ptr == &A && n != ncoefA) || (ptr == &B && n != ncoefA)) {
			if (ptr == &A)
				ncoefA = n;
			else
				ncoefB = n;
			if (*ptr) {
				delete[] (*ptr);
				(*ptr) = NULL;
			}
		}
		if (*ptr == NULL)
			(*ptr) = new double[n];
	}

	/*!
	 * @brief 参考平面坐标转换为世界坐标
	 * @param xi    参考平面X坐标, 量纲: 弧度
	 * @param eta   参考平面Y坐标, 量纲: 弧度
	 * @param newr  参考点(x0, y0)对应的新的赤经, 量纲: 弧度
	 * @param newd  参考点(x0, y0)对应的新的赤纬, 量纲: 弧度
	 * @param ra    世界坐标赤经, 量纲: 弧度
	 * @param dec   世界坐标赤纬, 量纲: 弧度
	 */
	void plane_to_wcs(double xi, double eta, double newr, double newd, double &ra, double &dec) {
		double fract = cos(newd) - eta * sin(newd);
		ra = cyclemod(newr + atan2(xi, fract), A2PI);
		dec = atan2(((eta * cos(newd) + sin(newd)) * cos(ra - newr)), fract);
	}

	double poly_val(double x, double y, double *coef, int order) {
		int i, j, k, m;
		double val(0.0), t, px(1.0), py;

		for (i = 0, k = 0; i <= order; ++i) {
			for (j = 0, py = 1.0, t = 0.0, m = order - i; j <= m; ++j, ++k) {
				t += coef[k] * py;
				py *= y;
			}

			val += t * px;
			px *= x;
		}

		return val;
	}

	void project_correct(double &x, double &y) {
		double dx(0.0), dy(0.0);
		dx = poly_val(x, y, A, orderA);
		dy = poly_val(x, y, B, orderB);
		x += dx;
		y += dy;
	}

public:
	bool load_wcs(const string &filepath) {
		fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
		char key[10];
		int status(0), ncoef, i, j, k, m;

		fits_open_file(&fitsptr, filepath.c_str(), 0, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRPIX1", &x0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRPIX2", &y0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRVAL1", &r0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRVAL2", &d0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD1_1", &cd[0][0], NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD1_2", &cd[0][1], NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD2_1", &cd[1][0], NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD2_2", &cd[1][1], NULL, &status);

		fits_read_key(fitsptr, TINT, "A_ORDER", &orderA, NULL, &status);
		if (status)
			return false;
		ncoef = term_count(orderA);
		alloc_coef(ncoef, &A);
		for (i = 0, k = 0; i <= orderA; ++i) {
			for (j = 0, m = orderA - i; j <= m; ++j, ++k) {
				sprintf(key, "A_%d_%d", i, j);
				fits_read_key(fitsptr, TDOUBLE, key, A + k, NULL, &status);
			}
		}

		fits_read_key(fitsptr, TINT, "B_ORDER", &orderB, NULL, &status);
		if (status)
			return false;
		ncoef = term_count(orderB);
		alloc_coef(ncoef, &B);
		for (i = 0, k = 0; i <= orderB; ++i) {
			for (j = 0, m = orderB - i; j <= m; ++j, ++k) {
				sprintf(key, "B_%d_%d", i, j);
				fits_read_key(fitsptr, TDOUBLE, key, B + k, NULL, &status);
			}
		}

		fits_close_file(fitsptr, &status);

		r0 *= D2R;
		d0 *= D2R;
		return !status;
	}

	/*!
	 * @brief 图像坐标转换为世界坐标
	 * @param x     图像X坐标, 量纲: 像素
	 * @param y     图像Y坐标, 量纲: 像素
	 * @param newr  参考点(x0, y0)对应的新的赤经, 量纲: 弧度
	 * @param newd  参考点(x0, y0)对应的新的赤纬, 量纲: 弧度
	 * @param ra    世界坐标赤经, 量纲: 弧度
	 * @param dec   世界坐标赤纬, 量纲: 弧度
	 */
	void image_to_wcs(double x, double y, double newr, double newd, double &ra, double &dec) {
		double xi, eta;
		x -= (x0 + x1);
		y -= (y0 + y1);
		project_correct(x, y);
		xi = (cd[0][0] * x + cd[0][1] * y) * D2R;
		eta = (cd[1][0] * x + cd[1][1] * y) * D2R;
		plane_to_wcs(xi, eta, newr, newd, ra, dec);
		ra *= R2D;
		dec *= R2D;
	}

	void image_to_wcs(double x, double y, double &ra, double &dec) {
		image_to_wcs(x, y, r0, d0, ra, dec);
	}
};

/*!
 * @struct OneFrame 单帧图像的特征信息
 */
struct OneFrame {
	int result;			//< 处理结果标志字
	/* FITS文件 */
	string filepath;	//< 文件路径
	string filename;	//< 文件名
	string tmobs;		//< 曝光起始时间, UTC. CCYY-MM-DDThh:mm:ss.ssssss
	string tmmid;		//< 曝光中间时间, UTC. CCYY-MM-DDThh:mm:ss.ssssss
	string imgtype;		//< 图像类型
	int wimg;			//< 图像宽度
	int himg;			//< 图像高度
	int fno;			//< 帧编号
	double expdur;		//< 曝光时间, 量纲: 秒
	double secofday;	//< 当日秒数
	double mjd;			//< 修正儒略日: 曝光中间时刻
	/* 网络标志 */
	string gid;		//< 组标志
	string uid;		//< 单元ID
	string cid;		//< 相机ID
	/* 处理结果 */
	double fwhm;		//< 中心区域统计FWHM, 量纲: 像素
	double rac, decc;	//< 中心视场指向, 量纲: 角度
	double azic, altc;	//< 中心视场指向, 量纲: 角度
	double airmass;		//< 大气质量: 中心指向
	int lastid;			//< 感兴趣目标的最后一个编号
	NFObjVec nfobjs;	//< 集合: 目标特征
	/*
	 * 仪器星等改正及大气消光
	 * @ -mag0参与拟合消光系数, 每天依据大气质量分布范围判定是否输出晓光系数
	 */
	double mag0;	//< 星等零点. 仪器星等0时对应的视星等
	double magk;	//< 拟合系数. magV=kmag*magInst+mag0

public:
	OneFrame() {
		result = SUCCESS_INIT;
		wimg = himg = 0;
		fno  = 0;
		secofday = 0;
		mjd  = 0;
		expdur = 0.0;
		fwhm = 0.0;
		rac  = decc = 0.0;
		azic = altc = 0.0;
		airmass  = 0.0;
		lastid   = 0;
		mag0     = 0.0;
		magk     = 0.0;
	}

	virtual ~OneFrame() {
		nfobjs.clear();
	}
};
typedef boost::shared_ptr<OneFrame> FramePtr;	//< 单帧图像特征信息访问指针
typedef std::deque<FramePtr> FrameQueue;		//< 图像特征存储队列

#endif

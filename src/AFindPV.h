/**
 * @class AFindPV 检测关联识别相邻图像中的运动目标
 * @version 0.1
 * @date 2019-11-12
 */

#ifndef AFINDPV_H_
#define AFINDPV_H_

#include <boost/smart_ptr.hpp>
#include <vector>
#include <string>
#include <cmath>
#include "airsdata.h"

namespace AstroUtil {
//////////////////////////////////////////////////////////////////////////////
typedef struct pv_point {// 单数据点
	// 图像帧
	string filename;	//< 文件名
	string tmmid;		//< 曝光中间时刻, 格式: CCYY-MM-DDThh:mm:ss<.sss<sss>>
	int fno;			//< 帧编号
	double secofday;	//< 曝光中间时间对应的当日秒数
	// 目标
	int id;				//< 在原始图像中的唯一编号
	int related;		//< 被关联次数
	int matched;		//< 恒星标志
	double x, y;		//< 星象质心在模板中的位置
	double ra, dc;		//< 赤道坐标, 量纲: 角度. 坐标系: J2000
	double mag;			//< 星等
	double magerr;		//< 星等误差
	double snr;			//< 积分信噪比

public:
	int inc_rel() {	// 增加一次关联次数
		return ++related;
	}

	int dec_rel() {	// 减少一次关联次数
		return --related;
	}
} PvPt;
typedef boost::shared_ptr<PvPt> PvPtPtr;
typedef std::vector<PvPtPtr> PvPtVec;

typedef struct pv_frame {// 单帧数据共性属性及数据点集合
	PvPtVec pts;	//< 数据点集合

public:
	virtual ~pv_frame() {
		pts.clear();
	}
} PvFrame;
typedef boost::shared_ptr<PvFrame> PvFrmPtr;

/*
 * pv_candidate使用流程:
 * 1. 构建对象
 * 2. xy_expect(): 评估输出的xy与数据点之间的偏差是否符合阈值
 * 3. add_point(): 将数据点加入候选体
 * 4. recheck_frame(): 在EndFrame()中评估当前帧数据点是否为候选体提供有效数据
 */
typedef struct pv_candidate {// 候选体
	PvPtVec pts;	//< 已确定数据点集合
	PvPtVec frmu;	//< 由当前帧加入的不确定数据点
	int mode;		//< 运动规律. 0: 初始; 1: 凝视(星象不动); 2: 穿越
	double vx, vy;	//< XY变化速度

public:
	PvPtPtr first_point() {// 构成候选体的第一个数据点
		return pts[0];
	}

	PvPtPtr last_point() {// 构成候选体的最后一个数据点
		return pts[pts.size() - 1];
	}

	void xy_expect(double secs, double &x, double &y) {	// 由候选体已知(加)速度计算其预测位置
		PvPtPtr pt = last_point();
		double t = secs - pt->secofday;
		x = pt->x + vx * t;
		y = pt->y + vy * t;
	}

	int track_mode(PvPtPtr pt1, PvPtPtr pt2) {
		double dt = pt2->secofday - pt1->secofday;
		double dr = pt2->ra - pt1->ra;
		double dd = pt2->dc - pt1->dc;
		double dx = pt2->x - pt1->x;
		double dy = pt2->y - pt1->y;
		double limit = 10.0 * AS2D * dt;

		if (dr < -180.0) dr += 360.0;
		else if (dr > 180.0) dr -= 360.0;
		if (fabs(dr) < limit && fabs(dd) < limit) return 0;
		if (fabs(dx) <= 2.0 && fabs(dy) <= 2.0) return 1;   // 2.0: "跟踪"误差
		return 2;
	}

	/*!
	 * @brief 将一个数据点加入候选体
	 * @return
	 * -1 : 模式不匹配
	 * -2 : 模式匹配但位置不匹配
	 *  0 : 模式错误
	 *  1 : 目标位置凝视模式
	 *  2 : 目标位置变化模式
	 */
	int add_point(PvPtPtr pt) {
		if (pts.size() >= 2) {
			PvPtPtr last = last_point();
			int mode_new = track_mode(last, pt);
			if (mode_new != mode) return -1;
			if (mode == 2) {// 计算期望位置与输入位置差异
				double x, y;
				xy_expect(pt->secofday, x, y);
				if (fabs(x - pt->x) > 2.0 || fabs(y - pt->y) > 2.0) return -2;
			}
			if (pts.size() == 2) last->inc_rel(); // 为初始构建的候选体的末尾添加关联标志
			pt->inc_rel();
			frmu.push_back(pt);
		}
		else {
			pts.push_back(pt);
			if (pts.size() == 2) {// 创建候选体
				PvPtPtr prev = first_point();
				if ((mode = track_mode(prev, pt)) == 2) {
					vx = (pt->x - prev->x) / (pt->secofday - prev->secofday);
					vy = (pt->y - prev->y) / (pt->secofday - prev->secofday);
				}
			}
		}
		return mode;
	}

	PvPtPtr update() {	// 检查/确认来自当前帧的数据点是否加入候选体已确定数据区
		PvPtPtr pt;
		int n(frmu.size());
		if (!n) return pt;
		if (n == 1) pt = frmu[0];
		else {
			double x, y; // 期望位置
			double dx, dy, dx2y2, dx2y2min(1E30);

			if (mode == 1) {
				x = last_point()->x;
				y = last_point()->y;
			}
			else {// mode == 2
				xy_expect(frmu[0]->secofday, x, y);
			}
			// 查找与期望位置最接近的点
			for (int i = 0; i < n; ++i) {
				if (!frmu[i]->matched) {
					dx = frmu[i]->x - x;
					dy = frmu[i]->y - y;
					dx2y2 = dx * dx + dy * dy;
					if (dx2y2 < dx2y2min) {
						dx2y2min = dx2y2;
						if (pt.use_count()) pt->dec_rel();
						pt = frmu[i];
					}
				}
			}
			if (pt.use_count() && mode == 2) {
				PvPtPtr last = last_point();
				double dt = pt->secofday - last->secofday;
				vx = (pt->x - last->x) / dt;
				vy = (pt->y - last->y) / dt;
			}
		}

		if (pt.use_count()) pts.push_back(pt);
		frmu.clear();
		return pt;
	}

	void complete() {
		PvPtPtr first = first_point();
		PvPtPtr last  = last_point();
		int mode_new = track_mode(first, last);
		if (mode_new == 0) pts.clear();
	}

	virtual ~pv_candidate() {
		pts.clear();
	}
} PvCan;
typedef boost::shared_ptr<PvCan> PvCanPtr;
typedef std::vector<PvCanPtr> PvCanVec;

typedef struct pv_object {	// PV目标
	PvPtVec pts;	//< 已确定数据点集合
} PvObj;
typedef boost::shared_ptr<PvObj> PvObjPtr;
typedef std::vector<PvObjPtr> PvObjVec;

class AFindPV {
public:
	AFindPV();
	virtual ~AFindPV();

protected:
	/* 成员变量 */
	int last_fno_;	//< 最后一个帧编号
	PvFrmPtr frmprev_;	//< 帧数据
	PvFrmPtr frmnow_;	//< 当前帧数据
	PvCanVec cans_;		//< 候选体集合
	PvObjVec objs_;		//< 目标集合

protected:
	/*!
	 * @brief 开始处理新的一帧数据
	 */
	void new_frame(int fno);
	/*!
	 * @brief 完成一帧数据处理
	 */
	void end_frame();
	/*!
	 * @brief 建立新的候选体
	 */
	void create_candidates();
	/*!
	 * @brief 尝试将当前帧数据加入候选体
	 */
	void append_candidates();
	/*!
	 * @brief 检查候选体, 确认其有效性
	 * @note
	 * 判据: 候选体时标与当前帧时标之差是否大于阈值
	 */
	void recheck_candidates();
	/*!
	 * @brief 处理所有候选体
	 */
	void complete_candidates();
	/*!
	 * @brief 将一个候选体转换为目标
	 */
	void candidate2object(PvCanPtr can);

public:
	/*!
	 * @brief 开始处理一段连续数据
	 */
	void NewSequence();
	/*!
	 * @brief 完成一段连续数据处理
	 */
	void EndSequence();
	/*!
	 * @brief 添加一个候选体
	 */
	void AddPoint(PvPtPtr pt);
	/*!
	 * @brief 查看被识别的目标
	 */
	PvObjVec& GetObject();
};
//////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */

#endif /* AFINDPV_H_ */

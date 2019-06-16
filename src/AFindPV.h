/*
 * @file AFindPV.h 查找位置变化的瞬变源
 */

#ifndef AFINDPV_H_
#define AFINDPV_H_

#include <string>
#include <string.h>
#include <boost/smart_ptr.hpp>
#include <boost/container/stable_vector.hpp>
#include <boost/container/deque.hpp>

using std::string;

typedef boost::shared_array<char> chararr;

namespace AstroUtil {
///////////////////////////////////////////////////////////////////////////////
struct param_pv { // 位置变源关联识别参数
	int nptmin;	//< 构成PV的最小数据点数量
	double dtmax;	//< 相邻关联数据点的最大时间间隔, 量纲: 天
	double stepmin;	//< 最小步长
	double stepmax;	//< 最大步长
	double dxymax;	//< XY坐标偏差的最大值, 量纲: 像素

public:
	param_pv() {
		nptmin = 5;
		dtmax = 2.0 / 86400.0;
		stepmin = 0.0;
		stepmax = 50.0;
		dxymax = 2.0;
	}
};

typedef struct pv_point { // 单数据点
	int related;	//< 被关联次数
	int fno;		//< 帧编号
	double mjd;		//< 曝光中间时间对应的修正儒略日
	double x, y;	//< 星象质心在模板中的位置
	double ra, dc;	//< 赤道坐标, 量纲: 角度. 坐标系: J2000
	double mag;		//< 星等

public:
	pv_point() {
		memset(this, 0, sizeof(pv_point));
	}

	int inc_rel() {		// 增加一次关联次数
		return ++related;
	}

	int dec_rel() {		// 减少一次关联次数
		return --related;
	}
} PVPT;
typedef boost::shared_ptr<PVPT> PPVPT;
typedef boost::container::stable_vector<PPVPT> PPVPTVEC;

typedef struct pv_frame {		// 单帧数据共性属性及数据点集合
	double mjd;		//< 曝光中间时间对应的修正儒略日
	PPVPTVEC pts;	//< 数据点集合

public:
	pv_frame() {
		mjd = 0.0;
	}

	pv_frame(double Mjd) {
		mjd = Mjd;
	}

	virtual ~pv_frame() {
		pts.clear();
	}
} PVFRM;
typedef boost::shared_ptr<PVFRM> PPVFRM;
typedef boost::container::deque<PPVFRM> PPVFRMDQ;

/*
 * pv_candidate使用流程:
 * 1. 构建对象
 * 2. xy_expect(): 评估输出的xy与数据点之间的偏差是否符合阈值
 * 3. add_point(): 将数据点加入候选体
 * 4. recheck_frame(): 在EndFrame()中评估当前帧数据点是否为候选体提供有效数据
 */
typedef struct pv_candidate {	// 候选体
	PPVPTVEC pts;	//< 已确定数据点集合
	PPVPTVEC frmu;	//< 由当前帧加入的不确定数据点
	double vx, vy;	//< XY变化速度
	double lastmjd;	//< 加入候选体的最后一个数据点对应的时间, 量纲: 天; 涵义: 修正儒略日
	int fx, fy;		//< 运动方向

public:
	PPVPT last_point() {	// 构成候选体的最后一个数据点
		return pts[pts.size() - 1];
	}

	void xy_expect(double mjd, double &x, double &y) {	// 由候选体已知(加)速度计算其预测位置
		PPVPT pt = last_point();
		double t = mjd - pt->mjd;
		x = pt->x + vx * t;
		y = pt->y + vy * t;
	}

	/*!
	 * @brief 将一个数据点加入候选体
	 */
	void add_point(PPVPT pt) {
		if (pts.size() >= 2) {	// 构成候选体的候选点
			pt->inc_rel();
			frmu.push_back(pt);
		} else {	// 构成候选体的初始2点
			pts.push_back(pt);
			if (pts.size() == 2) {
				PPVPT prev = pts[0];
				double t = (lastmjd = pt->mjd) - prev->mjd;
				vx = (pt->x - prev->x) / t;
				vy = (pt->y - prev->y) / t;
				// 评估运动方向
				fx = fabs(pt->x - prev->x) < 1.0 ?
						0 : (pt->x > prev->x ? 1 : -1);
				fy = fabs(pt->y - prev->y) < 1.0 ?
						0 : (pt->y > prev->y ? 1 : -1);
			}
		}
	}

	PPVPT update(int mode) {	// 检查/确认来自当前帧的数据点是否加入候选体已确定数据区
		PPVPT pt;
		double x, y;	// 期望位置
		double dx, dy, dx2y2, dx2y2max(1E30);
		if (!frmu.size())
			return pt;
		lastmjd = frmu[0]->mjd;
		/*
		 * 该算法解决: 多点加入一个候选体时带来的混淆
		 */
		if (mode == 0)
			xy_expect(lastmjd, x, y);
		else {
			x = last_point()->x;
			y = last_point()->y;
		}
		for (PPVPTVEC::iterator it = frmu.begin(); it != frmu.end(); ++it) {// 查找与候选体末端最接近的数据
			dx = (*it)->x - x;
			dy = (*it)->y - y;
			dx2y2 = dx * dx + dy * dy;
			if (dx2y2 < dx2y2max) {
				if (pt.use_count())
					pt->dec_rel();

				dx2y2max = dx2y2;
				pt = *it;
			}
		}

		if (mode == 0) {
			PPVPT last = last_point();
			double t = lastmjd - last->mjd;
			vx = (pt->x - last->x) / t;
			vy = (pt->y - last->y) / t;
		}
		pts.push_back(pt);
		frmu.clear();
		return pt;
	}

	virtual ~pv_candidate() {
		pts.clear();
	}
} PVCAN;
typedef boost::shared_ptr<PVCAN> PPVCAN;
typedef boost::container::stable_vector<PPVCAN> PPVCANVEC;

typedef struct pv_object {	// PV目标
	PPVPTVEC pts;	//< 已确定数据点集合
} PVOBJ;
typedef boost::shared_ptr<PVOBJ> PPVOBJ;
typedef boost::container::stable_vector<PPVOBJ> PPVOBJVEC;

class AFindPV {
public:
	AFindPV();
	virtual ~AFindPV();

protected:
	param_pv param_;	//< 数据处理参数
	int wmap_, hmap_;	//< 图像分辨率
	chararr lastmap_;	//< 上一帧的恒星标志位图
	chararr newmap_;	//< 当前图的恒星标志位图
	int rcross_;		//< 恒星交叉半径
	int track_mode_;	//< 观测模式. -1: Unknown; 0: (近似)跟踪恒星; 1: 跟踪目标
	int fno_;			//< 最新数据帧编号
	double frm_interval_;	//< 帧间隔时间
	PPVFRM frmprev_;	//< 前一数据帧
	PPVFRM frmlast_;	//< 最新数据帧
	PPVCANVEC cans_;	//< 候选体集合
	PPVOBJVEC objs_;	//< 目标集合

protected:
	/*!
	 * @brief 检查当前帧中数据点位置是否不变
	 */
	bool is_freeze(PPVPT pt);
	/*!
	 * @brief 检查前一帧中数据点位置是否不变
	 */
	bool prev_is_freeze(PPVPT pt);

public:
	/*!
	 * @brief 依据帧间隔设置候选体点数据之间最长时间间隔
	 */
	void UpdateFrameDelay(double dt);
	/*!
	 * @brief 设置图像分辨率
	 */
	void SetDimension(int wimg, int himg);
	/*!
	 * @brief 设置跟踪模式
	 * @param mode 跟踪模式. -1: Unknown; 0: Sidereal/Freeze; 1: Object
	 */
	void SetTrackMode(int mode);
	/*!
	 * @brief 开始新的查找序列
	 */
	void NewSequence();
	void EndSequence();
	void AddPoint(PPVPT pt);
	/*!
	 * @brief 查看候选体
	 */
	PPVCANVEC& GetCandidate();
	/*!
	 * @brief 查看被识别的目标数量
	 */
	int GetNumber();
	/*!
	 * @brief 查看被识别的目标
	 */
	PPVOBJVEC& GetObject();

protected:
	/*!
	 * @brief 以x/y为中心, 渲染星象位置
	 * @param x0 X坐标
	 * @param y0 Y坐标
	 */
	void render_star(double x0, double y0);
	void new_frame(double mjd);
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
	void candidate2object(PPVCAN can);
};
///////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */

#endif /* AFINDPV_H_ */

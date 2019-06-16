/*
 * @file AFindPV.cpp 查找位置变化的瞬变源
 */
#include <boost/make_shared.hpp>
#include "ADefine.h"
#include "AFindPV.h"

namespace AstroUtil {
///////////////////////////////////////////////////////////////////////////////
AFindPV::AFindPV() {
	wmap_ = hmap_ = 0;
	fno_ = -1;
	rcross_ = 2;
	track_mode_ = 0;
	frm_interval_ = 1.0;
}

AFindPV::~AFindPV() {
	objs_.clear();
	lastmap_.reset();
}

bool AFindPV::is_freeze(PPVPT pt) {
	if (!lastmap_.unique())
		return false;
	int x = int(pt->x + 0.5);
	int y = int(pt->y + 0.5);
	return lastmap_[y * wmap_ + x] > 0;
}

bool AFindPV::prev_is_freeze(PPVPT pt) {
	int x = int(pt->x + 0.5);
	int y = int(pt->y + 0.5);
	int w = 2 * rcross_ + 1;
	int x1(x - rcross_), x2(x + rcross_);
	int y1(y - rcross_), y2(y + rcross_);
	int k;
	if (x1 < 0)
		x1 = 0;
	if (x2 >= wmap_)
		x2 = wmap_ - 1;
	if (y1 < 0)
		y1 = 0;
	if (y2 >= hmap_)
		y2 = hmap_ - 1;
	for (y = y1, k = y1 * wmap_ + x1; y <= y2; ++y) {
		for (x = x1; x <= x2; ++x, ++k) {
			if (newmap_[k])
				return true;
		}
		k += (wmap_ - w);
	}
	return false;
}

void AFindPV::UpdateFrameDelay(double dt) {
	param_.dtmax = dt * 4;
}

void AFindPV::SetDimension(int wimg, int himg) {
	wmap_ = wimg;
	hmap_ = himg;
}

void AFindPV::SetTrackMode(int mode) {
	track_mode_ = mode;
}

void AFindPV::NewSequence() {
	fno_ = -1;
	frm_interval_ = 1.0;
	objs_.clear();
}

void AFindPV::EndSequence() {
	if (fno_ != -1) {
		recheck_candidates();	// 检查候选体的有效性
		append_candidates(); 	// 尝试将该帧数据加入候选体
		complete_candidates();	// 将所有候选体转换为目标
	}
	cans_.clear();
	frmprev_.reset();
	frmlast_.reset();
}

void AFindPV::AddPoint(PPVPT pt) {
	if (fno_ != pt->fno) {
		if (fno_ != -1)
			end_frame();
		new_frame(pt->mjd);
		fno_ = pt->fno;
	}
	render_star(pt->x, pt->y);
	// 依据运动模式
	bool isfreeze = is_freeze(pt);
	if ((track_mode_ == 0 && !isfreeze) || (track_mode_ == 1 && isfreeze))
		frmlast_->pts.push_back(pt);
}

/*!
 * @brief 查看候选体
 */
PPVCANVEC& AFindPV::GetCandidate() {
	return cans_;
}

/*!
 * @brief 查看被识别的目标数量
 */
int AFindPV::GetNumber() {
	return objs_.size();
}

/*!
 * @brief 查看被识别的目标
 */
PPVOBJVEC& AFindPV::GetObject() {
	return objs_;
}

void AFindPV::render_star(double x0, double y0) {
	int x = int(x0 + 0.5);
	int y = int(y0 + 0.5);
	int w = 2 * rcross_ + 1;
	int x1(x - rcross_), x2(x + rcross_);
	int y1(y - rcross_), y2(y + rcross_);
	int k;
	if (x1 < 0)
		x1 = 0;
	if (x2 >= wmap_)
		x2 = wmap_ - 1;
	if (y1 < 0)
		y1 = 0;
	if (y2 >= hmap_)
		y2 = hmap_ - 1;
	for (y = y1, k = y1 * wmap_ + x1; y <= y2; ++y) {
		for (x = x1; x <= x2; ++x, ++k)
			++newmap_[k];
		k += (wmap_ - w);
	}
}

void AFindPV::new_frame(double mjd) {
	frmprev_ = frmlast_;
	frmlast_ = boost::make_shared<PVFRM>(mjd);
	newmap_.reset(new char[wmap_ * hmap_]);
}

void AFindPV::end_frame() {
	recheck_candidates();	// 检查候选体的有效性, 释放无效候选体
	if (frmlast_->pts.size()) {
		append_candidates(); 	// 尝试将该帧数据加入候选体
		create_candidates();	// 为未关联数据建立新的候选体
	}
	lastmap_ = newmap_;
	newmap_.reset();
//	printf("candidates : %5d, frame : %5d\n", cans_.size(),
//			frmlast_->pts.size());
}

void AFindPV::create_candidates() {
	if (!(frmprev_.unique() && frmlast_.unique()))
		return;

	PPVPTVEC &pts1 = frmprev_->pts;
	PPVPTVEC &pts2 = frmlast_->pts;
	double stepmin = param_.stepmin;
	double stepmax = param_.stepmax;
	double dxymax = param_.dxymax;
	double dx, dy;
	// 由相邻帧未关联数据构建候选体
	if (track_mode_ == 0) {	// 候选目标位置变化
		for (PPVPTVEC::iterator it1 = pts1.begin(); it1 != pts1.end(); ++it1) {
			if (!prev_is_freeze(*it1)) {
				for (PPVPTVEC::iterator it2 = pts2.begin(); it2 != pts2.end();
						++it2) {
					dx = (*it2)->x - (*it1)->x;
					dy = (*it2)->y - (*it1)->y;
					if (dx < 0.0)
						dx = -dx;
					if (dy < 0.0)
						dy = -dy;

					if (stepmin <= dx && dx <= stepmax && stepmin <= dy
							&& dy <= stepmax) {
						PPVCAN can = boost::make_shared<PVCAN>();
						can->add_point(*it1);
						can->add_point(*it2);
						if (can->fx || can->fy)
							cans_.push_back(can);
					}
				}
			}
		}
	}
	else {	// 候选目标位置不变
		for (PPVPTVEC::iterator it1 = pts1.begin(); it1 != pts1.end(); ++it1) {
			for (PPVPTVEC::iterator it2 = pts2.begin(); it2 != pts2.end();
					++it2) {
				dx = (*it2)->x - (*it1)->x;
				dy = (*it2)->y - (*it1)->y;
				if (dx < 0.0)
					dx = -dx;
				if (dy < 0.0)
					dy = -dy;

				if (dx < dxymax && dy < dxymax) {
					PPVCAN can = boost::make_shared<PVCAN>();
					can->add_point(*it1);
					can->add_point(*it2);
					cans_.push_back(can);
				}
			}
		}
	}
}

void AFindPV::append_candidates() {
	if (!cans_.size())
		return; // 无候选体立即返回

	double stepmin = param_.stepmin;
	double stepmax = param_.stepmax;
	double dxy = param_.dxymax;
	double mjd = frmlast_->mjd;
	double x1, y1, x2, y2, dx, dy, fx, fy;
	PPVPTVEC &pts = frmlast_->pts;
	PPVPT pt;
	PPVCAN can;

	// 尝试将帧数据追加至候选体
	if (track_mode_ == 0) {
		for (PPVCANVEC::iterator itcan = cans_.begin(); itcan != cans_.end();
				++itcan) {
			can = *itcan;
			// 候选体最后一个点的坐标
			pt = can->last_point();
			x1 = pt->x;
			y1 = pt->y;
			can->xy_expect(mjd, x2, y2);
			for (PPVPTVEC::iterator itpt = pts.begin(); itpt != pts.end();
					++itpt) { // 与当前帧数据交叉比对
				pt = (*itpt);
				if (pt->related == 0) {
					fx = fabs(pt->x - x1) < 1.0 ? 0 : (pt->x > x1 ? 1 : -1);
					fy = fabs(pt->y - y1) < 1.0 ? 0 : (pt->y > y1 ? 1 : -1);
					if (can->fx == fx && can->fy == fy) {
						dx = pt->x - x2;
						dy = pt->y - y2;
						if (dx < 0.0)
							dx = -dx;
						if (dy < 0.0)
							dy = -dy;
						if (dx <= dxy && dy <= dxy)
							can->add_point(pt);
					}
				}
			}
			can->update(0);
		}
	}
	else {
		for (PPVCANVEC::iterator it = cans_.begin(); it != cans_.end(); ++it) {
			can = *it;
			// 候选体最后一个点的坐标
			pt = can->last_point();
			x1 = pt->x;
			y1 = pt->y;
			for (PPVPTVEC::iterator i = pts.begin(); i != pts.end(); ++i) { // 与当前帧数据交叉比对
				if ((*i)->related == 0) {
					dx = (*i)->x - x1;
					dy = (*i)->y - y1;
					if (dx < 0.0)
						dx = -dx;
					if (dy < 0.0)
						dy = -dy;
					if (dx < dxy && dy < dxy)
						can->add_point(*i);
				}
			}
			can->update(1);
		}
	}
	// 剔除已加入候选体的数据点
	for (PPVPTVEC::iterator it = pts.begin(); it != pts.end();) {
		if ((*it)->related)
			it = pts.erase(it);
		else
			++it;
	}
	if (!pts.size())
		frmlast_.reset();
}

/*
 * 检查候选体的有效性, 判据: 时标偏差
 * - 大于阈值, 1. 数据点多于阈值, 转换为目标; 2. 数据点小于阈值, 剔除
 * - 小于阈值, 保留
 */
void AFindPV::recheck_candidates() {
	int nptmin = param_.nptmin;
	double dtmax = param_.dtmax;
	double mjd = frmlast_->mjd;
	double dt;
	int n(0);

	for (PPVCANVEC::iterator it = cans_.begin(); it != cans_.end();) {
		dt = mjd - (*it)->lastmjd;
		if (0.0 <= dt && dt <= dtmax)
			++it; // 保留. dt > 0: 原始数据未按时间严格排序
		else { // 移出候选体集合
			if ((*it)->pts.size() >= nptmin)
				candidate2object(*it); // 转换为目标
			it = cans_.erase(it);
		}
	}
}

void AFindPV::complete_candidates() {
	int nptmin = param_.nptmin;

	for (PPVCANVEC::iterator it = cans_.begin(); it != cans_.end(); ++it) {
		if ((*it)->pts.size() >= nptmin)
			candidate2object(*it);
	}
}

void AFindPV::candidate2object(PPVCAN can) {
	PPVPTVEC & pts = can->pts;
	PPVOBJ obj = boost::make_shared<PVOBJ>();
	PPVPTVEC &npts = obj->pts;

	for (PPVPTVEC::iterator it = pts.begin(); it != pts.end(); ++it) {
		npts.push_back(*it);
	}
	objs_.push_back(obj);
}

///////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */

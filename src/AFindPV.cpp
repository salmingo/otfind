/**
 * @class AFindPV 检测关联识别相邻图像中的运动目标
 * @version 0.1
 * @date 2019-11-12
 */

#include <cstdio>
#include <boost/make_shared.hpp>
#include "AFindPV.h"

namespace AstroUtil {
//////////////////////////////////////////////////////////////////////////////
AFindPV::AFindPV() {
	last_fno_ = INT_MAX;
}

AFindPV::~AFindPV() {

}

void AFindPV::NewSequence() {
	last_fno_ = INT_MAX;
	objs_.clear();
}

void AFindPV::EndSequence() {
	if (last_fno_ != INT_MAX) {
		recheck_candidates();	// 检查候选体的有效性
		append_candidates(); 	// 尝试将该帧数据加入候选体
		complete_candidates();	// 将所有候选体转换为目标
	}
	frmprev_.reset();
	frmnow_.reset();
	cans_.clear();
}

void AFindPV::AddPoint(PvPtPtr pt) {
	if (pt->fno != last_fno_) {// 上一帧数据的结束操作
		end_frame();
		new_frame(pt->fno);
	}
	frmnow_->pts.push_back(pt);
}

PvObjVec& AFindPV::GetObject() {
	return objs_;
}

void AFindPV::new_frame(int fno) {
	last_fno_ = fno;
	frmprev_  = frmnow_;
	frmnow_   = boost::make_shared<PvFrame>();
}

void AFindPV::end_frame() {
	recheck_candidates();	// 检查候选体的有效性, 释放无效候选体
	if (frmnow_.unique() && frmnow_->pts.size()) {
		append_candidates(); 	// 尝试将该帧数据加入候选体
		create_candidates();	// 为未关联数据建立新的候选体
	}
}

void AFindPV::create_candidates() {
	if (!frmprev_.unique()) return;
	PvPtVec &prev = frmprev_->pts;
	PvPtVec &now  = frmnow_->pts;
	if (!(prev.size() && now.size())) return;

	// 使用相邻帧创建候选体
	PvPtVec::iterator itprev, itnow;
	int mode;
	double x1, x2, y1, y2, dx, dy;
	for (itprev = prev.begin(); itprev != prev.end() && now.size(); ++itprev) {
		if ((*itprev)->related) continue;
		mode = -1;
		x1 = (*itprev)->x;
		y1 = (*itprev)->y;
		for (itnow = now.begin(); mode && itnow != now.end(); ) {
			if ((*itnow)->related) {
				++itnow;
				continue;
			}
			x2 = (*itnow)->x;
			y2 = (*itnow)->y;
			dx = x2 - x1;
			dy = y2 - y1;
			if (fabs(x2 - x1) > 100 || fabs(y2 - y1) > 100) ++itnow;
			else {
				PvCanPtr can = boost::make_shared<PvCan>();
				can->add_point(*itprev);
				mode = can->add_point(*itnow);
				if (mode > 0) {
					cans_.push_back(can);
					++itnow;
				}
				else if (mode == 0) {
					(*itnow)->matched = 1;
					itnow = now.erase(itnow);
				}
			}
		}
	}

	// 回归检查新创建候选体的有效性
	for (PvCanVec::iterator it = cans_.begin(); it != cans_.end(); ) {
		if ((*it)->pts.size() > 2 || !(*it)->last_point()->matched) ++it;
		else it = cans_.erase(it);
	}

#ifdef NDEBUG
	printf("create_candidates(), candidates count = %d\n", cans_.size());
#endif
}

void AFindPV::append_candidates() {
	PvCanVec::iterator itcan;
	PvPtVec::iterator itpt;
	PvPtVec &pts = frmnow_->pts;

	for (itcan = cans_.begin(); itcan != cans_.end(); ++itcan) {
		for (itpt = pts.begin(); itpt != pts.end(); ++itpt) {
			(*itcan)->add_point(*itpt);
		}
		(*itcan)->update();
	}
}

void AFindPV::recheck_candidates() {
	for (PvCanVec::iterator it = cans_.begin(); it != cans_.end();) {
		if ((last_fno_ - (*it)->last_point()->fno) <= 5) ++it;
		else { // 移出候选体集合
			(*it)->complete();
			if ((*it)->pts.size() >= 5) candidate2object(*it); // 转换为目标
			it = cans_.erase(it);
		}
	}
#ifdef NDEBUG
	printf("recheck_candidates(), candidates count = %d\n", cans_.size());
#endif
}

void AFindPV::complete_candidates() {
	for (PvCanVec::iterator it = cans_.begin(); it != cans_.end(); ++it) {
		(*it)->complete();
		if ((*it)->pts.size() >= 5) candidate2object(*it);
	}
}

void AFindPV::candidate2object(PvCanPtr can) {
	PvPtVec & pts = can->pts;
	PvObjPtr obj = boost::make_shared<PvObj>();
	PvPtVec &npts = obj->pts;
	for (PvPtVec::iterator it = pts.begin(); it != pts.end(); ++it) npts.push_back(*it);
	objs_.push_back(obj);
}

//////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */

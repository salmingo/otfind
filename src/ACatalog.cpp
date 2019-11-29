/**
 * ACatalog.cpp
 *
 *  Created on: 2016年6月13日
 *      Author: lxm
 */

#include <cmath>
#include "ACatalog.h"

namespace AstroUtil {
///////////////////////////////////////////////////////////////////////////////
ACatalog::ACatalog() {
	memset(m_pathCat, 0, sizeof(m_pathCat));
	m_nstars= 0;
	m_max   = 0;
}

ACatalog::ACatalog(const char *pathdir) {
	strcpy(m_pathCat, pathdir);
	m_nstars= 0;
	m_max   = 0;
}

ACatalog::~ACatalog() {
}

double ACatalog::sphere_range(double a1, double d1, double a2, double d2) {
	double val = cos(d1) * cos(d2) * cos(a1 - a2) + sin(d1) * sin(d2);
	return acos(val);
}

void ACatalog::SetPathRoot(const char *pathdir) {
	strcpy(m_pathCat, pathdir);
}

bool ACatalog::FindStar(double ra0, double dec0, double radius) {
	if (ra0 < 0 || ra0 > 360 || dec0 < -90 || dec0 > 90 || radius < 0.001)
		return false;
	m_csb.new_seek(ra0, dec0, radius / 60.0);	// 计算搜索范围
	return true;
}
///////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */

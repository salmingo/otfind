/*!
 Name        : otfind.cpp
 Author      : Xiaomeng Lu
 Version     : 0.1
 Copyright   : SVOM@NAOC, CAS
 Description : 从前后关联的星象中提取位置变化的OT
 (1) 从.cat文件读取X/Y/FLUX信息
 (2) 从对应的.fits文件中读取文件名, 曝光中间时刻
 (3) 输出关联后OT信息, 包含:
 文件名
 曝光中间时刻
 X
 Y
 MAG(仪器星等/1秒曝光时间), 仪器0点: 25等

 @version 0.2
 @date 2019-11-12
 @note
 - 覆盖帧频<= 1Hz
 - 覆盖恒星跟踪模式
 - 覆盖静止模式
 - 覆盖轨迹跟踪模式
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <string>
#include <list>
#include <boost/make_shared.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "ADefine.h"
#include "AFindPV.h"
#include "ACatUCAC4.h"
#include "airsdata.h"

using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace boost::posix_time;
using namespace AstroUtil;

enum {
	NDX_MIN,
	NDX_X = 0,
	NDX_Y,
	NDX_FLUX,
	NDX_MAG,
	NDX_MAGERR,
	NDX_FWHM,
	NDX_ELLIP,
	NDX_BACK,
	NDX_MAX
};

ACatUCAC4 ucac4;
AFindPV finder;

//////////////////////////////////////////////////////////////////////////////
/* 临时以固定数组管理坏像素 */
int bad_col[] = {
	1380
};

int bad_pixel[][2] = {
	{4090,   79},
	{ 943,  179},
	{3568, 1069},
	{ 840, 1201},
	{3976, 1210},
	{2989, 1236},
	{3063, 1677},
	{2404, 2307},
	{2458, 2336},
	{1867, 2340},
	{2816, 2579},
	{3226, 2894},
	{3227, 2894},
	{3276, 2908},
	{3277, 2908},
	{3319, 2942},
	{3232, 3375},
	{3794, 3408},
	{4051, 3458},
	{4041, 3473},
	{3733, 3800},
	{1509, 3953}
};
//////////////////////////////////////////////////////////////////////////////
/*!
 * @brief 检查是否坏像素
 * @note
 * bad_pixel[][]先按[][1]排序, 若[1]相同, 则按[0]排序
 */
bool is_badpixel(double x, double y) {
	int x0 = int(x + 0.5);
	int y0 = int(y + 0.5);
	int n = sizeof(bad_pixel) / sizeof(int) / 2;
	int low(0), high(n - 1), now;

	if (y0 < bad_pixel[low][1] || y0 > bad_pixel[high][1]) return false;
	while (low < n && bad_pixel[low][1] < y0) ++low;
	while (high < n && bad_pixel[high][1] > y0) --high;
	if (low > high) return false;
	for (now = low; now <= high; ++now) {
		if (x0 == bad_pixel[now][0]) return true;
	}

	return false;
}

bool load_fits(const path &filepath, FramePtr frame) {
	fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
	int status(0);
	char dateobs[40], timeobs[30], timefull[40];
	bool datefull;

	fits_open_file(&fitsptr, filepath.c_str(), 0, &status);
	fits_read_key(fitsptr, TSTRING, "DATE-OBS", dateobs,  NULL, &status);
	if (!(datefull = NULL != strstr(dateobs, "T")))
		fits_read_key(fitsptr, TSTRING, "TIME-OBS", timeobs,  NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "EXPTIME",  &frame->expdur, NULL, &status);
	fits_read_key(fitsptr, TINT, "FRAMENO",  &frame->fno, NULL, &status);
	fits_close_file(fitsptr, &status);

	if (!status) {
		if (!datefull) sprintf(timefull, "%sT%s", dateobs, timeobs);
		frame->tmobs = datefull ? dateobs : timefull;
		ptime tmmid = from_iso_extended_string(frame->tmobs) + millisec(int(frame->expdur * 500.0));
		frame->tmmid = to_iso_extended_string(tmmid);
		frame->secofday = tmmid.time_of_day().total_milliseconds() * 0.001;
		frame->mjd = tmmid.date().modjulian_day() + frame->secofday / 86400.0;
		frame->filename = filepath.filename().string();
		frame->filepath = filepath.string();
	}

	return (!status);
}

bool match_catalog(FramePtr frame) {
	// 加载wcs
	path pathwcs(frame->filepath);
	pathwcs.replace_extension(path(".wcs"));

	wcsinfo wcs;
	if (!wcs.load_wcs(pathwcs.string())) {
		cout << "failed to load WCS: " << pathwcs.filename().string() << endl;
		return false;
	}
	// 交叉匹配
	FILE *fpcat;
	path pathcat(frame->filepath);
	char line[200], ch, *token;
	char seps[] = " ";
	int pos, n;
	double buff[NDX_MAX];
	double ra, dec;

	pathcat.replace_extension(path(".cat"));
	if (NULL == (fpcat = fopen(pathcat.c_str(), "r"))) {
		cout << "failed to load CAT: " << pathcat.filename().string() << endl;
		return false;
	}

	while (!feof(fpcat)) {
		// 从cat文件读取
		if (NULL == fgets(line, 200, fpcat) || line[0] == '#') continue;
		n = strlen(line);
		if ((ch = line[n - 1]) == '\r' || ch == '\n') line[n - 1] = 0;
		if ((ch = line[n - 2]) == '\r' || ch == '\n') line[n - 2] = 0;
		pos = 0;
		token = strtok(line, seps);
		while (token && pos < NDX_MAX) {
			buff[pos] = atof(token);
			token = strtok(NULL, seps);
			++pos;
		}
		if (abs(int(buff[NDX_X] + 0.5 - 1380)) < 2) continue; // 坏列
		if (is_badpixel(buff[NDX_X], buff[NDX_Y])) continue; // 坏像素
		if (buff[NDX_FLUX] < 1.0 || buff[NDX_FWHM] <= 0.5 || buff[NDX_BACK] > 50000.0) continue;
		// 与星表交叉验证
		wcs.image_to_wcs(buff[NDX_X], buff[NDX_Y], ra, dec);
		if (!ucac4.FindStar(ra, dec, 0.5)) {
			PvPtPtr pt = boost::make_shared<pv_point>();
			pt->filename = frame->filename;
			pt->tmmid    = frame->tmmid;
			pt->fno      = frame->fno;
			pt->secofday = frame->secofday;
			pt->id       = ++frame->lastid;
			pt->x        = buff[NDX_X];
			pt->y        = buff[NDX_Y];
			pt->mag      = buff[NDX_MAG];
			pt->magerr   = buff[NDX_MAGERR];
			pt->ra       = ra;
			pt->dc       = dec;
			finder.AddPoint(pt);
		}
	}
	fclose(fpcat);
	return true;
}

void scan_directory(const string &dirname, FrameQueue &allfrm) {
	path subpath(dirname);
	path pathfit, pathwcs, pathcat;
	string filename;

	for (directory_iterator x = directory_iterator(subpath); x != directory_iterator(); ++x) {
		if (x->path().extension().string().find(".fit") != string::npos) {
			pathfit = x->path().string();
			pathwcs = pathfit;
			pathcat = pathfit;
			pathwcs.replace_extension(path(".wcs"));
			pathcat.replace_extension(path(".cat"));

			if (exists(pathwcs) && exists(pathcat)) {
				FramePtr frame = boost::make_shared<OneFrame>();
				if (load_fits(pathfit, frame)) allfrm.push_back(frame);
			}
		}
	}

	sort(allfrm.begin(), allfrm.end(), [](const FramePtr &frm1, const FramePtr &frm2) {
		return frm1->mjd < frm2->mjd;
	});
}

/*!
 * @brief 输出已关联识别目标
 * @param pvrec  关联识别算法接口
 * @param dirDst 输出数据存储目录
 * @return
 * 导出目标的数量
 */
int OutputObjects(const string &dirname, int &fsn) {
	char filename[50];
	PvObjVec & objs = finder.GetObject();
	PvPtPtr pt;
	int iy, im, id, hh, mm;
	double ss;
	FILE *fpdst;
	path filepath;
	int count(0);

	for (PvObjVec::iterator it = objs.begin(); it != objs.end(); ++it) {
		PvObjPtr obj = *it;
		PvPtVec &pts = obj->pts;

		++count;
		// 生成文件路径
		ptime tmmid = from_iso_extended_string(pts[0]->tmmid);
		ptime::date_type datemid = tmmid.date();
		ptime::time_duration_type tdt = tmmid.time_of_day();
		iy = datemid.year();
		im = datemid.month();
		id = datemid.day();
		hh = tdt.hours();
		mm = tdt.minutes();
		ss = tdt.seconds() + tdt.fractional_seconds() * 1E-6;

		sprintf(filename, "%d%02d%02d_%04d_%02d%02d%03d.txt", iy, im, id, ++fsn,
				hh, mm, int(ss * 10));
		filepath = dirname;
		filepath /= filename;
		fpdst = fopen(filepath.c_str(), "w");
		printf(">>>> %s\n", filename);
		// 写入文件内容
		for (PvPtVec::iterator i = pts.begin(); i != pts.end(); ++i) {
			pt = *i;
			tmmid = from_iso_extended_string(pt->tmmid);
			datemid = tmmid.date();
			tdt = tmmid.time_of_day();
			iy = datemid.year();
			im = datemid.month();
			id = datemid.day();
			hh = tdt.hours();
			mm = tdt.minutes();
			ss = tdt.seconds() + tdt.fractional_seconds() * 1E-6;

			fprintf(fpdst, "%s %d %02d %02d %02d %02d %06.3f %5d %9.5f %9.5f ", pt->filename.c_str(),
					iy, im, id, hh, mm, ss, pt->fno, pt->ra, pt->dc);
			if (pt->mag > 20.0)
				fprintf(fpdst, "99.99 ");
			else
				fprintf(fpdst, "%5.2f ", pt->mag);
			fprintf(fpdst, "%7.2f %7.2f\r\n", pt->x, pt->y);
		}

		fclose(fpdst);
	}
	return count;
}

int main(int argc, char **argv) {
	path subpath;
	if (argc == 1) subpath = path("./");
	else subpath = path(argv[1]);
	// 预扫目录内所有有效文件, 并按照时间增量排序
	FrameQueue allfrm;
	scan_directory(subpath.string(), allfrm);
	if (allfrm.size() < 5) {
		cout << "not found enough image to match. valid file number is " << allfrm.size() << endl;
		return -1;
	}
	cout << allfrm.size() << " valid frames will attend cross matching" << endl;

	ucac4.SetPathRoot("/Users/lxm/Catalogue/UCAC4");

	// 识别剔除视场内的恒星, 未匹配目标参与关联识别
	path pathrslt = subpath;
	pathrslt /= "stat.txt";
	FILE *fpstat = fopen(pathrslt.c_str(), "w");
	FramePtr frame;
	int fno(INT_MAX), fsn(0), total(0), count(0);

	while (allfrm.size()) {// 遍历队列
		frame = allfrm.front();
		allfrm.pop_front();
		printf("%4d : %s\n", frame->fno, frame->filename.c_str());
		if (frame->fno < fno) {// 一段连续数据完成处理
			if (fno != INT_MAX) {
				finder.EndSequence();
				if ((count = OutputObjects(subpath.string(), fsn))) {
					cout << count << " objects were found" << endl;
					total += count;
				}
				fprintf (fpstat, "%3d  %5d\n", count, total);
				fflush (fpstat);
			}
			fprintf (fpstat, "%s \t ", frame->filename.c_str());
			finder.NewSequence();
		}
		fno = frame->fno;
		match_catalog(frame);
	}
	finder.EndSequence();
	if ((count = OutputObjects(subpath.string(), fsn))) {
		cout << count << " objects were found" << endl;
		total += count;
	}
	fprintf (fpstat, "%3d  %5d\n", count, total);
	fclose(fpstat);
	cout << "!!! Total !!! < " <<  total << " > objects were found" << endl;

	return 0;
}

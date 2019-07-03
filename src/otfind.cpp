/*
 Name        : otfind.cpp
 Author      : Xiaomeng Lu
 Version     :
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
 */

#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <longnam.h>
#include <fitsio.h>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/make_shared.hpp>
#include <boost/thread.hpp>
#include "ADefine.h"
#include "AFindPV.h"
#include "ATimeSpace.h"

using namespace std;
using namespace AstroUtil;
using namespace boost::filesystem;

struct param_wcs {
	double x0, y0;	//< XY参考点
	double r0, d0;	//< RA/DEC参考点, 量纲: 弧度
	double cd[2][2];	//< 转换矩阵
	int orderA, orderB;	//< SIP改正阶数
	int ncoefA, ncoefB;	//< SIP系数数量
	double *A, *B;	//< 线性改正系数

public:
	param_wcs() {
		x0 = y0 = 0.0;
		r0 = d0 = 0.0;
		orderA = orderB = 0;
		ncoefA = ncoefB = 0;
		A = B = NULL;
	}

	virtual ~param_wcs() {
		if (A) {
			delete[] A;
			A = NULL;
		}
		if (B) {
			delete[] B;
			B = NULL;
		}
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

	void plane_to_wcs(double xi, double eta, double &ra, double &dec) {
		double fract = cos(d0) - eta * sin(d0);
		ra = cyclemod(r0 + atan2(xi, fract), A2PI);
		dec = atan2(((eta * cos(d0) + sin(d0)) * cos(ra - r0)), fract);
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

	void image_to_wcs(double x, double y, double &ra, double &dec) {
		double xi, eta;
		x -= x0;
		y -= y0;
		project_correct(x, y);
		xi = (cd[0][0] * x + cd[0][1] * y) * D2R;
		eta = (cd[1][0] * x + cd[1][1] * y) * D2R;
		plane_to_wcs(xi, eta, ra, dec);
		ra *= R2D;
		dec *= R2D;
	}
};

bool resolve_date_obs(const string &filepath, ATimeSpace &ats, double &mjd,
		double &expt, int &wimg, int &himg) {
	fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
	int status(0);
	char dateobs[40], seps[] = "-T:";
	char *token;
	int year, month, day, hour, minute;
	double second;
	
	fits_open_file(&fitsptr, filepath.c_str(), 0, &status);
	fits_read_key(fitsptr, TINT, "NAXIS1", &wimg, NULL, &status);
	fits_read_key(fitsptr, TINT, "NAXIS2", &himg, NULL, &status);
	fits_read_key(fitsptr, TSTRING, "DATE-OBS", dateobs, NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "EXPTIME", &expt, NULL, &status);
	fits_close_file(fitsptr, &status);
	if (status)
		return false;
	
	token = strtok(dateobs, seps);
	year = atoi(token);
	token = strtok(NULL, seps);
	month = atoi(token);
	token = strtok(NULL, seps);
	day = atoi(token);
	token = strtok(NULL, seps);
	hour = atoi(token);
	token = strtok(NULL, seps);
	minute = atoi(token);
	token = strtok(NULL, seps);
	second = atof(token);

	ats.SetUTC(year, month, day,
			(hour + (minute + second / 60.0) / 60.0) / 24.0);
//	mjd = ats.ModifiedJulianDay() + expt * 0.5 / DAYSEC;
	mjd = ats.ModifiedJulianDay() - expt * 0.5 / DAYSEC;

	return true;
}

int resolve_file_sn(const char *filename) {
	int n = strlen(filename), sn;
	char *buff = new char[n + 1];
	char seps[] = "-_.";
	char *token;
	strcpy(buff, filename);
	token = strtok(buff, seps);
	token = strtok(NULL, seps);
	sn = atoi(token);
	delete[] buff;
	return sn;
}

void Days2HMS(double fd, int &hh, int &mm, double &ss) {
	hh = (int) fd;
	fd = (fd - hh) * 60.0;
	mm = (int) fd;
	ss = (fd - mm) * 60.0;
}

/*!
 * @brief 输出已关联识别目标
 * @param pvrec  关联识别算法接口
 * @param dirDst 输出数据存储目录
 * @return
 * 导出目标的数量
 */
int OutputObjects(AFindPV *finder, const char *dirDst, ATimeSpace &ats,
		int &sn) {
	namespace fs = boost::filesystem;
	char filename[50];
	PPVOBJVEC & objs = finder->GetObject();
	PPVPT pt;
	int iy, im, id, hh, mm;
	double ss, fd;
	FILE *fpdst;
	fs::path path;
	int count(0);

	for (PPVOBJVEC::iterator it = objs.begin(); it != objs.end(); ++it) {
		PPVOBJ obj = *it;
		PPVPTVEC &pts = obj->pts;

		++count;
		// 生成文件路径
		ats.Mjd2Cal(pts[0]->mjd, iy, im, id, fd);
		Days2HMS(fd * 24.0, hh, mm, ss);
		sprintf(filename, "%d%02d%02d_%04d_%02d%02d%03d.txt", iy, im, id, ++sn,
				hh, mm, int(ss * 10));
		path = dirDst;
		path /= filename;
		fpdst = fopen(path.c_str(), "w");
		printf(">>>> %s\n", filename);
		// 写入文件内容
		for (PPVPTVEC::iterator i = pts.begin(); i != pts.end(); ++i) {
			pt = *i;
			ats.Mjd2Cal(pt->mjd, iy, im, id, fd);
			Days2HMS(fd * 24.0, hh, mm, ss);
			fprintf(fpdst, "%d %02d %02d %02d %02d %06.3f %5d %9.5f %9.5f ", iy,
					im, id, hh, mm, ss, pt->fno, pt->ra, pt->dc);
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

struct file_tag {
	string filename; //< 文件名
	int fno;	//< 帧编号
	double mjd;	//< 曝光中间时刻对应的修正儒略日
	double expt;	//< 曝光时间
	int wimg, himg;	//< 图像宽度和高度
};

bool inc_fno(file_tag &x1, file_tag &x2) {
	return (x1.fno < x2.fno);
}

bool inc_mjd(file_tag &x1, file_tag &x2) {
	return (x1.mjd < x2.mjd);
}

typedef vector<file_tag> ftagvec;

/*!
 * @brief 扫描fits文件, 并按曝光时间增量排序
 * @param pathname 目录名
 * @return
 * 当文件数量过少时返回false
 */
bool scan_fits(const string &pathname, ftagvec &vec) {
	namespace fs = boost::filesystem;
	ATimeSpace ats;
	fs::directory_iterator itend = fs::directory_iterator();
	string extdef = ".fit", extname;

	for (fs::directory_iterator x = fs::directory_iterator(pathname);
			x != itend; ++x) {
		extname = x->path().filename().extension().string();
		if (extname.find(extdef) == 0) {
			file_tag tag;
			tag.filename = x->path().filename().string();
			tag.fno = resolve_file_sn(tag.filename.c_str());
			resolve_date_obs(x->path().string(), ats, tag.mjd, tag.expt,
					tag.wimg, tag.himg);
			vec.push_back(tag);
		}
	}
	sort(vec.begin(), vec.end(), inc_mjd);
	return vec.size() >= 5;
}

bool try_loadwcs(const string &filepath, param_wcs &wcs) {
	// 检查cat是否存在. 若不存在, 则延时1秒等待, 直至出现
	path pathcat = filepath;
	int i(0), j(0), k(0);
	pathcat.replace_extension(path(".cat"));
	cout << "waiting for " << pathcat.filename().string();
	while (!exists(pathcat) && ++k <= 0) {
		boost::this_thread::sleep_for(boost::chrono::seconds(1));
	}
	if (!exists(pathcat)) {
		cout << "\t\t !! Failed !!" << endl;
		return false;
	}
	else
		cout << "\t\t Found" << endl;
	// 检查wcs是否存在. 若不存在, 则延时1秒等待, 直至出现并且其大小不再变化
	path pathwcs = filepath;

	cout << "waiting for " << pathwcs.filename().string();
	if (!exists(pathwcs)) {// WCS不存在
		// 先检查axy
		path pathaxy = filepath;
		pathaxy.replace_extension(path(".axy"));
		if (!exists(pathaxy)) {
			while (!exists(pathaxy) && ++i <= 0)
				boost::this_thread::sleep_for(boost::chrono::seconds(1));
			if (!exists(pathaxy)) {
				cout << "\t\t !! Failed !!" << endl;
				return false;
			}
		}

		// 等待wcs
		while (!exists(pathwcs) && ++j <= 0) {
			boost::this_thread::sleep_for(boost::chrono::seconds(1));
		}
	}
	if (!exists(pathwcs)) {
		cout << "\t\t !! Failed !!" << endl;
		return false;
	}
	else {
		cout << "\t\t Found" << endl;
		return wcs.load_wcs(filepath);
	}
}

/*
 * 输入通配文件名, *.wcs
 * 文件输出结果: ot_xxxx.txt
 */
int main(int argc, char **argv) {
	int trkmode = argc > 1 ? atoi(argv[1]) : 0;
	if (trkmode != 0 && trkmode != 1) {
		cout << "wrong track mode" << endl;
		return -1;
	}

	namespace fs = boost::filesystem;
	ftagvec filevec; // 文件名集合
	fs::path pathroot = argc < 3 ? "./" : argv[2];
	AFindPV finder;
	ATimeSpace ats;
	int fno(0);

	if (!scan_fits(pathroot.string(), filevec)) {
		cout << "require more files to find objects" << endl;
		return -2;
	} else {
		int nfile = filevec.size();
		double lastmjd = filevec[0].mjd;
		double sum(0.0), dt;
		for (int i = 1; i < nfile; ++i) {
			dt = filevec[i].mjd - lastmjd;
			sum += dt;
			lastmjd = filevec[i].mjd;
		}
		finder.UpdateFrameDelay(sum / (nfile - 1));
	}
	// 遍历文件集合
	int sn(0);
	int wimg, himg;
	double expt, ra, dec, ra0(-10.0), dec0(
			-90.0), er, ed, mjd;
	double x, y, flux, fwhm, elon;
	double limit = 60.0 / 3600.0;
	double fluxmin, fwhmmin(1.5);
	param_wcs wcs;
	FILE *fp;
	char line[100];
	char seps[] = " \r\n";
	char *token;
	int nobj(0), npt;

	finder.SetTrackMode(trkmode);
	if (trkmode)
		finder.NewSequence();

	for (ftagvec::iterator it = filevec.begin(); it != filevec.end();
			++it) {
		cout << (*it).fno << " : " << (*it).filename << endl;
		fno = (*it).fno;
		mjd = (*it).mjd;
		expt = (*it).expt;
		fluxmin = 4000.0 * expt;
		finder.SetDimension(wimg = (*it).wimg, himg = (*it).himg);

		fs::path pathfit = pathroot;
		fs::path pathwcs, pathcat;
		pathfit /= (*it).filename;
		pathwcs = pathfit;
		pathcat = pathfit;
		pathwcs.replace_extension(fs::path(".wcs"));
		pathcat.replace_extension(fs::path(".cat"));

		// 2: 加载WCS信息
		if (!try_loadwcs(pathwcs.string(), wcs)) {
			cout << "failed to load WCS from file: " << pathwcs.string()
					<< endl;
			continue;
		}

		if (!trkmode) {	// 检验中心指向. trackmode == 0时适用
			wcs.image_to_wcs(wimg * 0.5, himg * 0.5, ra, dec);
			if ((er = ra - ra0) > API)
				er = A2PI - er;
			else if (er < -API)
				er = A2PI + er;
			ed = dec - dec0;
			ra0 = ra;
			dec0 = dec;
			if (fabs(er) > limit || fabs(ed) > limit) {
				// 输出OT关联结果至文件
				finder.EndSequence();
				if ((npt = OutputObjects(&finder, pathroot.string().c_str(),
						ats, sn)) > 0) {
					nobj += npt;
					cout << npt << " objects found." << endl;
				}
				// 启动新的关联序列
				finder.NewSequence();
			}
		}

		// 4: 将文件名、时间和文件中提取的候选OT赋值给类AFindPV
		fp = fopen(pathcat.string().c_str(), "r");
		if (!fp) {
			cout << "failed to open file : " << pathcat.filename().string()
					<< endl;
			continue;
		}

		while (!feof(fp)) {
			if (fgets(line, 200, fp) == NULL || line[0] == '#')
				continue;
			token = strtok(line, seps);
			x = atof(token);
			token = strtok(NULL, seps);
			y = atof(token);
			token = strtok(NULL, seps);
			flux = atof(token);
			token = strtok(NULL, seps);
			fwhm = atof(token);
			token = strtok(NULL, seps);
			elon = atof(token);
			if (flux >= fluxmin && fwhm >= fwhmmin) {
				PPVPT pt = boost::make_shared<PVPT>();
				pt->fno = fno;
				pt->mjd = mjd;
				pt->x = x;
				pt->y = y;
				wcs.image_to_wcs(x, y, pt->ra, pt->dc);
				pt->mag = 25.0 - 2.5 * log10(flux / expt);
				finder.AddPoint(pt);
			}
		}

		fclose(fp);
	}

	// 6: 结尾: 输出OT关联结果至文件
	finder.EndSequence();
	if ((npt = OutputObjects(&finder, pathroot.string().c_str(), ats, sn))
			> 0) {
		nobj += npt;
		cout << npt << " objects found." << endl;
	}
	cout << nobj << " objects totally found." << endl;

	return 0;
}


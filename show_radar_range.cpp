#include"show_radar_range.h"
#include <iostream>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 定义相控阵雷达参数
class PhasedArrayRadarParameters {

public:
	PhasedArrayRadarParameters() {
		k = 1.380649e-23;
	};
	double azimuth; // 方位角（以弧度表示）
	double elevation; // 俯仰角（以弧度表示）
	double R; // 雷达的最大探测范围（单位：米）
	double Gt; // 雷达天线主瓣增益(dB)
	double Pt; // 发射功率(W)
	double freq;	//工作频率(Hz)
	double wavelength; // 波长(m)
	double rcs;	//目标的雷达反射截面积(m2)
	int n; //脉冲积累数
	long double k;	//波尔兹曼常数 1.380649 × 10-23J/K
	double T0;	//	雷达接收机噪声温度(K)
	double Bn;	//接收机通频带宽度(Hz)
	double Fn;	//雷达接收机噪声系数(dB)
	double L;	//雷达损耗(dB)
	double SNR_omin;	//雷达接收机最小可检测信噪比(dB)
	double Smin;	//最小可检测信号功率(W)
	double Ni;	//接收机输入噪声功率(W)

	double J;	//雷达接收的干扰功率(W)
	double S;	//雷达接收的信号功率(W)
	double KJ;	//压制系数

	double PJ;	//干扰机的发射功率(W)
	double GJ;	//干扰机的发射增益(dB)
	double theta;//干扰机相对雷达的方位角
	double theta_point5;	//雷达天线在水平方向的半功率波束宽度（度）
	double K;	//雷达天线的方向性系数0.04-0.1，K 为增益修正系数．对于高增益、波束较窄的天线，K = 0.07~0.1 ；对于低增益、波束较宽的天线，K =0.04 ~0.06 ．
	double xi;	//干扰机在雷达平面极坐标系上的方位角（度）
	double delta;	//干扰机在雷达极坐标系下的俯仰角（度）
	double delta_point5;	//为雷达天线的垂直波束宽度
	double phi_z;	//雷达天线在俯仰角上的波束指向

	double Gt_theta;	//雷达天线波束扫描到theta方向时，在干扰机方向的接收增益(dB)
	double gammaJ;//雷达天线接收干扰机信号时的极化损耗(dB)
	double RJ;	//干扰机与雷达的距离(m)
	double BJ;	//干扰机进入雷达天线的信号带宽(Hz)
	
	double getNi();
	double GetDistance();
	double GetSmin();
	//计算雷达方向上的信号功率密度 （W/m2)
	double GetS1();
	//雷达天线波束扫描到theta方向时，在干扰机方向的接收增益
	double GetGt_theta(double theta);
	//计算雷达接收的干扰功率
	double GetJ(double theta);
	//计算雷达接收的信号功率
	double GetS();
	//雷达在干扰条件下的最大探测距离
	double GetJammedDistance(double theta);
	//雷达在垂直面上的天线方向图函数
	double Get_fy_phi(double phi);
};

double PhasedArrayRadarParameters::getNi() {
	return k * T0 * Bn;
}

double PhasedArrayRadarParameters::GetSmin() {
	return k * T0 * Bn * Fn * SNR_omin;
}

double PhasedArrayRadarParameters::GetDistance() {
	double log_Pt = 10 * log10(Pt);
	double log_wavelen = 10 * log10(wavelength);
	double log_rcs = 10 * log10(rcs);
	double log_kTB = 10 * log10(k*T0*Bn);
	double log_R4 = log_Pt + Gt + Gt + log_wavelen + log_wavelen + log_rcs + 10*log10(sqrt(n)) - 10 * log10(64 * M_PI*M_PI*M_PI) - log_kTB - Fn - L - SNR_omin;
	double R4 = pow(10.0, log_R4 / 10);
	R = sqrt(sqrt(R4));
	return R;
}

double PhasedArrayRadarParameters::GetGt_theta(double theta) {
	if (0 <= theta && theta<= (theta_point5 / 2)) {
		Gt_theta = Gt;
	}
	else if ((theta_point5 / 2) <= theta&& theta <= 90) {
		Gt_theta = K * (theta_point5 / theta)*(theta_point5 / theta)*Gt;
	}
	else if (theta >= 90) {
		Gt_theta = K * (theta_point5 / 90)*(theta_point5 / 90)*Gt;
	}
	return Gt_theta;
}

double PhasedArrayRadarParameters::GetJ(double theta) {
	J = PJ * GJ*GetGt_theta(theta)*wavelength*wavelength*Bn*gammaJ / (4 * M_PI*M_PI * RJ*RJ*BJ);
	return J;
}

double PhasedArrayRadarParameters::GetS() {
	J = Pt * Gt*Gt*wavelength*wavelength*rcs*sqrt(n) / (64 * M_PI*M_PI*M_PI * R*R*R*R);
	return J;
}

double PhasedArrayRadarParameters::Get_fy_phi(double phi) {
	double f_phi = sin(2 * M_PI*cos(phi / 90)) * sin((M_PI / 2) * cos(phi / 90));
	return f_phi;
}
// 计算相控阵雷达在给定方向上的最大探测距离
double PhasedArrayRadarParameters::GetJammedDistance(double theta) {
	double log_Pt = 10.0 * log10(Pt);
	double log_KJ = 10.0 * log10(KJ);
	double log_rcs = 10.0 * log10(rcs);
	double log_PJ = 10.0 * log10(PJ);

	double log_Bn = 10.0 * log10(Bn);
	double log_RJ = 10.0 * log10(RJ);
	double log_BJ = 10.0 * log10(BJ);

	double log_R4_theta = log_Pt + Gt + Gt + log_KJ + log_rcs + 10.0 * log10(sqrt(n)) - ( 10.0 * log10(4 * M_PI)+log_PJ +GJ+ GetGt_theta(theta) +log_Bn+gammaJ- log_RJ -log_RJ - log_BJ);
	double R4 = pow(10.0, (log_R4_theta / 10.0));
	R = sqrt(sqrt(R4));
	return R;
}


double PhasedArrayRadarParameters::GetS1() {
	// 计算雷达方向上的信号功率密度 （W/m2)，这可以考虑天线增益、发射功率等因素
	double signalPower = (Pt * Gt) / (4.0 * M_PI * M_PI * R * R);
	return signalPower;

}

//计算在指定方位角、俯仰角处雷达最大探测距离
double GetDistance(double dEleR, double dDirR, double &Percent) {
	// 设置相控阵雷达参数
	PhasedArrayRadarParameters params;
	params.azimuth = dDirR; // 45度方位角（弧度）
	params.elevation = dEleR; // 30度俯仰角（弧度）
	params.Gt = 40.0; // 天线增益
	params.Pt = 1500000.0; // 发射功率
	params.wavelength = 0.056; // 波长（假设为0.056米）
	params.n = 16;
	params.rcs = 3;
	params.Bn = 1600000;
	params.Fn = 10;
	params.SNR_omin = 20;
	params.T0 = 290;
	params.L = 4;

	params.delta = 10;
	params.theta_point5 = 20;
	params.delta_point5 = 5;
	params.phi_z = 10;

	params.PJ = 10;
	params.GJ = 30;
	params.BJ = 2000000;
	params.gammaJ = 0.5;
	params.KJ = 2;
	params.RJ = 50000;
	params.K = 0.08;

	params.xi = 180;

	bool Jammed = true;
	// 计算最大探测距离
	double maxDistance = 0;
	// 雷达天线在俯仰方向上的归一化增益函数
	double fy_phi = 0;

	if (Jammed) {
		if (abs(dDirR * 180 / M_PI - params.xi) <= 180) {
			maxDistance = params.GetJammedDistance(abs(dDirR * 180 / M_PI - params.xi));
		}
		else {
			maxDistance = params.GetJammedDistance(360 - abs(dDirR * 180 / M_PI - params.xi));
		}
		// 雷达天线在俯仰方向上的归一化增益函数
		//if (abs(params.delta - params.phi_z) <= params.delta_point5/2) {
		//	fy_phi = params.Get_fy_phi(dEleR * 180 / M_PI);
		//}
		//else {
		//	fy_phi = (params.delta / params.delta_point5) * (params.delta / params.delta_point5) * params.Get_fy_phi(dEleR * 180 / M_PI) / params.K;
		//}
		fy_phi = sqrt(sqrt(abs(sin(dEleR)))) * cos(dEleR) * cos(dEleR);
	}
	else {
		maxDistance = params.GetDistance();
		// 雷达天线在俯仰方向上的归一化增益函数
		fy_phi = sqrt(sqrt(abs(sin(dEleR)))) * cos(dEleR) * cos(dEleR);
		//fy_phi = params.Get_fy_phi(dEleR * 180 / M_PI);
	}
	
	double R_azimuth_elevation = maxDistance * fy_phi;

	// 输出结果
	std::cout << "在方位角 " << params.azimuth * 180.0 / M_PI << " 度和俯仰角 " << params.elevation * 180.0 / M_PI << " 度下，相控阵雷达的最大探测距离为 " << R_azimuth_elevation << " 米。" << std::endl;
	//归一化颜色色度
	Percent = R_azimuth_elevation/ params.R;
	return R_azimuth_elevation;
}

osg::Vec4 GetColor(double Parcent, double Transparency) {
	
	if (Parcent >= 0.6) {
		return osg::Vec4(1.0f, 0.0f, 0.0f, Transparency);
	}
	else if (Parcent >= 0.5) {
		return osg::Vec4(1.0f, 0.8f, 0.0f, Transparency);
	}
	else if (Parcent >= 0.4) {
		return osg::Vec4(0.3f, 1.0f, 0.0f, Transparency);
	}
	else if (Parcent >= 0.3) {
		return osg::Vec4(0.0f, 1.0f, 1.0f, Transparency);
	}
	else if (Parcent >= 0.2) {
		return osg::Vec4(0.0f, 0.5f, 1.0f, Transparency);
	}
	else if (Parcent >= 0.1) {
		return osg::Vec4(0.0f, 0.0f, 1.0f, Transparency);
	}
	else {
		return osg::Vec4(1.0f, 0.0f, 1.0f, Transparency);
	}
	
}

int show_radar_range() {

	osgEarth::initialize();
	osgEarth::ProfileOptions profileOpts;

	//地图配置：设置缓存目录
	osgEarth::Drivers::FileSystemCacheOptions cacheOpts;
	std::string cacheDir = "D:/OSG/TEMP";
	cacheOpts.rootPath() = cacheDir;

	//
	osgEarth::Map::Options mapOpts;
	mapOpts.cache() = cacheOpts;
	mapOpts.profile() = profileOpts;

	//创建地图节点
	osg::ref_ptr<osgEarth::Map> map = new osgEarth::Map(mapOpts);
	osg::ref_ptr<osgEarth::MapNode> m_mapNode = new osgEarth::MapNode(map);

	//osgEarth::GDAL::Options gdal;
	//gdal.url() = "D:\OSG\osgearth-3.1\data\world.tif";
	m_mapNode->addChild(osgDB::readNodeFile("D:/OSG/osgearth-3.1/tests/ArcgisImage.earth"));
	//osg::ref_ptr<osgEarth::GDALImageLayer> layer = new osgEarth::GDALImageLayer();
	//layer->setURL("D:\OSG\osgearth-3.1\data\world.tif");
	//map->addLayer(layer);

	//先在坐标原点绘制三维曲面，再通过GeoTransform将绘制结果转移到雷达的经纬坐标处
	osgEarth::GeoTransform* xform = new osgEarth::GeoTransform();
	xform->setPosition(osgEarth::GeoPoint(osgEarth::SpatialReference::get("wgs84"), 110.50, 33.04, 1000));
	osg::Geode* geode = new osg::Geode();
	//设置三维曲面采样步长
	double dStep = 5;
	//俯仰方向循环
	for (double dEle = 0; dEle < 90; dEle += dStep)
	{
		//三维曲面的半透明蒙皮，通过QUAD_STRIP图元合成
		osg::Geometry* polySkin = new osg::Geometry();
		//计算蒙皮所需要的采样点数量
		int numCoordsSkin = (360 / dStep + 1) * 2;
		//申请蒙皮采样点
		osg::Vec3* myCoordsSkin = new osg::Vec3[numCoordsSkin];
		osg::Vec4Array* colorsSkin = new osg::Vec4Array;
		//各采样点的法线方向，目前还没搞清楚这个法线方向有何作用
		osg::Vec3Array* normals = new osg::Vec3Array;
		normals->push_back(osg::Vec3(0.0f, -1.0f, 0.0f));
		//水平方向循环
		int index = 0;
		for (double dDir = 0; dDir <= 360; dDir += dStep)
		{
			double dEleR = osg::DegreesToRadians(dEle);
			double dDirR = osg::DegreesToRadians(dDir);
			double dEleR2 = osg::DegreesToRadians(dEle + dStep);
			double dPercent1, dPercent2;
			//GetDistance为计算在指定方位角、俯仰角处雷达最大探测距离的核心函数，dPercent为与最大探测距离相对应的颜色色度
			double dDis = GetDistance(dEleR, dDirR, dPercent1);
			double dDis2 = GetDistance(dEleR2, dDirR, dPercent2);
			//计算蒙皮的采样点坐标，同时计算对应的显示颜色
			myCoordsSkin[index * 2] = osg::Vec3(dDis*sin(dDirR)*cos(dEleR), dDis*cos(dDirR)*cos(dEleR), dDis*sin(dEleR));
			colorsSkin->push_back(GetColor(dPercent1, 0.3));//0.3为透明度，0为完全透明，1为完全不透明
			myCoordsSkin[index * 2 + 1] = osg::Vec3(dDis2*sin(dDirR)*cos(dEleR2), dDis2*cos(dDirR)*cos(dEleR2), dDis2*sin(dEleR2));
			colorsSkin->push_back(GetColor(dPercent2, 0.3));
			if (index < 360 / dStep)
			{
				//三维曲面的经络骨骼，通过LINE_LOOP图元进行合成
				osg::Geometry* polyBone = new osg::Geometry();
				osg::Vec4Array* colorsBone = new osg::Vec4Array;
				//与蒙皮不同，由LINE_LOOP的特点决定，经络骨骼的图像无法一笔画，必须一个一个拼凑，因此顶点数量为4
				int numCoordsBone = 4;
				osg::Vec3* myCoordsBone = new osg::Vec3[numCoordsBone];
				double dDirR2 = osg::DegreesToRadians(dDir + dStep);
				double dPercent3, dPercent4;
				double dDis3 = GetDistance(dEleR2, dDirR2, dPercent3);
				double dDis4 = GetDistance(dEleR, dDirR2, dPercent4);
				myCoordsBone[0] = osg::Vec3(dDis*sin(dDirR)*cos(dEleR), dDis*cos(dDirR)*cos(dEleR), dDis*sin(dEleR));
				colorsBone->push_back(GetColor(dPercent1, 0.3));
				myCoordsBone[1] = osg::Vec3(dDis2*sin(dDirR)*cos(dEleR2), dDis2*cos(dDirR)*cos(dEleR2), dDis2*sin(dEleR2));
				colorsBone->push_back(GetColor(dPercent2, 0.3));
				myCoordsBone[2] = osg::Vec3(dDis3*sin(dDirR2)*cos(dEleR2), dDis3*cos(dDirR2)*cos(dEleR2), dDis3*sin(dEleR2));
				colorsBone->push_back(GetColor(dPercent3, 0.3));
				myCoordsBone[3] = osg::Vec3(dDis4*sin(dDirR2)*cos(dEleR), dDis4*cos(dDirR2)*cos(dEleR), dDis4*sin(dEleR));
				colorsBone->push_back(GetColor(dPercent4, 0.3));
				//设置曲线的图元类型为LINE_LOOP
				polyBone->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, numCoordsBone));
				osg::Vec3Array* verticesBone = new osg::Vec3Array(numCoordsBone, myCoordsBone);

				polyBone->setVertexArray(verticesBone);
				//设置曲线的颜色为逐点设置
				polyBone->setColorArray(colorsBone, osg::Array::BIND_PER_VERTEX);
				polyBone->setNormalArray(normals, osg::Array::BIND_OVERALL);

				geode->addDrawable(polyBone);
			}
			index++;
		}
		osg::Vec3Array* verticesSkin = new osg::Vec3Array(numCoordsSkin, myCoordsSkin);
		polySkin->setVertexArray(verticesSkin);
		//设置曲面的颜色为逐点设置
		polySkin->setColorArray(colorsSkin, osg::Array::BIND_PER_VERTEX);
		polySkin->setNormalArray(normals, osg::Array::BIND_OVERALL);
		//设置曲面的图元类型为QUAD_STRIP
		polySkin->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUAD_STRIP, 0, numCoordsSkin));
		//设置曲面半透明效果的一条关键语句，少了该句会导致半透明显示错误，暂未搞清楚该句的原理
		polySkin->getOrCreateStateSet()->setAttribute(new osg::Depth(osg::Depth::LESS, 0.0, 1.0, false));
		//开启曲面的半透明效果
		polySkin->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::ON);
		polySkin->getOrCreateStateSet()->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);
		geode->addDrawable(polySkin);
	}

	xform->addChild(geode);
	m_mapNode->addChild(xform);

	osgViewer::Viewer viewer;
	osg::ref_ptr<osgEarth::Util::EarthManipulator> mainManipulator = new osgEarth::Util::EarthManipulator;
	viewer.setCameraManipulator(mainManipulator);
	viewer.setSceneData(m_mapNode);
	osgEarth::Util::Viewpoint vp("", 110.50, 33.04, 1000,
		mainManipulator->getViewpoint().heading()->getValue(),
		mainManipulator->getViewpoint().pitch()->getValue(),
		mainManipulator->getViewpoint().range()->getValue()
	);
	//焦距
	vp.range()->set(300000.0, osgEarth::Units::METERS);
	//垂直俯仰角
	vp.pitch()->set(-10, osgEarth::Units::DEGREES);
	//水平方位角
	vp.heading()->set(0, osgEarth::Units::DEGREES);

	mainManipulator->setHomeViewpoint(vp, 2);
	//viewer.setUpViewInWindow(100, 100, 1920, 1080);

	return viewer.run();
}

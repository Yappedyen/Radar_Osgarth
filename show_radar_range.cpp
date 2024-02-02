#include"show_radar_range.h"
#include <iostream>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ����������״����
class PhasedArrayRadarParameters {

public:
	PhasedArrayRadarParameters() {
		k = 1.380649e-23;
	};
	double azimuth; // ��λ�ǣ��Ի��ȱ�ʾ��
	double elevation; // �����ǣ��Ի��ȱ�ʾ��
	double R; // �״�����̽�ⷶΧ����λ���ף�
	double Gt; // �״�������������(dB)
	double Pt; // ���书��(W)
	double freq;	//����Ƶ��(Hz)
	double wavelength; // ����(m)
	double rcs;	//Ŀ����״ﷴ������(m2)
	int n; //���������
	long double k;	//������������ 1.380649 �� 10-23J/K
	double T0;	//	�״���ջ������¶�(K)
	double Bn;	//���ջ�ͨƵ�����(Hz)
	double Fn;	//�״���ջ�����ϵ��(dB)
	double L;	//�״����(dB)
	double SNR_omin;	//�״���ջ���С�ɼ�������(dB)
	double Smin;	//��С�ɼ���źŹ���(W)
	double Ni;	//���ջ�������������(W)

	double J;	//�״���յĸ��Ź���(W)
	double S;	//�״���յ��źŹ���(W)
	double KJ;	//ѹ��ϵ��

	double PJ;	//���Ż��ķ��书��(W)
	double GJ;	//���Ż��ķ�������(dB)
	double theta;//���Ż�����״�ķ�λ��
	double theta_point5;	//�״�������ˮƽ����İ빦�ʲ�����ȣ��ȣ�
	double K;	//�״����ߵķ�����ϵ��0.04-0.1��K Ϊ��������ϵ�������ڸ����桢������խ�����ߣ�K = 0.07~0.1 �����ڵ����桢�����Ͽ�����ߣ�K =0.04 ~0.06 ��
	double xi;	//���Ż����״�ƽ�漫����ϵ�ϵķ�λ�ǣ��ȣ�
	double delta;	//���Ż����״Ｋ����ϵ�µĸ����ǣ��ȣ�
	double delta_point5;	//Ϊ�״����ߵĴ�ֱ�������
	double phi_z;	//�״������ڸ������ϵĲ���ָ��

	double Gt_theta;	//�״����߲���ɨ�赽theta����ʱ���ڸ��Ż�����Ľ�������(dB)
	double gammaJ;//�״����߽��ո��Ż��ź�ʱ�ļ������(dB)
	double RJ;	//���Ż����״�ľ���(m)
	double BJ;	//���Ż������״����ߵ��źŴ���(Hz)
	
	double getNi();
	double GetDistance();
	double GetSmin();
	//�����״﷽���ϵ��źŹ����ܶ� ��W/m2)
	double GetS1();
	//�״����߲���ɨ�赽theta����ʱ���ڸ��Ż�����Ľ�������
	double GetGt_theta(double theta);
	//�����״���յĸ��Ź���
	double GetJ(double theta);
	//�����״���յ��źŹ���
	double GetS();
	//�״��ڸ��������µ����̽�����
	double GetJammedDistance(double theta);
	//�״��ڴ�ֱ���ϵ����߷���ͼ����
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
// ����������״��ڸ��������ϵ����̽�����
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
	// �����״﷽���ϵ��źŹ����ܶ� ��W/m2)������Կ����������桢���书�ʵ�����
	double signalPower = (Pt * Gt) / (4.0 * M_PI * M_PI * R * R);
	return signalPower;

}

//������ָ����λ�ǡ������Ǵ��״����̽�����
double GetDistance(double dEleR, double dDirR, double &Percent) {
	// ����������״����
	PhasedArrayRadarParameters params;
	params.azimuth = dDirR; // 45�ȷ�λ�ǣ����ȣ�
	params.elevation = dEleR; // 30�ȸ����ǣ����ȣ�
	params.Gt = 40.0; // ��������
	params.Pt = 1500000.0; // ���书��
	params.wavelength = 0.056; // ����������Ϊ0.056�ף�
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
	// �������̽�����
	double maxDistance = 0;
	// �״������ڸ��������ϵĹ�һ�����溯��
	double fy_phi = 0;

	if (Jammed) {
		if (abs(dDirR * 180 / M_PI - params.xi) <= 180) {
			maxDistance = params.GetJammedDistance(abs(dDirR * 180 / M_PI - params.xi));
		}
		else {
			maxDistance = params.GetJammedDistance(360 - abs(dDirR * 180 / M_PI - params.xi));
		}
		// �״������ڸ��������ϵĹ�һ�����溯��
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
		// �״������ڸ��������ϵĹ�һ�����溯��
		fy_phi = sqrt(sqrt(abs(sin(dEleR)))) * cos(dEleR) * cos(dEleR);
		//fy_phi = params.Get_fy_phi(dEleR * 180 / M_PI);
	}
	
	double R_azimuth_elevation = maxDistance * fy_phi;

	// ������
	std::cout << "�ڷ�λ�� " << params.azimuth * 180.0 / M_PI << " �Ⱥ͸����� " << params.elevation * 180.0 / M_PI << " ���£�������״�����̽�����Ϊ " << R_azimuth_elevation << " �ס�" << std::endl;
	//��һ����ɫɫ��
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

	//��ͼ���ã����û���Ŀ¼
	osgEarth::Drivers::FileSystemCacheOptions cacheOpts;
	std::string cacheDir = "D:/OSG/TEMP";
	cacheOpts.rootPath() = cacheDir;

	//
	osgEarth::Map::Options mapOpts;
	mapOpts.cache() = cacheOpts;
	mapOpts.profile() = profileOpts;

	//������ͼ�ڵ�
	osg::ref_ptr<osgEarth::Map> map = new osgEarth::Map(mapOpts);
	osg::ref_ptr<osgEarth::MapNode> m_mapNode = new osgEarth::MapNode(map);

	//osgEarth::GDAL::Options gdal;
	//gdal.url() = "D:\OSG\osgearth-3.1\data\world.tif";
	m_mapNode->addChild(osgDB::readNodeFile("D:/OSG/osgearth-3.1/tests/ArcgisImage.earth"));
	//osg::ref_ptr<osgEarth::GDALImageLayer> layer = new osgEarth::GDALImageLayer();
	//layer->setURL("D:\OSG\osgearth-3.1\data\world.tif");
	//map->addLayer(layer);

	//��������ԭ�������ά���棬��ͨ��GeoTransform�����ƽ��ת�Ƶ��״�ľ�γ���괦
	osgEarth::GeoTransform* xform = new osgEarth::GeoTransform();
	xform->setPosition(osgEarth::GeoPoint(osgEarth::SpatialReference::get("wgs84"), 110.50, 33.04, 1000));
	osg::Geode* geode = new osg::Geode();
	//������ά�����������
	double dStep = 5;
	//��������ѭ��
	for (double dEle = 0; dEle < 90; dEle += dStep)
	{
		//��ά����İ�͸����Ƥ��ͨ��QUAD_STRIPͼԪ�ϳ�
		osg::Geometry* polySkin = new osg::Geometry();
		//������Ƥ����Ҫ�Ĳ���������
		int numCoordsSkin = (360 / dStep + 1) * 2;
		//������Ƥ������
		osg::Vec3* myCoordsSkin = new osg::Vec3[numCoordsSkin];
		osg::Vec4Array* colorsSkin = new osg::Vec4Array;
		//��������ķ��߷���Ŀǰ��û�����������߷����к�����
		osg::Vec3Array* normals = new osg::Vec3Array;
		normals->push_back(osg::Vec3(0.0f, -1.0f, 0.0f));
		//ˮƽ����ѭ��
		int index = 0;
		for (double dDir = 0; dDir <= 360; dDir += dStep)
		{
			double dEleR = osg::DegreesToRadians(dEle);
			double dDirR = osg::DegreesToRadians(dDir);
			double dEleR2 = osg::DegreesToRadians(dEle + dStep);
			double dPercent1, dPercent2;
			//GetDistanceΪ������ָ����λ�ǡ������Ǵ��״����̽�����ĺ��ĺ�����dPercentΪ�����̽��������Ӧ����ɫɫ��
			double dDis = GetDistance(dEleR, dDirR, dPercent1);
			double dDis2 = GetDistance(dEleR2, dDirR, dPercent2);
			//������Ƥ�Ĳ��������꣬ͬʱ�����Ӧ����ʾ��ɫ
			myCoordsSkin[index * 2] = osg::Vec3(dDis*sin(dDirR)*cos(dEleR), dDis*cos(dDirR)*cos(dEleR), dDis*sin(dEleR));
			colorsSkin->push_back(GetColor(dPercent1, 0.3));//0.3Ϊ͸���ȣ�0Ϊ��ȫ͸����1Ϊ��ȫ��͸��
			myCoordsSkin[index * 2 + 1] = osg::Vec3(dDis2*sin(dDirR)*cos(dEleR2), dDis2*cos(dDirR)*cos(dEleR2), dDis2*sin(dEleR2));
			colorsSkin->push_back(GetColor(dPercent2, 0.3));
			if (index < 360 / dStep)
			{
				//��ά����ľ��������ͨ��LINE_LOOPͼԪ���кϳ�
				osg::Geometry* polyBone = new osg::Geometry();
				osg::Vec4Array* colorsBone = new osg::Vec4Array;
				//����Ƥ��ͬ����LINE_LOOP���ص���������������ͼ���޷�һ�ʻ�������һ��һ��ƴ�գ���˶�������Ϊ4
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
				//�������ߵ�ͼԪ����ΪLINE_LOOP
				polyBone->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, numCoordsBone));
				osg::Vec3Array* verticesBone = new osg::Vec3Array(numCoordsBone, myCoordsBone);

				polyBone->setVertexArray(verticesBone);
				//�������ߵ���ɫΪ�������
				polyBone->setColorArray(colorsBone, osg::Array::BIND_PER_VERTEX);
				polyBone->setNormalArray(normals, osg::Array::BIND_OVERALL);

				geode->addDrawable(polyBone);
			}
			index++;
		}
		osg::Vec3Array* verticesSkin = new osg::Vec3Array(numCoordsSkin, myCoordsSkin);
		polySkin->setVertexArray(verticesSkin);
		//�����������ɫΪ�������
		polySkin->setColorArray(colorsSkin, osg::Array::BIND_PER_VERTEX);
		polySkin->setNormalArray(normals, osg::Array::BIND_OVERALL);
		//���������ͼԪ����ΪQUAD_STRIP
		polySkin->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUAD_STRIP, 0, numCoordsSkin));
		//���������͸��Ч����һ���ؼ���䣬���˸þ�ᵼ�°�͸����ʾ������δ������þ��ԭ��
		polySkin->getOrCreateStateSet()->setAttribute(new osg::Depth(osg::Depth::LESS, 0.0, 1.0, false));
		//��������İ�͸��Ч��
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
	//����
	vp.range()->set(300000.0, osgEarth::Units::METERS);
	//��ֱ������
	vp.pitch()->set(-10, osgEarth::Units::DEGREES);
	//ˮƽ��λ��
	vp.heading()->set(0, osgEarth::Units::DEGREES);

	mainManipulator->setHomeViewpoint(vp, 2);
	//viewer.setUpViewInWindow(100, 100, 1920, 1080);

	return viewer.run();
}

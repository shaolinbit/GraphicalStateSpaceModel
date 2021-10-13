#include <stdio.h>
#include <fstream>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Marginals.h>
#include "gssmradarmodel.h"
#include "radarvar.h"


using namespace std;
using namespace gtsam;


int main() {

	ifstream file("radardata2.txt", ios::in);
	FILE* fopresult = fopen("opresultradar.txt","w+");
	fprintf(fopresult,"Timecount Esin Esind Ez\n");
	double timed;
	int timecount;
	double txk1,txk2,tsind,msin,Exk1,Exk2,Exk3,sqrtExk3;

	if (!file.is_open()) {
		return 0;
	}


	InitRadarVar_window(0);
	InitRadarPriori_window();

	while (!file.eof()) {
		file >> timed;//N
		file >> txk1;
		file >> txk2;
		file >> msin;
		file >> Exk1;
		file >> Exk2;
		file >> Exk3;
		file >> sqrtExk3;
		timecount = round(timed / CONSTANT_RADAR_TS);

		   AddRadarPreFactor(timecount);
		   AddRadarMeasureFactor(timecount, msin);
		   OptimizeGraph_Radar_window();
			if (RadarWindowIndex == SIN_WINDOW_LENGTH) {
				Update_Radar_window_Factors();
				std::cout << std::endl;
				for (auto fb : RADARvalues)
	{

		int keyi = fb.key;
		std::cout << keyi << std::endl;
		std::cout << std::endl;
		std::cout << RADARvalues.at(keyi).cast<Vector1>()<< std::endl;
		std::cout << std::endl;
	}
	std::cout << std::endl;

			}
			fprintf(fopresult,"%d %f %f %f %f\n",timecount, txk1,
			RADARvalues.at(timecount).cast<Vector1>()(0),
				RADARvalues.at(RADARXVKey).cast<Vector1>()(0),
				 RADARvalues.at(RADARHKey).cast<Vector1>()(0));
			if (timecount >= 400)
				break;
	}

	fclose(fopresult);
	file.close();
}

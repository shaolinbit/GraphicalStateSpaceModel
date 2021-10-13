#ifndef GSSMRADARMODEL_H
#define GSSMRADARMODEL_H
#include <gtsam/nonlinear/NonlinearFactorGraph.h>

namespace gtsam
{
void InitRadarVar_window(int sinkey);
void InitRadarPriori_window();
void OptimizeGraph_Radar_window();
void Update_Radar_window_Factors();
void AddRadarPreFactor(int sinkey);
void AddRadarMeasureFactor(int sinkey, double msin);
};

#endif //

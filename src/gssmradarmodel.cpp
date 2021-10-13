#include "gssmradarmodel.h"
#include "radarvar.h"
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/Marginals.h>
#include "kfpre.h"
#include "radarmeasure.h"
#include <iostream>

namespace gtsam
{

void InitRadarVar_window(int sinkey)
{
	RADARXKey = sinkey;
	RADARXVKey = RADAR_CONSTANTV_INTEGER + sinkey;
	RADARHKey= RADAR_CONSTANTH_INTEGER + sinkey;
	RADARKeys.push_back(RADARXKey);
}

void InitRadarPriori_window()
{

	Vector1 sininit;
	sininit(0) = -100.0;

    noiseModel::Diagonal::shared_ptr modelfactorsin = noiseModel::Diagonal::Sigmas(Vector1(7.0));



	RADARPriorigraph.emplace_shared<PriorFactor<Vector1>>(RADARXKey, sininit, modelfactorsin);
    RADARvalues.insert(RADARXKey, Vector1(-100.0));



	Vector1 omegainit;
	omegainit(0) = 200.0;

	noiseModel::Diagonal::shared_ptr modelfactoromega = noiseModel::Diagonal::Sigmas(Vector1(7.0));
    RADARPriorigraph.emplace_shared<PriorFactor<Vector1>>(RADARXVKey, omegainit, modelfactoromega);
    RADARvalues.insert(RADARXVKey, Vector1(200.0));



	Vector1 radarhinit;
	radarhinit(0) = 2000.0;

	noiseModel::Diagonal::shared_ptr modelfactorh = noiseModel::Diagonal::Sigmas(Vector1(7.0));
    RADARPriorigraph.emplace_shared<PriorFactor<Vector1>>(RADARHKey, radarhinit, modelfactorh);
    RADARvalues.insert(RADARHKey, Vector1(2000.0));




	RadarWindowIndex = 0;
}

void OptimizeGraph_Radar_window()
{
	NonlinearFactorGraph OptimizeFactors;
	for (auto& bf : RADARPriorigraph)
	{
		OptimizeFactors.push_back(bf);
	}
	for (auto& bf : RADARWindowgraph)
	{
		OptimizeFactors.push_back(bf);
	}


	for (auto fb : RADARvalues)
	{
		int keyi = fb.key;
		std::cout << keyi << std::endl;
		std::cout << std::endl;
		std::cout << RADARvalues.at(keyi).cast<Vector1>()<< std::endl;
		std::cout << std::endl;
	}
	std::cout << std::endl;


	for (int i=0;i<OptimizeFactors.size();i++)
	{
		for (int keyi : OptimizeFactors.at(i)->keys())
		{
		std::cout<<keyi<<std::endl;
		}
		std::cout << std::endl;


	}
	std::cout << std::endl;


	GaussNewtonParams GaussGNparams;
	GaussNewtonOptimizer optimizer(OptimizeFactors, RADARvalues,GaussGNparams);
	double error = optimizer.error();
	Values Optimize_result = optimizer.optimize();
	double nowerror = optimizer.error();

	RADARvalues.clear();
	for (auto fb : Optimize_result)
	{
	RADARvalues.insert(fb.key,Vector1(Optimize_result.at(fb.key).cast<Vector1>()));
	}

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


void Update_Radar_window_Factors()
{
	int sinkey1 = RADARKeys.at(1);
	int sinkey0 = RADARKeys.at(0);

    NonlinearFactorGraph OptimizeFactors;
	NonlinearFactorGraph::iterator factorsrase;
	for (auto& bf : RADARPriorigraph)
	{
		OptimizeFactors.push_back(bf);
	}
	for (auto& bf : RADARWindowgraph)
	{
		OptimizeFactors.push_back(bf);
	}
    boost::shared_ptr<GaussianFactorGraph> graph_=OptimizeFactors.linearize(RADARvalues);
    boost::shared_ptr<GaussianBayesTree> bayesTree_=
    graph_->eliminateMultifrontal(boost::none, EliminatePreferCholesky);
    Matrix sin_cov =bayesTree_->marginalFactor(sinkey1, EliminatePreferCholesky)->information();
    Matrix omega_cov = bayesTree_->marginalFactor(RADARXVKey, EliminatePreferCholesky)->information();
    Matrix h_cov =  bayesTree_->marginalFactor(RADARHKey, EliminatePreferCholesky)->information();

	factorsrase = RADARWindowgraph.begin();

	RADARWindowgraph.erase(factorsrase);

	if (RADARWindowgraph.at(0)->keys().size() == 2)
	{
		if (RADARWindowgraph.at(0)->keys()[0]== sinkey1)
		{
			factorsrase = RADARWindowgraph.begin();
			RADARWindowgraph.erase(factorsrase);
		}
	}


	RADARKeys.erase(RADARKeys.begin());

	RADARvalues.erase(sinkey0);


	RADARPriorigraph.erase(RADARPriorigraph.begin(),RADARPriorigraph.end());

	for (int i=0;i<RADARPriorigraph.size();i++)
	{
		for (int keyi : RADARPriorigraph.at(i)->keys())
		{
		std::cout<<keyi<<std::endl;
		}
		std::cout << std::endl;


	}


	noiseModel::Diagonal::shared_ptr modelfactorsin = noiseModel::Diagonal::Sigmas(
	Vector1(1.0 / sqrt(sin_cov(0,0))));
	Vector1 Rb0;
	Rb0(0)=RADARvalues.at(sinkey1).cast<Vector1>()(0);
    RADARPriorigraph.emplace_shared<PriorFactor<Vector1>>(sinkey1,
    Rb0,
     modelfactorsin);

    Vector1 Rb1;
	Rb1(0)=RADARvalues.at(RADARXVKey).cast<Vector1>()(0);

    noiseModel::Diagonal::shared_ptr modelfactoromega =
    noiseModel::Diagonal::Sigmas(Vector1(1.0 / sqrt(omega_cov(0,0))));
    RADARPriorigraph.emplace_shared<PriorFactor<Vector1>>(RADARXVKey,
    Rb1,
     modelfactoromega);

      Vector1 Rb2;
	Rb2(0)=RADARvalues.at(RADARHKey).cast<Vector1>()(0);

    noiseModel::Diagonal::shared_ptr modelfactorh =
    noiseModel::Diagonal::Sigmas(Vector1(1.0 / sqrt(h_cov(0,0))));
    RADARPriorigraph.emplace_shared<PriorFactor<Vector1>>(RADARHKey,
    Rb2,
    modelfactorh);
}
void AddRadarPreFactor(int sinkey)
{


	Vector thisinsvalue(1);

    thisinsvalue(0) = (RADARvalues.at(sinkey - 1)).cast<Vector1>()(0)+
		CONSTANT_RADAR_TS * (RADARvalues.at(RADARXVKey)).cast<Vector1>()(0);

    RADARvalues.insert(sinkey, Vector1(thisinsvalue(0)));


    noiseModel::Diagonal::shared_ptr modelpreposq = noiseModel::Diagonal::Sigmas(Vector1(0.1*
     CONSTANT_RADAR_TS));
    RADARWindowgraph.emplace_shared<KFPreFactor>(sinkey - 1, sinkey, RADARXVKey,
     CONSTANT_RADAR_TS, modelpreposq);
	RADARKeys.push_back(sinkey);

}

void AddRadarMeasureFactor(int sinkey, double mpos)
{

noiseModel::Diagonal::shared_ptr modelmeasureq = noiseModel::Diagonal::Sigmas(Vector1(5.0));

 RADARWindowgraph.emplace_shared<RadarMeasureFactor>(sinkey, RADARHKey,
     mpos, modelmeasureq);
	if (RadarWindowIndex < SIN_WINDOW_LENGTH) {
		RadarWindowIndex++;
	}
}
}

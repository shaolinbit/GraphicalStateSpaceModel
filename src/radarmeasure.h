#ifndef RADARMEASUREFACTOR_H
#define RADARMEASUREFACTOR_H

#pragma once
#include <gtsam/nonlinear/NonlinearFactor.h>
namespace gtsam
{

    class RadarMeasureFactor : public NoiseModelFactor2<Vector1, Vector1>
    {
    public:
        double measured_;
    public:

        /** default constructor - only use for serialization */
        RadarMeasureFactor() {}

        /** Constructor */
        RadarMeasureFactor(int key1,int key2, double measured,
            const SharedNoiseModel& model) :
            NoiseModelFactor2(model, key1,key2), measured_(measured)
        {
        }

        ~RadarMeasureFactor()
        {

        }

        virtual gtsam::NonlinearFactor::shared_ptr clone() const {
    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
        gtsam::NonlinearFactor::shared_ptr(new RadarMeasureFactor(*this))); }

 Vector evaluateError(const Vector1& p1, const Vector1& p2,
      boost::optional<Matrix&> H1 = boost::none,
      boost::optional<Matrix&> H2 = boost::none) const
    {
    double dnorm= sqrt(p1(0)*p1(0) + p2(0) * p2(0));
    if (H1)
    {
    *H1 = Matrix::Identity(1,1)* p1(0)/dnorm;
    }
    if(H2)
    {
    *H2 =  Matrix::Identity(1,1)* p2(0)/dnorm;
    }
    Vector1 resultm;
    resultm(0)=dnorm - measured_;
    return resultm;
  }

    };
};

#endif


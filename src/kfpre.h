#ifndef KFPREFACTOR_H
#define KFPREFACTOR_H

/**
 *  @file  KFPreFactor.h
 **/
#pragma once
#include <gtsam/nonlinear/NonlinearFactor.h>
namespace gtsam
{

    class KFPreFactor : public NoiseModelFactor3<Vector1, Vector1, Vector1>
    {
    public:
        double SampleT_;

    public:

        /** default constructor - only use for serialization */
        KFPreFactor() {}

        /** Constructor */
        KFPreFactor(Key key1, Key key2,Key keyvel, double SampleT,
            const SharedNoiseModel& model) :
            NoiseModelFactor3(model, key1,key2,keyvel), SampleT_(SampleT)
        {

        }

        ~KFPreFactor()
        {

        }

          virtual gtsam::NonlinearFactor::shared_ptr clone() const {
    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
        gtsam::NonlinearFactor::shared_ptr(new KFPreFactor(*this))); }

         Vector evaluateError(const Vector1& p1, const Vector1& p2, const Vector1& p3,
      boost::optional<Matrix&> H1 = boost::none,
      boost::optional<Matrix&> H2 = boost::none,
      boost::optional<Matrix&> H3 = boost::none) const
    {
    const size_t p = 1;
    if (H1) *H1 = -Matrix::Identity(p,p);
    if (H2) *H2 = Matrix::Identity(p,p);
    if (H3) *H3 = Matrix::Identity(p,p)*(-SampleT_);
    Vector1 resultm(p2(0)-p1(0));
    resultm(0)-=p3(0) * SampleT_;
    return resultm;
  }

    };
};

#endif

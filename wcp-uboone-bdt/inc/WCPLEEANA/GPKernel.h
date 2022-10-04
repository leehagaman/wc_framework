#pragma once

#include "GPPoint.h"

#include <cassert> 

#include "TVector.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/IFunction.h"

class GPKernel
{
public:
  GPKernel(std::vector<double> pars, std::vector<bool> log_scales={false,false,false,false,false}, double noise=0.)
    :  fPars(pars),
       doLogScales(log_scales),
       fNoise(noise)
  {
  }

  virtual double Mag(GPPoint pt1, GPPoint pt2) const = 0;

  virtual double Element(GPPoint pt1, GPPoint pt2, int dpar_idx = -1) const = 0;

  void SetParameters(std::vector<double> pars) { fPars = pars; }
  double NPar() const { return fPars.size(); }
  virtual void SetThetas(std::vector<double> thetas) = 0;

  std::vector<double> GetParameters() const { return fPars; }
  
  TMatrixDSym operator()(std::vector<GPPoint> pts, int dpar_idx = -1) const;
  
  TMatrixD KernelBlock(std::vector<GPPoint> pts1, std::vector<GPPoint> pts2) const;
  TVectorD KernelDiag(std::vector<GPPoint> pts) const;

protected:
  std::vector<double> fPars;
  std::vector<bool> doLogScales;
  double fNoise;
};

//----------------------------------------------------------------------

class RBFKernel: virtual public GPKernel
{
public:
  RBFKernel(std::vector<double> pars, std::vector<bool> log_scales={false,false,false,false,false}, double noise = 1.e-10) : 
   GPKernel(pars, log_scales, noise) { 
     assert(pars.size() == 6 && "Input hyperparameters need to be size 6 for RBF Kernel"); 
   }
  
  double Mag(GPPoint pt1, GPPoint pt2) const override;
  double Element(GPPoint pt1, GPPoint pt2, int dpar_idx = -1) const override;
  void SetThetas(std::vector<double> thetas) override;
};

//----------------------------------------------------------------------

class RationalQuadraticKernel: virtual public GPKernel
{
public:
  RationalQuadraticKernel(std::vector<double> pars, std::vector<bool> log_scales={false,false,false,false,false}, double noise = 1.e-10) : 
   GPKernel(pars, log_scales, noise) { 
     assert(pars.size() == 3 && "Input hyperparameters need to be size 3 for Rational Quadratic Kernel"); 
   }
  
  double Element(GPPoint pt1, GPPoint pt2, int dpar_idx = -1) const override;
  void SetThetas(std::vector<double> thetas) override;
};

#include "ClassHeaders.h"
#include "matplotlibcpp.h"
#include <cmath>
#include <vector>

namespace plt = matplotlibcpp;

int main()
{
  // Setup Problem
  int n = 0;
  double TargetManeuvers = 96.6;
  double HeadingErrorDeg = 20.; // initial deviation of the missile from the collision triangle
  double HeadingErrorRad = HeadingErrorDeg/57.3; //Heading error in radians
  double EffNavRatio = 5.; // Effective navigation ratio (gain) Larger ratio removes heading error faster
  double T = 0.; // Initiating time constant
  double S = 0.;
  double H = 0.01;
  double TargetAngularVelocity = 0;
  double GuidanceCommand = 0;

  AirVehicle Missile;
  AirVehicle Target;

  Missile.VelocityMagnitude = 3000.;
  Missile.LocationX = 0.;
  Missile.LocationY = 10000.;

  Target.VelocityMagnitude = 1000.;
  Target.LocationX = 40000.;
  Target.LocationY = 10000.;

  double TargetBeta = 135.; //Angle of target velocity wrt to reference frame
  Target.VelocityX = -Target.VelocityMagnitude*cos(TargetBeta);
  Target.VelocityY = Target.VelocityMagnitude*sin(TargetBeta);

  double MTRangeX = Target.LocationX - Missile.LocationX;
  double MTRangeY = Target.LocationY - Missile.LocationY;
  double RelativeSeparation = sqrt( MTRangeX*MTRangeX +
                                    MTRangeY*MTRangeY );

  double LineOfSight = atan2(MTRangeY,MTRangeX);
  double MissileLeadAngle = asin(Target.VelocityMagnitude *
      sin(TargetBeta+LineOfSight)/Missile.VelocityMagnitude);
  double Theta = MissileLeadAngle + LineOfSight;

  Missile.VelocityX = Missile.VelocityMagnitude * cos(Theta + HeadingErrorRad);
  Missile.VelocityY = Missile.VelocityMagnitude * sin(Theta + HeadingErrorRad);

  double MTRelativeVelocityX = Target.VelocityX - Missile.VelocityX;
  double MTRelativeVelocityY = Target.VelocityY - Missile.VelocityY;
  double ClosingVelocity = -(MTRangeX*MTRelativeVelocityX +
    MTRangeY*MTRelativeVelocityY) / RelativeSeparation;


  //Allocating storage for outputs
  std::vector<double> ArrayT;
  std::vector<double> ArrayTargetLocation1;
  std::vector<double> ArrayTargetLocation2;
  std::vector<double> ArrayMissileLocation1;
  std::vector<double> ArrayMissileLocation2;
  std::vector<double> ArrayGuidanceCommand;
  std::vector<double> ArrayRelativeSeparation;

  // Second-Order Runge Kutta Method
  while (ClosingVelocity >= 0)
  {
    // Decreases step size for lower RelativeSeparation to increase accuracy
    if (RelativeSeparation < 1000) { H = 0.0002;}
    else { H = 0.01;}

    double BetaOld = TargetBeta;
    double TOldX = Target.LocationX;
    double TOldY = Target.LocationY;
    double MOldX = Missile.LocationX;
    double MOldY = Missile.LocationY;
    double MVelOldX = Missile.VelocityX;
    double MVelOldY = Missile.VelocityY;
    int STEP = 1;
    int FLAG = 0;
    while (STEP <= 1)
    {
      if (FLAG == 1)
      {
        STEP = 2;
        TargetBeta = TargetBeta + H*TargetAngularVelocity;
        Target.LocationX += H*Target.VelocityX;
        Target.LocationY += H*Target.VelocityY;
        Missile.LocationX += H*Missile.VelocityX;
        Missile.LocationY += H*Missile.VelocityY;
        Missile.VelocityX += H*Missile.AccelerationX;
        Missile.VelocityY += H*Missile.AccelerationY;
        T += H;
      }

      MTRangeX = Target.LocationX - Missile.LocationX;
      MTRangeY = Target.LocationY - Missile.LocationY;
      RelativeSeparation = sqrt(MTRangeX*MTRangeX + MTRangeY*MTRangeY); //Calculate magnitude of current range vector
      MTRelativeVelocityX = Target.VelocityX - Missile.VelocityX;       //X-comp of relative velocity
      MTRelativeVelocityY = Target.VelocityY - Missile.VelocityY;       //Y-comp of relative velocity
      ClosingVelocity = -(MTRangeX*MTRelativeVelocityX+MTRangeY*MTRelativeVelocityY)/RelativeSeparation;          // Closing velocity
      LineOfSight = atan2(MTRangeY,MTRangeX);                                                                               // LOS angle
      double LineOfSightD = (MTRangeX*MTRelativeVelocityY-MTRangeY*MTRelativeVelocityX)/(RelativeSeparation*RelativeSeparation);  // Rate of change of LOS angle
      GuidanceCommand = EffNavRatio*ClosingVelocity*LineOfSightD;                       // Nc - magnitude of missile guidance command
      Missile.AccelerationX = -GuidanceCommand*sin(LineOfSight);                     // X-comp of missile acceleration
      Missile.AccelerationY = GuidanceCommand*cos(LineOfSight);                      // Y-comp of missile acceleration
      Target.VelocityX = -Target.VelocityMagnitude*cos(TargetBeta);                      // X-comp of target velocity
      Target.VelocityY = Target.VelocityMagnitude*sin(TargetBeta);                       // Y-comp of target velocity
      TargetAngularVelocity = TargetManeuvers/Target.VelocityMagnitude;                           // Angular velocity of target (Slope)
      FLAG = 1;
    }

    FLAG = 0;
    TargetBeta = .5*(BetaOld+TargetBeta+H*TargetAngularVelocity);
    Target.LocationX = .5*(TOldX+Target.LocationX+H*Target.VelocityX);
    Target.LocationY = .5*(TOldY+Target.LocationY+H*Target.VelocityY);
    Missile.LocationX = .5*(MOldX+Missile.LocationX+H*Missile.VelocityX);
    Missile.LocationY = .5*(MOldY+Missile.LocationY+H*Missile.VelocityY);
    Missile.VelocityX = .5*(MVelOldX+Missile.VelocityX+H*Missile.AccelerationX);
    Missile.VelocityY = .5*(MVelOldY+Missile.VelocityY+H*Missile.AccelerationY);
    S=S+H;

    if (S >= .09999)
    {
      S=0.;
      ArrayT.push_back(T);
      ArrayTargetLocation1.push_back(Target.LocationX);
      ArrayTargetLocation2.push_back(Target.LocationY);
      ArrayMissileLocation1.push_back(Missile.LocationX);
      ArrayMissileLocation2.push_back(Missile.LocationY);
      ArrayGuidanceCommand.push_back(GuidanceCommand/32.2);
      ArrayRelativeSeparation.push_back(RelativeSeparation);
    }
  }

  //Plot values using matplotlibcpp
  plt::plot(ArrayTargetLocation1,ArrayTargetLocation2,ArrayMissileLocation1,ArrayMissileLocation2);
  plt::title("Two-dimensional tactical missile-target engagement simulation");
  plt::show();
  return 0;
}

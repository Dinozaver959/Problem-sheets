using System;

namespace VasicekCalibrationAndPricing
{  
    public class VasicekModel
    {
        public const int numModelParams = 3;

        private const int kIndex = 0;
        private const int thetaIndex = 1;
        private const int sigmaIndex = 2;

        private double r0;
        private double k;
        private double theta;
        private double sigma;

        public VasicekModel() : this(0, 0, 0, 0)
        {
        }

        public VasicekModel(VasicekModel otherModel) : 
            this(otherModel.r0, otherModel.k, otherModel.theta, otherModel.sigma)
        {
        }

        public VasicekModel(double r0, double k, double theta, double sigma)
        {
            this.r0 = r0;
            this.k = k;
            this.theta = theta;
            this.sigma = sigma;
        }

        public VasicekModel(double r0, double[] paramsArray)
            : this(r0, paramsArray[kIndex], paramsArray[thetaIndex], paramsArray[sigmaIndex])
        {
        }

        public double GetK() { return k; }
        public double GetTheta() { return theta; }
        public double GetSigma() { return sigma; }

        private double A(double t, double T)
        {
            double B_t_T = B(t, T);
            double sigmaSq = sigma * sigma;
            double kSq = k * k;
            double insideExp = (theta - 0.5 * sigmaSq / kSq) * (B_t_T - T + t) - (0.25*sigmaSq / k) * B_t_T * B_t_T;
            return Math.Exp(insideExp);
        }

        private double B(double t, double T)
        {
            return (1.0 / k) * (1.0 - Math.Exp(-k * (T - t)));
        }

        public double PriceZCB(double bondMaturity)
        {
            double A_T = A(0, bondMaturity);
            double B_T = B(0, bondMaturity);
            return A_T * Math.Exp(-B_T * r0);
        }

        public double PriceZCBCall(double optionExercise, double strike, double bondMaturity)
        {
            double sigmaP = sigma * Math.Sqrt((1 - Math.Exp(-2 * k * optionExercise)) / (2 * k)) * B(optionExercise, bondMaturity);
            double h = (1.0 / sigmaP) * Math.Log(PriceZCB(bondMaturity) / (PriceZCB(optionExercise) * strike)) + sigmaP / 2.0;
            double N_h = MathNet.Numerics.Distributions.Normal.CDF(0,1,h);
            double N_hMinusSigmaP = MathNet.Numerics.Distributions.Normal.CDF(0, 1,h - sigmaP);
            return PriceZCB(bondMaturity) * N_h - strike * PriceZCB(optionExercise) * N_hMinusSigmaP;
        }

        // For use in calibration
        public double[] ConvertVasicekModeCalibrationParamsToArray()
        {
            double[] paramsArray = new double[VasicekModel.numModelParams];
            paramsArray[kIndex] = k;
            paramsArray[thetaIndex] = theta;
            paramsArray[sigmaIndex] = sigma;
            return paramsArray;
        }            
    }
}

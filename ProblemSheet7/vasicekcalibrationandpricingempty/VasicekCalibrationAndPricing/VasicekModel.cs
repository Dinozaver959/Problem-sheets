using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace VasicekCalibrationAndPricing
{
    

    public class VasicekModel
    {
        public const int numModelParams = 3;
        
        private double r0, k, theta, sigma;
        public VasicekModel() { r0 = 0; k = 0; theta = 0; sigma = 0; }
        public VasicekModel(VasicekModel otherModel) { r0 = otherModel.r0; k = otherModel.k; theta = otherModel.theta; sigma = otherModel.sigma; }
        public VasicekModel(double r0, double k, double theta, double sigma)
        {
            this.r0 = r0; this.k = k; this.theta = theta; this.sigma = sigma;
        }
               

        public double PriceZCB(double bondMaturity)
        {
            throw new NotImplementedException("Method not implemented... that's part of the exercise.");
        }

        public double PriceZCBCall(double optionExercise, double strike, double bondMaturity)
        {
            throw new NotImplementedException("Method not implemented... that's part of the exercise.");
        }


        

        public double GetK() { return k; }
        public double GetTheta() { return theta; }
        public double GetSigma() { return sigma; }
        
    }
}

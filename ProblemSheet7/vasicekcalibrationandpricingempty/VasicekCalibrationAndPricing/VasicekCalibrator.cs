using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace VasicekCalibrationAndPricing
{
    public class CalibrationFailedException : Exception
    {
        public CalibrationFailedException()
		{
		}
        public CalibrationFailedException(string message)
            : base(message)
		{
		}
    }
    
    public enum CalibrationOutcome { NotStarted, FinishedOK, FailedMaxItReached, FailedOtherReason };
    
    
    public class VasicekCalibrator
    {
        // Set whatever member variables you see fit

        
        public VasicekCalibrator(double r0, double accuracy, int maxIterations)
        {
            // You should modify this as you see fit
            throw new NotImplementedException("Method not implemented... that's part of the exercise.");
        }

        public void SetGuessParameters(double k, double theta, double sigma)
        {
            throw new NotImplementedException("Method not implemented... that's part of the exercise.");
        }
              
        // Used by Alglib minimisation algorithm, this is where you will need to 
        // calculate the mean square error between market and model prices
        public void CalibrationObjectiveFunction(double[] paramsArray, ref double func, object obj)
        {
            throw new NotImplementedException("Method not implemented... that's part of the exercise.");
        }


        // This should use Alglib BFGS method for minimising the mean square error and thus choosing the parameters
        public void Calibrate()
        {
            throw new NotImplementedException("Method not implemented... that's part of the exercise.");
        }

        // This should be able to report how that calibration went and what the mean square error is
        public void GetCalibrationStatus(ref CalibrationOutcome calibOutcome, ref double pricingError)
        {
            throw new NotImplementedException("Method not implemented... that's part of the exercise.");
        }

        // This should return a calibrated Vasicek model
        public VasicekModel GetCalibratedModel()
        {
            throw new NotImplementedException("Method not implemented... that's part of the exercise.");
        }


    }
}

using System;

namespace VasicekCalibrationAndPricing
{
    public static class Vasicek 
    {
        public static double VasicekZCBPrice(double maturity, double r0, double k, double theta, double sigma)
        {
            VasicekModel m = new VasicekModel(r0, k, theta, sigma);
            return m.PriceZCB(maturity);            
        }

        public static double VasicekZCBCallPrice(double bondMaturity, double optionExercise, double strike, double r0, double k, double theta, double sigma)
        {
            VasicekModel m = new VasicekModel(r0, k, theta, sigma);
            return m.PriceZCBCall(optionExercise, strike, bondMaturity);
        }
        
        public static object[,] CalibrateVasicekParameters(
            double kGuess,
            double thetaGuess,
            double sigmaGuess,
            double r0,
            double[] strikesArray,
            double[] bondMaturitiesArray,
            double[] optionExerciseTimesArray,
            double[] observedPricesArray,
            double accuracy,
            int maxIterations)
        {
           
            if (strikesArray.Length != bondMaturitiesArray.Length
                || bondMaturitiesArray.Length != optionExerciseTimesArray.Length
                || optionExerciseTimesArray.Length != observedPricesArray.Length)
            {
                throw new ArgumentException("CalibrateVasicekParameters: strikes, optionExerciseTimes, bondMaturities and observedPrices must be of same length.");
            }         

            VasicekCalibrator vasicekCalibrator = new VasicekCalibrator(r0, accuracy, maxIterations);
            vasicekCalibrator.SetGuessParameters(kGuess, thetaGuess, sigmaGuess);

            int numObservedOptions = strikesArray.Length;
            for (int optionIdx = 0; optionIdx < numObservedOptions; ++optionIdx)
            {
                vasicekCalibrator.AddObservedOption(bondMaturitiesArray[optionIdx],
                                                        optionExerciseTimesArray[optionIdx],
                                                        strikesArray[optionIdx],
                                                        observedPricesArray[optionIdx]);
            }
            vasicekCalibrator.Calibrate();
            CalibrationOutcome outcome = CalibrationOutcome.NotStarted;
            double calibrationError = 0;
            vasicekCalibrator.GetCalibrationStatus(ref outcome, ref calibrationError);

            VasicekModel calibratedModel = vasicekCalibrator.GetCalibratedModel();

            // for output
            const int numCols = 2;
            const int numRows = 5;
            object[,] output = new object[numRows, numCols];
            output[0, 0] = "k"; output[0, 1] = calibratedModel.GetK();
            output[1, 0] = "theta"; output[1, 1] = calibratedModel.GetTheta();
            output[2, 0] = "sigma"; output[2, 1] = calibratedModel.GetSigma();
            output[3, 0] = "Minimizer Status";
            if (outcome == CalibrationOutcome.FinishedOK)
            {
                output[3, 1] = "OK";
            }
            else if (outcome == CalibrationOutcome.FailedMaxItReached)
            {
                output[3, 1] = "Reached max. num. iterations.";
            }
            else if (outcome == CalibrationOutcome.FailedOtherReason)
            {
                output[3, 1] = "Failed.";
            }
            else
            {
                output[3, 1] = "Unknown outcome.";
            }
            {
                output[4, 0] = "Pricing error"; output[4, 1] = Math.Sqrt(calibrationError);
            }
            return output;
        }
    }    
}

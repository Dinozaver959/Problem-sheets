using System;
using System.Collections.Generic;

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
    
    public enum CalibrationOutcome
    {
        NotStarted,
        FinishedOK,
        FailedMaxItReached,
        FailedOtherReason
    };
    
    public struct ZCBCallOptionMarketData
    {
        public double bondMaturity;
        public double optionExercise;
        public double strike;
        public double marketMidPrice;
    }
    
    public class VasicekCalibrator
    {
        private const double defaultAccuracy = 10e-3;
        private const int defaultMaxIterations = 500;
        private double accuracy;
        private int maxIterations;

        private LinkedList<ZCBCallOptionMarketData> marketOptionsList;
        private double r0; // initial interest rate, this is observed, no need to calibrate to options

        private CalibrationOutcome outcome;
        
        private double[] calibratedParams;

        public VasicekCalibrator()
        {
            accuracy = defaultAccuracy;
            maxIterations = defaultMaxIterations;
            marketOptionsList = new LinkedList<ZCBCallOptionMarketData>();
            r0 = 0;
            calibratedParams = new double[] { 0.1, 0.05, 0.05 };
        }

        public VasicekCalibrator(double r0, double accuracy, int maxIterations)
        {
            this.r0 = r0;
            this.accuracy = accuracy;
            this.maxIterations = maxIterations;
            marketOptionsList = new LinkedList<ZCBCallOptionMarketData>();
            calibratedParams = new double[] { 0.1, 0.1, 0.5 };
        }

        public void SetGuessParameters(double k, double theta, double sigma)
        {
            VasicekModel m = new VasicekModel(r0, k, theta, sigma);
            calibratedParams = m.ConvertVasicekModeCalibrationParamsToArray();
        }

        public void AddObservedOption(double bondMaturity, double optionExercise, double strike, double mktMidPrice)
        {
            ZCBCallOptionMarketData observedOption;
            observedOption.bondMaturity = bondMaturity; 
            observedOption.optionExercise = optionExercise; 
            observedOption.strike = strike;
            observedOption.marketMidPrice = mktMidPrice;
            marketOptionsList.AddLast(observedOption);
        }

        // Calculate difference between observed and model prices
        public double CalcMeanSquareErrorBetweenModelAndMarket(VasicekModel m)
        {
            double meanSqErr = 0;
            foreach (ZCBCallOptionMarketData option in marketOptionsList)
            {
                double bondMaturity = option.bondMaturity;
                double optionExercise = option.optionExercise;
                double strike = option.strike;
                double modelPrice = m.PriceZCBCall(optionExercise, strike, bondMaturity);

                double difference = modelPrice - option.marketMidPrice;
                meanSqErr += difference * difference;
            }
            return meanSqErr;
        }
        
        // Used by Alglib minimisation algorithm
        public void CalibrationObjectiveFunction(double[] paramsArray, ref double func, object obj)
        {
            VasicekModel m = new VasicekModel(r0, paramsArray);
            func = CalcMeanSquareErrorBetweenModelAndMarket(m);
        }
        
        public void Calibrate()
        {
            outcome = CalibrationOutcome.NotStarted;

            double[] initialParams = new double[VasicekModel.numModelParams];
            calibratedParams.CopyTo(initialParams, 0);  // a reasonable starting guees
            double epsg = accuracy;
            double epsf = accuracy; //1e-4;
            double epsx = accuracy;
            double diffstep = 1.0e-6;
            int maxits = maxIterations;
            double stpmax = 0.05;

            

            alglib.minlbfgsstate state;
            alglib.minlbfgsreport rep;
            alglib.minlbfgscreatef(1, initialParams, diffstep, out state);
            alglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits);
            alglib.minlbfgssetstpmax(state, stpmax);

            // this will do the work
            alglib.minlbfgsoptimize(state, CalibrationObjectiveFunction, null, null);
            double[] resultParams = new double[VasicekModel.numModelParams];
            alglib.minlbfgsresults(state, out resultParams, out rep);

            System.Console.WriteLine("Termination type: {0}", rep.terminationtype);
            System.Console.WriteLine("Num iterations {0}", rep.iterationscount);
            System.Console.WriteLine("{0}", alglib.ap.format(resultParams, 5));

            if (rep.terminationtype == 1			// relative function improvement is no more than EpsF.
                || rep.terminationtype == 2			// relative step is no more than EpsX.
                || rep.terminationtype == 4)
            {    	// gradient norm is no more than EpsG
                outcome = CalibrationOutcome.FinishedOK;
                // we update the ''inital parameters''
                calibratedParams = resultParams;
            }
            else if (rep.terminationtype == 5)
            {	// MaxIts steps was taken
                outcome = CalibrationOutcome.FailedMaxItReached;
                // we update the ''inital parameters'' even in this case
                calibratedParams = resultParams;

            }
            else
            {
                outcome = CalibrationOutcome.FailedOtherReason;
                throw new CalibrationFailedException("Vasicek model calibration failed badly.");
            }
        }

        public void GetCalibrationStatus(ref CalibrationOutcome calibOutcome, ref double pricingError)
        {
            calibOutcome = outcome;
            VasicekModel m = new VasicekModel(r0, calibratedParams);
            pricingError = CalcMeanSquareErrorBetweenModelAndMarket(m);
        }

        public VasicekModel GetCalibratedModel()
        {
            VasicekModel m = new VasicekModel(r0, calibratedParams);
            return m;
        }
    }
}

using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;

using VasicekCalibrationAndPricing;

namespace VasicekCmdLine
{
    public class Examples
    {
        public static void TestVasicekPricing()
        {
            double r0 = 0.01;
            double k = 0.5;
            double theta = 0.1;
            double sigma = 0.1;

            VasicekModel model = new VasicekModel(r0, k, theta, sigma);
            double[] bondMaturities = new double[] { 2, 10 };
            Console.WriteLine("Bond prices");
            for (int i = 0; i < bondMaturities.Length; ++i)
            {
                double T = bondMaturities[i];
                double p = model.PriceZCB(T);
                Console.WriteLine("T={0}, p={1}", T, p);
            }
            
            double[] optionExerciseTimes = new double[] {1, 2 };
            double[] optionStrikes = new double[] { 0.9, 0.9 };

            Console.WriteLine("European call options on ZCBs, prices");
            for (int i = 0; i < bondMaturities.Length; ++i)
            {
                double T = bondMaturities[i], S = optionExerciseTimes[i], K = optionStrikes[i];
                double v = model.PriceZCBCall(S, K, T); 
                Console.WriteLine("T={0}, S={1}, K={2}, v={3}", T, S, K, v);
            }
        }

        public static void TestVasicekCalibration()
        {
            double r0 = 0.01;
            double k = 0.5;
            double theta = 0.01;
            double sigma = 0.2;

            VasicekModel model = new VasicekModel(r0, k, theta, sigma);
            double[] bondMaturities = new double[] { 1, 2, 10 };
            double[] optionExerciseTimes = new double[] {0.5, 1, 2 };
            double[] optionStrikes = new double[] { 0.9, 0.9, 0.9 };
            double[] prices = new double[bondMaturities.Length];
            //Console.WriteLine("European call options on ZCBs, prices");
            for (int i = 0; i < bondMaturities.Length; ++i)
            {
                double T = bondMaturities[i], S = optionExerciseTimes[i], K = optionStrikes[i];
                //prices[i] = model.PriceZCBCall(S, K, T);
                prices[i] = 0.05;
            }

            VasicekCalibrator calibrator = new VasicekCalibrator(r0, 1e-15, 2000);
            for (int i = 0; i < bondMaturities.Length; ++i)
            {
                calibrator.AddObservedOption(bondMaturities[i], optionExerciseTimes[i], optionStrikes[i], prices[i]);
            }
            calibrator.Calibrate();
            double error = 0;
            CalibrationOutcome outcome = CalibrationOutcome.NotStarted;
            calibrator.GetCalibrationStatus(ref outcome, ref error);
            Console.WriteLine("Calibration outcome: {0} and error: {1}", outcome, error);

        }

        
    }
}

using System;

namespace MinimisationWithAlglib
{
    class Example
    {
        
        // Alglib needs the function to minimize to be implemented in a method with 
        // the following 'signature'
        private static void FuncToMinimizeForAlglib(double[] xIn, ref double func, object obj)
        {
            // f(x,y) = x^4 + y^5
            double x = xIn[0]; double y = xIn[1];
            func = Math.Pow(x,4) + Math.Pow(y,5);
        }
        
        // This just runs the example minimisation and prints output to console
        public static void MinimisationExample()
        {
            int dim = 2;
            double[] xZeroAlglib = new double[dim];
            xZeroAlglib[0] = -100; xZeroAlglib[1] = 100; 


            double[] xFinalAlglib = new double[dim];

            // these are parameters for various stopping criteria going into the LBFGS algorithm
            // you can read about the precise meaning in optimization.cs
            double epsg = 1e-20;
            double epsf = 1e-20; 
            double epsx = 1e-20;
            int maxits = 500;

            // how big delta to take when approximating derivatives
            double diffstep = 1.0e-6;
            
            // What is the maximum 'step' the algorithm is allowed to take from iteration to iteration.
            // Small value will lead to slow convergence. Large value may lead to 'overshoot'
            double stpmax = 0.5;

            
            alglib.minlbfgsstate state;
            alglib.minlbfgsreport rep;
            alglib.minlbfgscreatef(1, xZeroAlglib, diffstep, out state);
            alglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits);
            alglib.minlbfgssetstpmax(state, stpmax);

            // this will do the work
            alglib.minlbfgsoptimize(state, FuncToMinimizeForAlglib, null, null);
            alglib.minlbfgsresults(state, out xFinalAlglib, out rep);

            // You can read about the meaning of different termination types in optimization.cs
            System.Console.WriteLine("Termination type: {0}", rep.terminationtype);
            System.Console.WriteLine("Num iterations {0}", rep.iterationscount);
            System.Console.WriteLine("{0}", alglib.ap.format(xFinalAlglib, 5)); 
        }
    }
}

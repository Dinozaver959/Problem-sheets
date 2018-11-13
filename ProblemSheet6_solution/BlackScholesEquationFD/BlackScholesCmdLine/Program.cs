using System;

class MainClass
{
	public static void Main (string[] args)
	{
		// Model params
		double r = 0.05;
		double sigma = 0.1;
		double K = 100;
		double T = 1;
        double S0 = 100;
        BlackScholesFormula bsFormula = new BlackScholesFormula(new BlackScholesModelParams(r, sigma));
        double bsPrice = bsFormula.PutPrice(new EuropeanCallPutParams(T, S0, K));
        Func<double, double> putPayoff = (S) => Math.Max(K - S, 0);
        BlackScholesFiniteDifferenceSolver solverFD =
                new BlackScholesFiniteDifferenceSolver(T, putPayoff, r, sigma, 3 * K, 10, 10);
        double bsPriceFD = solverFD.Price(100);



        uint N, M;

        int numberRefinments = 5;

        // test convergence w.r.t. number of partitions of space interval
        N = 200; M = 100;
        for (int refinementIndex = 0; refinementIndex < numberRefinments; ++refinementIndex,M *= 2) 
        {
            BlackScholesFiniteDifferenceSolver solverForThisLevelOfRefinement = 
                new BlackScholesFiniteDifferenceSolver (T, putPayoff, r, sigma, 5 * K, N, M);
            double error = Math.Abs (bsPrice - solverForThisLevelOfRefinement.Price (S0));
            Console.WriteLine ("Space partitions: {0}, time steps: {1}, error: {2}", M, N, error);
        }

        // test convergence w.r.t. number of time steps
        N = 10; M = 8001;
        for (int refinementIndex = 0; refinementIndex < numberRefinments; ++refinementIndex, N*=2) 
        {
            BlackScholesFiniteDifferenceSolver solverForThisLevelOfRefinement = 
                new BlackScholesFiniteDifferenceSolver (T, putPayoff, r, sigma, 5 * K, N, M);
            double error = Math.Abs (bsPrice - solverForThisLevelOfRefinement.Price (S0));
            Console.WriteLine ("Space partitions: {0}, time steps: {1}, error: {2}", M, N, error);
        }
        Console.WriteLine ("Finished. Press any key.");
		//Console.ReadKey ();
	}
}


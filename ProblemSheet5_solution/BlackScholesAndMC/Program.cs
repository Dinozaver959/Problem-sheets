using System;

class TestBlackScholeMC
{
	private double sigma = 0.1;
	private double S = 100;
	private double K = 100;
	private double r = 0.05;
	private double T = 1;

	public TestBlackScholeMC (double sigma = 0.1, //volatilty
	                          double S = 100, //underlying price
	                          double K = 100, //strike price
	                          double r = 0.05, //risk free rate
	                          double T = 1) //time to maturity in years
	{
		this.sigma = sigma;
		this.S = S;
		this.K = K;
		this.r = r;
		this.T = T;
	}

	public double CalcAvgMonteCaroError (int nForAverage, int N)
	{
		double callOptionPriceBS = BlackScholesFormula.CalculateCallOptionPrice (sigma, S, K, r, T);
		double avgErr = 0;
		for (int i = 0; i < nForAverage; ++i) {
			BlackScholesMonteCarlo mc = new BlackScholesMonteCarlo (sigma, r, N);
			double mcPrice = mc.CalculateEuropeanCallOptionPrice (S, K, T);
			avgErr += Math.Abs (mcPrice - callOptionPriceBS);
		}
		return avgErr / nForAverage;
	}
}

class Program
{
    
    static void TestImpliedVol(double sigma = 0.1, //volatilty
		double S = 100, //underlying price
		double K = 100, //strike price
		double r = 0.05, //risk free rate
		double T = 1) //time to maturity in years
	{
		double CallOptionPrice = BlackScholesFormula.CalculateCallOptionPrice (sigma, S, K, r, T);
		double PutOptionPrice = BlackScholesFormula.CalculatePutOptionPrice (sigma, S, K, r, T);

		System.Console.WriteLine ("Price of a call option using Black Scholes formula is: {0}", BlackScholesFormula.CalculateCallOptionPrice (sigma, S, K, r, T));
		System.Console.WriteLine ("Price of a put option using Black Scholes formula is: {0}", BlackScholesFormula.CalculatePutOptionPrice (sigma, S, K, r, T));

		System.Console.WriteLine ("Implied volatiliy of a call option given the price from Black Scholes is: {0}", BlackScholesImpliedVolEuropean.CalculateImpliedVolCall (10, S, K, r, T));
		System.Console.WriteLine ("Implied volatiliy of a put option given the price from Black Scholes is: {0}", BlackScholesImpliedVolEuropean.CalculateImpliedVolPut (3, S, K, r, T));

	}

    static void TestMonteCarloBS(double sigma = 0.1, //volatilty
        double S = 100, //underlying price
        double K = 100, //strike price
        double r = 0.05, //risk free rate
        double T = 1) //time to maturity in years
	{
        TestBlackScholeMC test = new TestBlackScholeMC(sigma, S, K, r, T);

		double numberOfScalings = 5;
		int N = 100;
		for (int i = 0; i < numberOfScalings; ++i) 
		{
			Console.WriteLine("N = {0}, err={1:0.00}", N, test.CalcAvgMonteCaroError(1000,N));
			N *= 4;
		}
		Console.WriteLine ("We should be seeing the error roughly halving every time we quadruple the number of samples.");
	}

	static void TestMonteCarloBSAsian (double sigma, double S, double K, double r)
	{
		int N = 100000;
		BlackScholesMonteCarlo mc = new BlackScholesMonteCarlo (sigma, r, N);	

		double[] monitorDates = new double[]{ 0.25, 0.5, 1 };
		double exerciseT = 1;
		Console.WriteLine("Asian: {0}", mc.CalculateAsianCallOptionPrice (S, K, monitorDates, exerciseT));
	}

	static void Main (string[] args)
	{

		//test our option pricing methods
		double sigma = 0.1; //volatilty
		double S = 100; //underlying price
		double K = 100; //strike price
		double r = 0.05; //risk free rate
        double T = 1; // time to maturity as year fraction

        try
        {
            TestImpliedVol(sigma, S, K, r, T);
            TestMonteCarloBS(sigma, S, K, r, T);
            TestMonteCarloBSAsian(sigma, S, K, r);
        }
        catch (Exception e)
        {
            System.Console.WriteLine("Error: " + e.Message);
        }

		System.Console.ReadKey ();

           
	}
}


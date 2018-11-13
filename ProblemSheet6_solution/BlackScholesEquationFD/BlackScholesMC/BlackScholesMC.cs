using System;

using MathNet.Numerics;

public class BlackScholesMC
{
	// problem parameters
	private double T;
	private Func<double, double> payoff;
	private double r; // risk-free rate
	private double sigma; // Black--Scholes vol

	private int N; // number of samples

	public BlackScholesMC ()
	{
		T = 0;
		payoff = (x) => 0;
		r = 0; 
		sigma = 0;
	}

	public BlackScholesMC(BlackScholesMC other)
	{
		T = other.T;
		payoff = new Func<double, double> (other.payoff);
		r = other.r;
		sigma = other.sigma;
	}

	public BlackScholesMC(double maturity, 
		Func<double, double> payoffFunction,  
		double riskFreeRate, 
		double sigma, 
		int numMCSamples)
	{
		T = maturity;
		payoff = payoffFunction;
		r = riskFreeRate;
		this.sigma = sigma;
		this.N = numMCSamples;
	}

	public double Price(double S)
	{
		double[] normals = new double[N];
		MathNet.Numerics.Distributions.Normal.Samples(normals, 0,1);
		double sumOfSamplePayoffs = 0;

		for (int i = 0; i < N; ++i) {
			double S_T = S * Math.Exp ((r - 0.5 * sigma * sigma) * T + sigma * Math.Sqrt (T) * normals [i]);
			double samplePayoff = payoff (S_T);
			sumOfSamplePayoffs += Math.Exp(-r*T)*samplePayoff; // don't forget discounting
		}
		return sumOfSamplePayoffs / N;
	}
}



using System;
using MathNet.Numerics.Distributions;

class BlackScholesMonteCarlo
{

	private double sigma;
	private double r; 
	private int N;

	public BlackScholesMonteCarlo (double sigma, double r, int N) 
	{
        if (sigma <= 0 || N <= 0)
            throw new System.ArgumentException("Need sigma >0, N > 0.");
        this.sigma = sigma;
		this.r = r;
		this.N = N;
	}

	//Returns the price of a European Call Option using a Monte Carlo method with N iterations, volatility sigma, underlying price S,
	//strike price K, risk free rate r and time to maturity in years T.
	public double CalculateEuropeanCallOptionPrice (double S, double K, double T)
	{
        if (S <= 0 || K <= 0 || T <= 0)
            throw new System.ArgumentException("Need S, K, T > 0");   
		double[] Z = new double[N];
		Normal.Samples (Z, 0, 1);

		double C_n = 0;
		double sqrtOfT = Math.Sqrt (T);
		for (int i = 0; i < N; i++) {
			C_n += Math.Exp (-r * (T)) * Math.Max (S * Math.Exp ((r - 0.5 * sigma * sigma) * (T) + sigma * sqrtOfT * Z [i]) - K, 0);
		}

		return C_n / N;
	}

	public double CalculatePutOptionPrice (double S, double K, double T)
	{
        if (S <= 0 || K <= 0 || T <= 0)
            throw new System.ArgumentException("Need S, K, T > 0"); 

        return CalculateEuropeanCallOptionPrice (S, K, T) - S + K * Math.Exp (-r * (T));
	}

	private double SolveGBM(double S0, double T, double dW)
	{
		return S0 * Math.Exp ((r - 0.5 * sigma * sigma) * T + sigma * dW);
	}

    private void CheckAsianOptionInputs(double[] T, double exerciseT)
    {
        if (T.Length == 0)
            throw new System.ArgumentException("Need at least one monitoring date for Asian option.");

        if (T[0] <= 0)
            throw new System.ArgumentException("The first monitoring date must be positive.");


        for (int i = 1; i < T.Length; ++i)
        {
            if (T[i - 1] >= T[i])
                throw new System.ArgumentException("Monitoring dates must be increasing.");
        }

        if (T[T.Length - 1] > exerciseT)
            throw new System.ArgumentException("Last monitoring dates must not be greater than the exercise time.");

    }

	public double CalculateAsianCallOptionPrice (double S, double K, double[] T, double exerciseT)
	{
        if (S <= 0 || K <= 0 || exerciseT <= 0)
            throw new System.ArgumentException("Need S, K, T > 0");

        CheckAsianOptionInputs(T, exerciseT);

        Func<double, double> callPayoff = (x) => Math.Max (x - K, 0);
		return CalculateAsianOptionPrice (S, K, T, exerciseT, callPayoff);
	}

	private double  CalculateAsianOptionPrice (double S, double K, double[] T, double exerciseT, 
	                                        Func<double, double> payoffFn)
	{
		// Note that here we don't need to use an Euler scheme as we have the exact solution
		int M = T.Length;
		double[,] pathsOfS = new double[N,M];

		// We generate the 'path' for the T_i that we need
		double[] Z = new double[N];
		for (int m = 0; m < M; ++m)
		{
			Normal.Samples (Z, 0, 1);
			double deltaT = T[0];
			if(m>0)
				deltaT = T[m]-T[m-1];
			
			double sqrtOfDeltaT = Math.Sqrt (deltaT);
			for (int i = 0; i < N; ++i)
			{
				double s_i_m_minus_one = S;
				if (m>0)
					s_i_m_minus_one = pathsOfS[i,m-1];
				double dW = sqrtOfDeltaT*Z[i];
				pathsOfS[i,m] = SolveGBM(s_i_m_minus_one, deltaT, dW);
			}
		}

		// the is the average for Asian option price
		double averageAlongPath;

		// this is the MC average
		double averagePayoff = 0;
		for (int i = 0; i < N; ++i)
		{
			averageAlongPath = 0;
			for (int m = 0; m < M; ++m)
			{
				averageAlongPath += pathsOfS[i,m];
			}
			averageAlongPath *= 1.0/M;
			averagePayoff += payoffFn(averageAlongPath);
		}
		averagePayoff *= 1.0/N;

		return Math.Exp(-r*exerciseT)*averagePayoff;
	}
}


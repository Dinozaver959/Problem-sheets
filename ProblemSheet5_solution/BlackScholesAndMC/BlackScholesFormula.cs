using System;

class BlackScholesFormula
{
	public static double CalculateCallOptionPrice (double sigma, double S, double K, double r, double T)
	{
		if (sigma <= 0 || T <= 0 || K <= 0 || S <= 0)
            throw new System.ArgumentException("Need sigma > 0, T > 0, K > 0 and S > 0.");

		double d1 = (Math.Log (S / K) + (r + (sigma * sigma) / 2) * (T)) / (sigma * Math.Sqrt (T));
		double d2 = d1 - sigma * Math.Sqrt (T);
        return S * MathNet.Numerics.Distributions.Normal.CDF(0, 1, d1) - K * Math.Exp(-r * T) * MathNet.Numerics.Distributions.Normal.CDF(0, 1, d2);
    }

    public static double GetCallFromPutPrice(double S, double K, double r, double T, double putPrice)
    {
        return putPrice - Math.Exp(-r * T) * K + S;
    }

    public static double GetPutFromCallPrice(double S, double K, double r, double T, double callPrice)
    {
        return callPrice - S + K * Math.Exp(-r * T);
    }

	public static double CalculatePutOptionPrice (double sigma, double S, double K, double r, double T)
	{
        if (sigma <= 0 || T <= 0 || K <= 0 || S <= 0)
            throw new System.ArgumentException("Need sigma > 0, T > 0, K > 0 and S > 0.");

        return GetPutFromCallPrice(S, K, r, T, CalculateCallOptionPrice(sigma, S, K, r, T));
	}

    
}






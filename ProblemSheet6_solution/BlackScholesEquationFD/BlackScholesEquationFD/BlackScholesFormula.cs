using System;


public struct EuropeanCallPutParams
{

	public double T; // time to expiry
	public double strike; // strike
	public double s0; // initial asset price

	public EuropeanCallPutParams(double timeToExpiry, double initialAsstePrice, double strike)
	{
		this.T = timeToExpiry;
		this.s0 = initialAsstePrice;
		this.strike = strike;
	}
}


public struct BlackScholesModelParams
{
    public double sigma; // vol 
    public double r; // risk-free rate

    public BlackScholesModelParams(double riskFreeRate, double sigma)
    {
        this.r = riskFreeRate;
        this.sigma = sigma; 
    }
}

public class BlackScholesFormula
{
    private BlackScholesModelParams p;
    

    public BlackScholesFormula(BlackScholesModelParams modelParameters)
    {
        p = modelParameters;    
    }


    public double CallPrice(EuropeanCallPutParams optionParameters)
    {
        double r = p.r; // Interest rate
        double sig = p.sigma; // Volatility
        double K = optionParameters.strike; // Strike price
        double T = optionParameters.T; // Expiry date
        double S = optionParameters.s0;
        
        double tmp = sig * Math.Sqrt(T);
        double d1 = (Math.Log(S / K) + (r + 0.5 * sig * sig) * T) / tmp;
        double d2 = d1 - tmp;
        return S * SpecialFunctions.N(d1) - (K * Math.Exp(-r * T) * SpecialFunctions.N(d2));
    }

	public double PutPrice(EuropeanCallPutParams optionParameters)
	{
		double callPrice = CallPrice (optionParameters);
		// put-call parity
		return callPrice - optionParameters.s0 + optionParameters.strike*Math.Exp(-p.r * optionParameters.T);
	}
    
}


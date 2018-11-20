using System;

namespace Problem_Sheet_5_Assignment_EwenGillespie
{
	//Class for various applications of the Black Scholes Formula for option pricing
	public class BlackScholesModel
	{
		// Key B-S model parameters
		double S;
		double K;
		double r;
		double sigma; 
		double T; 

		// Class constructor instance for specified parameters
		public BlackScholesModel(double currentAssetPrice,double strikeAssetPrice,double riskFreeRate, double volatility, double contractTerm)
		{
			S = currentAssetPrice;
			K = strikeAssetPrice;
			r = riskFreeRate;
			sigma = volatility;
			T = contractTerm;
		}



		//Returns price of call option given specified model parameters using B-S Formula
		public double CallOptionPrice()
		{
			double d1 = (Math.Log( S / K) + (r + sigma*sigma / 2) * T) / (sigma * Math.Sqrt(T));

			double d2 = d1 - sigma * Math.Sqrt(T);

			return S * N(d1)- K * Math.Exp(-r * T) *N(d2);
		}

		//Returns price of put option determined by put-call parity and price of call option under B-S Formula
		public double PutOptionPrice()
		{			
			return CallOptionPrice() + Math.Exp(-r * T) * K - S;
		}
	}
}


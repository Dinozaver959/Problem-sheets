
using System;

class BlackScholesImpliedVolEuropean
{
	public static double CalculateImpliedVolCall (double C, double S, double K, double r, double T, double x0 = 0.5, double maxErr = 1e-6, int N = 10000)
	{
        if (C <= 0 || T <= 0 || K <= 0 || S <= 0 || maxErr <= 0 || N <= 0)
            throw new System.ArgumentException("Need C > 0, T > 0, K > 0, S > 0, maxErr > 0, N > 0.");

		Func<double, double> F = (x) => {
			return C - BlackScholesFormula.CalculateCallOptionPrice (x, S, K, r, T);
		};
		NewtonSolver s = new NewtonSolver (maxErr, N);
		return s.Solve (F, null, x0);

	}

	public static double CalculateImpliedVolPut (double P, double S, double K, double r, double T, double x0 = 0.5, double maxErr = 1e-6, int N = 10000)
	{
        if (P <= 0 || T <= 0 || K <= 0 || S <= 0 || maxErr <= 0 || N <= 0)
            throw new System.ArgumentException("Need P > 0, T > 0, K > 0, S > 0, maxErr > 0, N > 0.");

        double callPrice = BlackScholesFormula.GetCallFromPutPrice(S, K, r, T, P);
        if (callPrice < 0)
            throw new System.ArgumentException("Input arguments violate put/call parity.");

        return CalculateImpliedVolCall(callPrice, S, K, r, T, x0, maxErr, N);

	}
}


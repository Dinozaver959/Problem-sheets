using System;

public class NewtonSolver
{
	const double h = 1e-6;
	private double maxError;
	private int maxIter;

	public NewtonSolver(double maxError, int maxIter) { this.maxError = maxError; this.maxIter = maxIter; }

	//Returns an approximate solution s to f(x)=0 using Newtons method starting at x0 and such that |f(x)| < maxError
	//where f is a function R->R and fPrime is f' (if null then approximated using central difference).
	//Throws exception if derivative is too small or number of iterations exceeds maxIter.
	public double Solve(Func<double, double> f, Func<double, double> fPrime, double x0)
	{
		if (fPrime == null)
        {
			fPrime = (x) => (f (x + h) - f (x-h)) / (2*h);
		}

        for (int i = 0; i < maxIter; i++)
        {
            double f_x0 = f(x0);
            if (Math.Abs(f(x0)) < maxError)
                return x0;
            double fPrime_x0 = fPrime(x0);
            if (Math.Abs(fPrime_x0) < 1e-16)
                throw new SystemException("Newton's Method failed - derivative too small");
            x0 = x0 - f(x0) / fPrime(x0);
        }
        throw new SystemException("Newton's Method did not converge");       
	}
}


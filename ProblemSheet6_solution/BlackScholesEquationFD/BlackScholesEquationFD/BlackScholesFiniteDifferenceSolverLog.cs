using System;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Solvers;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;




public class BlackScholesFiniteDifferenceSolverLog
{
	// problem parameters
	private double T;
	private Func<double, double> payoff;
	private double r; // risk-free rate
	private double sigma; // Black--Scholes vol

	// discretization parameters
	private double R1,R2; // we are solving on (0,T) x (0,R).
	private int N; // number of steps in space
	private int M; // number of steps in time
	private double h;

	// needed for work
	// matrix A represents FD discretization of the `spatial' part of the PDE
	private Matrix<double> A; 

	// matrix S represents the system needed to be solved for every step of implicit scheme: S u^{J,K}_n = u^{J,K}_{n-1}.
	private Matrix<double> S; 

	// for solving linear system iteratively;
	Iterator<double> monitor;
	BiCgStab solver;

	// check whether we've computed the solution
	bool solved;
	Vector<double> u; // the solution.

	// default constructor
	public BlackScholesFiniteDifferenceSolverLog()
	{
		T = 0; payoff = null; r = 0; sigma = 0; 
		R1 = 0; R2 = 0; N = 0; M = 0; h = 0;

		A = Matrix<double>.Build.Sparse(M, M);
		S = Matrix<double>.Build.Sparse(M, M);
		u = Vector<double>.Build.Dense (M);

		solved = false;
	}

	// copy constructor
	public BlackScholesFiniteDifferenceSolverLog(BlackScholesFiniteDifferenceSolverLog other)
	{
		T = other.T; r = other.r; sigma = other.sigma;
		payoff = new Func<double, double> (other.payoff); // make a copy

		R1 = other.R1; R2 = other.R2; N = other.N; M = other.N; h = other.h;

		A = Matrix<double>.Build.SparseOfMatrix(other.A); // make a copy
		S = Matrix<double>.Build.SparseOfMatrix(other.S); // make a copy

		solved = other.solved;
		u = Vector<double>.Build.DenseOfVector (other.u); // make a copy
	}

	// useful constructor
	public BlackScholesFiniteDifferenceSolverLog(double maturity, 
		Func<double, double> payoffFunction,  
		double riskFreeRate, 
		double sigma, 
		double R1, double R2,
		uint numTimeSteps,
		uint numSpacePartitions)
	{
		T = maturity; 
		payoff = payoffFunction;
		r = riskFreeRate;
		this.sigma = sigma; 
		this.R1 = R1; this.R2 = R2; N = (int)numTimeSteps; M = (int)numSpacePartitions; h = (R2-R1) / (M - 1);

		SetUpFiniteDifferenceMatrix();
		SetUpSolverMatrix();

		u = Vector<double>.Build.Dense (M);
	}

	private void SetUpFiniteDifferenceMatrix()
	{
		double rMinusHalfSigmaSq = r - 0.5 * sigma * sigma;
		double rMinusHalfSigmaSqPlus = Math.Max (rMinusHalfSigmaSq, 0); 
		double rMinusHalfSigmaSqMinus = -Math.Min (rMinusHalfSigmaSq, 0);

		double c = -0.5 * sigma * sigma / (h * h) - rMinusHalfSigmaSqPlus / h; 
		double b = sigma * sigma / (h * h) + rMinusHalfSigmaSqPlus / h + rMinusHalfSigmaSqMinus / h + r;
		double a = -0.5 * sigma * sigma / (h * h) - rMinusHalfSigmaSqMinus / h;

		// i is the row index, j is the column index
		Func<int, int, double> matrixEntry = (i, j) =>
		{
			if (i == j && (i == 0 || i == M-1)) // diagonal entry first and last row for bdry conds.
				return 1.0; // for Dirichle bdry at the left end of interval
			else if (i == j) // diagonal entry
			{
				return b;
			}
			else if (i == M - 1 && j == i - 1) //below diagonal, last row
			{
				return 1.0;
			}
			else if (i > 0 && i < M - 1 &&  j == i + 1) // above diagonal
			{
				return c; 
			}
			else if (i < M - 1 && j == i - 1) //below diagonal
			{
				return a; 
			}
			else 
				return 0;
		};
		A = Matrix<double>.Build.Sparse(M, M, matrixEntry);
		//Console.WriteLine (A.ToString ());
	}

	private void SetUpSolverMatrix()
	{
		double tau = T / N;
		Matrix<double> I = Matrix<double>.Build.SparseIdentity(M); 
		S = (I + tau*A);
		//Console.WriteLine ("Determinant of S: {0}", S.Determinant ());
		System.Diagnostics.Debug.Assert (Math.Abs(S.Determinant ()) >  1e-7);
		Console.WriteLine("Set up matrix for solver");
	}

	private void SetUpSolver()
	{
		// Stop calculation if 1000 iterations reached during calculation
		IterationCountStopCriterion<double> iterationCountStopCriterion = new IterationCountStopCriterion<double>(100);

		// Stop calculation if residuals are below 1E-10 --> the calculation is considered converged
		ResidualStopCriterion<double> residualStopCriterion = new ResidualStopCriterion<double>(1e-8);

		// Create monitor with defined stop criteria
		monitor = new Iterator<double>(iterationCountStopCriterion, residualStopCriterion);

		// Load all suitable solvers from current assembly. Below in this example, there is user-defined solver
		// "class UserBiCgStab : IIterativeSolverSetup<double>" which uses regular BiCgStab solver. But user may create any other solver
		// and solver setup classes which implement IIterativeSolverSetup<T> and pass assembly to next function:
		solver = new BiCgStab();
	}

	public Vector<double> ApproxInitialCondition()
	{
		Vector<double> u0 = Vector<double>.Build.Dense(M);
		for (int m = 0; m < M; ++m) {
			double S = Math.Exp (R1 + m * h);
			u0 [m] = payoff (S); 
		}
		//Console.WriteLine (u0.ToString ());
		return u0;
	}


	public Vector<double> Solve()
	{
		double tau = T / N;
		SetUpSolver();

		Vector<double> uOld = ApproxInitialCondition();
		Vector<double> uNew = Vector<double>.Build.Dense(M);
		double gamma = uOld [M - 1] - uOld [M - 2]; // for neumann bdry at x=R2
		for (int k = 0; k < N; k++)
		{
			double uOldZero = Math.Exp(-r*tau*k)*payoff(Math.Exp(R1));
			uOld [0] = uOldZero*(1+tau);
			uOld [M - 1] = gamma*(1+tau); // enforce boundary

			// Must solve ( I - tau * A ) * uNew =  uOld i.e. S * uNew = uOld
			uNew = S.SolveIterative(uOld, solver, monitor);
			//Console.WriteLine ("Error at 0 is: {0}", Math.Abs((S * uNew) [0] - uOld[0]));
			uOld = uNew;
			//Console.Write("Step {0}, ", k);
		}
		//Console.WriteLine();
		u = uNew;
		solved = true;
		return u;
	}
	public double Price(double S)
	{
		double x = Math.Log (S);
		if (x + h >= R2)
			throw new Exception ("S too large for solver domain.");
		if (x-h <= R1)
			throw new Exception ("S must be > 0.");
		if (!solved)
			Solve ();
		int indexBelow = (int)Math.Floor ((x-R1) / h);
		int indexAbove = (int)Math.Ceiling ((x-R1) / h); 

		return u [indexBelow] + (x-(R1+indexBelow*h))*(u [indexAbove] - u [indexBelow]) / h; 
	}
}


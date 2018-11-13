using System;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Solvers;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;

public class BlackScholesFiniteDifferenceSolver
{
	// problem parameters
	private double T;
	private Func<double, double> payoff;
	private double r; // risk-free rate
	private double sigma; // Black--Scholes vol

	// discretization parameters
	private double R; // we are solving on (0,T) x (0,R).
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
	public BlackScholesFiniteDifferenceSolver()
	{
		T = 0;  r = 0; sigma = 0; payoff = null;
		R = 0; N = 0; M = 0; h = 0;

		A = Matrix<double>.Build.Sparse(M, M);
		S = Matrix<double>.Build.Sparse(M, M);
		u = Vector<double>.Build.Dense (M);

		solved = false;
	}

	// copy constructor
	public BlackScholesFiniteDifferenceSolver(BlackScholesFiniteDifferenceSolver other)
	{
		T = other.T; r = other.r; sigma = other.sigma; 
		payoff = new Func<double, double> (other.payoff); // make a copy

		R = other.R; N = other.N; M = other.N; h = other.h;

		A = Matrix<double>.Build.SparseOfMatrix(other.A); // make a copy
		S = Matrix<double>.Build.SparseOfMatrix(other.S); // make a copy

		solved = other.solved;
		u = Vector<double>.Build.DenseOfVector (other.u); // make a copy
	}

	// useful constructor
	public BlackScholesFiniteDifferenceSolver(double maturity, 
		Func<double, double> payoffFunction, 
		double riskFreeRate, 
		double sigma, 
		double R,
		uint numTimeSteps,
		uint numSpacePartitions)
	{
        if (maturity <= 0 || sigma <= 0 || R <= 0 || numTimeSteps <= 0 || numSpacePartitions <= 1)
            throw new System.ArgumentException("Cannot set-up BlackScholesFiniteDifferenceSolver - invalid parameters.");
        if (payoffFunction == null)
            throw new System.ArgumentException("Cannot set-up BlackScholesFiniteDifferenceSolver - payoff function musnt' be NULL.");
        
        T = maturity; this.payoff = payoffFunction; this.r = riskFreeRate; this.sigma = sigma; 
		this.R = R; N = (int)numTimeSteps; M = (int)numSpacePartitions; h = R / (M - 1);

		SetUpFiniteDifferenceMatrix();
		SetUpSolverMatrix();

		u = Vector<double>.Build.Dense (M);
	}

	private void SetUpFiniteDifferenceMatrix()
	{
		double rPlus = Math.Max (r, 0); 
		double rMinus = -Math.Min (r, 0);

		// i is the row index, j is the column index
		Func<int, int, double> matrixEntry = (i, j) =>
		{
            if (i == j && (i == 0 || i == M - 1))
                return 1.0;
            else if (i == j) // diagonal entry
            {
                double b = sigma * sigma * i * i + rPlus * i + rMinus * i + r;
                return b;
            }
            //else if (i == 1 && j == 0) // below diagonal, first column
            //    return 0.0;
            else if (i == M - 1 && j == i - 1) //below diagonal, last row
            {
                return 1.0;
            }

            else if (i > 0 && i < M - 1 && j == i + 1)
            {
                double a = -0.5 * sigma * sigma * i * i - rPlus * i;
                return a;
            }
            else if (i < M - 1 && j == i - 1)
            {
                double c = -0.5 * sigma * sigma * i * i - rMinus * i;
                return c;
            }
            else
                return 0;
		};
		A = Matrix<double>.Build.Sparse(M, M, matrixEntry);
//		Console.WriteLine (A.ToString ());
	}

	private void SetUpSolverMatrix()
	{
		double tau = T / N;
		Matrix<double> I = Matrix<double>.Build.SparseIdentity(M);
		//Console.WriteLine (A.ToString());
		S = I + tau*A;
		//Console.WriteLine("Set up matrix for solver");
	}

	private void SetUpSolver()
	{
		// Stop calculation if 1000 iterations reached during calculation
		IterationCountStopCriterion<double> iterationCountStopCriterion = new IterationCountStopCriterion<double>(20000);

		// Stop calculation if residuals are below 1E-? --> the calculation is considered converged
		ResidualStopCriterion<double> residualStopCriterion = new ResidualStopCriterion<double>(1e-3);

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
		for (int m = 0; m < M; ++m)
			u0[m] = payoff( m * h);

		return u0;
	}


	public Vector<double> Solve()
	{
		double tau = T / N;
		SetUpSolver();

		Vector<double> uOld = ApproxInitialCondition();
		Vector<double> uNew = Vector<double>.Build.Dense(M);

		double gamma = uOld [M - 1] - uOld [M - 2]; // for neumann bdry at S=R
		for (int k = 0; k < N; k++)
		{
            uOld[0] = (1 + tau) * Math.Exp(-r * tau * k) * payoff(0); // enforce boundary at 0
            uOld[M - 1] = gamma; // enforce boundary

			// Must solve ( I + tau * A ) * uNew =  uOld i.e. S * uNew = uOld
			uNew = S.SolveIterative(uOld, solver, monitor);

            if (monitor.Status != IterationStatus.Converged)
                throw new System.ArithmeticException("Iterative solver failed to converge at time step " + k);
            

			uOld = uNew;
		}
		Console.WriteLine();
		u = uNew;
		solved = true;
		return u;
	}
	public double Price(double S)
	{
		if (S + h > R)
			throw new Exception ("S too large for solver domain.");
		if (S <= 0)
			throw new Exception ("S must be > 0.");
		if (!solved)
			Solve ();


		int indexBelow = (int)Math.Floor (S / h);
		int indexAbove = (int)Math.Ceiling (S / h); 

		return u [indexBelow] + (S-indexBelow*h)*(u [indexAbove] - u [indexBelow]) / h; 
	}
}

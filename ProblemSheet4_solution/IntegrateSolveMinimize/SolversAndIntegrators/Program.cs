using System;
using MathNet.Numerics.LinearAlgebra;

using MathNet.Numerics.LinearAlgebra.Solvers;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;


namespace SolversAndIntegrators
{
	/// <summary>
    /// Test the integrator, solver and minimizer code.
    /// </summary>
    class MainClass
	{
		public static void Main (string[] args)
		{
			Func<double, double> normalDensity
                = (x) => Math.Exp(-x * x / 2.0) / Math.Sqrt(2 * Math.PI);
            CompositeIntegrator integrator = new CompositeIntegrator(4);
            double integral = integrator.Integrate(normalDensity, -4, 4, 100);
            Console.WriteLine("Integral is approximately {0}", integral);

            NewtonSolver newtonSolver = new NewtonSolver(1e-16, 100);

			// we want a function F(x,y)_1 = x^2 + y^2 - 2xy, F(x,y)_2 = x^2 - y^2
			// you can check that F(x,x)_1 = 0 and F(x,x)_2 = 0 i.e. uncountably many solutions
			Func<Vector<double>, Vector<double>> F = (x) => {
				Vector<double> y = Vector<double>.Build.Dense(x.Count);
				y[0] = x[0]*x[0] + x[1]*x[1] - 2*x[0]*x[1];
				y[1] = x[0]*x[0] - x[1]*x[1];
				return y;
			};

			Func<Vector<double>,Matrix<double>> J_F = (x) => {
				int d = x.Count;
				Matrix<double> J_F_vals = Matrix<double>.Build.Dense (d, d);
				J_F_vals[0,0] = 2*x[0] -2*x[1]; J_F_vals[0,1] = 2*x[1]-2*x[0];
				J_F_vals[1,0] = 2*x[0]; J_F_vals[1,1] = -2*x[1];
				return J_F_vals;
			};

			Vector<double> x0 = Vector<double>.Build.Dense (2);
			x0[0] = 1; x0[1] = -1;
            Vector<double> x1;
            x1 = newtonSolver.Solve(F,J_F,x0);
			Console.WriteLine (x1.ToString ());



			// we want a function F(x,y)_1 = x^2 + y^2 - 2xy - 1, F(x,y)_2 = x^2 - y^2 - 7
			// you can check that [F(4,3), F(4,3)] = [0,0] 
			Func<Vector<double>, Vector<double>> F2 = (x) => {
				Vector<double> y = Vector<double>.Build.Dense(x.Count);
				y[0] = x[0]*x[0] + x[1]*x[1] - 2*x[0]*x[1] - 1;
				y[1] = x[0]*x[0] - x[1]*x[1] - 7;
				return y;
			};

			Func<Vector<double>,Matrix<double>> J_F2 = (x) => {
				int d = x.Count;
				Matrix<double> J_F_vals = Matrix<double>.Build.Dense (d, d);
				J_F_vals[0,0] = 2*x[0] -2*x[1]; J_F_vals[0,1] = 2*x[1]-2*x[0];
				J_F_vals[1,0] = 2*x[0]; J_F_vals[1,1] = -2*x[1];
				return J_F_vals;
			};


			x1 = newtonSolver.Solve(F2,J_F2,x0);
			Console.WriteLine (x1.ToString ());

			NewtonMnimizer minimizer = new NewtonMnimizer (1e-5, 100);
			Func<Vector<double>, double> f = (x) => x [0] * x [0] + x [1] * x [1];

			Vector<double> startPt = Vector<double>.Build.Dense (2);
			startPt[0] = 1; startPt[1] = -1;

			Console.WriteLine ("With approximate grad of f and hessian of f using f.d.");
			Console.WriteLine(minimizer.Minimize(f, null, null, startPt));

			Console.WriteLine ("With exact grad_f and approximate hessian f using f.d.");
			Func<Vector<double>, Vector<double>> grad_f = (x) => {
				Vector<double> grad = Vector<double>.Build.Dense(x.Count);
				grad[0] = 2*x[0]; grad[1] = 2*x[1];
				return grad;
			};
			Console.WriteLine(minimizer.Minimize(f, grad_f, null, startPt));

			Console.WriteLine ("With exact grad of f and hessian of f.");
			Func<Vector<double>, Matrix<double>> hessian_f = (x) => {
				Matrix<double> hessian = Matrix<double>.Build.Dense(x.Count,x.Count);
				hessian[0,0] = 2; hessian[0,1] = 0;
				hessian[1,0] = 0; hessian[1,1] = 2;
				return hessian;
			};
			Console.WriteLine(minimizer.Minimize(f, grad_f, hessian_f, startPt));
            Console.ReadKey ();
		}
	}
}

using System;
using MathNet.Numerics.LinearAlgebra;

using MathNet.Numerics.LinearAlgebra.Solvers;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;


namespace SolversAndIntegrators
{
    /// <summary>
    ///  NewtonMinimizer finds a local mimimum of a given function.
    /// </summary>
    public class NewtonMnimizer
    {
        private const double delta = 1e-7; // for approximating grad
        private NewtonSolver solver;

        /// <summary>
        /// Initializes a new instance of the <see cref="T:NewtonSolverMultiDim.NewtonMnimizer"/> class.
        /// </summary>
        /// <param name="tol">Tol.</param>
        /// <param name="maxIt">Max it.</param>
        public NewtonMnimizer(double tol, int maxIt)
        {
            solver = new NewtonSolver(tol, maxIt);
        }

        /// <summary>
        /// Approximates the gradient of a given function at a given point.
        /// </summary>
        /// <returns>Approximation of the gradient at the point.</returns>
        /// <param name="f">Function f:R^d -> R.</param>
        /// <param name="x">The point x at which the gradient is to be approximated.</param>
        private Vector<double> ApproximateGrad(Func<Vector<double>, double> f, Vector<double> x)
        {
            int d = x.Count;
            Vector<double> grad = Vector<double>.Build.Dense(d);
            for (int i = 0; i < d; ++i)
            {
                Vector<double> xPlus = Vector<double>.Build.DenseOfVector(x);
                xPlus[i] = x[i] + delta;

                Vector<double> xMinus = Vector<double>.Build.DenseOfVector(x);
                xMinus[i] = x[i] - delta;

                double partial_f_parital_x_i = (f(xPlus) - f(xMinus)) / (2 * delta);
                grad[i] = partial_f_parital_x_i;

            }
            return grad;
        }

        /// <summary>
        /// Approximates the Hessian matrix of a given function at a given point.
        /// </summary>
        /// <returns>Approximation of the Hessian matrix.</returns>
        /// <param name="f">Function f:R^d -> R.</param>
        /// <param name="x">The point x at which the Hessian is to be approximated.</param>
        private Matrix<double> ApproximateHessian(Func<Vector<double>, double> f,
                                                    Vector<double> x)
        {
            int d = x.Count;
            Matrix<double> hessian = Matrix<double>.Build.Dense(d, d);

            for (int i = 0; i < d; ++i)
            {
                for (int j = 0; j < d; ++j)
                {

                    Vector<double> x_i_plus_j_plus = Vector<double>.Build.DenseOfVector(x);
                    x_i_plus_j_plus[i] = x[i] + delta;
                    x_i_plus_j_plus[j] = x_i_plus_j_plus[j] + delta;
                    double f_plus_plus = f(x_i_plus_j_plus);

                    Vector<double> x_i_plus_j_minus = Vector<double>.Build.DenseOfVector(x);
                    x_i_plus_j_minus[i] = x[i] + delta;
                    x_i_plus_j_minus[j] = x_i_plus_j_minus[j] - delta;
                    double f_plus_minus = f(x_i_plus_j_minus);

                    Vector<double> x_i_minus_j_plus = Vector<double>.Build.DenseOfVector(x);
                    x_i_minus_j_plus[i] = x[i] - delta;
                    x_i_minus_j_plus[j] = x_i_minus_j_plus[j] + delta;
                    double f_minus_plus = f(x_i_minus_j_plus);

                    Vector<double> x_i_minus_j_minus = Vector<double>.Build.DenseOfVector(x);
                    x_i_minus_j_minus[i] = x[i] - delta;
                    x_i_minus_j_minus[j] = x_i_minus_j_minus[j] - delta;
                    double f_minus_minus = f(x_i_minus_j_minus);

                    double partial_f_sq_parital_x_i_partial_x_j
                    = (1 / (2 * delta)) * ((1 / (2 * delta)) * (f_plus_plus - f_plus_minus) - (1 / (2 * delta)) * (f_minus_plus - f_minus_minus));

                    hessian[i, j] = partial_f_sq_parital_x_i_partial_x_j;
                }
            }

            return hessian;
        }


        /// <summary>
        /// Try to find a (local) minimum of a given function from a given starting point. 
        /// </summary>
        /// <returns>Approximate local minimum.</returns>
        /// <param name="f">Function f:R^d -> R to be minimized.</param>
        /// <param name="grad_f">Exact gradient of f.</param>
        /// <param name="hessian_f">Exact Hessian of f.</param>
        /// <param name="x_0">X 0.</param>
        public Vector<double> Minimize(Func<Vector<double>, double> f,
            Func<Vector<double>, Vector<double>> grad_f,
            Func<Vector<double>, Matrix<double>> hessian_f, Vector<double> x_0)
        {
            Vector<double> x = Vector<double>.Build.Dense(x_0.Count);
            // need to create the functions for the solver
            if (grad_f != null)
            {
                x = solver.Solve(grad_f, hessian_f, x_0);
            }
            else {
                Func<Vector<double>, Vector<double>> grad_approx = (point) => ApproximateGrad(f, point);
                Func<Vector<double>, Matrix<double>> hessian_approx = (point) => ApproximateHessian(f, point);
                x = solver.Solve(grad_approx, hessian_approx, x_0);
            }
            return x;
        }
    }
}

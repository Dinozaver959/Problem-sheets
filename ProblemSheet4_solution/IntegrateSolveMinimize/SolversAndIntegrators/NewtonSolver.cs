using System;
using MathNet.Numerics.LinearAlgebra;

using MathNet.Numerics.LinearAlgebra.Solvers;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;


namespace SolversAndIntegrators
{

    /// <summary>
    /// Our own exception class to let the user know in case 
    /// the solver fails to converge
    /// </summary>
    public class NewtonSolverFailedException : Exception
    {
        public NewtonSolverFailedException()
        {
        }
        public NewtonSolverFailedException(string message) : base(message)
        {
        }
    }


    /// <summary>
    /// A multidimensional Newton solver for solving nonlinear systems of the 
    /// form F(x) = 0 given in an initial guess. 
    /// </summary>
    public class NewtonSolver
    {
        private const double delta = 1e-10; // for approximating partial derivatives

        private double tol;
        private int maxIt;

        /// <summary>
        /// Initializes a new instance of the <see cref="T:NewtonSolverMultiDim.NewtonSolver"/> class.
        /// </summary>
        /// <param name="tolerance">Tolerance.</param>
        /// <param name="maximumIterations">Maximum iterations.</param>
        public NewtonSolver(double tolerance, int maximumIterations)
        {
            tol = tolerance;
            maxIt = maximumIterations;
        }

        /// <summary>
        /// Calculates an approximation for Jacobian matrix for a function F:R^d \to R^d at the point x. 
        /// </summary>
        /// <returns>The approximation.</returns>
        /// <param name="F">F.</param>
        /// <param name="x">The x coordinate.</param>
        public Matrix<double> ApproximateJacobian(Func<Vector<double>, Vector<double>> F, Vector<double> x)
        {
            int d = x.Count;
            Matrix<double> J_F = Matrix<double>.Build.Dense(d, d);

            for (int j = 0; j < d; ++j)
            {
                Vector<double> xPlus = Vector<double>.Build.DenseOfVector(x);
                xPlus[j] = x[j] + delta;
                Vector<double> xMinus = Vector<double>.Build.DenseOfVector(x);
                xMinus[j] = x[j] - delta;
                Vector<double> partial_F_parital_x_j = ((F(xPlus) - F(xMinus)) / (2 * delta));
                J_F.SetColumn(j, partial_F_parital_x_j);
            }


            return J_F;
        }

        /// <summary>
        /// Approximates solution to a nonlinear equation F(x) = 0.
        /// Stops if the l_2 error of the function evaluated at an approximation point is smaller than tolerance.
        /// Throws an exception if either the determinant of the Jacobian is very small 
        /// for any point x in the approximating sequence
        /// or if the maximum number of iterations is reached.
        /// </summary>
        /// <returns>The approximation solution to F(x) = 0.</returns>
        /// <param name="F">The function F:R^d -> R^d.</param>
        /// <param name="J_F">Analytiaclly computed Jacobian matrix. May be null in which case approximation is used.</param>
        /// <param name="x_0">The initial guess.</param>
        public Vector<double> Solve(Func<Vector<double>, Vector<double>> F, Func<Vector<double>, Matrix<double>> J_F, Vector<double> x_0)
        {
            Vector<double> F_of_x_0;
            Matrix<double> J_F_of_x_0;
            Vector<double> error; // this will effectively be x_1 - x_0
            int iterationsCount = 0;
            while ((F_of_x_0 = F(x_0)).L2Norm() > tol)
            {
                ++iterationsCount;

                if (J_F != null)
                    J_F_of_x_0 = J_F(x_0);
                else
                    J_F_of_x_0 = ApproximateJacobian(F, x_0);


                F_of_x_0 = F(x_0);
                if (Math.Abs(J_F_of_x_0.Determinant()) < 1e-25)
                {
                    throw new NewtonSolverFailedException("Jacobian matrix has almost zero determinant."
                        + "matrix: \n" + J_F_of_x_0.ToString());
                }
                error = J_F_of_x_0.Solve(F_of_x_0);
                x_0 = error - x_0; // this is now the new approximation 

                if (iterationsCount > maxIt)
                    throw new NewtonSolverFailedException("Solver failed to converge.");
            }

            return x_0;
        }
    }
}

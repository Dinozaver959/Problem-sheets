using System;

namespace SolversAndIntegrators
{
    // Note that this is beyond what was asked in the problem sheet.


    /// <summary>
    /// An adaptive extension of the <see cref="T:CompositeIntegrator"/> class.
    /// </summary>
    public class AdaptiveCompositeIntegrator : CompositeIntegrator
    {
        private const double MaxPrecision = 1e-14;
        private const int defaultN = 10;

        public AdaptiveCompositeIntegrator() : base() { }
        public AdaptiveCompositeIntegrator(AdaptiveCompositeIntegrator integrator) : base(integrator) { }
        public AdaptiveCompositeIntegrator(int newtonCotesOrder) : base(newtonCotesOrder) { }


        /// <summary>
        /// Approximate the integral \int_a^b f(x) dx to given tol using a recursive adaptive algorithm.
        /// </summary>
        /// <param name="f">The integrand.</param>
        /// <param name="a">The lower bound.</param>
        /// <param name="b">The upper bound.</param>
        /// <param name="tol">The desired tolerance.</param>
        public double Integrate(Func<double, double> f, double a, double b, double tol)
        {
            double integralRough = base.Integrate(f, a, b, defaultN);
            double integralBetter = base.Integrate(f, a, b, 2 * defaultN);
            double err = Math.Abs(integralRough - integralBetter);
            if (err <= tol || tol < MaxPrecision)
                return integralBetter;
            else {
                double mid = a + (b - a) / 2.0;
                double leftInt = Integrate(f, a, mid, tol / 2);
                double rightInt = Integrate(f, mid, b, tol / 2);
            return leftInt + rightInt;
        }
    }

}


}
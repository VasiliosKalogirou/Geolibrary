using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoLibrary
{
    public class LambertConformalConic
    {
        #region Members

        #region StaticMembers
        public static readonly LambertConformalConic ETRS89 = new LambertConformalConic(10.0, 52.0, 35.0, 65.0, 4000000.0, 2800000.0, Ellipsoid.GRS80);
        public static readonly LambertConformalConic NAD27 = new LambertConformalConic(27.83333333, -99.0, 28.383333333, 30.2833333333, 609601.2, 0.0, Ellipsoid.Clarke1866);

        #endregion


        /// <summary>
        /// The latitude of the true origin.
        /// </summary>
        private double phi0;

        /// <summary>
        /// The longitude of the true origin.
        /// </summary>
        private double lambda0;

        /// <summary>
        /// The sine of the latitude of the true origin.
        /// </summary>
        private double sinphi0;

        /// <summary>
        /// The latitude of the upper parallel.
        /// </summary>
        private double phiu;

        /// <summary>
        /// The latitude of the lower parallel.
        /// </summary>
        private double phil;

        /// <summary>
        /// The latitude of the (false) grid origin.
        /// </summary>
        private double phib;

        /// <summary>
        /// False eastings to be added to all coordinates, equivalent to the eastings at the true origin.
        /// </summary>
        private double FE;

        /// <summary>
        /// False northings to be added to all coordinates, equivalent to the northings at the true origin.
        /// </summary>
        private double FN;

        /// <summary>
        /// The geodetic system on which the projected system is based.
        /// </summary>
        private GeoLibrary.Ellipsoid ellipsoid = new GeoLibrary.Ellipsoid(1.0, 1.0);

        /// <summary>
        /// Grid scale factor calculated at a general point.
        /// </summary>
        private double k;

        /// <summary>
        /// The value of the convergence angle.
        /// </summary>
        private double gamma;

        /// <summary>
        /// The isometric latitude of the upper parallel.
        /// </summary>
        private double Qu;

        /// <summary>
        /// The isometric latitude of the lower parallel.
        /// </summary>
        private double Ql;

        /// <summary>
        /// The isometric latitude of the parallel of the grid origin.
        /// </summary>
        private double Qb;

        /// <summary>
        /// The Wu parameter used in the mapping equations.
        /// </summary>
        private double Wu;

        /// <summary>
        /// The parameter Wl used in the mapping equations.
        /// </summary>
        private double Wl;

        /// <summary>
        /// The mapping radius at the equator.
        /// </summary>
        private double K;

        /// <summary>
        /// The mapping radius at latitude phi0;
        /// </summary>
        private double R0;

        /// <summary>
        /// The sine of the latitude of the lower standard parallel.
        /// </summary>
        private double sinphil;

        /// <summary>
        /// The sine of the latitude of the upper standard parallel.
        /// </summary>
        private double sinphiu;

        /// <summary>
        /// The sine of the latitude of the parallel of the grid origin.
        /// </summary>
        private double sinphib;

        /// <summary>
        /// The flattening of the ellipdoid.
        /// </summary>
        private double f;

        /// <summary>
        /// The first eccentricity of the ellipsoid.
        /// </summary>
        private double ea;

        /// <summary>
        /// The square of the first eccentricity of the ellipsoid.
        /// </summary>
        private double ea2;
        #endregion

        #region Construction

        /// <summary>
        /// Creates a new instance of the Lambert Conformal Conic class.
        /// </summary>
        /// <param name="lambda0">The longitude of the true origin.</</param>
        /// <param name="phib">The latitude of the (false) grid origin.</param>
        /// <param name="phiu">The latitude of the upper parallel.</param>
        /// <param name="phil">The latitude of the lower parallel.</param>
        /// <param name="FE">False eastings to be added to all coordinates, equivalent to the eastings at the true origin.</param>
        /// <param name="FN">False northings to be added to all coordinates, equivalent to the northings at the true origin.</param>
        public LambertConformalConic(double lambda0, double phib, double phil, double phiu, double FE, double FN, GeoLibrary.Ellipsoid ellipsoid)
        {
            this.lambda0 = lambda0;
            this.phib = phib;
            this.phiu = phiu;
            this.phil = phil;
            this.FE = FE;
            this.FN = FN;
            this.ellipsoid = ellipsoid;

            this.f = 1 / ellipsoid.InverseFlattening;
            this.ea = Math.Sqrt((2 * f) - (f * f));
            this.ea2 = ea * ea;
            this.sinphil = Math.Sin(phil);
            this.sinphiu = Math.Sin(phiu);
            this.sinphib = Math.Sin(phib);

            //Calculations on "Q" Parameters
            this.Qu = 0.5 * (Math.Log((1 + sinphiu) / (1 - sinphiu)) - ea * Math.Log((1 + ea * sinphiu) / (1 - ea * sinphiu)));
            this.Ql = 0.5 * (Math.Log((1 + sinphil) / (1 - sinphil)) - ea * Math.Log((1 + ea * sinphil) / (1 - ea * sinphil)));
            this.Qb = 0.5 * (Math.Log((1 + sinphib) / (1 - sinphib)) - ea * Math.Log((1 + ea * sinphib) / (1 - ea * sinphib)));

            //Calculations on "W" Parameters
            this.Wu = Math.Sqrt(1 - ea2 * Math.Pow(sinphiu, 2));
            this.Wl = Math.Sqrt(1 - ea2 * Math.Pow(sinphil, 2));

            double cosphil = Math.Cos(phil);
            double cosphiu = Math.Cos(phiu);
            double a = ellipsoid.SemimajorAxis;

            this.sinphi0 = (Math.Log((Wu * cosphil) / (Wl * cosphiu))) / (Qu - Ql);
            this.phi0 = Math.Asin(sinphi0);
            this.K = a * cosphiu * Math.Exp(Qu * sinphi0) / (Wu * sinphi0);
            this.R0 = K / (Math.Exp(Qb * sinphi0));
        }
        #endregion

        #region Properties

        public double LambdaOrigin
        {
            get { return this.lambda0; }
        }

        public double SinePhiOrigin
        {
            get { return this.sinphi0; }
        }
        public double UpperStandardParallel
        {
            get { return this.phiu; }
        }
        public double LowerStandardParallel
        {
            get { return this.phil; }
        }
        public double FalseEasting
        {
            get { return this.FE; }
        }
        public double FalseNorthing
        {
            get { return this.FN; }
        }

        public double Kappa
        {
            get { return this.K; }
        }

        public double Gamma
        {
            get { return this.gamma; }
        }

        public double Rzero
        {
            get { return this.R0; }
        }

        #endregion

        #region Methods

        public void GeodeticToProjected(double phi, double lambda, out double E, out double N)
        {
            double sinphi = Math.Sin(phi);
            double Q = 0.5 * (Math.Log((1 + sinphi) / (1 - sinphi)) - ea * Math.Log((1 + ea * sinphi) / (1 - ea * sinphi)));
            double R = K / (Math.Exp(Q * sinphi0));
            gamma = (lambda - lambda0) * sinphi0;

            //Calculation of E, N, k.
            E = FE + R * Math.Sin(gamma);
            N = R0 + FN - R * Math.Cos(gamma);
            k = Math.Sqrt(1 - ea2 * sinphi * sinphi);

            E = Math.Round(E, 3);
            N = Math.Round(N, 3);
        }

        public void ProjectedToGeodetic(double E, double N, out double phi, out double lambda)
        {
            double Rb = R0 - N + FN;
            double Rb2 = Rb * Rb;
            double Eb = FE - E;
            double Eb2 = Eb * Eb;
            gamma = Math.Atan(Eb / Rb);

            lambda = lambda0 - (gamma / sinphi0);

            double R = Math.Sqrt(Rb2 + Eb2);
            double Q = (Math.Log(K / R)) / sinphi0;

            double sinphitemp;
            double sinphi = (Math.Exp(2 * Q) - 1) / (Math.Exp(2 * Q) + 1);
            do
            {
                sinphitemp = sinphi;
                double f1 = (Math.Log((1 + sinphi) / (1 - sinphi)) - ea * Math.Log((1 + ea * sinphi) / (1 - ea * sinphi))) / 2 - Q;
                double f2 = 1 / (1 - sinphi * sinphi) - ea2 / (1 - ea2 * sinphi * sinphi);
                sinphi = sinphi - (f1 / f2);

            } while (Math.Abs(sinphi - sinphitemp) > 1e-9);
            phi = Math.Asin(sinphi);

            phi *= 180 / Math.PI;
            lambda *= 180 / Math.PI;

            phi = Math.Round(phi, 6);
            lambda = Math.Round(lambda, 6);
        }

        #endregion
    }
}

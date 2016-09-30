using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoLibrary
{
    public class LambertConformalConic_2
    {
        #region Members

        #region StaticMembers

        public static readonly LambertConformalConic_2 ETRS89 = new LambertConformalConic_2(10.0, 52.0, 35.0, 65.0, 4000000.0, 2800000.0, Ellipsoid.GRS80);
        public static readonly LambertConformalConic_2 NAD27 = new LambertConformalConic_2(27.83333333, -99.0, 28.383333333, 30.2833333333, 609601.2, 0.0, Ellipsoid.Clarke1866);

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
        /// The latitude of the upper parallel.
        /// </summary>
        private double phi1;

        /// <summary>
        /// The latitude of the lower parallel.
        /// </summary>
        private double phi2;

        /// <summary>
        /// False eastings to be added to all coordinates, equivalent to the eastings at the true origin.
        /// </summary>
        private double FE;

        /// <summary>
        /// False northings to be added to all coordinates, equivalent to the northings at the true origin.
        /// </summary>
        private double FN;

        /// <summary>
        /// The semi-major axis of the ellipsoid.
        /// </summary>
        private double a;

        /// <summary>
        /// The inverse of the flattening constant of the ellipsoid.
        /// </summary>
        private double invf;

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

        /// <summary>
        /// The parameter m1 used in calculations.
        /// </summary>
        private double m1;

        /// <summary>
        /// The parameter m2 used in calculations.
        /// </summary>
        private double m2;

        /// <summary>
        /// The parameter t0 used in calculations.
        /// </summary>
        private double t0;

        /// <summary>
        /// The parameter t1 used in calculations.
        /// </summary>
        private double t1;

        /// <summary>
        /// The parameter t2 used in calculations.
        /// </summary>
        private double t2;

        /// <summary>
        /// The parameter n used in calculations.
        /// </summary>
        private double n;

        /// <summary>
        /// The parameter F used in calculations.
        /// </summary>
        private double F;
        /// <summary>
        /// The rectifying rotation of the map grid for oblique projections.
        /// </summary>
        private double gamma;

        /// <summary>
        /// The parameter R0 used in calculations.
        /// </summary>
        private double R0;

        /// <summary>
        /// The geodetic system on which the projected system is based.
        /// </summary>
        private GeoLibrary.Ellipsoid ellipsoid = new GeoLibrary.Ellipsoid(1.0, 1.0);

        #endregion

        #region Construction

        /// <summary>
        /// Creates a new instance of the LambertConformalConic_2 class.
        /// </summary>
        /// <param name="phi0">The latitude of the true origin.</param>
        /// <param name="lambda0">The longitude of the true origin.</param>
        /// <param name="phi1">The latitude of the upper parallel.</param>
        /// <param name="phi2">The latitude of the lower parallel.</param>
        /// <param name="FE">False eastings to be added to all coordinates, equivalent to the eastings at the true origin.</param>
        /// <param name="FN">False northings to be added to all coordinates, equivalent to the northings at the true origin.</param>
        /// <param name="ellipsoid">The geodetic system on which the projected system is based.</param>
        public LambertConformalConic_2(double phi0, double lambda0, double phi1, double phi2, double FE, double FN, GeoLibrary.Ellipsoid ellipsoid)
        {
            this.phi0 = phi0 * Math.PI / 180;
            this.lambda0 = lambda0 * Math.PI / 180;
            this.phi1 = phi1 * Math.PI / 180;
            this.phi2 = phi2 * Math.PI / 180;
            this.FE = FE;
            this.FN = FN;
            this.ellipsoid = ellipsoid;

            this.a = ellipsoid.SemimajorAxis;
            this.invf = ellipsoid.InverseFlattening;
            this.f = 1 / invf;
            ea2 = 2 * f - f * f;
            ea = Math.Sqrt(ea2);

            double cosphi1 = Math.Cos(phi1);
            double cosphi2 = Math.Cos(phi2);
            double sinphi1 = Math.Sin(phi1);
            double sinphi2 = Math.Sin(phi2);
            double sinphi0 = Math.Sin(phi0);
            this.m1 = (cosphi1) / Math.Sqrt(1 - ea2 * sinphi1 * sinphi1);
            this.m2 = (cosphi2) / Math.Sqrt(1 - ea2 * sinphi2 * sinphi2);
            this.t1 = Math.Tan(Math.PI / 4 - phi1 / 2) / Math.Pow((1 - ea * sinphi1) / (1 + ea * sinphi1), ea / 2);
            this.t2 = Math.Tan(Math.PI / 4 - phi2 / 2) / Math.Pow((1 - ea * sinphi2) / (1 + ea * sinphi2), ea / 2);
            this.t0 = Math.Tan(Math.PI / 4 - phi0 / 2) / Math.Pow((1 - ea * sinphi0) / (1 + ea * sinphi0), ea / 2);
            this.n = (Math.Log(m1) - Math.Log(m2)) / (Math.Log(t1) - Math.Log(t2));
            this.F = m1 / (n * Math.Pow(t1, n));
            this.R0 = a * F * Math.Pow(t0, n);
        }

        #endregion

        #region Properties

        public double PhiOrigin
        {
            get { return this.phi0; }
        }

        public double LambdaOrigin
        {
            get { return this.lambda0; }
        }

        public double PhiUpper
        {
            get { return this.phi1; }
        }

        public double PhiLower
        {
            get { return this.phi2; }
        }

        public double FalseEasting
        {
            get { return this.FE; }
        }

        public double FalseNorthing
        {
            get { return this.FN; }
        }

        public double SMAxis
        {
            get { return this.a; }
        }

        public double FParameter
        {
            get { return this.F; }
        }

        public double nParameter
        {
            get { return this.n; }
        }

        public double FirstEccentricity
        {
            get { return this.ea; }
        }

        #endregion

        #region Methods

        public void GeodeticToProjected(double phi, double lambda, out double E, out double N)
        {
            double sinphi = Math.Sin(phi);
            double t = Math.Tan(Math.PI / 4 - phi / 2) / Math.Pow((1 - ea * sinphi) / (1 + ea * sinphi), ea / 2);
            double R = a * F * Math.Pow(t, n);
            gamma = n * (lambda - lambda0);
            N = FN + R0 - R * Math.Cos(gamma);
            E = FE + R * Math.Sin(gamma);
            E = Math.Round(E, 3);
            N = Math.Round(N, 3);

        }

        public void ProjectedToGeodetic(double E, double N, out double phi, out double lambda)
        {
            double Rb;
            if (n > 0)
            {
                Rb = +Math.Pow(Math.Pow((E - FE), 2) + Math.Pow((R0 - (N - FN)), 2), 0.5);
            }
            else
            {
                Rb = -Math.Pow(Math.Pow((E - FE), 2) + Math.Pow((R0 - (N - FN)), 2), 0.5);

            }
            double tb = Math.Pow(Rb / (a * F), 1 / n);
            double gammab = Math.Atan((E - FE) / (R0 - (N - FN)));
            double phitemp;
            phi = Math.PI / 2 - 2 * Math.Atan(tb);

            do
            {
                phitemp = phi;
                double sinphi = Math.Sin(phi);
                phi = (Math.PI / 2) - 2 * Math.Atan(tb * Math.Pow(((1 - ea * sinphi) / (1 + ea * sinphi)), ea / 2));

            } while (Math.Abs(phi - phitemp) > 1e-9);

            lambda = (gammab / n) + lambda0;

            phi = phi * 180 / Math.PI;
            lambda = lambda * 180 / Math.PI;

            phi = Math.Round(phi, 6);
            lambda = Math.Round(lambda, 6);
        }

        #endregion
    }
}

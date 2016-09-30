using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoLibrary
{
    public class TransverseMercator
    {
        #region Members

        #region StaticMembers

        public static readonly TransverseMercator TM07 = new TransverseMercator(0.0, 24.0, 500000.0, -2000000.0, 0.9996);
        public static readonly TransverseMercator TM87 = new TransverseMercator(0.0, 24.0, 500000.0, 0.0, 0.9996);
        public static readonly TransverseMercator LTM93 = new TransverseMercator(0.0, 33.0, 200000.0, -3500000.0, 0.99995);

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
        /// False eastings to be added to all coordinates, equivalent to the eastings at the true origin.
        /// </summary>
        private double FE;

        /// <summary>
        /// False northings to be added to all coordinates, equivalent to the northings at the true origin.
        /// </summary>
        private double FN;

        /// <summary>
        /// Overall scaling factor to be applied. For the TM it may be referred to as the central meridian scale factor. 
        /// </summary>
        private double k;

        /// <summary>
        /// The geodetic system on which the projected system is based.
        /// </summary>
        private GeoLibrary.Ellipsoid ellipsoid = new GeoLibrary.Ellipsoid(1.0, 1.0);

        #endregion

        #region Construction
        /// <summary>
        /// Creates a new instance of the Transverse Mercator class.
        /// </summary>
        /// <param name="phi0">The latitude of the true origin.</param>
        /// <param name="lambda0">The longitude of the true origin.</param>
        /// <param name="FE">False eastings to be added to all coordinates, equivalent to the eastings at the true origin.</param>
        /// <param name="FN">False northings to be added to all coordinates, equivalent to the northings at the true origin.</param>
        /// <param name="k">Overall scaling factor to be applied. For the TM it may be referred to as the central meridian scale factor. </param>
        public TransverseMercator(double phi0, double lambda0, double FE, double FN, double k)
        {
            this.phi0 = phi0 * Math.PI / 180;
            this.lambda0 = lambda0 * Math.PI / 180;
            this.FE = FE;
            this.FN = FN;
            this.k = k;
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

        public double FalseEasting
        {
            get { return this.FE; }
        }

        public double FalseNorthing
        {
            get { return this.FN; }
        }

        public double ScaleFactor
        {
            get { return this.k; }
        }

        #endregion

        #region Methods
        public void GeodeticToProjected(double phi, double lambda, GeoLibrary.Ellipsoid ellipsoid, out double E, out double N)
        {
            // Trigonometric Powers
            double sinphi = Math.Sin(phi); // CD: Optimization
            double sinphi2 = sinphi * sinphi;
            double sinphi4 = sinphi2 * sinphi2;
            double cosphi = Math.Cos(phi); // CD: Optimization
            double cosphi2 = cosphi * cosphi; //
            double cosphi3 = cosphi2 * cosphi; //
            double cosphi5 = cosphi3 * cosphi2; //
            double cosphi7 = cosphi5 * cosphi2; //


            // Difference on longitude
            double DL = (lambda - lambda0);
            double DL2 = DL * DL;
            double DL3 = DL2 * DL;
            double DL4 = DL2 * DL2;
            double DL5 = DL2 * DL3;
            double DL6 = DL2 * DL4;
            double DL7 = DL2 * DL5;
            double DL8 = DL4 * DL4;

            // Calculations on flattening and eccentricity
            double invf = ellipsoid.InverseFlattening;
            double a = ellipsoid.SemimajorAxis;
            double f = 1 / invf;
            double ea = Math.Sqrt(2 * f - Math.Pow(f, 2));
            double ea2 = Math.Pow(ea, 2);
            double eb = Math.Sqrt(ea2 / (1 - ea2));
            double Ni = a / Math.Sqrt(1 - (ea2 * sinphi2));  // radius of curvature of the prime vertical           

            // t & h - Calculations.
            double t = Math.Tan(phi);
            double tsq = t * t;   //VK: tsquare = t*t
            double h2 = eb * eb * Math.Cos(phi) * Math.Cos(phi);
            double h4 = h2 * h2;
            double h6 = h2 * h4;
            double h8 = h4 * h4;
            //double Dl = (lamvda - lo);

            // Calculation of the MeridianArc - Sf.
            double Sf = ellipsoid.MeridianArc(phi);

            // "T" Parameters used to calculate N and E.
            double T1 = k * Sf;
            double T2 = (sinphi * cosphi) / 2;
            double T3 = ((sinphi * cosphi3) / 24) * (5 - tsq + (9 * h2) + (4 * h4));
            double T4 = ((sinphi * cosphi5) / 720) * (61 - (58 * tsq) + (tsq * tsq) + (270 * h2) - (330 * tsq * h2) + (445 * h4) + (324 * h6) - (680 * tsq * h4) + (88 * h8) - (660 * tsq * h6) - (192 * tsq * h8));
            double T5 = (sinphi * cosphi7 / 40320) * (1385 - (3111 * tsq) + (543 * tsq * tsq) - (tsq * tsq * tsq));
            double T6 = cosphi;
            double T7 = (cosphi3 / 6) * (1 - tsq + h2);
            double T8 = (cosphi5 / 120) * (5 - (18 * tsq) + (tsq * tsq) + (14 * h2) - (58 * tsq * h2) + (13 * h4) + (4 * h6) - (64 * tsq * h4) - (24 * tsq * h6));
            double T9 = (cosphi7 / 5040) * (61 - (479 * tsq) + (179 * tsq * tsq) - (tsq * tsq * tsq));

            // Calculation of Easting & Northing.
            N = FN + T1 + (k * Ni) * ((T2 * DL2) + (T3 * DL4) + (T4 * DL6) + (T5 * DL8));
            E = FE + (k * Ni) * ((T6 * DL) + (T7 * DL3) + (T8 * DL5) + (T9 * DL7));
            E = Math.Round(E, 3);
            N = Math.Round(N, 3);
        }

        public void ProjectedToGeodetic(double E, double N, GeoLibrary.Ellipsoid ellipsoid, out double phi, out double lambda)
        {
            // Calculations on flattening and eccentricity       
            double a = ellipsoid.SemimajorAxis;
            double invf = ellipsoid.InverseFlattening;
            double f = 1 / invf;
            double ea = Math.Sqrt(2 * f - Math.Pow(f, 2));
            double ea2 = Math.Pow(ea, 2);
            double eb = Math.Sqrt(ea2 / (1 - ea2));

            double Sfb = (N - this.FN) / this.k;
            double phib = ellipsoid.FootprintLatitude(Sfb);
            double sinphib = Math.Sin(phib);
            double cosphib = Math.Cos(phib);
            double sinphib2 = sinphib * sinphib;
            double Nb = a / Math.Sqrt(1 - (ea2 * sinphib2));  // radius of curvature of the prime vertical for the phib = aktina kamoylotitas tis protis kathetis tomis. 
            double Mb = (a * (1 - ea2)) / Math.Sqrt(1 - (ea2 * sinphib2)); // meridian radius of curvature = aktina kampylotitas tis mesimvrinis tomis.
            double Nb3 = Math.Pow(Nb, 3);
            double Nb5 = Math.Pow(Nb, 5);
            double Nb7 = Math.Pow(Nb, 7);

            double tb = Math.Tan(phib);
            double tb2 = tb * tb;
            double tb4 = tb2 * tb2;
            double tb6 = tb2 * tb4;
            double hb2 = eb * eb * Math.Cos(phib) * Math.Cos(phib);
            double hb4 = hb2 * hb2;
            double hb6 = hb4 * hb2;
            double hb8 = hb2 * hb6;
            double Eb = E - FE;
            double Eb2 = Eb * Eb;
            double Eb4 = Eb2 * Eb2;
            double Eb6 = Eb4 * Eb2;
            double Eb8 = Eb4 * Eb4;
            double Q = Eb / (k * Nb);
            double Q2 = Q * Q;
            double Q3 = Q2 * Q;
            double Q5 = Q3 * Q2;
            double Q7 = Q5 * Q2;
            double k2 = k * k;
            double k4 = k2 * k2;
            double k6 = k4 * k2;
            double k8 = k6 * k2;

            // "T" Parameters used to calculate phi and lambda. 
            double T10 = tb / (2 * k2 * Mb * Nb);
            double T11 = (tb) * (5 + (3 * tb * tb) + hb2 - (4 * hb4) - (9 * tb2 * hb2)) / (24 * k4 * Mb * Nb3);
            double T12 = (tb) * (61 + (90 * tb2) + (45 * tb4) + (46 * hb2) - (252 * tb2 * hb2) - (3 * hb4) + (100 * hb6) - (66 * tb2 * hb4) - (90 * tb4 * hb2) + (88 * hb8) + (225 * tb4 * hb4) + (84 * tb2 * hb6) - (192 * tb2 * hb8)) / (720 * k6 * Mb * Nb5);
            double T13 = (tb) * (1385 + (3633 * tb2) + (4095 * tb4) + (1575 * tb6)) / (40320 * k8 * Mb * Nb7);
            double T14 = 1 / cosphib;
            double T15 = (1 + (2 * tb2) + hb2) / (6 * cosphib);
            double T16 = (5 + (6 * hb2) + (28 * tb2) - (3 * hb4) + (8 * tb2 * hb2) + (24 * tb4) - (4 * hb6) + (4 * tb2 * hb4) + (24 * tb2 * hb6)) / (120 * cosphib);
            double T17 = (61 + (66 * tb2) + (1320 * tb4) + (720 * tb6)) / (5040);

            // Phi and Lambda Calculations.
            phi = phib - (T10 * Eb2) + (T11 * Eb4) - (T12 * Eb6) + (T13 * Eb8);
            lambda = lambda0 + (T14 * Q) - (T15 * Q3) + (T16 * Q5) - (T17 * Q7);

            phi *= 180 / Math.PI;
            lambda *= 180 / Math.PI;

            phi = Math.Round(phi, 6);
            lambda = Math.Round(lambda, 6);
        }
        #endregion
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoLibrary
{
    public class Ellipsoid
    {
        #region Members

        #region Static Members

        public static readonly Ellipsoid WGS84 = new Ellipsoid(6378137.0, 298.257223563);
        public static readonly Ellipsoid GRS80 = new Ellipsoid(6378137.0, 298.257222101);
        public static readonly Ellipsoid Clarke1858 = new Ellipsoid(6378293.645, 294.26068);
        public static readonly Ellipsoid Clarke1866 = new Ellipsoid(6378206.400, 294.97870);
        public static readonly Ellipsoid WGS72 = new Ellipsoid(6378135.0, 298.2599997755319);
        public static readonly Ellipsoid Bessel1841 = new Ellipsoid(6377397.155, 299.15281535132334);

        #endregion

        /// <summary>
        /// The semi-major axis of the ellipsoid.
        /// </summary>
        private readonly double a;

        /// <summary>
        /// The semi-minor axis of the ellipsoid.
        /// </summary>
        private readonly double b;

        /// <summary>
        /// The inverse of the flattening constant of the ellipsoid.
        /// </summary>
        private readonly double invf;

        /// <summary>
        /// The flattening of the ellipdoid.
        /// </summary>
        private readonly double f;

        /// <summary>
        /// The first eccentricity of the ellipsoid.
        /// </summary>
        private readonly double ea;

        /// <summary>
        /// The square of the first eccentricity of the ellipsoid.
        /// </summary>
        private readonly double ea2;

        /// <summary>
        /// The first eccentricity of the ellipsoid to the fourth.
        /// </summary>
        private readonly double ea4;

        /// <summary>
        /// The first eccentricity of the ellipsoid to the sixth.
        /// </summary>
        private readonly double ea6;

        /// <summary>
        /// The first eccentricity of the ellipsoid to the eighth.
        /// </summary>
        private readonly double ea8;

        /// <summary>
        /// The second eccentricity of the ellipsoid.
        /// </summary>
        private readonly double eb;

        /// <summary>
        /// The square of the second eccentricity of the ellipsoid.
        /// </summary>
        private readonly double eb2;

        /// <summary>
        /// The Meridian Arc constant A0. It is depended only on the first eccentricity ea.
        /// </summary>
        private readonly double a0;

        /// <summary>
        /// The Meridian Arc constant A2. It is depended only on the first eccentricity ea.
        /// </summary>
        private readonly double a2;

        /// <summary>
        /// The Meridian Arc constant A4. It is depended only on the first eccentricity ea.
        /// </summary>
        private readonly double a4;

        /// <summary>
        /// The Meridian Arc constant A6. It is depended only on the first eccentricity ea.
        /// </summary>
        private readonly double a6;

        /// <summary>
        /// The Meridian Arc constant A8. It is depended only on the first eccentricity ea.
        /// </summary>
        private readonly double a8;

        #endregion

        #region Construction

        /// <summary>
        /// Creates a new instance of the Ellipsoid class.
        /// </summary>
        /// <param name="a">The semi-major axis of the ellipsoid.</param>
        /// <param name="invf">The inverse of the flattening constant.</param>
        public Ellipsoid(double a, double invf)
        {
            this.a = a;
            this.invf = invf;
            this.f = 1 / invf;
            this.b = a - (a * f);

            // Calculations on flattening and eccentricity.
            this.ea = Math.Sqrt((2 * f) - (f * f));
            this.ea2 = ea * ea;
            this.eb2 = ((a * a) - (b * b)) / (b * b);
            this.ea4 = ea2 * ea2;         //Math.Pow(e,4);
            this.ea6 = ea4 * ea2;         //Math.Pow(e,6);
            this.ea8 = ea6 * ea2;         // Math.Pow(e, 8);

            // "A" Parameters used in calculating the arc of a meridian Sf.
            this.a0 = 1.0 - (1.0 / 4.0) * ea2 - (3.0 / 64.0) * ea4 - (5.0 / 256.0) * ea6 - (175.0 / 16384.0) * ea8; // CD: Parentheses
            this.a2 = (3.0 / 8.0) * ea2 * (1.0 + ((1.0 / 4.0) * ea2) + ((15.0 / 128.0) * ea4) + ((35.0 / 512.0) * ea6));
            this.a4 = (15.0 / 256.0) * ea4 * (1.0 + ((3.0 / 4.0) * ea2) + ((35.0 / 64.0) * ea4));
            this.a6 = (35.0 / 3072.0) * ea6 * (1.0 + ((5.0 / 4.0) * ea2));
            this.a8 = (315.0 / 131072.0) * ea8;
        }

        #endregion

        #region Properties

        public double InverseFlattening
        {
            get { return this.invf; } 
        }

        public double SemimajorAxis
        {  
            get {return this.a; }
        }

        public double A0
        {
            get {return this.a0;}
        }

        public double A2
        {
            get { return this.a2; }
        }

        public double A4
        {
            get { return this.a4; }
        }

        public double A6
        {
            get { return this.a6; }
        }

        public double A8
        {
            get { return this.a8; }
        }
       
        #endregion

        #region Methods

        public double MeridianArc(double phi)
        {
            // Calculations on the sine of phi.
            double sinphi = Math.Sin(phi);
            double sin2phi = Math.Sin(2 * phi);
            double sin4phi = Math.Sin(4 * phi);
            double sin6phi = Math.Sin(6 * phi);          
            
            // VK double Sf = a * (a0*phi - a2*Math.Sin(2*phi) + a4*Math.Sin(4*phi) - a6 * Math.Sin(6*phi) + a8*Math.Sin(8*phi));
            double Sf = a * ((a0 * phi) - (a2 * sin2phi) + (a4 * sin4phi) - (a6 * sin6phi));
            return Sf;
        }
        
        public double FootprintLatitude(double s)
        {
            double phitemp;
            double phi = s / (this.a * this.a0 );
            do 
            {
                phitemp = phi;
                double sinphi = Math.Sin(phi);
                double sin2phi = Math.Sin(2 * phi);
                double sin4phi = Math.Sin(4 * phi);
                double sin6phi = Math.Sin(6 * phi);          

                phi = s/(this.a*this.a0) + this.a2/this.a0*sin2phi - this.a4/this.a0*sin4phi + this.a6/this.a0*sin6phi;

            } while (Math.Abs(phi - phitemp) > 1e-9);

            return phi;
        }
        
        public void GeodeticToGeocentric(double phi, double lambda, double h, out double X, out double Y, out double Z)
        {
            double sinphi = Math.Sin(phi);
            double sinphi2 = sinphi * sinphi;
            double Ni = a / (Math.Sqrt(1 - (ea2 * sinphi2)));

            // Geocentric coordinates calculations
            X = (Ni + h) * Math.Cos(phi) * Math.Cos(lambda);
            Y = (Ni + h) * Math.Cos(phi) * Math.Sin(lambda);
            Z = (((1 - ea2) * Ni) + h) * Math.Sin(phi);
            X = Math.Round(X, 3);
            Y = Math.Round(Y, 3);
            Z = Math.Round(Z, 3);
        }   

        public void GeocentricToGeodetic(double X, double Y, double Z, out double phi, out double lambda, out double h)
        {
            double P = Math.Sqrt((X * X) + (Y * Y));
            double sinphi;
            double sinphi2;
            double Ni;
            
            //Calculation of lambda
            lambda = Math.Atan2(Y, X);            
            
            //Calculation of phi
            phi = Math.Atan((Z*(1+eb2))/P);
            double phitemp = 0;
            while ((phi - phitemp) > 1.0e-9)
	        {
                sinphi = Math.Sin(phi);
                sinphi2 = sinphi * sinphi;
                Ni = a / (Math.Sqrt(1 - (ea2 * sinphi2))); 
                phitemp = phi;
                phi = Math.Atan( (Z +(ea2 * Ni * sinphi)) / P);              
            }
            
            //Calculation of h
            sinphi = Math.Sin(phi);
            sinphi2 = sinphi * sinphi;
            Ni = a / (Math.Sqrt(1 - (ea2 * sinphi2))); 
            h = (Z / sinphi) - ((1 - ea2) * Ni);
            
            //Conversion of lamvda and phi from radians to degrees
            lambda = lambda * 180 / Math.PI;
            phi = phi * 180 / Math.PI;

            phi = Math.Round(phi, 6);
            lambda = Math.Round(lambda, 6);
            h = Math.Round(h, 3);        
        }

        #endregion
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoLibrary
{
    public class GeoPoint
    {
        #region Members

        /// <summary>
        /// The X coordinate of the point.
        /// </summary>
        double X;

        /// <summary>
        /// The Y coordinate of the point.
        /// </summary>
        double Y;

        /// <summary>
        /// The Z coordinate of the point.
        /// </summary>
        double Z;

        /// <summary>
        /// The id number of the Coordinate Reference System.
        /// </summary>
        int CRSid;

        #endregion

        #region Construction

        /// <summary>
        /// Creates a new instance of a 3-dimensional point.
        /// </summary>
        /// <param name="X">The X coordinate of the point.</param>
        /// <param name="Y">The Y coordinate of the point.</param>
        /// <param name="Z">The Z coordinate of the point.</param>
        public GeoPoint(double X, double Y, double Z)
        {
            this.X = X;
            this.Y = Y;
            this.Z = Z;
        }

        /// <summary>
        /// Creates a new instance of a 2-dimensional point.
        /// </summary>
        /// <param name="X">The X coordinate of the point.</param>
        /// <param name="Y">The Y coordinate of the point.</param>
        public GeoPoint(double X, double Y)
        {
            this.X = X;
            this.Y = Y;
        }

        #endregion

        #region Properties
        public double A
        {
            get { return this.X; }
        }
        public double B
        {
            get { return this.Y; }
        }
        public double C
        {
            get { return this.Z; }
        }
        #endregion

        #region Methods

        public string ToWkt()
        {
            return string.Format("POINT({0} {1} {2})", X, Y, Z);
        }

        #endregion
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoLibrary
{
    public class BursaWolfTransformation
    {
        #region Members

        #region StaticMembers
        public static readonly BursaWolfTransformation TM07toTM87 = new BursaWolfTransformation(203.437, -73.461, -243.594, -0.170, -0.060, -0.151, -0.294);
        public static readonly BursaWolfTransformation TM87toTM07 = new BursaWolfTransformation(-203.437, 73.461, 243.594, 0.170, 0.060, 0.151, 0.294);
        public static readonly BursaWolfTransformation WGS72toWGS84 = new BursaWolfTransformation(0.0, 0.0, 4.5, 0.0, 0.0, -0.554, 0.219);
        #endregion

        private double tX;

        private double tY;

        private double tZ;

        private double rX;

        private double rY;

        private double rZ;

        private double dS;

        #endregion

        #region Construction

        public BursaWolfTransformation(double DX, double DY, double DZ, double ax, double ay, double az, double m)
        {
            this.tX = DX;
            this.tY = DY;
            this.tZ = DZ;
            this.rX = ax * Math.PI / 180 / 3600; //arcseconds to rad
            this.rY = ay * Math.PI / 180 / 3600; //arcseconds to rad
            this.rZ = az * Math.PI / 180 / 3600; //arcseconds to rad
            this.dS = m;
        }
        #endregion

        #region Properties

        public double DXtranslation
        {
            get { return this.tX; }
        }

        public double DYtranslation
        {
            get { return this.tY; }
        }

        public double DZtranslation
        {
            get { return this.tZ; }
        }

        public double xrotation
        {
            get { return this.rX; }
        }

        public double yrotation
        {
            get { return this.rY; }
        }

        public double zrotation
        {
            get { return this.rZ; }
        }

        public double scalefactor
        {
            get { return this.dS; }
        }

        #endregion

        #region Methods

        public void DirectBursaWolfTransformation(double Xsource, double Ysource, double Zsource, out double Xtarget, out double Ytarget, out double Ztarget)
        {
            double S = 1 + dS / 1000000;
            double[,] R = new double[3, 3] { { 1, rZ, -rY }, { -rZ, 1, -rX }, { rY, -rX, 1 } };
            Xtarget = S * (Xsource + Ysource * rZ + Zsource * (-rY)) + tX;
            Ytarget = S * (Xsource * (-rZ) + Ysource + Zsource * rX) + tY;
            Ztarget = S * (Xsource * rY + Ysource * (-rY) + Zsource) + tZ;

            Xtarget = Math.Round(Xtarget, 3);
            Ytarget = Math.Round(Ytarget, 3);
            Ztarget = Math.Round(Ztarget, 3);
        }
        #endregion
    }
}

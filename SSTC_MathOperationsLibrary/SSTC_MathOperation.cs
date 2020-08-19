using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;

namespace SSTC_MathOperationsLibrary
{
    /// <summary>
    /// Class that contains methods to calculate sag-tension related parameters.
    /// </summary>
    public static class SSTC_MathOperation
    {
        /// <summary>
        /// Calculates bare conductor weight per unit length - mC1g [N/m]
        /// </summary>
        /// <param name="conductorWeightPerUnitLength">mC1 [kg/m]</param>
        /// <param name="gravitionalAcceleration">g [m/s^2]</param>
        /// <returns>A double type value.</returns>
        public static double ConductorBareLoadPerUnitLength(double conductorWeightPerUnitLength, double gravitionalAcceleration)
        {
            return conductorWeightPerUnitLength * gravitionalAcceleration;
        }
        /// <summary>
        /// Calculates total conductor weight per unit length - mC2g [N/m]
        /// </summary>
        /// <param name="bareConductorLoad">mC1g [N/m]</param>
        /// <param name="conductorRelativeWindLoad">ww [N/m]</param>
        /// <param name="conductorRelativeIceLoad">wi [N/m]</param>
        /// <returns>A double type value.</returns>
        public static double ConductorTotalLoadPerUnitLenght(double bareConductorLoad, double conductorRelativeWindLoad, double conductorRelativeIceLoad)
        {
            double mC1g = bareConductorLoad;
            double ww = conductorRelativeWindLoad;
            double wi = conductorRelativeIceLoad;

            return Math.Sqrt(Math.Pow(mC1g + wi, 2) + Math.Pow(ww, 2));
        }
        /// <summary>
        /// Calculates vertical difference in span between neighboring CONDUCTOR attachment points - hc [m]. Negative value indicates that left end of conductor is attached higher than the right one.
        /// </summary>
        /// <param name="tower1Ordinate">YT1 [m]</param>
        /// <param name="Tower2Ordinate">YT2 [m]</param>
        /// <param name="insulatorAttachmentPiont1Heigth">hp1 [m]</param>
        /// <param name="insulatorAttachmentPoint2Heigth">hp2 [m]</param>
        /// <param name="insulator1ArmLength">LINS_Length [m].0 if insulator set is of strain type.</param>
        /// <param name="insulator2ArmLength">RINS_Length [m].0 if insulator set is of strain type.</param>
        /// <param name="insulator1VerticalOffset">LeINS [m]</param>
        /// <param name="insulator2VerticalOffset">ReINS [m]</param>
        /// <returns>A double type value.</returns>
        public static double SpanVerticalDifference(double tower1Ordinate, double Tower2Ordinate, double insulatorAttachmentPiont1Heigth, double insulatorAttachmentPoint2Heigth, double insulator1ArmLength, double insulator2ArmLength, double insulator1VerticalOffset, double insulator2VerticalOffset)
        {
            double YT1 = tower1Ordinate;
            double YT2 = Tower2Ordinate;
            double hp1 = insulatorAttachmentPiont1Heigth;
            double hp2 = insulatorAttachmentPoint2Heigth;
            double LINS_Length = insulator1ArmLength;
            double RINS_Length = insulator2ArmLength;
            double LeINS = insulator1VerticalOffset;
            double ReINS = insulator2VerticalOffset;

            return (YT2 + hp2 - RINS_Length) - (YT1 + hp1 - LINS_Length) + (ReINS - LeINS);
        }
        /// <summary>
        /// Calculates span length measured between neighboring CONDUCTOR attachment points - ac [m].
        /// </summary>
        /// <param name="finalHorizontalInsulatorAttachmentPointPosition1">XT1 [m]</param>
        /// <param name="finalHorizontalInsulatorAttachmentPointPosition2">XT2 [m]</param>
        /// <param name="insulator1HorizontalOffset">LdINS [m]</param>
        /// <param name="insulator2HorizontalOffset">RdINS [m]</param>
        /// <returns>A double type value.</returns>
        public static double SpanHorizontalLength(double finalHorizontalInsulatorAttachmentPointPosition1, double finalHorizontalInsulatorAttachmentPointPosition2, double insulator1HorizontalOffset, double insulator2HorizontalOffset)
        {
            double XT1 = finalHorizontalInsulatorAttachmentPointPosition1;
            double XT2 = finalHorizontalInsulatorAttachmentPointPosition2;
            double LdINS = insulator1HorizontalOffset;
            double RdINS = insulator2HorizontalOffset;

            return Math.Abs(XT2 - XT1) + (RdINS - LdINS);
        }
        /// <summary>
        /// Calculates inclined span length - L0 [m].
        /// </summary>
        /// <param name="spanVerticalDifference">h [m]</param>
        /// <param name="spanHorizontalLength">a [m]</param>
        /// <returns>A double type value.</returns>
        public static double SpanInclinedLength(double spanVerticalDifference, double spanHorizontalLength)
        {
            double h = spanVerticalDifference;
            double a = spanHorizontalLength;

            return Math.Sqrt(Math.Pow(a, 2) + Math.Pow(h, 2));
        }
        /// <summary>
        /// Calculates span catenary length - L [m].
        /// </summary>
        /// <param name="horizontalTension">H [N]</param>
        /// <param name="conductorLoadPerUnitLength">mCg [N/m]</param>
        /// <param name="spanVerticalDifference">h [m]</param>
        /// <param name="spanHorizontalLength">a [m]</param>
        /// <returns>A double type value.</returns>
        public static double SpanCatenaryLength(double horizontalTension, double conductorLoadPerUnitLength, double spanVerticalDifference, double spanHorizontalLength)
        {
            double H = horizontalTension;
            double mCg = conductorLoadPerUnitLength;
            double h = spanVerticalDifference;
            double a = spanHorizontalLength;

            return Math.Sqrt(Math.Pow(h, 2) + Math.Pow(((2 * H) / mCg) * Math.Sinh(mCg * (a / (2 * H))), 2));
        }
        // _L_ [-] = f(L0, L):[m]
        /// <summary>
        /// Calculates abbreviated span length - _L_ [-].
        /// </summary>
        /// <param name="spanInclinedLength">L0 [m]</param>
        /// <param name="spanCatenaryLength">L [m]</param>
        /// <returns>A double type value.</returns>
        public static double SpanAbbreviateLength(double spanInclinedLength, double spanCatenaryLength)
        {
            double L0 = spanInclinedLength;
            double L = spanCatenaryLength;

            return (L / L0) - 1;
        }
        /// <summary>
        /// Calculates span vertex - xV [m]. Which is represented by distance measured from vertex do left support point.
        /// </summary>
        /// <param name="horizontalTension">H [N]</param>
        /// <param name="conductorLoadPerUnitLength">mCg [N/m]</param>
        /// <param name="catenaryLength">L [m]</param>
        /// <param name="spanVerticalDifference">h [m]</param>
        /// <param name="spanHorizontalLength">a[m]</param>
        /// <returns>A double type value.</returns>
        public static double SpanVertex(double horizontalTension, double conductorLoadPerUnitLength, double catenaryLength, double spanVerticalDifference, double spanHorizontalLength)
        {
            double H = horizontalTension;
            double mCg = conductorLoadPerUnitLength;
            double L = catenaryLength;
            double h = spanVerticalDifference;
            double a = spanHorizontalLength;

            return (H / mCg) * Math.Log((H / (mCg * (L - h))) * (1 - Math.Exp(-mCg * (a / H))));
        }
        /// <summary>
        /// Calculates sag - f [m]. For x = 0 calculates sag in lowest point of catenary curve (vertex point).
        /// </summary>
        /// <param name="horizontalTension">H [N]</param>
        /// <param name="conductorLoad">mCg [N/m]</param>
        /// <param name="spanVertex">xV [m]</param>
        /// <param name="spanVerticalDifference">h [m]</param>
        /// <param name="spanHorizontalLength">a [m]</param>
        /// <param name="xRelativePosition">x [m]. Sugested range: 0 - 1*. For example 0,5 means sag in middle of span. *Exceeding sugested range may result negative sag values.</param>
        /// <returns>A double type value.</returns>
        public static double SpanSag(double horizontalTension, double conductorLoad, double spanVertex, double spanVerticalDifference, double spanHorizontalLength, double xRelativePosition)
        {
            double H = horizontalTension;
            double mCg = conductorLoad;
            double xV = spanVertex;
            double h = spanVerticalDifference;
            double a = spanHorizontalLength;
            double x = xV + (xRelativePosition * a);

            return (h / a) * (x - xV) + (H / mCg) * (Math.Cosh(mCg * xV / H) - Math.Cosh(mCg * x / H));
        }
        /// <summary>
        /// Calculates sag correction for spans containing strain attachment set or sets - p [m]. 
        /// </summary>
        /// <param name="spanLenght">a [m]</param>
        /// <param name="spanVerticalDifference">h [m]</param>
        /// <param name="leftInsulatorArmLength">LINS_Length [m]</param>
        /// <param name="leftInsulatorHorizontalOffset">LdINS [m]</param>
        /// <param name="leftInsulatorVerticalOffset">LeINS [m]</param>
        /// <param name="rightInsulatorArmLength">RINS_Length [m]</param>
        /// <param name="rightInsulatorHorizontalOffset">RdINS [m]</param>
        /// <param name="rightInsulatorVerticalOffset">ReINS [m]</param>
        /// <param name="spanConductorHorizontalLength">ac [m]</param>
        /// <param name="spanConductorRelativePosition">xc [m]. Sugested range: 0 - 1*. For example 0,5 means correction for sag in middle of span. *Exceeding sugested range results in reaching one of range border.</param>
        /// <returns>A double type value.</returns>
        public static double SpanSagCorrectionFromAttachmentSet(double spanLenght, double spanVerticalDifference, double leftInsulatorArmLength, double leftInsulatorHorizontalOffset, double leftInsulatorVerticalOffset, double rightInsulatorArmLength, double rightInsulatorHorizontalOffset, double rightInsulatorVerticalOffset, double spanConductorHorizontalLength, double spanConductorRelativePosition)
        {
            double a = spanLenght;
            double h = spanVerticalDifference;
            double L = leftInsulatorArmLength;
            double Ld = leftInsulatorHorizontalOffset;
            double Le = leftInsulatorVerticalOffset;
            double R = rightInsulatorArmLength;
            double Rd = rightInsulatorHorizontalOffset;
            double Re = rightInsulatorVerticalOffset;
            double ac = spanConductorHorizontalLength;
            double xc;

            if (spanConductorRelativePosition < 0) xc = 0;
            else if (spanConductorRelativePosition > 1) xc = 1;
            else 
            {
                xc = spanConductorRelativePosition * ac;
            }
            double ComponentA = (L - Le) + (Ld * (h / a));
            double ComponentB = (R - Re) + (Rd * (h / a));

            return (ComponentA * (ac - xc) + (ComponentB * xc)) * (1 / ac);
        }
        /// <summary>
        /// Calculates sag adjustment to sag in the middle of complementary span length - yf [m].
        /// </summary>
        /// <param name="spanVertex">xV [m]</param>
        /// <param name="spanVerticalDifference">h [m]</param>
        /// <param name="spanHorizontalLength">a [m]</param>
        /// <returns>A double type value.</returns>
        public static double SpanSagToComplementarySagAdjustment(double spanVertex, double spanVerticalDifference, double spanHorizontalLength)
        {
            double xV = spanVertex;
            double h = spanVerticalDifference;
            double a = spanHorizontalLength;

            return xV * (h / a);
        }
        /// <summary>
        /// Calculates horizontal resultant force.
        /// </summary>
        /// <param name="precedingHorizontalTension">Hprv [N]</param>
        /// <param name="followingHorizontalTension">Hnxt [N] </param>
        /// <param name="insulatorTotalVerticalForce">G [N]</param>
        /// <returns>A double type value.</returns>
        public static double SupportResultantForce(double precedingHorizontalTension, double followingHorizontalTension, double insulatorTotalVerticalForce = 0)
        {
            double H1 = precedingHorizontalTension;
            double H2 = followingHorizontalTension;
            double G = insulatorTotalVerticalForce;

            if (G == 0) return H2 - H1;
            else return Math.Sqrt(Math.Pow(H2 - H1, 2) + Math.Pow(G, 2));
        }
        /// <summary>
        /// Calculates insulator arm dead weigth - JK [N].
        /// </summary>
        /// <param name="insulatorWeight">mINS [kg]</param>
        /// <param name="gravitationalAcceleration">g [m/s^2]</param>
        /// <param name="insulatorIceLoad">wINS [N]. 0 if calculated for non ice conditions.</param>
        /// <returns>A double type value.</returns>
        public static double InsulatorArmDeadWeight(double insulatorWeight, double gravitationalAcceleration = 9.80665, double insulatorIceLoad = 0)
        {
            double mINS = insulatorWeight;
            double wINS = insulatorIceLoad;
            double g = gravitationalAcceleration;

            return (mINS * g) + wINS;
        }
        // GK [N] = f(JK [N], H [N], mCg0 [N/m], mCg1 [N/m], xV0 [m], xV1 [m], a0 [m], X0 [-], X1 [-], alphaINS [deg])
        /// <summary>
        /// Calculates total vertical force which is acting on corresponding tower - GK [N].
        /// </summary>
        /// <param name="insulatorDeadWeight">JK [N]</param>
        /// <param name="initialHorizontalTension">H [N]</param>
        /// <param name="precedingSpanConductorLoad">mCg0 [N/m]</param>
        /// <param name="followingSpanConductorLoad">mCg1 [N/m]</param>
        /// <param name="precedingSpanVertex">xV0 [m]</param>
        /// <param name="followingSpanVertex">xV1 [m]</param>
        /// <param name="precedingSpanHorizontalLength">a0 [m]</param>
        /// <param name="precedingSpanFinalTensionModifier">X0 [-]. 1 if GK is calculated for initial state.</param>
        /// <param name="followingSpanFinalTensionModifier">X1 [-]. 1 if GK is calculated for initial state.</param>
        /// <param name="insulatorInvertedVOpeningAngle">alphaINS [deg]. 0 for long rod insulator's sets.</param>
        /// <returns>A double type value.</returns>
        public static double InsulatorTotalVerticalForce(double insulatorDeadWeight, double initialHorizontalTension, double precedingSpanConductorLoad, double followingSpanConductorLoad, double precedingSpanVertex, double followingSpanVertex, double precedingSpanHorizontalLength, double precedingSpanFinalTensionModifier = 1, double followingSpanFinalTensionModifier = 1, double insulatorInvertedVOpeningAngle = 0)
        {
            double mINSg = insulatorDeadWeight;
            double H = initialHorizontalTension;
            double mCg0 = precedingSpanConductorLoad;
            double mCg1 = followingSpanConductorLoad;
            double xV0 = precedingSpanVertex;
            double xV1 = followingSpanVertex;
            double a0 = precedingSpanHorizontalLength;
            double X0 = precedingSpanFinalTensionModifier;
            double X1 = followingSpanFinalTensionModifier;
            double alphaINS = insulatorInvertedVOpeningAngle;

            double armFactor;

            double Component_M = mINSg / 2;                                             // Force component from insulator;
            double Component_0 = X0 * H * Math.Sinh((mCg0 * (a0 + xV0)) / (X0 * H));    // Force component from preceding span conductor;
            double Component_1 = X1 * H * Math.Sinh((mCg1 * xV1) / (X1 * H));           // Force component from following span conductor;

            if (alphaINS > 0) armFactor = 2;
            else armFactor = 1;

            return (armFactor * Component_M) + Component_0 - Component_1;
        }
        // dINS [m] = f(Hprv [N], Hnxt [N], LK [m], GK [N], alphaINS [deg]) - if calculated for invertedV insulator set (alphaINS != 0), returns offset for right arm!
        /// <summary>
        /// Calculates insulator horizontal offset - dINS [m].
        /// </summary>
        /// <param name="precedingHorizontalTension">Hprv [N]</param>
        /// <param name="followingHorizontalTension">Hnxt [N]</param>
        /// <param name="insulatorArmLength">LK [m]</param>
        /// <param name="insulatorTotalVerticalForce">GK [N]</param>
        /// <param name="insulatorInvertedVOpeningAngle">alphaINS [deg]. If calculated for inverted-V insulator set (alphaINS != 0), method returns offset for right arm!</param>
        /// <returns>A double type value.</returns>
        public static double InsulatorHorizontalOffset(double precedingHorizontalTension, double followingHorizontalTension, double insulatorArmLength, double insulatorTotalVerticalForce, double insulatorInvertedVOpeningAngle = 0)
        {
            double H0 = precedingHorizontalTension;
            double H1 = followingHorizontalTension;
            double LK = insulatorArmLength;
            double GK = insulatorTotalVerticalForce;

            double alphaINS = insulatorInvertedVOpeningAngle;
            double alphaINS_RAD, psiINS_RAD;

            if (alphaINS > 0)
            {
                alphaINS_RAD = alphaINS * (Math.PI / 180);
                psiINS_RAD = InsulatorInvertedVOffsetAngle(H0, H1, alphaINS, GK) * (Math.PI / 180);

                return LK * (Math.Sin(alphaINS_RAD / 2 + psiINS_RAD) - Math.Sin(alphaINS_RAD / 2));
            }

            return ((H1 - H0) * LK) / Math.Sqrt(Math.Pow(GK, 2) + Math.Pow(H1 - H0, 2));
        }
        /// <summary>
        /// Calculates insulator vertical offset - eINS [m].
        /// </summary>
        /// <param name="precedingHorizontalTension">Hprv [N]</param>
        /// <param name="followingHorizontalTension">Hnxt [N]</param>
        /// <param name="insulatorLength">LK [m]</param>
        /// <param name="insulatorTotalVerticalForce">GK [N]</param>
        /// <param name="insulatorInvertedVOpeningAngle">alphaINS [deg]. If calculated for inverted-V insulator set (alphaINS != 0), method returns offset for right arm!</param>
        /// <returns>A double type value.</returns>
        public static double InsulatorVerticalOffset(double precedingHorizontalTension, double followingHorizontalTension, double insulatorLength, double insulatorTotalVerticalForce, double insulatorInvertedVOpeningAngle = 0)
        {
            double H0 = precedingHorizontalTension;
            double H1 = followingHorizontalTension;
            double LK = insulatorLength;
            double GK = insulatorTotalVerticalForce;

            double alphaINS = insulatorInvertedVOpeningAngle;
            double alphaINS_RAD, psiINS_RAD;

            if (alphaINS > 0)
            {
                alphaINS_RAD = alphaINS * (Math.PI / 180);
                psiINS_RAD = InsulatorInvertedVOffsetAngle(H0, H1, alphaINS, GK) * (Math.PI / 180);

                return LK * (Math.Cos(alphaINS_RAD / 2) - Math.Cos(alphaINS_RAD / 2 + psiINS_RAD));
            }

            return LK * (1 - (GK / (Math.Sqrt(Math.Pow(GK, 2) + Math.Pow(H1 - H0, 2)))));
        }
        /// <summary>
        /// Calculates (generally) final offset angle for inverted-V instulator set - psiINS [deg].
        /// </summary>
        /// <param name="precedingHorizontalTension">H2prv [N]</param>
        /// <param name="followingHorizontalTension">H2nxt [N]</param>
        /// <param name="insulatorInvertedVOpeningAngle">alphaINS [deg]</param>
        /// <param name="insulatorInvertedVFinalTotalVerticalForce">GK [N]</param>
        /// <returns>A double type value.</returns>
        public static double InsulatorInvertedVOffsetAngle(double precedingHorizontalTension, double followingHorizontalTension, double insulatorInvertedVOpeningAngle, double insulatorInvertedVFinalTotalVerticalForce)
        {
            double H0 = precedingHorizontalTension;
            double H1 = followingHorizontalTension;
            double alphaINS_RAD = insulatorInvertedVOpeningAngle * (Math.PI / 180);
            double GK = insulatorInvertedVFinalTotalVerticalForce;

            double NUMERATOR = (H1 - H0) * Math.Cos(alphaINS_RAD / 2) + GK * Math.Sin(alphaINS_RAD / 2);
            double DENOMINATOR = (H1 + H0) * Math.Sin(alphaINS_RAD / 2) + GK * Math.Cos(alphaINS_RAD / 2);

            return NUMERATOR / DENOMINATOR;
        }
        /// <summary>
        /// Calculates inverted-V insulator set's bridge horizontal span length - aINS [m].
        /// </summary>
        /// <param name="insulatorLength">LK [m]</param>
        /// <param name="insulatorInvertedVOpeningAngle">alphaINS [deg]</param>
        /// <param name="insulatorInvertedVVerticalDifference">hINS [m]</param>
        /// <returns>A double type value.</returns>
        public static double InsulatorInvertedVBridgeHorizontalSpan(double insulatorLength, double insulatorInvertedVOpeningAngle, double insulatorInvertedVVerticalDifference)
        {
            double LK = insulatorLength;
            double alphaINS_RAD = insulatorInvertedVOpeningAngle * (Math.PI / 180);
            double h = insulatorInvertedVVerticalDifference;

            double pointToPointSpan = 2 * Math.Pow(LK, 2) * (1 - Math.Cos(alphaINS_RAD));
            return Math.Sqrt(Math.Pow(pointToPointSpan, 2) - Math.Pow(h, 2));
        }
        /// <summary>
        /// Calculates inverted-V insulator set's bridge vertical difference length - hINS [m]. 
        /// </summary>
        /// <param name="rightArmHorizontalOffset">dINS</param>
        /// <param name="rightArmVerticalOffset">eINS</param>
        /// <returns>A double type value.</returns>
        public static double InsulatorInvertedVBridgeVerticalDifference(double rightArmHorizontalOffset, double rightArmVerticalOffset)
        {
            double LdINS = rightArmHorizontalOffset;
            double LeINS = rightArmVerticalOffset;

            return Math.Abs(LdINS) + Math.Abs(LeINS);
        }
    }
}

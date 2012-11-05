/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.iso2d;

import java.util.LinkedList;
import java.util.List;
import net.epsilony.iso.IsoFunction;
import net.epsilony.iso.SegmentLengthTooSmall;

/**
 *
 * @author epsilon
 */
public class Polygonizer {

    double minThetaCos, minLen, alpha, initLen, error;
    IsoFunction fun;

    public Polygonizer(double initLen, IsoFunction fun) {
        this.initLen = initLen;
        this.fun = fun;
        minLen=initLen/10;
        alpha=0.8;
        error=1e-10;
        minThetaCos=Math.cos(Math.PI/3);
    }

    public void setMinLen(double minLen) {
        this.minLen = minLen;
    }

    public void setInitLen(double initLen) {
        this.initLen = initLen;
    }

    public void setError(double error) {
        this.error = error;
    }

    public void setFun(IsoFunction fun) {
        this.fun = fun;
    }
    
    /**
     * Determine a candidate point for calculating accurate isoline point. The candidate
     * point should be near to the iso line. Given a point <i>P</i> on the isoline and the
     * gradient on point <i>P</i>, the candidate point <i>Q</i> will be on the tagent line
     * across <i>P</i>. In order to get a unique solution, the distance betwean <i>P</i> 
     * and <i>Q</i> must be supplied, and the marching direction (counterclockwise or not) 
     * of the isoline construction must be also given.
     * @param funValOnP The function value of the given {@link IsoFunction} at point <i> P </i>, which contains the function value and gradient.
     * If {@code funValOnP} is {@code} null. The paramenter {@fun} may not be {@code null}. Otherwise, if {@code funValOnP} is not {@null} the given {@code fun}
     * will be ignored.
     * @param fun The field function of which isoline is going to be determined. Only when {@code funValOp} is {@code null} this parameter is available.
     * @param p The coordinate of point <i>P</i> which is on the isoline.
     * @param countercw Wheather the isoline is marching counterclockwisely.
     * @param dis The distance betwean the result point <i>Q</i> and given isoline point <i>P</i>
     * @param result The result containter for return. Can be null.
     * @return {@code result} or a new {@code double[]} which contains <i>Q</i>
     */
    public static double[] candidatePoint(double[] funValOnP, IsoFunction fun, double[] p, boolean countercw, double dis, double[] result) {
        final int gradIndex = 1;
        if (null == funValOnP) {
            funValOnP = fun.value(p, 1, null);
        }
        double pX = p[0], pY = p[1], nX = funValOnP[gradIndex], nY = funValOnP[gradIndex + 1];
        double qX, qY;
        double r = dis / Math.sqrt(nX * nX + nY * nY);
        if (countercw) {
            qX = pX - r * nY;
            qY = pY + r * nX;
        } else {
            qX = pX + r * nY;
            qY = pY - r * nX;
        }
        if (null == result) {
            return new double[]{qX, qY};
        } else {
            result[0] = qX;
            result[1] = qY;
            return result;
        }
    }

    /**
     * Determine a point <i>P</i> on the isoline near the given candidate point <i>Q</i>
     * @param fun
     * @param q
     * @param err
     * @param funResult The result of {@link IsoFunction#value(double[], int, double[]) } at final point <i>P</i>.
     * @return the coordinate of point <i>P</i>
     */
    public static double[] pointOnCurve(IsoFunction fun, double[] q, double err, double[] funResult) {
        double[] result = new double[2];
        result[0] = q[0];
        result[1] = q[1];
        if (null == funResult || funResult.length < 3) {
            funResult = new double[3];
        }
        fun.value(q, 1, funResult);
        while (Math.abs(funResult[0]) > err) {
            double f = funResult[0];
            double nX = funResult[1];
            double nY = funResult[2];
            double c = f / (nX * nX + nY * nY);
            result[0] = result[0] - c * nX;
            result[1] = result[1] - c * nY;
            fun.value(result, 1, funResult);
        }
        return result;
    }

    /**
     * Determine a new candidate point distance from current activate point. Which is
     * designed to be called when the {@code formalDist} lead to a new isoline point which
     * cannot pass the {@link #checkNewPoint(double[], double[], double[], double[]) checkNewPoint} test after
     * calling {@link #candidatePoint(double[], net.epsilony.isomesh.IsoFunction, double[], boolean, double, double[]) candidatePoint}
     * and {@link #pointOnCurve(net.epsilony.isomesh.IsoFunction, double[], double, double[])  pointOnCurve}.
     * @param formalDist
     * @return 
     */
    double newDistance(double formalDist) {
        return formalDist * alpha;
    }

    /**
     * Using a marching method to construct a isoline made by line segments.
     * @param init Initial length of line segment or the distance betwean nearby vertes. 
     * This parameter affacts the biggest length of the line segments which constructs the isoline.
     * @return A list of the vertex coordinates of isoline, which is counterclockwise.
     * @throws SegmentLengthTooSmall 
     */
    List<double[]> mesh(double[] init) throws SegmentLengthTooSmall {
        double[] funValOnEnd = new double[3];
        double[] funValOnStart = new double[3];
        double[] start = pointOnCurve(fun, init, error, funValOnStart);
        double[] end = pointOnCurveWithCheck(start, funValOnStart, initLen, funValOnEnd);
        LinkedList<double[]> polygon = new LinkedList<>();
        polygon.add(start);
        polygon.add(end);
        double[] funValOnNewEnd=new double[3];
        do{
            double[] newEnd=pointOnCurveWithCheck(end, funValOnEnd, initLen, funValOnNewEnd);
            polygon.add(newEnd);
            double[] t=funValOnEnd;
            funValOnEnd=funValOnNewEnd;
            funValOnNewEnd=t;
            
            t=end;
            end=newEnd;
            newEnd=t;
            if(checkNewPoint(end, start, funValOnEnd, funValOnStart)){
                break;
            }
        }while(true);
        return polygon;
    }

    /**
     * Get a new point, which passes the {@link #checkNewPoint(double[], double[], double[], double[]) checkNewPoint}.
     * The new point is on isoline from activate end point <i>P</i>. With a given estimate distance.
     * @param p
     * @param funValOnP
     * @param initDis 
     * @param funValOnNewPt
     * @return
     * @throws SegmentLengthTooSmall 
     */
    double[] pointOnCurveWithCheck(double[] p, double[] funValOnP,double initDis,  double[] funValOnNewPt) throws SegmentLengthTooSmall {
        double len = initDis;
        boolean checkPass;
        double[] q = new double[2];
        double[] newPt;
        do {
            candidatePoint(funValOnP, fun, p, true, len, q);
            newPt = pointOnCurve(fun, q, error, funValOnNewPt);
            checkPass = checkNewPoint(p, funValOnP,newPt,  funValOnNewPt);
            if (checkPass) {
                break;
            } else {
                double t1 = p[0] - newPt[0];
                double t2 = p[1] - newPt[1];
                len = newDistance(Math.sqrt(t1 * t1 + t2 * t2));
                if (!checkLength(len)) {
                    throw new SegmentLengthTooSmall();
                }
            }
        } while (true);
        return newPt;
    }

    /**
     * Only check if the given segment lenght is too small now.
     * @param len
     * @return 
     */
    boolean checkLength(double len) {
        if (len < minLen) {
            return false;
        } else {
            return true;
        }
    }

    /**
     * Check weather the line segment <i>P-Q</i> fit the condition. Where <i>P</i> is the
     * start point and <i>Q</i> is the end point. The condition is that the cross angle betwean 
     * <i>grad_<sub>P</sub></i> and <i>n</i> and betwean <i>grad_<sub>Q</sub></i> and <i>n</i> must
     * be smaller that a given parameter. Where <i>grad_<sub>Q</sub></i> and <i>grad_<sub>Q</sub></i>
     * are the gradient at <i>P</i> and <i>Q</i> and <i>n</i> is the out normal of segment <i>P-Q</i>.
     * @param start
     * @param funValOnStart
     * @param end
     * @param funValOnEnd
     * @return 
     */
    boolean checkNewPoint(double[] start, double[] funValOnStart,double[] end,  double[] funValOnEnd) {
        final int gradIndex=1;
        double startX = start[0], startY = start[1],
                endX = end[0], endY = end[1],
                startGradX = funValOnStart[gradIndex], startGradY = funValOnStart[gradIndex + 1],
                endGradX = funValOnEnd[gradIndex], endGradY = funValOnEnd[gradIndex + 1];
        double vX = endX - startX, vY = endY - startY;
        double vLen = Math.sqrt(vX * vX + vY * vY);
        double segNormalX = -vY / vLen, segNormY = vX / vLen;
        double startNormalLen = Math.sqrt(startGradX * startGradX + startGradY * startGradY);
        double endNormalLen = Math.sqrt(endGradX * endGradX + endGradY * endGradY);
        double startCos = (startX * segNormalX + startY * segNormY) / startNormalLen;
        double endCos = (endX * segNormalX + endY * segNormY) / endNormalLen;
        if (startCos < minThetaCos || endCos < minThetaCos) {
            return false;
        } else {
            return true;
        }

    }
}

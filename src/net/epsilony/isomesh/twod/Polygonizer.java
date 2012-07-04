/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.isomesh.twod;

import java.util.LinkedList;
import java.util.List;
import net.epsilony.isomesh.IsoFunction;
import net.epsilony.isomesh.util.SegmentLengthTooSmall;

/**
 *
 * @author epsilon
 */
public class Polygonizer {

    double minThetaCos, minLen, alpha, initLen, error;
    IsoFunction fun;

    public static double[] candidatePoint(IsoFunction fun, double[] p, boolean pRot, double len, double[] funValOnP, double[] result) {
        final int gradIndex = 1;
        if (null == funValOnP) {
            funValOnP = fun.value(p, 1, null);

        }
        double pX = p[0], pY = p[1], nX = funValOnP[gradIndex], nY = funValOnP[gradIndex + 1];
        double qX, qY;
        double r = len / Math.sqrt(nX * nX + nY * nY);
        if (pRot) {
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

    double newLenght(double len) {
        return len * alpha;
    }

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

    double[] pointOnCurveWithCheck(double[] p, double[] funValOnP,double initLen,  double[] funValOnNewPt) throws SegmentLengthTooSmall {
        double len = initLen;
        boolean checkPass;
        double[] q = new double[2];
        double[] newPt;
        do {
            candidatePoint(fun, p, true, len, funValOnP, q);
            newPt = pointOnCurve(fun, q, error, funValOnNewPt);
            checkPass = checkNewPoint(p, funValOnP,newPt,  funValOnNewPt);
            if (checkPass) {
                break;
            } else {
                double t1 = p[0] - newPt[0];
                double t2 = p[1] - newPt[1];
                len = newLenght(Math.sqrt(t1 * t1 + t2 * t2));
                if (!checkLength(len)) {
                    throw new SegmentLengthTooSmall();
                }
            }
        } while (true);
        return newPt;
    }

    boolean checkLength(double len) {
        if (len < minLen) {
            return false;
        } else {
            return true;
        }
    }

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

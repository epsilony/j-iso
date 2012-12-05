/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.functions;

import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.spfun.CommonUtils;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class EllipseLSF implements CoordinatePartDiffFunction {

    private int diffOrder;
    double xc, yc, a, b;

    public EllipseLSF(double xc, double yc, double a, double b) {
        this.xc = xc;
        this.yc = yc;
        this.a = a;
        this.b = b;

    }
    private final static double LOG2 = Math.log(2);

    @Override
    public double[] values(Coordinate coord, double[] results) {
        if (null == results) {
            results = new double[CommonUtils.len2DBase(diffOrder)];
        }
        double x = coord.x - xc;
        double y = coord.y - yc;
        double delta = x * x / (a * a)+y * y / (b * b);
        double t1=Math.exp(-delta * LOG2);
        results[0] = 2 * t1 - 1;
        if (diffOrder > 0) {
            double t2=-4*LOG2*t1;
            results[1] = t2*x/(a*a);
            results[2] = t2*y/(b*b);
        }
        return results;
    }

    @Override
    public void setDiffOrder(int order) {
        this.diffOrder = order;
    }

    @Override
    public int getDiffOrder() {
        return diffOrder;
    }

    @Override
    public int getDim() {
        return 2;
    }
}

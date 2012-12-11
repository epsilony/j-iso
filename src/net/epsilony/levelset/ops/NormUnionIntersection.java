/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.ops;

import static java.lang.Math.*;
import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class NormUnionIntersection extends UnionIntersection {

    double b;
    double x1, x2;
    double a1, a2;
    private final double B_SHRINK = 0.98;

    private void deterCoefs() {
        x1 = isUnion ? -2 + sqrt(2) : -2 - sqrt(2);
        x2 = isUnion ? 2 + sqrt(2) : 2 - sqrt(2);

        x1 *= pow(2, continueM / 2.0);
        x2 *= pow(2, continueM / 2.0);

        b = Math.min(-2 / x1, 2 / x2) * B_SHRINK;

        a1 = (-1 - b * x1) / (x1 * x1);
        a2 = (1 - b * x2) / (x2 * x2);
    }

    public NormUnionIntersection(boolean isUnion, int continueM, boolean throwSingular, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
        super(isUnion, continueM, throwSingular, fun1, fun2, dim);
        deterCoefs();
    }

    @Override
    public double[] values(Coordinate coord, double[] results) {
        double[] res = super.values(coord, results);

        double f = res[0];
        double a = f <= 0 ? a1 : a2;
        res[0] = a * f * f + b * f;
        double r_f = 2 * a * f + b;
        if (diffOrder > 0) {
            for (int i = 0; i < dim; i++) {
                res[1 + i] *= r_f;
            }
        }

        return res;
    }
}

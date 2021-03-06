/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.functions;

import net.epsilony.levelset.ops.AtomOperations;
import net.epsilony.levelset.ops.LogisticNorm;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class LogistNormHalfPlane extends HalfPlaneLinearLSF {

    double k;

    public LogistNormHalfPlane(Coordinate p1, Coordinate p2, double k) {
        super(p1, p2, 1);
        this.k = k;
    }

    @Override
    public double[] values(Coordinate coord, double[] results) {
        results = super.values(coord, results);
        LogisticNorm.valueOfLogistic(results, results, k, diffOrder, getDim());
        return results;
    }
}

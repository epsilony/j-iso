/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.ops;

import net.epsilony.math.CoordinatePartDiffFunction;

/**
 *
 * @author epsilon
 */
public class AtomOperations {

    public static CoordinatePartDiffFunction scale(CoordinatePartDiffFunction fun, double scale) {
        return new Scale(fun, scale);
    }

    public static CoordinatePartDiffFunction logisticNorm(CoordinatePartDiffFunction fun, double k) {
        return new LogisticNorm(fun, k);
    }

    public static CoordinatePartDiffFunction union_intersection(boolean isUnion, int m, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
        return new UnionIntersection(isUnion, m, false, fun1, fun2, dim);
    }

    public static CoordinatePartDiffFunction union_intersection(boolean isUnion, int m, boolean throwZeroDivider, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
        return new UnionIntersection(isUnion, m, throwZeroDivider, fun1, fun2, dim);
    }
    
    public static CoordinatePartDiffFunction normed_union_intersection(boolean isUnion, int m, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
        return new NormUnionIntersection(isUnion, m, false, fun1, fun2, dim);
    }

    public static CoordinatePartDiffFunction normed_union_intersection(boolean isUnion, int m, boolean throwZeroDivider, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
        return new NormUnionIntersection(isUnion, m, throwZeroDivider, fun1, fun2, dim);
    }

    public static CoordinatePartDiffFunction complement(CoordinatePartDiffFunction fun) {
        return new Complement(fun);
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset;

import java.util.Arrays;
import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.spfun.CommonUtils;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class AtomOperations {

    public static CoordinatePartDiffFunction scale(CoordinatePartDiffFunction fun, double scale) {
        return new Scale(fun, scale);
    }

    public static CoordinatePartDiffFunction logistic(CoordinatePartDiffFunction fun, double k) {
        return new Logistic(fun, k);
    }

    public static CoordinatePartDiffFunction union_intersection(boolean isUnion,int m, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
        return new UnionIntersection(isUnion, m, false, fun1, fun2, dim);
    }

    public static CoordinatePartDiffFunction union_intersection(boolean isUnion,int m, boolean throwZeroDivider, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
        return new UnionIntersection(isUnion, m, throwZeroDivider, fun1, fun2, dim);
    }

    public static CoordinatePartDiffFunction complement(CoordinatePartDiffFunction fun) {
        return new Complement(fun);
    }

    public static void valueOfComplement(double[] ori, double[] results) {
        for (int i = 0; i < ori.length; i++) {
            results[i] = -ori[i];
        }
    }

    /**
     *
     * @param ori
     * @param results
     * @param scale
     * @param k
     * @param partDiffOrder
     * @param dim
     */
    public static void valueOfScale(double[] ori, double[] results, double k, int partDiffOrder, int dim) {
        System.arraycopy(ori, 0, results, 0, CommonUtils.lenBase(dim, partDiffOrder));
    }

    /**
     * $\frac{2s}{1+e^{-kt}}-s$
     *
     * @param ori
     * @param results
     * @param scale
     * @param k
     * @param partDiffOrder
     * @param dim
     */
    public static void valueOfLogistic(double[] ori, double[] results, double k, int partDiffOrder, int dim) {
        double t = ori[0];
        double e = Math.exp(-k * t);
        double b = (1 + e);
        double f = 2 / b - 1;
        results[0] = f;
        if (partDiffOrder > 0) {
            double f_t = e * k * 2 / (b * b);
            for (int i = 0; i < dim; i++) {
                results[i + 1] = ori[i + 1] * f_t;
            }
        }
    }

    private static void checkSupport(int dim, int partDiffOrder) {
        if (dim > 3 || dim < 2 || partDiffOrder < 0 || partDiffOrder > 1) {
            throw new UnsupportedOperationException("Unsupported dimention or partial differential order");
        }
    }

    public static void valueOfUnionIntersection(boolean isUnion, int m, double[] v1, double[] v2, int dim, int partDiffOrder, boolean throwSingular, double[] results) {
        checkSupport(dim, partDiffOrder);
        double f1 = v1[0];
        double f2 = v2[0];
        double sqr = f1 * f1 + f2 * f2;
        double sqrtv = Math.sqrt(sqr);

        double sn = 0;
        if (isUnion) {
            sn = 1;
        } else {
            sn = -1;
        }
        double f_m0 = f1 + f2 + sn * sqrtv;
        double f = 0;
        switch (m) {
            case 0:
                f = f_m0;
                break;
            case 1:
                f = f_m0 * sqrtv;
                break;
            case 2:
                f = f_m0 * sqr;
                break;
            default:
                throw new IllegalArgumentException("m should be 0,1 or 2 but given" + m);
        }
        results[0] = f;
        if (partDiffOrder > 0) {
            if (f1 == 0 && f2 == 0) {
                if (m == 0 && throwSingular) {
                    throw new ArithmeticException("Singular at the both zero value point!");
                } else {
                    for (int i = 1; i <= dim; i++) {
                        results[i] = 0;
                    }
                }
            } else {
                double f_f1 = 0, f_f2 = 0;

                switch (m) {
                    case 0:
                        f_f1 = 1 + sn * f1 / sqrtv;
                        f_f2 = 1 + sn * f2 / sqrtv;
                        break;
                    case 1:
                        f_f1 = f1 * f_m0 / sqrtv + sqrtv * (1 + sn * f1 / sqrtv);
                        f_f2 = f2 * f_m0 / sqrtv + sqrtv * (1 + sn * f2 / sqrtv);
                        break;
                    case 2:
                        f_f1 = 2 * f1 * f_m0 + sqr * (1 + sn * f1 / sqrtv);
                        f_f2 = 2 * f2 * f_m0 + sqr * (1 + sn * f2 / sqrtv);
                        break;
                    default:
                        throw new IllegalArgumentException("m should be 0,1 or 2 but given" + m);

                }

                for (int i = 1; i <= dim; i++) {
                    results[i] = v1[i] * f_f1 + v2[i] * f_f2;
                }
            }
        }
    }

    static class UnionIntersection implements CoordinatePartDiffFunction {

        private final boolean isUnion;
        protected final CoordinatePartDiffFunction fun1;
        protected final CoordinatePartDiffFunction fun2;
        protected final int dim;
        protected int diffOrder;
        protected final boolean throwWhenSingular;
        protected int continueM;

        public UnionIntersection(boolean isUnion, int continueM, boolean throwSingular, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
            this.fun1 = fun1;
            this.fun2 = fun2;
            this.dim = dim;
            this.diffOrder = 0;
            this.continueM = continueM;
            this.throwWhenSingular = throwSingular;
            fun1.setDiffOrder(0);
            fun2.setDiffOrder(0);
            this.isUnion = isUnion;
        }

        @Override
        public double[] values(Coordinate coord, double[] results) {
            double[] res1 = fun1.values(coord, null);
            double[] res2 = fun2.values(coord, null);
            results = initAndCheckSingular(results, res1, res2);
            valueOfUnionIntersection(isUnion, continueM, res1, res2, dim, diffOrder, throwWhenSingular, results);
            return results;
        }

        public boolean isThrowWhenSingular() {
            return throwWhenSingular;
        }

        @Override
        public void setDiffOrder(int order) {
            if (order > 1) {
                throw new UnsupportedOperationException("Only 0 and 1 differential order is supported.");
            }
            this.diffOrder = order;
            fun1.setDiffOrder(order);
            fun2.setDiffOrder(order);
        }

        @Override
        public int getDiffOrder() {
            return diffOrder;
        }

        @Override
        public int getDim() {
            return dim;
        }

        protected double[] initAndCheckSingular(double[] results, double[] res1, double[] res2) throws ArithmeticException {
            if (null == results) {
                results = new double[CommonUtils.lenBase(dim, diffOrder)];
            }
            if (res1[0] == 0 && res2[0] == 0) {
                if (throwWhenSingular) {
                    throw new ArithmeticException("At sharp corner!");
                } else {
                    Arrays.fill(results, 0);
                }
            }
            return results;
        }
    }

    public static class Complement implements CoordinatePartDiffFunction {

        private final CoordinatePartDiffFunction fun;

        public Complement(CoordinatePartDiffFunction fun) {
            this.fun = fun;
        }

        @Override
        public double[] values(Coordinate coord, double[] results) {
            fun.values(coord, results);
            for (int i = 0; i < results.length; i++) {
                results[i] *= -1;
            }
            return results;
        }

        @Override
        public void setDiffOrder(int order) {
            fun.setDiffOrder(order);
        }

        @Override
        public int getDiffOrder() {
            return fun.getDiffOrder();
        }

        @Override
        public int getDim() {
            return fun.getDim();
        }
    }

    public static class Scale implements CoordinatePartDiffFunction {

        CoordinatePartDiffFunction oriFun;
        double scale;

        public Scale(CoordinatePartDiffFunction oriFun, double scale) {
            this.oriFun = oriFun;
            this.scale = scale;
        }

        @Override
        public double[] values(Coordinate coord, double[] results) {
            results = oriFun.values(coord, results);
            valueOfScale(results, results, scale, getDiffOrder(), getDim());
            return results;
        }

        @Override
        public void setDiffOrder(int order) {
            if (order > 1 || order < 0) {
                throw new UnsupportedOperationException("Only supports order 0 or 1.");
            }
            oriFun.setDiffOrder(order);
        }

        @Override
        public int getDiffOrder() {
            return oriFun.getDiffOrder();
        }

        @Override
        public int getDim() {
            return oriFun.getDim();
        }
    }

    public static class Logistic implements CoordinatePartDiffFunction {

        CoordinatePartDiffFunction oriFun;
        double k;

        public Logistic(CoordinatePartDiffFunction oriFun, double k) {
            this.oriFun = oriFun;
            this.k = k;
        }

        @Override
        public double[] values(Coordinate coord, double[] results) {
            results = oriFun.values(coord, results);
            valueOfLogistic(results, results, k, getDiffOrder(), getDim());
            return results;
        }

        @Override
        public void setDiffOrder(int order) {
            if (order > 1 || order < 0) {
                throw new UnsupportedOperationException("Only supports order 0 or 1.");
            }
            oriFun.setDiffOrder(order);
        }

        @Override
        public int getDiffOrder() {
            return oriFun.getDiffOrder();
        }

        @Override
        public int getDim() {
            return oriFun.getDim();
        }
    }
}

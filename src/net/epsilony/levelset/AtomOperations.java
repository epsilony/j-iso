/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset;

import java.util.Arrays;
import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.spfun.CommonUtils;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class AtomOperations {

    public static CoordinatePartDiffFunction logistic(CoordinatePartDiffFunction fun, double scale, double k) {
        return new Logistic(fun, scale, k);
    }

    public static CoordinatePartDiffFunction union(int m, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
        return new UnionIntersection(true, m, false, fun1, 1.0, fun2, 1.0, dim);
    }

    public static CoordinatePartDiffFunction intersection(int m, CoordinatePartDiffFunction fun1, CoordinatePartDiffFunction fun2, int dim) {
        return new UnionIntersection(false, m, false, fun1, 1.0, fun2, 1.0, dim);
    }

    public static CoordinatePartDiffFunction union(int m, CoordinatePartDiffFunction fun1, double scale1, CoordinatePartDiffFunction fun2, double scale2, int dim) {
        return new UnionIntersection(true, m, false, fun1, scale1, fun2, scale2, dim);
    }

    public static CoordinatePartDiffFunction intersection(int m, CoordinatePartDiffFunction fun1, double scale1, CoordinatePartDiffFunction fun2, double scale2, int dim) {
        return new UnionIntersection(false, m, false, fun1, scale1, fun2, scale2, dim);
    }

    public static CoordinatePartDiffFunction union(boolean throwSinglur, CoordinatePartDiffFunction fun1, double scale1, CoordinatePartDiffFunction fun2, double scale2, int dim) {
        return new UnionIntersection(true, 0, throwSinglur, fun1, scale1, fun2, scale2, dim);
    }

    public static CoordinatePartDiffFunction intersection(boolean throwSinglur, CoordinatePartDiffFunction fun1, double scale1, CoordinatePartDiffFunction fun2, double scale2, int dim) {
        return new UnionIntersection(false, 0, throwSinglur, fun1, scale1, fun2, scale2, dim);
    }

    public static CoordinatePartDiffFunction complement(CoordinatePartDiffFunction fun) {
        return new Complement(fun);
    }

    public static void valueOfComplement(double[] ori, double[] results) {
        for (int i = 0; i < ori.length; i++) {
            results[i] = -ori[i];
        }
    }

    public static void valueOfLogistic(double[] ori, double[] results, double scale, double k, int partDiffOrder, int dim) {
        double t = ori[0];
        double e = Math.exp(-k * t);
        double b = (1 + e);
        double f = scale * 2 / b - scale;
        results[0] = f;
        if (partDiffOrder > 0) {
            double f_t = e * k * scale * 2 / (b * b);
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

    public static void valueOfUnion(double[] v1, double scale1, double[] v2, double scale2, int dim, int partDiffOrder, boolean throwSingular, double[] results) {
        checkSupport(dim, partDiffOrder);
        double f1 = v1[0] * scale1;
        double f2 = v2[0] * scale2;
        double sqr = f1 * f1 + f2 * f2;
        double sqrtv = Math.sqrt(sqr);
        double f = f1 + f2 + sqrtv;
        results[0] = f;
        if (partDiffOrder > 0) {
            if (f1 == 0 && f2 == 0) {
                if (throwSingular) {
                    throw new ArithmeticException("Singular at the both zero value point!");
                } else {
                    for (int i = 1; i <= dim; i++) {
                        results[i] = 0;
                    }
                }
            } else {
                double f_f1 = 1 + f1 / sqrtv;
                double f_f2 = 1 + f2 / sqrtv;
                for (int i = 1; i <= dim; i++) {
                    results[i] = v1[i] * f_f1 * scale1 + v2[i] * f_f2 * scale2;
                }
            }
        }
    }

    public static void valueOfUnionC1(double[] v1, double scale1, double[] v2, double scale2, int dim, int partDiffOrder, double[] results) {
        checkSupport(dim, partDiffOrder);
        double f1 = v1[0] * scale1;
        double f2 = v2[0] * scale2;
        double sqr = f1 * f1 + f2 * f2;
        double sqrtv = Math.sqrt(sqr);
        results[0] = (f1 + f2) * sqrtv + sqr;


        if (partDiffOrder > 0) {
            if (0 == sqrtv) {
                for (int i = 1; i < dim + 1; i++) {
                    results[i] = 0;
                }
            } else {
                double f_f1 = sqrtv + 2 * f1 + (f1 + f2) * f1 / sqrtv;
                double f_f2 = sqrtv + 2 * f2 + (f1 + f2) * f2 / sqrtv;
                for (int i = 1; i < dim + 1; i++) {
                    results[i] = f_f1 * v1[i] * scale1 + f_f2 * v2[i] * scale2;
                }
            }
        }
    }

    public static void valueOfUnionC2(double[] v1, double scale1, double[] v2, double scale2, int dim, int partDiffOrder, double[] results) {
        checkSupport(dim, partDiffOrder);
        double f1 = v1[0] * scale1;
        double f2 = v2[0] * scale2;
        double sqr = f1 * f1 + f2 * f2;
        double sqrtv = Math.sqrt(sqr);
        results[0] = (f1 + f2 + sqrtv) * sqr;


        if (partDiffOrder > 0) {
            double f_f1 = 3 * f1 * f1 + 2 * f1 * f2 + f2 * f2 + 3 * f1 * sqrtv;
            double f_f2 = f1 * f1 + 2 * f1 * f2 + 3 * f2 * f2 + 3 * f2 * sqrtv;
            for (int i = 1; i < dim + 1; i++) {
                results[i] = f_f1 * v1[i] * scale1 + f_f2 * v2[i] * scale2;
            }
        }
    }

    public static void valueOfIntersection(double[] v1, double scale1, double[] v2, double scale2, int dim, int partDiffOrder, boolean throwSingular, double[] results) {
        checkSupport(dim, partDiffOrder);
        double f1 = v1[0] * scale1;
        double f2 = v2[0] * scale2;
        double sqr = f1 * f1 + f2 * f2;
        double sqrtv = Math.sqrt(sqr);
        results[0] = f1 + f2 - sqrtv;

        if (partDiffOrder > 0) {
            if (f1 == 0 && f2 == 0) {
                if (throwSingular) {
                    throw new ArithmeticException("Singular at the both zero value point!");
                } else {
                    for (int i = 1; i <= dim; i++) {
                        results[i] = 0;
                    }
                }
            } else {
                double f_f1 = 1 - f1 / sqrtv;
                double f_f2 = 1 - f2 / sqrtv;
                for (int i = 1; i < dim + 1; i++) {
                    results[i] = f_f1 * v1[i] * scale1 + f_f2 * v2[i] * scale2;
                }
            }
        }
    }

    public static void valueOfIntersectionC1(double[] v1, double scale1, double[] v2, double scale2, int dim, int partDiffOrder, double[] results) {
        checkSupport(dim, partDiffOrder);
        double f1 = v1[0] * scale1;
        double f2 = v2[0] * scale2;
        double sqr = f1 * f1 + f2 * f2;
        double sqrtv = Math.sqrt(sqr);
        results[0] = (f1 + f2) * sqrtv - sqr;


        if (partDiffOrder > 0) {
            if (0 == sqrtv) {
                for (int i = 1; i < dim + 1; i++) {
                    results[i] = 0;
                }
            } else {
                double f_f1 = sqrtv - 2 * f1 + (f1 + f2) * f1 / sqrtv;
                double f_f2 = sqrtv - 2 * f2 + (f1 + f2) * f2 / sqrtv;
                for (int i = 1; i < dim + 1; i++) {
                    results[i] = f_f1 * v1[i] * scale1 + f_f2 * v2[i] * scale2;
                }
            }
        }
    }

    public static void valueOfIntersectionC2(double[] v1, double scale1, double[] v2, double scale2, int dim, int partDiffOrder, double[] results) {
        checkSupport(dim, partDiffOrder);
        double f1 = v1[0] * scale1;
        double f2 = v2[0] * scale2;
        double sqr = f1 * f1 + f2 * f2;
        double sqrtv = Math.sqrt(sqr);
        results[0] = (f1 + f2 - sqrtv) * sqr;


        if (partDiffOrder > 0) {
            double f_f1 = 3 * f1 * f1 + 2 * f1 * f2 + f2 * f2 - 3 * f1 * sqrtv;
            double f_f2 = f1 * f1 + 2 * f1 * f2 + 3 * f2 * f2 - 3 * f2 * sqrtv;
            for (int i = 1; i < dim + 1; i++) {
                results[i] = f_f1 * v1[i] * scale1 + f_f2 * v2[i] * scale2;
            }
        }
    }

    static abstract class CombineAdapter implements CoordinatePartDiffFunction {

        protected final CoordinatePartDiffFunction fun1, fun2;
        protected final double scale1, scale2;
        protected final int dim;
        protected int diffOrder;
        protected final boolean throwWhenSingular;
        protected int continueM;

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

        protected CombineAdapter(int continueM, boolean throwSingular, CoordinatePartDiffFunction fun1, double scale1, CoordinatePartDiffFunction fun2, double scale2, int dim) {
            this.fun1 = fun1;
            this.scale1 = scale1;
            this.fun2 = fun2;
            this.scale2 = scale2;
            this.dim = dim;
            this.diffOrder = 0;
            this.continueM = continueM;
            this.throwWhenSingular = throwSingular;
            fun1.setDiffOrder(0);
            fun2.setDiffOrder(0);
        }
    }

    static class UnionIntersection extends CombineAdapter {

        private final boolean isUnion;

        public UnionIntersection(boolean isUnion, int continueM, boolean throwSingular, CoordinatePartDiffFunction fun1, double scale1, CoordinatePartDiffFunction fun2, double scale2, int dim) {
            super(continueM, throwSingular, fun1, scale1, fun2, scale2, dim);
            this.isUnion = isUnion;
        }

        @Override
        public double[] values(Coordinate coord, double[] results) {
            double[] res1 = fun1.values(coord, null);
            double[] res2 = fun2.values(coord, null);
            results = initAndCheckSingular(results, res1, res2);
            if (isUnion) {
                switch (continueM) {
                    case 0:
                        valueOfUnion(res1, scale1, res2, scale2, dim, diffOrder, throwWhenSingular, results);
                        break;
                    case 1:
                        valueOfUnionC1(res1, scale1, res2, scale2, dim, diffOrder, results);
                        break;
                    case 2:
                        valueOfUnionC2(res1, scale1, res2, scale2, dim, diffOrder, results);
                        break;
                    default:
                        throw new UnsupportedOperationException("Only supports m=0,1 or 2");
                }
            } else {
                switch (continueM) {
                    case 0:
                        valueOfIntersection(res1, scale1, res2, scale2, dim, diffOrder, throwWhenSingular, results);
                        break;
                    case 1:
                        valueOfIntersectionC1(res1, scale1, res2, scale2, dim, diffOrder, results);
                        break;
                    case 2:
                        valueOfIntersectionC2(res1, scale1, res2, scale2, dim, diffOrder, results);
                        break;
                    default:
                        throw new UnsupportedOperationException("Only supports m=0,1 or 2");
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

    public static class Logistic implements CoordinatePartDiffFunction {

        CoordinatePartDiffFunction oriFun;
        double k;
        double scale;

        public Logistic(CoordinatePartDiffFunction oriFun, double scale, double k) {
            this.oriFun = oriFun;
            this.k = k;
            this.scale = scale;
        }

        @Override
        public double[] values(Coordinate coord, double[] results) {
            results = oriFun.values(coord, results);
            valueOfLogistic(results, results, scale, k, getDiffOrder(), getDim());
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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset;

import net.epsilony.levelset.ops.AtomOperations;
import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.spfun.CommonUtils;
import net.epsilony.utils.geom.Coordinate;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author epsilon
 */
public class AtomOperationsTest {

    public AtomOperationsTest() {
    }

    @Test
    public void testLogistic() {
        CoordinatePartDiffFunction fun = AtomOperations.logisticNorm(fun1, 3);
        double[] samples = new double[]{0, 0, 2, 0, -18, 0, -8, 5, -8, -5};
        fun.setDiffOrder(0);
        double[][] exps = new double[][]{{(3.28457352668571 - 2.2)/2.2},
            {0},
            {0},
            {0},
            {0}};
        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }

        fun.setDiffOrder(1);
        samples = new double[]{1.5, 2.2, -1.5, 2.3};
        exps = new double[][]{{1.88504846666306 - 2.2, -0.614149807509037, -0.568896663797845},
            {3.29924719020606 - 2.2, -0.321896747677589, -0.455607704405203},};
        for (int i=0;i<exps.length;i++){
            for(int j=0;j<3;j++){
                exps[i][j]/=2.2;
            }
        }
        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }
    }

    @Test
    public void testIntersection() {
        CoordinatePartDiffFunction fun = AtomOperations.union_intersection(false,0, fun1, fun2, 2);
        fun.setDiffOrder(0);

        double[] samples = new double[]{2, 0, -10, 0, -4, 4.58257569495584, -4, -4.58257569495584, -4, 0, 10, 0};
        double[][] exps = new double[][]{
            {0},
            {0},
            {0},
            {0},
            {0.4920606076066001},
            {-4.48}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }

        fun.setDiffOrder(1);
        samples = new double[]{1.5, 2.2, -1.5, 2.3, -4, 4.58257569495584, -4, -4.58257569495584,};
        exps = new double[][]{{-0.10196858635932082, -0.213342353582030, -0.198723708576420},
            {0.282985756481431, -0.0710301751086696, -0.122656331947563},
            {0, 0, 0},
            {0, 0, 0}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }
    }

    @Test
    public void testUnion() {
        CoordinatePartDiffFunction fun = AtomOperations.union_intersection(true,0, fun1, fun2, 2);
        fun.setDiffOrder(0);

        double[] samples = new double[]{2, 0, -10, 0, -4, 4.58257569495584, -4, -4.58257569495584, -18, -2, 0.1, 0.2};
        double[][] exps = new double[][]{
            {1.92},
            {1.92},
            {0},
            {0},
            {-0.154672579460335},
            {2.39595405433437}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }

        fun.setDiffOrder(1);
        samples = new double[]{1.5, 2.2, -1.5, 2.3, -4, 4.58257569495584, -4, -4.58257569495584,};
        exps = new double[][]{{1.47756858635932, -0.226657646417970, -0.505276291423580},
            {1.98061424351857, -0.128969824891330, -0.613343668052437},
            {0, 0, 0},
            {0, 0, 0}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }
    }

    @Test
    public void testUnionC1() {
        CoordinatePartDiffFunction fun = AtomOperations.union_intersection(true,1, fun1, fun2, 2);
        fun.setDiffOrder(0);

        double[] samples = new double[]{2, 0, -10, 0, -4, 4.58257569495584, -4, -4.58257569495584, -18, -2, 0.1, 0.2};
        double[][] exps = new double[][]{
            {1.8432},
            {1.8432},
            {0},
            {0},
            {-0.372038196581543},
            {2.52857982524066}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }

        fun.setDiffOrder(1);
        samples = new double[]{1.5, 2.2, -1.5, 2.3, -4, 4.58257569495584, -4, -4.58257569495584,};
        exps = new double[][]{{1.16693725369794, -0.188844218205330, -0.625527575639619},
            {1.68117358081432, -0.166849472163865, -1.00654600511658},
            {0, 0, 0},
            {0, 0, 0}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }
    }

    @Test
    public void testIntersectionC1() {
        CoordinatePartDiffFunction fun = AtomOperations.union_intersection(false,1, fun1, fun2, 2);
        fun.setDiffOrder(0);

        double[] samples = new double[]{2, 0, -10, 0, -4, 4.58257569495584, -4, -4.58257569495584, -4, 0, 10, 0};
        double[][] exps = new double[][]{
            {0},
            {0},
            {0},
            {0},
            {0.584538179220912},
            {-10.0352000000000}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }

        fun.setDiffOrder(1);
        samples = new double[]{1.5, 2.2, -1.5, 2.3, -4, 4.58257569495584, -4, -4.58257569495584,};
        exps = new double[][]{{-0.0805315863020591, -0.167812218205330, -0.141316375639619},
            {0.240202340814316, -0.0684894721638646, -0.173541205116581},
            {0, 0, 0},
            {0, 0, 0}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }
    }

    @Test
    public void testIntersectionC2() {
        CoordinatePartDiffFunction fun = AtomOperations.union_intersection(false,2, fun1, fun2, 2);
        fun.setDiffOrder(0);

        double[] samples = new double[]{2, 0, -10, 0, -4, 4.58257569495584, -4, -4.58257569495584, -4, 0, 10, 0};
        double[][] exps = new double[][]{
            {0},
            {0},
            {0},
            {0},
            {0.694395929454434},
            {-22.478848000000006}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }

        fun.setDiffOrder(1);
        samples = new double[]{1.5, 2.2, -1.5, 2.3, -4, 4.58257569495584, -4, -4.58257569495584,};
        exps = new double[][]{{-0.0636013170710509, -0.131996667518768, -0.0992636513274873},
            {0.203887168209693, -0.0650934592556352, -0.206236370110498},
            {0, 0, 0},
            {0, 0, 0}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }
    }

    @Test
    public void testUnionC2() {
        CoordinatePartDiffFunction fun = AtomOperations.union_intersection(true,2, fun1, fun2, 2);
        fun.setDiffOrder(0);

        double[] samples = new double[]{10, 0, -18, 0, 0, 5, 0, -5, -8, 5, -8, -5, -4, 4.58257569495584, -4, -4.58257569495584, -18, -2, 0.1, 0.2};
        double[][] exps = new double[][]{
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {0},
            {-0.894873675725712},
            {2.66854697027583}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }

        fun.setDiffOrder(1);
        samples = new double[]{1.5, 2.2, -1.5, 2.3};
        exps = new double[][]{{0.921610385223051, -0.156912286881232, -0.672885843712513},
            {1.42700408122231, -0.190327512744365, -1.26683587884950}};

        for (int i = 0; i < exps.length; i++) {
            double[] acts = fun.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }
    }
    CoordinatePartDiffFunction fun1 = new SampleEllipse(-8, 0, 10, 5);
    CoordinatePartDiffFunction fun2 = new SampleEllipse(0, 0, 10, 5);

    public static class SampleEllipse implements CoordinatePartDiffFunction {

        final double xc, yc, a, b;
        private int partOrder;

        public SampleEllipse(double xc, double yc, double a, double b) {
            this.xc = xc;
            this.yc = yc;
            this.a = a;
            this.b = b;
        }

        @Override
        public double[] values(Coordinate coord, double[] results) {
            if (null == results) {
                results = new double[CommonUtils.len2DBase(partOrder)];
            }
            double x = coord.x - xc;
            double y = coord.y - yc;
            double f = 1 - x * x / (a * a) - y * y / (b * b);
            results[0] = f;
            if (partOrder > 0) {
                results[1] = -2 * x / (a * a);
                results[2] = -2 * y / (b * b);
            }
            return results;
        }

        @Override
        public void setDiffOrder(int order) {
            this.partOrder = order;
        }

        @Override
        public int getDiffOrder() {
            return partOrder;
        }

        @Override
        public int getDim() {
            return 2;
        }
    }
}

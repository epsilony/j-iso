/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.ops;

import net.epsilony.levelset.functions.QuadrangleLSF;
import net.epsilony.utils.geom.Coordinate;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author epsilon
 */
public class QuadrangleLSFTest {

    public QuadrangleLSFTest() {
    }

    @Test
    public void testOnlySign() {

        double[] vertesxys = new double[]{-3, -4, 6, -1, 4, 7, -5, 3.5};
        Coordinate[] coords = new Coordinate[4];
        for (int i = 0; i < 4; i++) {
            coords[i] = new Coordinate(vertesxys[i * 2], vertesxys[i * 2 + 1]);
        }

        QuadrangleLSF q = new QuadrangleLSF(coords, 0.3);

        double[] samples = new double[]{1.5, -2.5, 5, 3, -0.5, 5.25};
        double[][] exps = new double[][]{{0},{0},{0}};

        for (int i = 0; i < exps.length; i++) {
            double[] act = q.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], act, 1e-10);
        }

        samples = new double[]{3, -4, -3.1, -3, 3, 7};
        exps = new double[][]{{-1}, {1}, {-1}};

        for (int i = 0; i < exps.length; i++) {
            double[] act = q.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertEquals(Math.signum(exps[i][0]), Math.signum(act[0]), 1e-10);
        }

    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.ops;

import net.epsilony.levelset.functions.HalfPlaneLinearLSF;
import net.epsilony.utils.geom.Coordinate;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author epsilon
 */
public class HalfPlaneLinearLSFTest {

    public HalfPlaneLinearLSFTest() {
    }

    @Test
    public void testSample() {
        double ramp = 0.3;
        double x1 = 1.5, y1 = 0.5, x2 = -4.5, y2 = -8;
        HalfPlaneLinearLSF hp = new HalfPlaneLinearLSF(new Coordinate(x1, y1), new Coordinate(x2, y2), ramp);

        double[] samples = new double[]{-10.5, -16.5, -3, 1, 10, 0.5};
        double exps1[] = new double[]{0, -1.1894090950472263, 2.0832680513251418};
        double exps2[] = new double[]{0.81696786326476156*ramp};
        double exps3[] = new double[]{-0.57668319759865527*ramp};

        hp.setDiffOrder(1);

        for (int i = 0; i < samples.length / 2; i++) {
            double[] acts = hp.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertEquals(exps1[i], acts[0], 1e-10);
            assertEquals(exps2[0], acts[1], 1e-10);
            assertEquals(exps3[0], acts[2], 1e-10);
        }
    }
}

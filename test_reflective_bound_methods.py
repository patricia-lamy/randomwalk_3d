import unittest
import numpy as np
from reflective_bc_functions import sph2cart, reflected_point_method, cart2sph
from reflective_bc_functions import collision_point_method
from reflective_bc_functions import distance_between_2points


class TestingReflectiveBoundaryMethods(unittest.TestCase):

    def test_valid_coordinates(self):
        # Add test to validate that new coordinates of points are within
        # expected range
        pass

    def test_total_distance(self):
        np.random.seed(41)
        r = 1.02
        step_length = 0.1
        p1_sph = (1.0199, 0, 0)
        p1_cart = sph2cart(p1_sph[0], p1_sph[1], p1_sph[2])
        d_theta = (np.random.uniform() * (2 * np.pi)) + 1e-14
        d_phi = (np.random.uniform() * (np.pi)) + 1e-14
        p2_cart = (p1_cart[0] + step_length*np.sin(d_phi)*np.cos(d_theta),
                   p1_cart[1] + step_length*np.sin(d_phi)*np.sin(d_theta),
                   p1_cart[2] + step_length*np.cos(d_phi)
                   )
        p2_sph = cart2sph(p2_cart[0], p2_cart[1], p2_cart[2])
        self.assertTrue(p2_sph[0] > r)
        pc = collision_point_method(p1_cart, p2_cart, r)
        p3_cart = reflected_point_method(p1_cart, p2_cart, r)[0]
        dist = (distance_between_2points(p1_cart, pc)
                + distance_between_2points(pc, p3_cart))
        self.assertAlmostEqual(dist, step_length)


if __name__ == '__main__':
    unittest.main()

import unittest
import numpy as np
from reflective_bc_functions import sph2cart, reflected_point_method_sph, cart2sph
from reflective_bc_functions import collision_point_method


class TestingReflectiveBoundaryMethods(unittest.TestCase):
    def total_distance(self):
        np.random.seed(42)
        r = 1.02
        p1 = (1.0199, 0, 0)
        d_theta = (np.random.uniform() * (2 * np.pi)) + 1e-14
        d_phi = (np.random.uniform() * (np.pi)) + 1e-14
        p2 = (p1[0]*np.sin(p1[2])*np.cos(p1[1]), p1[0]*np.sin(p1[2])*np.sin(p1[1]), p1[0]*np.cos(p1[2]))
        print(p2)
        p2 = cart2sph(p1[0]*np.sin(p1[2])*np.cos(p1[1]), p1[0]*np.sin(p1[2])*np.sin(p1[1]), p1[0]*np.cos(p1[2]))
        print(p2)
        pc = collision_point_method(p1, p2, r)
        p3 = reflected_point_method_sph(p1, p2, r)



if __name__ == '__main__':
    unittest.main()

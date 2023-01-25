from reflective_bc_functions import reflected_point_method_sph
import sys


def reflect(initial_pt, end_of_step_point, radius):
    """Returns the reflected point, using the initial point and the end of
    step point that crosses the sphere radius.

    Args:
        initial_pt (tuple): Spherical coordinates of initial point, in the
            form (rho, theta, phi).
        end_of_step_pt (tuple): Spherical coordinates of the end of step point,
            in the form (rho, theta, phi).

    Returns:
        Spherical coordinates of the reflected point, with the certainty that
        the reflected point is inside the sphere.
    """
    [rho1, theta1, phi1] = [initial_pt[0], initial_pt[1], initial_pt[2]]
    [rho2, theta2, phi2] = [end_of_step_point[0], end_of_step_point[1],
                            end_of_step_point[2]]
    reflected_pt = reflected_point_method_sph((rho1, theta1, phi1),
                                              (rho2, theta2, phi2), radius)[0]
    [rho3, theta3, phi3] = [reflected_pt[0], reflected_pt[1], reflected_pt[2]]

    # If reflected point not inside the boundary (the sphere), we place it on
    # sphere
    if rho3 >= radius:
        if rho3 > 1.03:
            print('Initial pt:', [rho1, theta1, phi1])
            print('End of step pt:', [rho2, theta2, phi2])
            print('Reflected point:', [rho3, theta3, phi3])
            sys.exit('Reflected point farther than what it should')
        else:
            rho3 = radius
    return [rho3, theta3, phi3]

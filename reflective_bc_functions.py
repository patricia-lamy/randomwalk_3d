import numpy as np
import sys


def sph2cart(rho, theta, phi):
    """Converts a point from spherical to cartesian coordinates.

    Args:
        rho (float): Radial distance (>= 0)
        theta (float): Azimuthal angle (in [0, 2pi])
        phi (float): Polar angle (in [0, pi])

    Returns:
        The cartesian coordinates (x, y, z)
    """
    x = rho*np.sin(phi)*np.cos(theta)
    y = rho*np.sin(phi)*np.sin(theta)
    z = rho*np.cos(phi)
    return (x, y, z)


def cart2sph(x, y, z):
    """Converts a point from cartesian to spherical coordinates.

    Args:
        The cartesian coordinates (x, y, z)

    Returns:
        The point (rho, theta, phi), where
            rho (float): Radial distance (>= 0)
            theta (float): Azimuthal angle (in [0, 2pi])
            phi (float): Polar angle (in [0, pi])
    """
    rho = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(y, (x+1e-16)) + np.pi
    phi = np.arccos((z + 1e-16)/(rho + 1e-16))
    if rho < 0 or theta < 0 or phi < 0:
        print('Rho:', rho)
        print('Theta:', theta)
        print('Phi:', phi)
        sys.exit("Negative value of rho or theta or phi")
    return (rho, theta, phi)


def collision_point_method(initial_pt, end_of_step_pt, sphere_radius):
    """Finds the collision point of a particle on the sphere, knowing that the
    end of the step of the random walk would have brought it to outside of the
    boundary.

    Args:
        inital_pt (tuple): Contains cartesian coordinates of the initial point
            of the step, in the form (x, y, z).
        end_of_step_pt (tuple): Contains cartesian coordinates of the point
            where the particle would have landed after the step, if it had not
            been blocked by the reflective sphere, in the form (x, y, z).

    Returns:
        Collision point of the particle, in cartesian coordinates, of the form
        (x, y, z)
    """
    (x1, y1, z1) = (initial_pt[0], initial_pt[1], initial_pt[2])
    (x2, y2, z2) = (end_of_step_pt[0], end_of_step_pt[1], end_of_step_pt[2])
    (rho1, theta1, phi1) = cart2sph(x1, y1, z1)
    alpha = 2*(x1*(x2 - x1) + y1*(y2 - y1) + z1*(z2 - z1))
    beta = (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2
    potential_d1 = ((-alpha + np.sqrt(alpha**2 -
                     4*beta*(rho1**2 - sphere_radius**2))) / (2*beta))
    potential_d2 = ((-alpha - np.sqrt(alpha**2 -
                     4*beta*(rho1**2 - sphere_radius**2))) / (2*beta))
    if np.abs(potential_d1) < np.abs(potential_d2):
        d = potential_d1
    else:
        d = potential_d2
    x_c = x1 + d*(x2 - x1)
    y_c = y1 + d*(y2 - y1)
    z_c = z1 + d*(z2 - z1)
    return (x_c, y_c, z_c)


def projected_point_method(initial_pt, collision_pt, sphere_radius):
    """Let v be the vector formed by the initial position point p1 and the
    collision point on the sphere pc. Then the function projects the vector v
    onto the tangent plane of the sphere at the collision point pc.

    Args:
        initial_pt (tuple): Contains cartesian coordinates of initial point in
        the form (x, y, z)
        collision_pt (tuple): Contains cartesian coordinates of collision
        point, in the form (x, y, z)

    Returns:
        p3 (tuple): Point at the end of the projected vector starting at pc,
        which is on the tangent plane of the sphere at pc, in cartesian
        coordinates (of the form (x,y,z))
    """
    (x1, y1, z1) = (initial_pt[0], initial_pt[1], initial_pt[2])
    (xc, yc, zc) = (collision_pt[0], collision_pt[1], collision_pt[2])
    # We first determine the normal of the tangent plane (normalised)
    pc_norm = sphere_radius
    n = (xc/pc_norm, yc/pc_norm, zc/pc_norm)
    # We calculate the scalar product of the vector p1pc with the normal n
    vec_p1pc = (xc - x1, yc - y1, zc - z1)
    p1pc_dot_n = np.dot(vec_p1pc, n)
    # We now determine the point p3 forming of the projected vector v with pc
    v = (vec_p1pc[0] - n[0]*p1pc_dot_n, vec_p1pc[1] - n[1]*p1pc_dot_n,
         vec_p1pc[2] - n[2]*p1pc_dot_n)
    p3 = (xc - v[0], yc - v[1], zc - v[2])
    return p3


def distance_between_2points(p1, p2):
    """Calculates the distance between points p1 and p2.

    Args:
        p1 (tuple): Contains cartesian coordinates of p1, in the form (x, y, z)
        p2 (tuple): Contains cartesian coordinates of p2, in the form (x, y, z)

    Returns:
        d (float): Distance between p1 and p2.
    """
    d = np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)
    return d


def initial_point_image(initial_pt, projected_pt):
    """Calculates the coordinates of the image of the initial point across the
    tangent plane of the sphere at the collision point.

    Args:
        initial_pt (tuple): Contains cartesian coordinates of the initial
            point, in the form (x, y, z).
        projected_pt (tuple): Contains cartesian coordinates of the initial
            point projected onto the tangent plane, in the form (x, y, z).

    Returns:
        image_pt (tuple): Contains cartesian coordinates of the image of the
            initial point across the tangent plane, in the form (x, y, z).
    """
    t = 2
    image_pt_x = initial_pt[0] + t*(projected_pt[0] - initial_pt[0])
    image_pt_y = initial_pt[1] + t*(projected_pt[1] - initial_pt[1])
    image_pt_z = initial_pt[2] + t*(projected_pt[2] - initial_pt[2])
    return (image_pt_x, image_pt_y, image_pt_z)


def reflected_point_method(initial_pt, end_of_step_pt, sphere_radius):
    """Method that calculates the reflected point for a particle that bounces
    of the sphere.

    Args:
        inital_pt (tuple): Contains cartesian coordinates of the initial point
            of the step, in the form (x, y, z).
        end_of_step_pt (tuple): Contains cartesian coordinates of the point
            where the particle would have landed after the step, if it had not
            been blocked by the reflective sphere, in the form (x, y, z).
        sphere_radius (float): Sphere radius defining the boundary condition.

    Returns: an array containing
        reflected_pt (tuple): Contains the cartesian coordinates of the point
            where the reflected particle stops, in the form (x, y, z).
        collision_pt (tuple): Contains the cartesian coordinates of the
            collision point.

    """
    step_length = distance_between_2points(initial_pt, end_of_step_pt)
    collision_pt = collision_point_method(initial_pt,
                                          end_of_step_pt,
                                          sphere_radius)
    projected_pt = projected_point_method(initial_pt,
                                          collision_pt,
                                          sphere_radius)
    image_pt = initial_point_image(initial_pt, projected_pt)
    dist_imagept_collpt = distance_between_2points(image_pt, collision_pt)
    translation_factor = step_length/(dist_imagept_collpt+1e-16)
    ref_x = image_pt[0] + translation_factor*(collision_pt[0] - image_pt[0])
    ref_y = image_pt[1] + translation_factor*(collision_pt[1] - image_pt[1])
    ref_z = image_pt[2] + translation_factor*(collision_pt[2] - image_pt[2])
    reflected_pt = (ref_x, ref_y, ref_z)
    return [reflected_pt, collision_pt]


def reflected_point_method_sph(initial_pt, end_of_step_pt, sphere_radius):
    """Method that yields the spherical coordinates point of the reflected
    particle, using reflected_point_method.

    Args:
        inital_pt (tuple): Contains spherical coordinates of the initial point
            of the step, in the form (rho, theta, phi).
        end_of_step_pt (tuple): Contains spherical coordinates of the point
            where the particle would have landed after the step, if it had not
            been blocked by the reflective sphere, in the form (rho, theta,
            phi).
        sphere_radius (float): Sphere radius defining the boundary condition.

    Returns: an array containing
        reflected_pt (tuple): Contains the spherical coordinates of the point
            where the reflected particle stops, in the form (rho, theta, phi).
        collision_pt (tuple): Contains spherical coordinates of collision
            point that was part of the reflection.

    """
    initial_pt_cart = sph2cart(initial_pt[0], initial_pt[1], initial_pt[2])
    end_of_step_pt_cart = sph2cart(end_of_step_pt[0], end_of_step_pt[1],
                                   end_of_step_pt[2])
    [ref_pt_cart, coll_pt_cart] = reflected_point_method(initial_pt_cart,
                                                         end_of_step_pt_cart,
                                                         sphere_radius)
    ref_pt_sph = cart2sph(ref_pt_cart[0], ref_pt_cart[1], ref_pt_cart[2])
    coll_pt_sph = cart2sph(coll_pt_cart[0], coll_pt_cart[1], coll_pt_cart[2])
    return [ref_pt_sph, coll_pt_sph]

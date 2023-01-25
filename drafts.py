        
# Part to do second reflection:
        coll_pt = reflected_point_method_sph((rho1, theta1, phi1),
                                             (rho2, theta2, phi2), radius)[1]
        [rho_c, theta_c, phi_c] = [coll_pt[0], coll_pt[1], coll_pt[2]]
        reflected_pt2 = reflected_point_method_sph((rho_c, theta_c, phi_c),
                                                   (rho3, theta3, phi3),
                                                   radius)[0]
        [rho4, theta4, phi4] = [reflected_pt2[0], reflected_pt2[1],
                                reflected_pt2[2]]
        if rho4 >= radius:
            print('Position after second reflection:', rho4)
        [rho_r, theta_r, phi_r] = [radius, theta4, phi4]
    else:
        [rho_r, theta_r, phi_r] = [rho3, theta3, phi3]
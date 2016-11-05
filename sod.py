import numpy as np
import matplotlib.pyplot as plt
import copy


def kernel_cubic_spline( r ) :
    q = abs( r ) / h
    if q < 1.0 :
        w = (2.0 / (3.0 * h)) * (1.0 - 1.5 * q**2 * (1.0 - q / 2.0))
    elif q < 2.0 :
        w = (1.0 / (6.0 * h)) * (2.0 - q)**3
    else :
        w = 0.0
    return w
    
def derivative_kernel_cubic_spline( r ) :
    q = abs( r ) / h
    if 0.0 < q < 1.0 :
        dw = (2.0 / h) * (-q + 0.75 * q**2) * abs( r ) / r
    elif 1.0 <= q < 2.0 :
        dw = (-1.0 / (2.0 * h)) * (2.0 - q)**2 * abs( r ) / r
    else : 
        dw = 0.0
    return dw / h

def initialise() :
    x_l = np.linspace( -0.5, 0.0, 320, endpoint = False )
    p_l = np.full( 320, 1.0 )
    rho_l = np.full( 320, 1.0 )
    x_r = np.linspace( 0.0, 0.5, 41)
    x_r = x_r[ 1 : ]
    p_r = np.full( 40, 0.1 )
    rho_r = np.full( 40, 0.125 )
    x = np.concatenate( [ x_l, x_r ] )
    p = np.concatenate( [ p_l, p_r ] )
    rho = np.concatenate( [ rho_l, rho_r ] )
    velocity = np.full( 360, 0.0 )
    energy = p / ( rho * ( gamma - 1.0 ) )
    return x, rho, velocity, p, energy

def get_neighbourhood( x, index ) :
    low = max( x[index] - 2.0000000001 * h, min( x ) )
    high = min( x[index] + 2.0000000001 * h, max( x ) )
    neighbourhood = []
    index1 = copy.deepcopy( index )
    while True:
        neighbourhood.append( x[index] )
        if index == 0 :
            break
        index -= 1
        if x[index] < low :
            break
    neighbourhood = neighbourhood[ : : -1 ]
    while True:
        if index1 == len( x ) - 1 :
            break
        index1 += 1
        if x[index1] > high :
            break
        neighbourhood.append( x[index1] )
    delta_x_l = x[1] - x[0]
    return np.array( neighbourhood ), np.where( x == neighbourhood[0] )[0][0]

def artificial_viscosity( p_a, p_b, rho_a, rho_b, x_a, x_b, v_a, v_b ) :
    v_ab = v_a - v_b
    x_ab = x_a - x_b
    approach = v_ab * x_ab
    if approach >= 0.0 :
        return 0.0
    else :
        mu_ab = h * approach / ( x_ab * x_ab + 0.01 * h * h )
        c_a = np.sqrt( gamma * p_a / rho_a )
        c_b = np.sqrt( gamma * p_b / rho_b )
        c_ab = (c_a + c_b ) / 2.0
        rho_ab = ( rho_a + rho_b ) / 2.0
        pi_ab = ( -alpha * c_ab * mu_ab + beta * mu_ab * mu_ab ) / rho_ab
        return pi_ab

def get_derivatives( p, rho, velocity, x ) :
    rho_rate = np.zeros_like( x )
    v_rate = np.zeros_like( x )
    e_rate = np.zeros_like( x )
    x_sph = np.zeros_like( x )
    for i, x_i in enumerate( x ) :
        x_near, index = get_neighbourhood( x, i )
        for j, x_j in enumerate( x_near ) :
            r = x_i - x_j
            j_ind = j + index
            v_ij = velocity[i] - velocity[j_ind]
            rho_rate[i] += rho[i] * mass * v_ij * derivative_kernel_cubic_spline( r ) / rho[j_ind]
            pi_ab = artificial_viscosity( p[i], p[j_ind], rho[i], rho[j_ind], 
                x[i], x[j_ind], velocity[i], velocity[j_ind] )
            v_rate[i] -= mass * ( p[i] / ( rho[i] * rho[i] ) + p[j_ind] / 
                ( rho[j_ind] * rho[j_ind] ) + pi_ab ) * derivative_kernel_cubic_spline( r )
            e_rate[i] += 0.5 * mass * ( p[i] / ( rho[i] * rho[i] ) + p[j_ind] / 
                ( rho[j_ind] * rho[j_ind] ) + pi_ab ) * v_ij * derivative_kernel_cubic_spline( r )
            
            rho_av = ( rho[i] + rho[j_ind] ) / 2.0
            x_sph[i] -= 0.5 * v_ij * kernel_cubic_spline( r ) * mass / rho_av
    position_rate = x_sph + velocity
    return v_rate, e_rate, position_rate

def get_density( x ) :
    rho = np.zeros_like( x )
    for i, x_i in enumerate( x ) :
        x_near, index = get_neighbourhood( x, i )
        for j, x_j in enumerate( x_near ) :
            r = x_i - x_j
            rho[i] += mass * kernel_cubic_spline( r )
    return rho

def sod() :
    x, rho, v, p, e = initialise()
    t = 0.0
    while t <= total_t :
        print t
        v_rate, e_rate, x_rate = get_derivatives( p, rho, v, x )
        rho = get_density( x )
        v += v_rate * dt
        e += e_rate * dt
        x += x_rate * dt
        p = ( gamma - 1.0 ) * rho * e
        t += dt
    return rho, v, p, e, x

def exact_solution( p1, p5, rho1, rho5 ) :
    p1 = 1.0
    p5 = 0.1
    rho1 = 1.0
    rho5 = 0.125
    gamma = 1.4
    sigma = ( gamma - 1.0 ) / ( 1.0 + gamma )
    b = ( gamma - 1.0 ) / ( 2.0 * gamma )
    residue = []
    p3_check = np.linspace( p5, p1, 10001 )
    for p3 in p3_check :
        u4 = ( p3 - p5 ) * np.sqrt( ( 1 - sigma ) / ( rho5 * ( p3 + sigma * p5 ) ) )
        u2 = ( p1**b - p3**b ) * np.sqrt( ( 1.0 - sigma * sigma ) * p1**(1.0/gamma) / rho1 ) / sigma
        residue.append( abs( u2 - u4 ) )
    index = np.argmin( residue )
    p3 = p3_check[index]
    u5 = np.sqrt( gamma * p5 / rho5 )
    u1 = np.sqrt( gamma * p1 / rho1 )
    u3 = u5 + ( p3 - p5 ) / ( np.sqrt( rho5 * 0.5 * ( ( gamma + 1.0 ) * p3 + ( gamma - 1.0 ) * p5 ) ) )
    u4 = u3
    rho3 = rho1 * ( p3 / p1 )**(1.0/gamma)
    p4 = p3
    rho4 = rho5 * ( p4 + sigma * p5 ) / ( p5 + sigma * p4 )
    x = np.array( [ -0.5, -0.25, 0.0, 0.199, 0.2, 0.399, 0.4, 0.5 ] )
    p = np.array( [ p1, p1, p3, p3, p4, p4, p5, p5 ] )
    v = np.array( [ u1, u1, u3, u3, u4, u4, u5, u5 ] )
    rho = np.array( [ rho1, rho1, rho3, rho3, rho4, rho4, rho5, rho5 ] )
    e = p / ( ( gamma - 1.0 ) * rho )
    return rho, v, p, e, x

def plot() :
    rho, v, p, e, x = sod()
    rho_exact, v_exact, p_exact, e_exact, x_exact = exact_solution( 1.0, 0.1, 1.0, 0.125 )
    plt.figure( figsize = ( 17.0, 10.0 ) )
    plt.plot( x, rho, label = 'SPH' )
    plt.plot( x_exact, rho_exact, 'r', label = 'Exact solution' )
    plt.legend()
    plt.title( "Density variation across SOD shock tube", fontsize = 28, fontweight = 'bold' )
    plt.ylabel( "Density, $\\rho$ in $\\frac{kg}{m^3}$", fontweight = 'bold', fontsize = 24 )
    plt.xlabel( "x", fontweight = 'bold', fontsize = 24 )
    plt.savefig( "density.png" )
    plt.close()
    plt.clf()
    plt.figure( figsize = ( 17.0, 10.0 ) )
    plt.plot( x, p, label = 'SPH' )
    plt.plot( x_exact, p_exact, 'r', label = 'Exact solution' )
    plt.legend()
    plt.title( "Pressure variation across SOD shock tube", fontsize = 28, fontweight = 'bold' )
    plt.ylabel( "Pressure, P in $\\frac{N}{m^2}$", fontweight = 'bold', fontsize = 24 )
    plt.xlabel( "x", fontweight = 'bold', fontsize = 24 )
    plt.savefig( "pressure.png" )
    plt.close()
    plt.clf()
    plt.figure( figsize = ( 17.0, 10.0 ) )
    plt.plot( x, e, label = 'SPH' )
    plt.plot( x_exact, e_exact, 'r', label = 'Exact solution' )
    plt.legend()
    plt.title( "Energy variation across SOD shock tube", fontsize = 28, fontweight = 'bold' )
    plt.ylabel( "Energy, E in $\\frac{m^2}{s^2}$", fontweight = 'bold', fontsize = 24 )
    plt.xlabel( "x", fontweight = 'bold', fontsize = 24 )
    plt.savefig( "energy.png" )
    plt.close()
    plt.clf()    
    plt.figure( figsize = ( 17.0, 10.0 ) )
    plt.plot( x, v, label = 'SPH' )
    plt.plot( x_exact, v_exact, 'r', label = 'Exact solution' )
    plt.legend()
    plt.title( "Velocity variation across SOD shock tube", fontsize = 28, fontweight = 'bold' )
    plt.ylabel( "Velocity, U in $\\frac{m}{s}$", fontweight = 'bold', fontsize = 24 )
    plt.xlabel( "x", fontweight = 'bold', fontsize = 24 )
    plt.savefig( "velocity.png" )
    plt.close()
    plt.clf()

if __name__ == '__main__':
    mass = 1.0 * 0.0015625 # Mass = density_right * delta_x_right = density_left * delta_x_left
    h = 2.0 * 0.0125
    gamma = 1.4
    dt = 1.0e-4
    total_t = 0.20000001
    alpha = 1.0
    beta = 1.0
    plt.ioff()
    plot()
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
    if approach > 0.0 :
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
    return rho_rate, v_rate, e_rate, position_rate

def sod() :
    x, rho, v, p, e = initialise()
    t = 0.0
    ind = 0
    while t <= total_t :
        print t
        rho_rate, v_rate, e_rate, x_rate = get_derivatives( p, rho, v, x )
        rho += rho_rate * dt
        v += v_rate * dt
        e += e_rate * dt
        x += x_rate * dt
        p = ( gamma - 1.0 ) * rho * e
        if np.isnan(np.prod(x)) :
            print 'x is out of bound'
        if np.isnan(np.prod(v)) :
            print 'v is out of bound'
        if np.prod(rho) == 0.0 :
            print 'density is 0'
        t += dt
        if ind%100 == 0 :
            plt.figure( figsize = ( 17.0, 10.0 ) )
            plt.plot( x, rho, '*' )
            plt.savefig( str( ind ) + '.png' )
            plt.close()
        ind += 1
    return p, x, rho, v, e

        

if __name__ == '__main__':
    mass = 1.0 * 0.0015625 # Mass = density_right * delta_x_right = density_left * delta_x_left
    h = 2.0 * 0.0125
    gamma = 1.4
    dt = 1.0e-4
    total_t = 0.200001
    alpha = 1.0
    beta = 1.0
    plt.ioff()
    plt.clf()
    p, x, rho, vel, e = sod()
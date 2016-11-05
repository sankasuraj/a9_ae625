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
    low = max( x[index] - 2.00000001 * h, min( x ) )
    high = min( x[index] + 2.00000001 * h, max( x ) )
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
    if index != 0 :
        index += 1
    return np.array( neighbourhood ), index

def get_density_rate( rho, velocity, x ) :
    density_rate = np.zeros_like( x )
    for i, x_i in enumerate( x ) :
        x_near, index = get_neighbourhood( x, i )
        for j, x_j in enumerate( x_near ) :
            r = x_i - x_j
            v_ij = ( velocity[i] + velocity[j + index] ) / 2.0
            density_rate[i] += rho[i] * v_ij * derivative_kernel_cubic_spline( r ) * mass / rho[j + index]
    return rho

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

def get_velocity_rate( rho, p, velocity, x ) :
    velocity_rate = np.zeros_like( x )
    for i, x_i in enumerate( x ) :
        x_near, index = get_neighbourhood( x, i )
        for j, x_j in enumerate( x_near ) :
            r = x_i - x_j
            pi = artificial_viscosity( p[i], p[j + index], rho[i], rho[j + index], x_i, x_j, velocity[i], velocity[j + index] )
            velocity_rate[i] -= mass * ( p[i] / ( rho[i] * rho[i] ) + p[j + index] / ( rho[j + index] * rho[j + index] ) + pi ) * derivative_kernel_cubic_spline( r )
    return velocity_rate

def get_energy_rate( rho, p, velocity, x ) :
    energy_rate = np.zeros_like( x )
    for i, x_i in enumerate( x ) :
        x_near, index = get_neighbourhood( x, i )
        for j, x_j in enumerate( x_near ) :
            r = x_i - x_j
            pi = artificial_viscosity( p[i], p[j + index], rho[i], rho[j + index], x_i, x_j, velocity[i], velocity[j + index] )
            v_ij = ( velocity[i] + velocity[j] ) / 2.0
            energy_rate[i] += 0.5 * mass * v_ij * ( p[i] / ( rho[i] * rho[i] ) + p[j + index] / ( rho[j + index] * rho[j + index] ) + pi ) * derivative_kernel_cubic_spline( r )    
    return energy_rate

def get_position_rate( rho, velocity, x ) :
    position_rate = np.zeros_like( x )
    for i, x_i in enumerate( x ) :
        x_near, index = get_neighbourhood( x, i )
        for j, x_j in enumerate( x_near ) :
            r = x_i - x_j
            xsph = 0.5 * ( velocity[j + index] - velocity[i] ) * kernel_cubic_spline( r ) * mass / rho[j + index]
            position_rate[i] = position_rate[i] + xsph + velocity[j + index] * kernel_cubic_spline( r ) * mass / rho[j + index]
    return position_rate

def euler_integrator( quantity, rate, time_step ) :
    new_quantity = quantity + rate * time_step
    return new_quantity

def sod() :
    x_old, rho_old, v_old, p_old, e_old = initialise()
    t = 0.0
    ind = 0
    while t <= total_t :
        print t
        rho_rate = get_density_rate( rho_old, v_old, x_old )
        rho_new = euler_integrator( rho_old, rho_rate, dt )
        v_rate = get_velocity_rate( rho_old, p_old, v_old, x_old )
        v_new = euler_integrator( v_old, v_rate, dt )
        e_rate = get_energy_rate( rho_old, p_old, v_old, x_old )
        e_new = euler_integrator( e_old, e_rate, dt )
        position_rate = get_position_rate( rho_old, v_old, x_old )
        x_new = euler_integrator( x_old, position_rate, dt )
        p_new = ( gamma - 1.0 ) * rho_new * e_new
        x_old = copy.deepcopy( x_new )
        rho_old = copy.deepcopy( rho_new )
        p_old = copy.deepcopy( p_new )
        e_old = copy.deepcopy( e_new )
        v_old = copy.deepcopy( v_new )
        if np.isnan(np.prod(x_old)) :
            print 'x is out of bound'
        if np.isnan(np.prod(v_old)) :
            print 'v is out of bound'
        if np.prod(rho_old) == 0.0 :
            print 'density is 0'
        t += dt
        if ind%200 == 0 :
            plt.figure( figsize = ( 17.0, 10.0 ) )
            plt.plot( x_old, rho_old )
            plt.savefig( str( t ) + '.png' )
            plt.close()
        ind += 1
    return p_old, x_old, rho_old, v_old, e_old
    

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
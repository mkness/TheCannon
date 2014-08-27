def main():
    nobservations = 4
    a, b, c = 3.0, 2.0, 1.0
    f, x, y, z = generate_data(nobservations, a, b, c)

#    print 'Linear results (should be {}, {}, {}):'.format(a, b, c)
#    print linear_invert(f, x, y, z)
#
#    print 'Non-linear results (should be {}, {}, {}):'.format(a, b, c)
#    print nonlinear_invert(f, x, y, z)

def generate_data(nobservations, a, b, c, noise_level=0.01):
    x, y, z = np.random.random((3, nobservations))
    noise = noise_level * np.random.normal(0, noise_level, nobservations)
    f = func(x, y, z, a, b, c) + noise
    return f, x, y, z

def func(x, y, z, a, b, c):
    f = (a**2
         + x * b**2
         + y * a * b * np.cos(c)
         + z * a * b * np.sin(c))
    return f

def nonlinear_invert(f, x, y, z):
    def wrapped_func(observation_points, a, b, c):
        x, y, z = observation_points
        return func(x, y, z, a, b, c)

    xdata = np.vstack([x, y, z])
    model, cov = opt.curve_fit(wrapped_func, xdata, f)
    return model

main() 
#nonlinear_invert(0.98,1,2,3) 

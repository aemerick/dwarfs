import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt
from astropy import units as u
from astropy import constants as const




class parabolic:

    def __init__(self, M_v = 2.0E10*u.Msun, M_p = 5.32E10*u.Msun,
                       b   = 20.0*u.kpc, t_o = -1250*u.Myr,
                       alpha = 0.0, dt = 2.0*u.Myr,
                       t_final = None):

        self.M_v = M_v
        self.M_p = M_p
        self.b = b
        self.t_o = t_o
        self.alpha = alpha
        self.dt  = dt

        # To Do: compute dynamical time and just set
        #        t_o to some fixed factor of t_dyn

        if t_final is None:
            self.t_final = -1.0 * self.t_o
        else:
            self.t_final = t_final

        return

    def pos(self, t, cartesian = True):

        xi  = self.xi(t)
        phi = 2.0 * np.arctan(xi)

        r = (self.b * (1.0 + xi**2.0)).to(u.kpc)
        
        if not cartesian:
            return np.array([r, phi])
        
        # cartesian positions in orbital plane
        x_o = r * np.cos(phi)
        y_o = r * np.sin(phi)
        z_o = 0.0

        # now in "victim" frame, assuming y-axis in orbital
        # "victim" frame are aligned, therefore that euler angle theta
        # is zero, and that euler angle phi and euler angle psi are 
        # both equal to the inclination angle (alpha)
 
        alpha_n = self.alpha #+ np.pi / 2.0  # renormalize angle to be correct

        x =  x_o * np.cos(alpha_n) + y_o * np.sin(alpha_n)
        y = -x_o * np.sin(alpha_n) * np.cos(alpha_n) + y_o * np.cos(alpha_n) * np.cos(alpha_n) + z_o * np.sin(alpha_n)
        z =  x_o * np.sin(alpha_n) * np.sin(alpha_n) - y_o * np.sin(alpha_n) * np.cos(alpha_n) + z_o * np.cos(alpha_n)

        
        # x = r * np.cos(phi) * np.sin(self.alpha + np.pi/2.0)
        # y = r * np.sin(phi) * np.sin(self.alpha + np.pi/2.0)
        # z = r * np.cos(self.alpha + np.pi/2.0)


        return np.array([x.value,y.value,z.value]) * r.unit

    def vmag(self, t = None, r = None):

        if r is None: 
            r, phi = self.pos(t, cartesian = False)

            return (self.v_o() * (self.b / r)**(0.5)).to(u.km/u.s)
        elif t is None:
            return (self.v_o() * (self.b / r)**(0.5)).to(u.km/u.s)
        else:
            print "Need to supply a time or position to compute velocity"
            raise ValueError

        return

    def v_o(self):

        return ((2.0 * const.G * (self.M_p + self.M_v) / self.b)**(0.5)).to(u.km/u.s)

    def xi(self, t):

        t = (t.to(u.Myr)).value
        b = ((self.b).to(u.kpc)).value
        v_o = ((self.v_o()).to(u.kpc/u.Myr)).value
        

        f = lambda x, t: (x + x*x*x/3.0)*(2.0*b/v_o) - t

        return (opt.brentq(f, -10000, 10000, args=(t,)))


    def save_orbit(self):

    # compute the orbit in spacings of dt from t_o to t_fi

        f = open('parabolic_orbit.txt', 'w')

        f.write("# Orbital Properties:\n")

        m_v = self.M_v.to('Msun').value
        m_p = self.M_p.to('Msun').value
        b   = self.b.to('kpc').value
        alpha = self.alpha * 180 / (np.pi)
        v_o   = (self.v_o()).to('km/s').value

        f.write("# M_v = %4.4E Msun  -  M_p = %4.4E Msun\n"%(m_v, m_p))
        f.write("# b = %4.4E kpc  -  alpha = %4.4E deg-  v_o = %4.4E km/s \n"%(b, alpha, v_o))
        f.write("# --------------------------------------------------------------------- \n")
        f.write("# t x y z\n")

        t = self.t_o

        norm = -1.0 * np.sign(self.t_o.value) * np.abs(self.t_o)

        while t < self.t_final:
            x, y, z = self.pos(t, cartesian=True)
            x = x.value; y = y.value; z = z.value

            t_write = t + norm

            f.write("%5.5E %8.8E %8.8E %8.8E\n"%(t_write.value, x, y, z))
            t = t + self.dt

        f.close()


def plot_orbit():

    

    return


if __name__ == "__main__":
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt


    print 'running test orbit'

    orbit = parabolic()

    orbit.save_orbit()

    data = np.genfromtxt('parabolic_orbit.txt', names = True, skip_header = 4)


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    ax.scatter(data['x'], data['y'], data['z'],  c = plt.cm.viridis(data['t']/np.max(data['t'])))
    
    npoints = 1.0E3
    r = np.random.rand(npoints) * 10
    theta = np.random.rand(npoints) * 2.0 * np.pi

    ax.scatter(r*np.cos(theta), r * np.sin(theta), np.zeros(npoints), c = 'black')
    ax.view_init(elev = 0.0, azim = 0.0)

    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    fig.savefig('parabolic_orbit.png')

    i = 0
    for x in [90, 180, 270]:
        ax.view_init(elev = 0.0, azim = x)
        fig.savefig('parabolic_orbit_%i.png'%(i))
        i = i + 1

    for x in [45, 90, 135, 180]:
        ax.view_init(elev = x, azim = 0.0)
        fig.savefig('parabolic_orbit_%i.png'%(i))
        i = i + 1

	


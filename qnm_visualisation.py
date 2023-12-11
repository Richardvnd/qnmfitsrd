
import numpy as np 
import matplotlib.pyplot as plt 
import qnmfitsrd.qnm as qnm
from matplotlib.animation import FuncAnimation
import quaternionic
import spherical
import qnmfitsrd as qnmfits

class qnm_viz:
    def __init__(self, sim, l_max=9, precomp_sYlm = False):
        lon = np.linspace(-np.pi, np.pi, 200)
        lat = np.linspace(-np.pi/2, np.pi/2, 200)
        self.Lon, self.Lat = np.meshgrid(lon, lat)

        theta_vals = np.linspace(0, np.pi, 100)
        phi_vals = np.linspace(0, 2*np.pi, 100)
        self.theta, self.phi = np.meshgrid(theta_vals, phi_vals) 
        
        self.x = np.sin(self.theta) * np.cos(self.phi)
        self.y = np.sin(self.theta) * np.sin(self.phi)
        self.z = np.cos(self.theta)

        self.index = np.where(sim.times == 0)[0][0]

        self.wigner = spherical.Wigner(l_max)

        if precomp_sYlm:
            self.all_mesh_sYlm = np.empty((l_max+1, 2*l_max+1, self.theta.shape[0], self.theta.shape[1]), dtype=complex)
            for l in range(l_max+1):
                for m in range(-l, l+1):
                    self.all_mesh_sYlm[l, m] = self.mesh_sYlm(l, m, self.theta, self.phi)


    def latlon(self):
        return self.Lon, self.Lat


    def mesh_sYlm(self, l, m, theta, phi, s=-2):
        Y = np.empty_like(theta, dtype=complex)
        for i in range(theta.shape[0]):
            for j in range(theta.shape[1]):
                R = quaternionic.array.from_spherical_coordinates(theta[i, j], phi[i, j])
                Y[i, j] = self.wigner.sYlm(s, R)[self.wigner.Yindex(l, m)]
        return Y


    def plot_mapping_projection(self, mapping, mode_mapping, expected): 

        fig, axs = plt.subplots(nrows=2, ncols=2, 
                            subplot_kw={'projection': 'mollweide'}, 
                            figsize=(12,8))

        Lon, Lat = self.Lon, self.Lat
        F = mode_mapping 
        G = expected 

        axs[0,0].title.set_text(str(mapping)+'_real')
        axs[0,0].pcolormesh(Lon, Lat, np.real(F), cmap=plt.cm.jet)

        axs[1,0].title.set_text(str(mapping)+'_imag')
        axs[1,0].pcolormesh(Lon, Lat, np.imag(F), cmap=plt.cm.jet)

        axs[0,1].title.set_text('S_real')
        axs[0,1].pcolormesh(Lon, Lat, np.real(G), cmap=plt.cm.jet)

        axs[1,1].title.set_text('S_imag')
        axs[1,1].pcolormesh(Lon, Lat, np.imag(G), cmap=plt.cm.jet)

        return fig, axs 
    

    def plot_mapping_sphere(self, mapping, mode_mapping, expected):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.plot_surface(self.x, self.y, self.z, facecolors=plt.cm.jet(np.real(mode_mapping)), rstride=2, cstride=2, vmin=0, vmax=1)

        return ax


    def animate_sYlm_sphere(self, sim, min_t0 = 0, max_t0 = 100, step = 1, modes=None, save = False):

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        if modes is None:
            spherical_modes = list(sim.h.keys())
        else:
            spherical_modes = modes 

        def update(step):
            index = self.index + step
            t = sim.times[index] 
            h = np.sum([sim.h[mode][index] * self.all_mesh_sYlm[mode[0], mode[1]] for mode in spherical_modes], axis=0)
            ax.clear() 
            ax.plot_surface(self.x, self.y, self.z, facecolors=plt.cm.viridis(h.real), rstride=2, cstride=2,)
            ax.set_title(f'Time: {t*sim.Mf}')  
            return ax 

        ani = FuncAnimation(fig, update, frames=range(min_t0, max_t0, step), interval=1)

        if save:
            ani.save('animation.mp4', writer='ffmpeg')

        return ani 

    def animate_qnm_sphere(self, sim, model_modes, min_t0 = 0, max_t0 = 100, step = 1, spherical_modes=None, save = False):

        # NEEDS FIXING / MAKING WORK 

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        if spherical_modes is None:
            spherical_modes = list(sim.h.keys())
        else:
            spherical_modes = spherical_modes 

        best_fit = qnmfits.multimode_ringdown_fit(
            sim.times,
            sim.h,
            model_modes,
            Mf=sim.Mf,
            chif=sim.chif_mag,
            t0=0,
            spherical_modes = spherical_modes
            )

        def update(step):
            
            for loop in range(len(best_fit['C'])):
                A = best_fit['C'][loop]
                ans += A * self.mesh_sYlm(self, l, m, theta, phi, s=-2)
                i += 1
            ans /= np.max(np.abs(ans)) # normalise peak value

            ax.clear() 
            ax.plot_surface(self.x, self.y, self.z, facecolors=plt.cm.viridis(h.real), rstride=2, cstride=2)
            ax.set_title(f'Time: {t*sim.Mf}')  
            return ax 

        ani = FuncAnimation(fig, update, frames=range(min_t0, max_t0, step), interval=1)

        if save:
            ani.save('animation.mp4', writer='ffmpeg')

        return ani 

import numpy as np
from numpy import pi
from quakeio.core import QuakeComponent, QuakeSeries

try:
    import jax.numpy as jnp
except:
    import numpy as jnp


from pathlib import Path
try:
    import matplotlib.pyplot as plt
    try:
        style_file = Path(__file__).parents[1]/"assets/brace2.mplstyle"
        plt.style.use([style_file])
        import tempfile
        with tempfile.NamedTemporaryFile() as tf:
            plt.plot([0.0, 1.0], [0.0, 0.1], label=r"$\sigma")
            plt.savefig(tf.name)
    except (FileNotFoundError, RuntimeError, Exception) as e:
        import matplotlib
        matplotlib.rcParams.update({
            "text.usetex": False,
            "text.latex.preamble": ""
        })

    finally:
        plt.close()

except ImportError:
    pass


def plot_grid(series, label=None):
    fig, ax = plt.subplots(len(series),1,sharex=True)
    cycler = ax[0]._get_lines.prop_cycler
    lbl = (lambda x: x[label]) if label is not None else lambda x: None
    if label is not None and not isinstance(label,(tuple,list)):
        label = [s[label] for s in series]
    for i,s in enumerate(series):
        s.plot(ax=ax[i],label=label[i],color=next(cycler)['color'])

    if label is not None:
        fig.legend()


def _plot_func(plot):

    def __call__(self, *args, **kwds):
        if "ax" in kwds:
            self.ax = kwds.pop("ax")
            self.fig = self.ax.figure
        else:
            self.fig, self.ax = plt.subplots()

        if "title" in kwds:
            self.ax.set_title(kwds.pop("title"))

        plot(self, *args, **kwds)
        return self.ax
    return __call__


class Spectrum:
    def __init__(self, accel, *args, **kwds):
        self._accel = accel
        self.kwds = kwds

    def accel(self):
        pass

    def time(self):
        pass

    def spect(self,accel=None,dt=None,damping=None,per=None,gamma=1/2,beta=1/4,interp=None):
        if interp is None:
            from scipy.interpolate import interp1d as interp
        if damping is None:
            damping = ("damping" in self.kwds and self.kwds.pop("damping")) or 0.0

        if accel is None:
            accel = self._accel
        if isinstance(accel, QuakeComponent):
            dt = accel.accel["time_step"]
            accel = accel.accel.data
        elif isinstance(accel, QuakeSeries):
            dt = accel["time_step"]
            accel = accel.data

        if isinstance(damping, float):
            damping = [damping]

        if per is None:
            per = self.kwds.get("per",None)
        if per is None:
            per = np.arange(0.02, 1.0, 0.01)
        elif isinstance(per,tuple):
            per = np.linspace(*per)

        return _accel_spectrum(
            accel=accel,
            dt=dt,
            per=per,
            gamma=gamma,
            beta=beta,
            damping=damping,
            interp=interp
        )

    @_plot_func
    def plot(self, **kwds):
        #sa = spectrum(self._accel, **kwds)
        tsa = self.spect(self._accel, **kwds)
        if len(tsa) > 2:
            for sa in tsa[1:]:
                self.ax.plot(tsa[0], sa)
        else:
            self.ax.plot(tsa[0],tsa[1])
        self.ax.set_xlabel(f"Period, (sec.)")
        chn = self._accel["channel"]
        self.ax.set_title(f"Response spectrum (Chn. {chn})")


class rstf(Spectrum):
    def __init__(self, pairs, **kwds):
        self._pairs = pairs
        self.kwds = kwds

    def tf(self):
        self.tf = tsa = transfer_function(self._pairs, **self.kwds)
        return tsa


    @_plot_func
    def plot(self, **kwds):
        self.tf = tsa = transfer_function(self._pairs, **self.kwds)
        if len(tsa) >= 2:
            for sa in tsa[1:]:
                self.ax.plot(tsa[0], sa)
        self.ax.set_xlabel(f"Period, (sec.)")
        chn1 = self._pairs[0]["channel"]
        chn2 = self._pairs[1]["channel"]
        self.ax.set_title(f"Transfer function (Chn. {chn1} vs. {chn2})")


def transfer_function(pairs, *args, **kwds):
    s1,s2 = Spectrum(pairs[0], **kwds).spect(), Spectrum(pairs[1], **kwds).spect()
    t = s2/s1
    t[0,:] = s1[0,:]
    return t

def spectrum(*args, **kwds):
    return _accel_spectrum(*args, **kwds)

def _accel_spectrum(accel,dt,damping,per,gamma=1/2,beta=1/4,interp=None):
    numper = len(per)
    m = 1.0
    numdata = len(accel)
    t = np.arange(0,(numdata)*dt,dt)

    u0,v0 = 0.0, 0.0
    SA = np.zeros((1+len(damping), numper))
    SA[0,:] = per[:]

    for di,dmp in enumerate(damping):

        for i in range(numper):
            if dt/per[i]>0.02:
                dtp = per[i]*0.02
                dtpx = np.arange(0,max(t),dtp)
                dtpx = dtpx
                accfrni = interp(t, accel)(dtpx)
                accfrn = accfrni[1:len(accfrni)-1]
                numdatan = len(accfrn)
                p  =  -m*accfrn
            else:
                dtp = dt
                accfrn = accel
                p = -m*accfrn
                numdatan = numdata

            u = np.zeros((numdatan))
            v = np.zeros((numdatan))
            a = np.zeros((numdatan))

            k = 4*pi**2*m/per[i]**2
            c = 2*dmp*np.sqrt(k/m)
            kstar = k + gamma*c/(beta*dtp) + m/(beta*dtp**2.0)
            acons = m/(beta*dtp) + gamma*c/beta
            bcons = m/(2*beta) + dtp*(gamma/(2*beta)-1)*c
            u[0] = u0
            v[0] = v0
            a[0] = (p[0]-c*v[0]-k*u[0])/m
            for j in range(1,numdatan):
                deltap = p[j]-p[j-1]
                deltaph = deltap+acons*v[j-1]+bcons*a[j-1]
                deltau = deltaph/kstar
                deltav = gamma*deltau/(beta*dtp) - gamma*v[j-1]/beta+dtp*(1-gamma/(2*beta))*a[j-1]
                deltaa = deltau/(beta*dtp**2)-v[j-1]/(beta*dtp)-a[j-1]/(2*beta)
                u[j] = u[j-1] + deltau
                v[j] = v[j-1] + deltav
                a[j] = a[j-1] + deltaa
            atot = a + accfrn
            SA[1+di,i] = abs(max(atot, key=abs))
    return SA


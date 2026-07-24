"""
Interactive Plotly visualizations of complex modes in
non-classically damped oscillators.

Migrated from the companion code for C. Chern and O. M. O'Reilly,
"On Complex Modal Inverse Analysis: The Cautionary Tale of the
Single Surviving Mode," ASME Journal of Applied Mechanics, 2026.
doi:10.1115/1.4072378
"""
import warnings
import numpy as np
import scipy.integrate
import plotly.graph_objects as go
from plotly.subplots import make_subplots

class TwoStoryFrame:
    """
    The dynamics of a two-story shear-frame structural system,
    idealized as two stacked damped harmonic oscillators.

      ________
      |__m2__| --> u2(t)
      |      |
      |      | k2, c2
      |______|
      |__m1__| --> u1(t)
      |      |
      |      | k1, c1
    __|______|__
    //////////// -> f(t)
    """

    def __init__(self,
                 m1=1,
                 m2=1.2,
                 damping='proportional',
                 Lambda_c=[0.75,0.5],
                 c=None,
                 k1=100,
                 k2=200,
                 L = 3,
                 w = 2,
                 h = 1,
                 u0 = None,
                 up0 = None,
                 m = None, # optional: directly set the mass matrix
                 k = None, # optional: directly set the stiffness matrix
                 i = 0
                 ):

        if m is not None:
            self.minv = np.linalg.inv(m)   # full matrix provided: invert it
        else:
            self.minv = np.array([          # build from scalar args m1, m2
                [1/m1, 0],
                [0, 1/m2]
            ])

        if k is not None:
            self.k = k                      # full matrix provided: use directly
        else:
            self.k = np.array([             # build from scalar args k1, k2
                [k1+k2, -k2],
                [-k2,    k2]
            ])

        self.minvk = self.minv @ self.k
        Lambda_k,Phi = np.linalg.eig(self.minvk)

        # Get the damping matrix c
        if damping=='proportional':
            if (Lambda_c is None) or (not isinstance(Lambda_c,list) and not isinstance(Lambda_c,np.ndarray)):
                raise ValueError("Proportional damping selected. Please provide valid proportional damping factors Lambda_c (list of length 2 or array of shape (2,))")
            elif isinstance(Lambda_c,list) and not len(Lambda_c)==2:
                raise ValueError("Lambda_c must be length 2.")
            elif isinstance(Lambda_c,np.ndarray) and not Lambda_c.shape==(2,):
                raise ValueError("Lambda_c must be shape (2,)")
            else:
                self.c = Phi@np.diag(Lambda_c)@np.linalg.inv(Phi)
        else:
            if not isinstance(c,np.ndarray):
                raise ValueError("Non-proportional damping selected. Please provide valid damping matrix c as an array.")    
            else:
                self.c = c

        self.minvc = self.minv@self.c

        self.A = np.array([
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [-self.minvk[0,0], -self.minvk[0,1], -self.minvc[0,0], -self.minvc[0,1]],
            [-self.minvk[1,0], -self.minvk[1,1], -self.minvc[1,0], -self.minvc[1,1]],
        ])

        self.freqs,self.shapes = np.linalg.eig(self.A) # TODO: Figure out how to get mode shape from A's eigendecomposition

        self.L = L
        self.w = w
        self.h = h
        
        if i is None:
            if u0 is None or up0 is None:
                raise ValueError("Provide initial displacements and velocities, or provide a mode index (0 or 1).")
            else:
                self.u0 = u0
                self.up0 = up0
        elif i in [0,1]:
            if u0 is not None or up0 is not None:
                warnings.warn("Both initial value and mode index are provided. Mode index takes precedence, so initial values will be thrown away.")
            self.u0  = self.shapes[:2,i] # TODO
            self.up0 = self.shapes[2:,i] # TODO
        else:
            if u0 is not None and up0 is not None:
                warnings.warn("Invalid mode index provided. Using specified initial values.")
                self.u0 = u0
                self.up0 = up0
            else:
                raise ValueError("Provide valid initial displacements and velocities, or provide valid mode index (0 or 1).")


    def us(self,ts):
        """
        Compute the displacements of each floor at each time sample
        in ts.

        ts:         time samples
        us:         array of u values corresponding to each t in ts
        uprimes:    array of u' values corresponding to each t in ts
        """

        def dydt(t,y):
            """
            Returns the derivative, with respect to t,
            of the state vector, (u1,u2,u1',u2'),
            where u1' and u2' are the derivatives of u1 and u2 with respect to t.
            Equation of motion:
                (u1,u2)'' = -m^(-1)@k@(u1,u2) -m^(-1)@c@(u1,u2)'

            t:   time
            y:   (u1,u2,u1',u2') vector
            """
            
            return self.A@y
        
        ode_sol = scipy.integrate.solve_ivp(fun=dydt,
                                            t_span=(ts[0],ts[-1]),
                                            y0=np.array([*self.u0,*self.up0]),
                                            t_eval=ts,)
        y = ode_sol.y
        us = y[:2]
        uprimes = y[2:]

        return us,uprimes


    def get_x_frame_floors(self,u):
        """
        Get the coordinates of the centroids of both floors.
        """
        x11 = u[0]
        x12 = self.L
        x21 = u[1]
        x22 = 2*self.L
        return np.array([
            [x11, x12], # Floor 1 centroid
            [x21, x22]  # Floor 2 centroid
        ])

    def get_border(self,x_floor):
        """
        Get the rectangular border of a certain floor or mass,
        based on its centroid x_floor.
        """
        w = self.w
        h = self.h
        x1 = x_floor[0]
        x2 = x_floor[1]

        border = np.array([
            [x1-w/2, x2-h/2],
            [x1+w/2, x2-h/2],
            [x1+w/2, x2+h/2],
            [x1-w/2, x2+h/2],
            [x1-w/2, x2-h/2],
        ])

        return border


    def column_shape(self,ui,uj,L):
        """
        Get the deformed column shape of a certain floor. Requires
        the transverse end displacements and the height of the floor.
        Fits a cubic polynomial with end displacements of ui and uj
        and end rotations of 0.

        ax^3 + bx^2 + cx + d = y

        3ax^2 + 2bx + c = y'

        |  x  |  y  |  y' |
        |-----|-----|-----|
        |  0  |  ui |  0  |
        |  L  |  uj |  0  |
        """

        van = np.array([
            [     0,    0, 0, 1],
            [  L**3, L**2, L, 1],
            [     0,    0, 1, 0],
            [3*L**2,  2*L, 1, 0]
        ])

        bdc = np.array([
            ui,
            uj,
            0,
            0
        ])

        coeffs = np.linalg.solve(van,bdc)

        return lambda x: np.dot(coeffs, np.array([x**3,x**2,x,x**0]))
    
    
    def get_columns(self,x):
        """
        Get the first and second floor columns
        """
        w = self.w
        h = self.h

        bases = np.array([
            [     0-w/2,          0],
            [     0+w/2,          0],
            [x[0,0]-w/2, x[0,1]+h/2],
            [x[0,0]+w/2, x[0,1]+h/2]
        ])
        tops = np.array([
            [x[0,0]-w/2, x[0,1]-h/2],
            [x[0,0]+w/2, x[0,1]-h/2],
            [x[1,0]-w/2, x[1,1]-h/2],
            [x[1,0]+w/2, x[1,1]-h/2]
        ])

        # columns: 4 (columns) x 100 (points) x 2 (coordinates)
        columns = np.array([
            [(self.column_shape(bases[i,0],tops[i,0],tops[i,1]-bases[i,1]))(np.linspace(0,tops[i,1]-bases[i,1],100)), np.linspace(0,tops[i,1]-bases[i,1],100)+bases[i,1]]
            for i in range(4)
        ])

        return columns
    
    def plot_as_frame(self,ts,fig,scale=5):
        us,uprimes = self.us(ts)
        us = scale*us.T
        xs = [np.real(self.get_x_frame_floors(u)) for u in us] # TODO: check out real casting

        # borders: nt x 2 (floors) x 5 (points) x 2 (coordinates)
        borders = np.array([
            [self.get_border(x[0]), self.get_border(x[1])]
            for x in xs
            ])
        
        # columns: nt x 4 (columns) x 100 (points) x 2 (coordinates)
        columns = np.array([
            self.get_columns(x)
            for x in xs
        ])

        xm = np.min(borders[:,:,:,0])
        xM = np.max(borders[:,:,:,0])
        ym = -1
        yM = 2*self.L+1

        # Initial postions
        fig.add_trace(
            go.Scatter(x=borders[0,0,:,0],y=borders[0,0,:,1],
                        name="Floor 1",
                        mode="lines",
                        # line=dict(color="#1f77b4"), # muted blue
                        ),
            # row=1, col=2
        )
        fig.add_trace(
            go.Scatter(x=borders[0,1,:,0],y=borders[0,1,:,1],
                        name="Floor 2",
                        mode="lines",
                        # line=dict(color="#ff7f0e"), # safety orange
                        ),
            # row=1, col=2
        )
        fig.add_traces(
            [
                go.Scatter(x=columns[0,i,0,:],y=columns[0,i,1,:],
                           mode="lines",
                           name=None,
                           line=dict(color="gray"),
                           showlegend=False)
                for i in range(4)
            ],
            # row=1, col=2
        )
        # Layout
        fig.update_layout(
                xaxis=dict(range=[xm,xM],zeroline=False),
                yaxis=dict(range=[ym,yM]),
                updatemenus=[dict(
                                type="buttons",
                                buttons=[dict(label="Play",
                                            method="animate",
                                            args=[
                                                None,
                                                {"frame": {"duration": 1000*(ts[1]-ts[0]), "redraw": False},
                                                "fromcurrent": True,
                                                "transition": {"duration": 100*(ts[1]-ts[0]), "mode":"immediate"}}]
                                            )],
                                x=-0.1
                                )],

                width=500, height=400,
                xaxis_showgrid=False, yaxis_showgrid=False,
                margin=dict(l=60, r=20, t=40, b=50)
                ),
        # Frames
        fig.frames = [go.Frame(
                data = [
                    # Update the Floor 1 location
                    go.Scatter(
                        x=borders[k,0,:,0],
                        y=borders[k,0,:,1],
                        mode="lines",
                    ),
                    # Update the Floor 2 location
                    go.Scatter(
                        x=borders[k,1,:,0],
                        y=borders[k,1,:,1],
                        mode="lines",
                    ),
                    *[
                    # Update the column deformations
                    go.Scatter(
                        x=columns[k,i,0,:],
                        y=columns[k,i,1,:],
                        mode="lines",
                        line=dict(color="gray"),
                        showlegend=False
                    )
                    for i in range(4)
                    ]
                       ],
                traces=[0,1,2,3,4,5], # The Floor 1 frames replace trace 1;
                                      # The Floor 2 frames replace trace 2;
                                      # The columns frames replace traces 3-6
                )
                for k in range(len(ts))]

        return fig


class SimplestOscillator:
    """
    The simplest non-classically damped oscillator.

                  ^ y(t)
               ___|___
    /|___vvv___|  m  | --> x(t)
    /|         |_____|
    /|       /    |
    /|      z     z
    /|     z      z
    /|____/_______|___
    /////////////////// -> f(t)=0 (homogeneous)

    Dimensionless equations of motion:
        x'' + 2ζ₁ x'            + x  = 0
        y'' + 2ζ₃ x' + 2ζ₂ y'   + Δy = 0

    State vector: z = [x, y, x', y']
    State-space form (paper Eq. 2.4):

        A = [[ 0,   0,   1,    0  ],
             [ 0,   0,   0,    1  ],
             [-1,   0, -2ζ₁,   0  ],
             [ 0,  -Δ, -2ζ₃, -2ζ₂]]

    The ζ₃ term couples x-velocity into the y-equation (non-classical /
    non-proportional damping). Setting ζ₃=0 gives a classically damped system.

    Eigenvalues of A:
        λ = -ζ₁ ± i√(1-ζ₁²)   (x-mode)
        λ = -ζ₂ ± i√(Δ-ζ₂²)   (y-mode)
    Note: ζ₃ does not appear in the eigenvalues.

    Parameters
    ----------
    zeta1 : float   Damping ratio of x-mode (ζ₁)
    zeta2 : float   Damping ratio of y-mode (ζ₂)
    zeta3 : float   Non-classical coupling coefficient (ζ₃); 0 ↔ classical
    Delta : float   Stiffness ratio Δ = ωy²/ωx²; 1 ↔ equal stiffness
    u0    : array   Initial displacements [x₀, y₀]   (default [1, 1])
    up0   : array   Initial velocities    [x'₀, y'₀] (default [1, 1])
    """

    def __init__(self, zeta1=0.05, zeta2=0.04, zeta3=0.02, Delta=1.5,
                 u0=None, up0=None, L=3, w=2, h=1):
        self.zeta1 = zeta1
        self.zeta2 = zeta2
        self.zeta3 = zeta3
        self.Delta  = Delta
        self.u0  = np.array(u0  if u0  is not None else [1.0, 1.0])
        self.up0 = np.array(up0 if up0 is not None else [1.0, 1.0])
        self.A = self._build_A(zeta1, zeta2, zeta3, Delta)
        self.L = L
        self.w = w
        self.h = h

    @staticmethod
    def _build_A(z1, z2, z3, D):
        return np.array([
            [ 0,   0,     1,     0],
            [ 0,   0,     0,     1],
            [-1,   0, -2*z1,     0],
            [ 0,  -D, -2*z3, -2*z2],
        ])

    def solve(self, ts):
        """
        Integrate dz/dτ = A z from t=ts[0] to t=ts[-1].

        Returns
        -------
        y : ndarray, shape (4, len(ts))
            Rows are [x(τ), y(τ), x'(τ), y'(τ)].
        """
        y0 = np.concatenate([self.u0, self.up0])
        sol = scipy.integrate.solve_ivp(
            lambda t, y: self.A @ y,
            (ts[0], ts[-1]),
            y0,
            t_eval=ts,
        )
        return np.real(sol.y)

    # Geometry
    def get_x_simple_mass(self, u):
        """
        Get the centroid of the mass from displacement vector u = [x, y].

        Returns
        -------
        ndarray, shape (1, 2) : [[cx, cy]] where cx = L+x, cy = L+y.
        """
        return np.array([[self.L + u[0], self.L + u[1]]])

    def get_border(self, x_mass):
        """
        Get the rectangular border of the mass from its centroid x_mass = [cx, cy].

        Returns
        -------
        ndarray, shape (5, 2) : closed polygon (5th point repeats 1st).
        """
        cx, cy = x_mass[0], x_mass[1]
        w, h   = self.w, self.h
        return np.array([
            [cx - w/2, cy - h/2],
            [cx + w/2, cy - h/2],
            [cx + w/2, cy + h/2],
            [cx - w/2, cy + h/2],
            [cx - w/2, cy - h/2],
        ])

    def get_springs(self, x):
        """
        Get the three springs attached to the mass.

        Parameters
        ----------
        x : ndarray, shape (1, 2)
            Current mass centroid [[cx, cy]], as returned by get_x_simple_mass.

                          ^ y
            /|            |
            /|--horiz-----+-------> x
            /|           /|
            /|     diag / |  vert
            /|         /  |
            /|________/___|___
            ///////////////////

        Returns
        -------
        springs : ndarray, shape (3, 2, 2)
            springs[i, 0, :] = x-coordinates [base_x, top_x] of spring i
            springs[i, 1, :] = y-coordinates [base_y, top_y] of spring i
        """
        w  = self.w
        h  = self.h
        cx = x[0, 0]
        cy = x[0, 1]

        bases = np.array([
            [self.L, 0      ],  # vertical spring:   anchored at (L, 0) on ground
            [0,      self.L ],  # horizontal spring: anchored at (0, L) on left wall
            [0,      0      ],  # diagonal spring:   anchored at corner (0, 0)
        ])
        tops = np.array([
            [cx,       cy - h/2],  # bottom edge of mass
            [cx - w/2, cy      ],  # left edge of mass
            [cx - w/2, cy - h/2],  # bottom-left corner of mass
        ])

        springs = np.array([
            [
                [bases[i, 0], tops[i, 0]],  # x: [base_x, top_x]
                [bases[i, 1], tops[i, 1]],  # y: [base_y, top_y]
            ]
            for i in range(3)
        ])
        return springs

    def plot_as_simple_mass(self, ts, fig, scale=1):
        """
        Animate the mass moving on a smooth plane under free vibration.

        Parameters
        ----------
        ts    : 1-D array   Time samples
        fig   : go.Figure   Figure to add traces/frames to
        scale : float       Amplification factor applied to displacements

        Returns
        -------
        fig : go.Figure (modified in place and returned)
        """
        state = self.solve(ts)
        us    = scale * state[:2].T   # shape (nt, 2)

        xs_geom = [self.get_x_simple_mass(u) for u in us]

        borders = np.array([self.get_border(x[0]) for x in xs_geom])  # (nt, 5, 2)
        springs = np.array([self.get_springs(x)   for x in xs_geom])  # (nt, 3, 2, 2)

        xm = np.min(borders[:, :, 0]) - 0.5
        xM = np.max(borders[:, :, 0]) + 0.5
        ym = np.min(borders[:, :, 1]) - 0.5
        yM = np.max(borders[:, :, 1]) + 0.5

        fig.add_trace(
            go.Scatter(x=borders[0, :, 0], y=borders[0, :, 1],
                       name='simple mass', mode='lines'),
        )
        fig.add_traces([
            go.Scatter(x=springs[0, i, 0, :], y=springs[0, i, 1, :],
                       mode='lines', name=None,
                       line=dict(color='gray'), showlegend=False)
            for i in range(3)
        ])
        fig.update_layout(
            xaxis=dict(range=[xm, xM], zeroline=False),
            yaxis=dict(range=[ym, yM]),
            updatemenus=[dict(
                type='buttons',
                buttons=[dict(label='Play', method='animate',
                              args=[None, {'frame': {'duration': 1000*(ts[1]-ts[0]), 'redraw': False},
                                          'fromcurrent': True,
                                          'transition': {'duration': 100*(ts[1]-ts[0]), 'mode': 'immediate'}}])],
                x=-0.1,
            )],
            width=500, height=400,
            xaxis_showgrid=False, yaxis_showgrid=False,
            margin=dict(l=60, r=20, t=40, b=50),
        )
        fig.frames = [
            go.Frame(
                data=[
                    go.Scatter(x=borders[k, :, 0], y=borders[k, :, 1], mode='lines'),
                    *[go.Scatter(x=springs[k, i, 0, :], y=springs[k, i, 1, :],
                                 mode='lines', line=dict(color='gray'), showlegend=False)
                      for i in range(3)],
                ],
                traces=[0, 1, 2, 3],
            )
            for k in range(len(ts))
        ]
        return fig

    def mode_ics(self, k):
        """
        Compute initial conditions that excite only mode k (1 or 2).

        Uses the paper's closed-form eigenvectors (Eq. 2.4 and following):

            λ₁ = -ζ₂ + i√(Δ-ζ₂²),   λ₂ = -ζ₁ + i√(1-ζ₁²)
            Γ  = (ζ₁ - ζ₂) / ζ₃      (requires ζ₃ ≠ 0)

            ψ₁ = [0,  λ₁*/Δ,  0,  1]
            ψ₂ = [λ₂*(α + Γ),  λ₂*,  α + Γ,  1]
                 where α = (1-Δ)λ₂* / (2ζ₃)

        Pure-mode initial state: state₀ = 2·Re(ψₖ),
        so u0 = state₀[:2], up0 = state₀[2:].

        For ζ₃ = 0 (classical): modes decouple and are orthogonal

        Parameters
        ----------
        k : int   1 or 2

        Returns
        -------
        u0  : ndarray, shape (2,)
        up0 : ndarray, shape (2,)
        """
        z1, z2, z3, D = self.zeta1, self.zeta2, self.zeta3, self.Delta
        lam1 = -z2 + 1j * np.sqrt(max(D  - z2**2, 0.0))
        lam2 = -z1 + 1j * np.sqrt(max(1.0 - z1**2, 0.0))

        if k == 1:
            # Mode 1: y-mode.  ψ₁ = [0, λ₁*/Δ, 0, 1]  (valid for all ζ₃)
            psi = np.array([0, np.conj(lam1) / D, 0, 1], dtype=complex)
        elif k == 2:
            if z3 == 0:
                # Classical limit: mode 2 is pure x-oscillation.
                psi = np.array([1, 0, lam2, 0], dtype=complex)
            else:
                Gamma = (z1 - z2) / z3
                alpha = (1 - D) * np.conj(lam2) / (2 * z3)
                psi   = np.array([np.conj(lam2) * (alpha + Gamma),
                                  np.conj(lam2),
                                  alpha + Gamma,
                                  1], dtype=complex)
        else:
            raise ValueError("k must be 1 or 2")

        # Apply monophase rotation: rotate by -arg(φ[j]) where j is the
        # index of the y-component.
        phi   = psi[:2]
        # j     = int(np.argmax(np.abs(phi)))
        j     = len(phi) - 1  # use the last component (y) for rotation
        phase = -np.angle(phi[j])
        psi   = np.exp(1j * phase) * psi

        state0 = np.real(psi)
        u0, up0 = state0[:2], state0[2:]
        scale = np.linalg.norm(u0) # unitize the vector
        return u0 / scale, up0 / scale

    # Physics checks
    def check_eigenvalues(self, tol=1e-10):
        """
        Verify eigenvalues of A match the closed-form expressions and that
        ζ₃ does not appear in them.

        Returns
        -------
        ok      : bool    True if eigenvalues match within tol
        computed: ndarray Eigenvalues of A (sorted by imaginary part)
        expected: ndarray Closed-form eigenvalues (sorted to match)
        """
        z1, z2, D = self.zeta1, self.zeta2, self.Delta
        lx = np.array([-z1 + 1j*np.sqrt(max(1 - z1**2, 0)),
                       -z1 - 1j*np.sqrt(max(1 - z1**2, 0))])
        ly = np.array([-z2 + 1j*np.sqrt(max(D  - z2**2, 0)),
                       -z2 - 1j*np.sqrt(max(D  - z2**2, 0))])
        expected = np.sort_complex(np.concatenate([lx, ly]))
        computed = np.sort_complex(np.linalg.eigvals(self.A))
        ok = bool(np.allclose(computed, expected, atol=tol))
        return ok, computed, expected

# The three examples used throughout the paper
PRESETS = [
    {'label': 'Ex 1 — ζ₃=0 (classical)',      'z1': 0.05, 'z2': 0.02, 'z3': 0.00, 'D': 1.5, 'shortname': '(Classical)'},
    {'label': 'Ex 2 — Δ=1 (monophase)',       'z1': 0.05, 'z2': 0.02, 'z3': 0.10, 'D': 1.0, 'shortname': '(Monophase Non-Classical)'},
    {'label': 'Ex 3 — Δ=1.5 (non-monophase)', 'z1': 0.05, 'z2': 0.02, 'z3': 0.10, 'D': 1.5, 'shortname': '(Non-Monophase Non-Classical)'},
]


def plot_simplest_oscillator(tf=50.0, nt=200, speed_up=False, springs=False):
    """
    Build and return a Plotly figure showing Mode 1 (left) and Mode 2 (right)
    side-by-side for the simplest non-classically damped oscillator.

    Mode ICs are computed per-preset from the paper's closed-form eigenvectors
    (monophase-rotated, unit-magnitude displacement).  IC subtitles update when
    switching presets via the example buttons.

    Preset buttons (below figure):
        Ex 1 — ζ₃=0  (classically damped → orthogonal straight-line modes)
        Ex 2 — Δ=1   (both modes monophase → straight-line trajectories)
        Ex 3 — Δ=1.5 (mode 2 non-monophase → curved trajectory)

    Parameters
    ----------
    tf       : float   Final dimensionless time τ
    nt       : int     Number of animation frames
    speed_up : bool    If True, play at 10x real time

    Returns
    -------
    fig : go.Figure
    """

    presets = PRESETS

    ts           = np.linspace(0, tf, nt)
    frame_dur_ms = round(1000 * tf / nt, 1)
    if speed_up:
        frame_dur_ms /= 10

    L, w, h = 1.5, 0.9, 0.6
    geom = SimplestOscillator(L=L, w=w, h=h)

    def geometry_at(x, y):
        u        = np.array([x, y])
        centroid = geom.get_x_simple_mass(u)
        border   = geom.get_border(centroid[0])
        springs  = geom.get_springs(centroid)
        return border[:, 0], border[:, 1], springs

    def fmt_ic(u0, up0):
        """Format IC as subtitle string."""
        u  = ', '.join(f'{v:.3f}' for v in u0)
        up = ', '.join(f'{v:.3f}' for v in up0)
        return f'u₀=[{u}], u\'₀=[{up}]'

    # Solve both modes for every preset
    preset_data = []
    for p in presets:
        osc = SimplestOscillator(p['z1'], p['z2'], p['z3'], p['D'])
        u0_1, up0_1 = osc.mode_ics(1)
        u0_2, up0_2 = osc.mode_ics(2)

        sol1 = SimplestOscillator(p['z1'], p['z2'], p['z3'], p['D'],
                                  list(u0_1), list(up0_1)).solve(ts)
        sol2 = SimplestOscillator(p['z1'], p['z2'], p['z3'], p['D'],
                                  list(u0_2), list(up0_2)).solve(ts)

        preset_data.append(dict(
            u0_1=u0_1, up0_1=up0_1,
            u0_2=u0_2, up0_2=up0_2,
            xs1=np.real(sol1[0]), ys1=np.real(sol1[1]),
            xs2=np.real(sol2[0]), ys2=np.real(sol2[1]),
            z1=p['z1'], z2=p['z2'], z3=p['z3'], D=p['D'],
            shortname=p['shortname'],
            phi1=_displacement_eigvec(osc, 1),
            phi2=_displacement_eigvec(osc, 2),
        ))

    # Axis ranges
    pad = 0.4

    def axis_range(all_disp_x, all_disp_y):
        xl = min(-pad, L + all_disp_x.min() - w/2 - pad)
        xh = L + all_disp_x.max() + w/2 + pad
        yl = min(-pad, L + all_disp_y.min() - h/2 - pad)
        yh = L + all_disp_y.max() + h/2 + pad
        # square it up
        xspan, yspan = xh - xl, yh - yl
        if xspan > yspan:
            yl -= (xspan - yspan) / 2
            yh += (xspan - yspan) / 2
        else:
            xl -= (yspan - xspan) / 2
            xh += (yspan - xspan) / 2
        return xl, xh, yl, yh

    all_x = np.concatenate([d['xs1'] for d in preset_data]+[d['xs2'] for d in preset_data])
    all_y = np.concatenate([d['ys1'] for d in preset_data]+[d['ys2'] for d in preset_data])

    x_lo, x_hi, y_lo, y_hi = axis_range(all_x, all_y)

    # Build figure: row 1 = animation subplots, row 2 = phase normalization subplots
    fig = make_subplots(rows=2, cols=2,
                        horizontal_spacing=0.02,
                        vertical_spacing=0.05,
                        row_heights=[3, 3])

    def add_subplot_traces(d, xs, ys, col, springs=False):
        """Add the 6 initial traces for one subplot (one mode)."""
        cpx = L + xs
        cpy = L + ys
        bx, by, sp = geometry_at(xs[0], ys[0])
        fig.add_trace(go.Scatter(x=cpx, y=cpy, mode='lines',
                                line=dict(color='lightgray', width=1),
                                showlegend=False, hoverinfo='skip'),
                    row=1, col=col)
        fig.add_trace(go.Scatter(x=[cpx[0]], y=[cpy[0]], mode='lines',
                                 line=dict(color='cornflowerblue', width=1.5),
                                 showlegend=False, hoverinfo='skip'),
                      row=1, col=col)
        fig.add_trace(go.Scatter(x=bx, y=by, mode='lines', fill='toself',
                                 fillcolor='rgba(100,149,237,0.35)',
                                 line=dict(color='steelblue', width=2),
                                 showlegend=False),
                      row=1, col=col)
        if springs:
            for j in range(3):
                fig.add_trace(go.Scatter(x=sp[j, 0, :], y=sp[j, 1, :], mode='lines',
                                        line=dict(color='dimgray', width=1.5),
                                        showlegend=False),
                            row=1, col=col)

    # Phase trace helper
    _PHASE_COLORS  = ['blue', 'red']
    _PHASE_SYMBOLS = ['x', 'circle']
    _PHASE_SIZES   = [10, 10]
    _PHASE_LABELS  = ['x', 'y']
    _THETA         = np.linspace(0, 2 * np.pi, 120)

    def phase_traces_data(phi):
        """Return 4 go.Scatter objects for a phase-normalization plot of phi (2D complex)."""
        traces = []
        for c, color, sym, sz, lbl in zip(
                phi, _PHASE_COLORS, _PHASE_SYMBOLS, _PHASE_SIZES, _PHASE_LABELS):
            r = abs(c)
            # Vector from origin to component, marker + label at tip only
            traces.append(go.Scatter(
                x=[0, float(c.real)], y=[0, float(c.imag)],
                mode='lines+markers+text',
                line=dict(color=color, width=2.5),
                marker=dict(size=[0, sz], symbol=sym, color=color),
                text=['', lbl],
                textposition='top right',
                textfont=dict(size=16, color=color),
                showlegend=False, hoverinfo='skip',
            ))
            # Circle at radius r
            traces.append(go.Scatter(
                x=(r * np.cos(_THETA)).tolist(),
                y=(r * np.sin(_THETA)).tolist(),
                mode='lines',
                line=dict(color='lightgray', width=0.8),
                showlegend=False, hoverinfo='skip',
            ))
        return traces  # [x-vec, x-circle, y-vec, y-circle]

    # Initial display: Ex 3 (index 2)
    default_idx = 2
    d0 = preset_data[default_idx]
    n_per = 6 if springs else 3
    add_subplot_traces(d0, d0['xs1'], d0['ys1'], col=1, springs=springs)  # traces 0–(n_per-1)
    add_subplot_traces(d0, d0['xs2'], d0['ys2'], col=2, springs=springs)  # traces n_per–(2*n_per-1)

    # Phase traces: 4 per mode × 2 modes = 8 traces
    # Mode 1 phase: traces 2*n_per … 2*n_per+3  (row 2, col 1)
    # Mode 2 phase: traces 2*n_per+4 … 2*n_per+7  (row 2, col 2)
    for tr in phase_traces_data(d0['phi1']):
        fig.add_trace(tr, row=2, col=1)
    for tr in phase_traces_data(d0['phi2']):
        fig.add_trace(tr, row=2, col=2)

    # Crosshairs for phase subplots (static shapes, never touched by frames)
    for xref, yref in [('x3', 'y3'), ('x4', 'y4')]:
        for x0, y0, x1, y1 in [(-1.3, 0, 1.3, 0), (0, -1.3, 0, 1.3)]:
            fig.add_shape(type='line', x0=x0, y0=y0, x1=x1, y1=y1,
                          xref=xref, yref=yref,
                          line=dict(color='black', width=1))

    # Static ground + wall lines for each subplot
    for xref, yref, xlo, xhi, ylo, yhi in [
        ('x',  'y',  x_lo, x_hi, y_lo, y_hi),
        ('x2', 'y2', x_lo, x_hi, y_lo, y_hi),
    ]:
        fig.add_shape(type='line', x0=xlo, y0=0, x1=xhi, y1=0,
                      line=dict(color='black', width=1), xref=xref, yref=yref)
        fig.add_shape(type='line', x0=0, y0=ylo, x1=0, y1=yhi,
                      line=dict(color='black', width=1), xref=xref, yref=yref)

    # Annotations
    def fmt_params(d):
        return (f'ζ₁={d["z1"]:g},  ζ₂={d["z2"]:g},  '
                f'ζ₃={d["z3"]:g},  Δ={d["D"]:g}'   )

    def make_annotations(d,phase_axis_lim = 1.3):
        return [
            dict(text='<b>Mode 1</b>',
                 x=0.0, y=1.16, xref='paper', yref='paper',
                 xanchor='left', showarrow=False, font=dict(size=22)),
            dict(text='<b>Mode 2</b>',
                 x=0.51, y=1.16, xref='paper', yref='paper',
                 xanchor='left', showarrow=False, font=dict(size=22)),
            dict(text=fmt_ic(d['u0_1'], d['up0_1']),
                 x=0.0, y=1.08, xref='paper', yref='paper',
                 xanchor='left', showarrow=False, font=dict(size=16)),
            dict(text=fmt_ic(d['u0_2'], d['up0_2']),
                 x=0.51, y=1.08, xref='paper', yref='paper',
                 xanchor='left', showarrow=False, font=dict(size=16)),
            dict(text=f'{fmt_params(d)} {d["shortname"]}',
                 x=0.5, y=1.25, xref='paper', yref='paper',
                 xanchor='center', yanchor='top',
                 showarrow=False, font=dict(size=18)),
            # Axis labels pinned to ends of the axis lines (data coords)
            dict(text='<i>x</i>', x=x_hi, y=0, xref='x', yref='y',
                 xanchor='right', yanchor='top', showarrow=False, font=dict(size=16)),
            dict(text='<i>y</i>', x=-0.1, y=y_hi, xref='x', yref='y',
                 xanchor='right', yanchor='top', showarrow=False, font=dict(size=16)),
            dict(text='<i>x</i>', x=x_hi, y=0, xref='x2', yref='y2',
                 xanchor='right', yanchor='top', showarrow=False, font=dict(size=16)),
            dict(text='<i>y</i>', x=-0.1, y=y_hi, xref='x2', yref='y2',
                 xanchor='right', yanchor='top', showarrow=False, font=dict(size=16)),
            dict(text='Re', x=phase_axis_lim, y=0, xref='x3', yref='y3',
                 xanchor='right', yanchor='top', showarrow=False, font=dict(size=16)),
            dict(text='Im', x=-0.1, y=phase_axis_lim, xref='x3', yref='y3',
                 xanchor='right', yanchor='top', showarrow=False, font=dict(size=16)),
            dict(text='Re', x=phase_axis_lim, y=0, xref='x4', yref='y4',
                 xanchor='right', yanchor='top', showarrow=False, font=dict(size=16)),
            dict(text='Im', x=-0.1, y=phase_axis_lim, xref='x4', yref='y4',
                 xanchor='right', yanchor='top', showarrow=False, font=dict(size=16)),
        ]

    # Pre-compute phase traces per preset (constant within a preset's frames)
    for d in preset_data:
        d['_ph1'] = phase_traces_data(d['phi1'])
        d['_ph2'] = phase_traces_data(d['phi2'])

    # Animation Frames
    frames = []
    phase_axis_lim = 1.3
    for i, d in enumerate(preset_data):
        cpx1 = L + d['xs1'];  cpy1 = L + d['ys1']
        cpx2 = L + d['xs2'];  cpy2 = L + d['ys2']
        anns = make_annotations(d,phase_axis_lim)
        for k in range(nt):
            bx1, by1, sp1 = geometry_at(d['xs1'][k], d['ys1'][k])
            bx2, by2, sp2 = geometry_at(d['xs2'][k], d['ys2'][k])
            mode1_traces = [
                go.Scatter(x=cpx1, y=cpy1),
                go.Scatter(x=cpx1[:k+1], y=cpy1[:k+1]),
                go.Scatter(x=bx1, y=by1, fill='toself',
                           fillcolor='rgba(100,149,237,0.35)'),
            ]
            mode2_traces = [
                go.Scatter(x=cpx2, y=cpy2),
                go.Scatter(x=cpx2[:k+1], y=cpy2[:k+1]),
                go.Scatter(x=bx2, y=by2, fill='toself',
                           fillcolor='rgba(100,149,237,0.35)'),
            ]
            if springs:
                mode1_traces += [go.Scatter(x=sp1[j, 0, :], y=sp1[j, 1, :])
                                  for j in range(3)]
                mode2_traces += [go.Scatter(x=sp2[j, 0, :], y=sp2[j, 1, :])
                                  for j in range(3)]
            frames.append(go.Frame(
                data=mode1_traces + mode2_traces + d['_ph1'] + d['_ph2'],
                traces=list(range(2 * n_per + 8)),
                layout=dict(annotations=anns),
                name=f'ex{i}_t{k}',
            ))
    fig.frames = frames

    # Preset buttons
    preset_buttons = [
        dict(
            label=p['label'],
            method='animate',
            args=[
                [f'ex{i}_t{k}' for k in range(nt)],
                dict(frame=dict(duration=frame_dur_ms, redraw=True),
                     transition=dict(duration=0),
                     fromcurrent=False, mode='immediate'),
            ],
        )
        for i, p in enumerate(presets)
    ]

    fig.update_layout(
        title=dict(text='Simplest Damped Oscillator',
                   font=dict(size=25), x=0.5, xanchor='center',
                   y=0.99, yanchor='top',
                   pad=dict(b=20, t=0)),
        font=dict(size=16),
        annotations=make_annotations(d0),
        xaxis=dict(range=[x_lo, x_hi], showgrid=False,
                   zeroline=False, showticklabels=False,
                   constrain='domain', constraintoward='left'),
        yaxis=dict(range=[y_lo, y_hi], showgrid=False,
                   zeroline=False, showticklabels=False,
                   scaleanchor='x', scaleratio=1, constraintoward='bottom'),
        xaxis2=dict(range=[x_lo, x_hi], showgrid=False,
                    zeroline=False, showticklabels=False,
                    constrain='domain', constraintoward='left'),
        yaxis2=dict(range=[y_lo, y_hi], showgrid=False,
                    zeroline=False, showticklabels=False,
                    scaleanchor='x2', scaleratio=1, constraintoward='bottom'),
        # Phase normalization subplots (row 2)
        xaxis3=dict(range=[-phase_axis_lim,phase_axis_lim], showgrid=False, zeroline=False,
                    showticklabels=False,
                    # tickmode='array', tickvals=[-1, 1], showticklabels=True,
                    constrain='domain', constraintoward='left'),
        yaxis3=dict(range=[-phase_axis_lim,phase_axis_lim], showgrid=False, zeroline=False,
                    showticklabels=False,
                    # tickmode='array', tickvals=[-1, 1], showticklabels=True,
                    scaleanchor='x3', scaleratio=1, constraintoward='bottom'),
        xaxis4=dict(range=[-phase_axis_lim,phase_axis_lim], showgrid=False, zeroline=False,
                    showticklabels=False,
                    # tickmode='array', tickvals=[-1, 1], showticklabels=True,
                    constrain='domain', constraintoward='left'),
        yaxis4=dict(range=[-phase_axis_lim,phase_axis_lim], showgrid=False, zeroline=False,
                    showticklabels=False,
                    # tickmode='array', tickvals=[-1, 1], showticklabels=True,
                    scaleanchor='x4', scaleratio=1, constraintoward='bottom'),
        width=820, height=840,
        plot_bgcolor='white',
        margin=dict(l=20, r=20, t=170, b=160),
        updatemenus=[
            dict(type='buttons', direction='left',
                 x=0.0, y=-0.07, xanchor='left', yanchor='top',
                 showactive=False,
                 buttons=[
                     dict(label='▶ Play', method='animate',
                          args=[[f'ex{default_idx}_t{k}' for k in range(nt)],
                                dict(frame=dict(duration=frame_dur_ms, redraw=True),
                                     fromcurrent=True,
                                     transition=dict(duration=0))]),
                     dict(label='⏹ Stop', method='animate',
                          args=[[None], dict(frame=dict(duration=0, redraw=False),
                                             mode='immediate',
                                             transition=dict(duration=0))]),
                 ]),
            dict(type='buttons', direction='left',
                 x=0.0, y=-0.20, xanchor='left', yanchor='top',
                 showactive=True, active=default_idx,
                 buttons=preset_buttons),
        ],
    )

    # --- Scale the whole figure down uniformly so it fits the page width ---
    scale = 0.75
    for _fr in fig.frames:
        for _a in (_fr.layout.annotations or ()):
            if _a.font is not None and _a.font.size:
                _a.font.size = _a.font.size * scale
    for _a in fig.layout.annotations:
        if _a.font is not None and _a.font.size:
            _a.font.size = _a.font.size * scale
    if fig.layout.title.font.size:
        fig.layout.title.font.size = fig.layout.title.font.size * scale
    fig.layout.font.size = (fig.layout.font.size or 16) * scale
    fig.update_layout(width=round(820 * scale), height=round(840 * scale),
                      margin=dict(l=15, r=15,
                                  t=round(170 * scale), b=round(210 * scale)))

    return fig


def main_frame_figure():
    m1 = 2
    m2 = 3
    k1 = 100
    k2 = 200
    m = np.array([[2,0],[0,3]])
    k = np.array([[300,-200],[-200,200]])

    ## Custom Proportional Damping
    cprop = 0.1*m+0.01*k

    prop_system = TwoStoryFrame(
        m1=m1,
        m2=m2,
        damping='nonproportional',
        c=cprop,
        k1=k1,
        k2=k2,
    )

    nt = 1000
    tf = 10
    ts = np.linspace(0,tf,nt)

    fig = go.Figure()
    fig = prop_system.plot_as_frame(ts,fig)
    fig.update_layout({'title': 'Proportional Damping'})
    fig.show()


    ## Custom Non-Proportional Damping
    cnprop = np.array([[3,-1],[-1,5]])

    nprop_system = TwoStoryFrame(
        m1=m1,
        m2=m2,
        damping='nonproportional',
        c=cnprop,
        k1=k1,
        k2=k2,
    )

    fig_np = go.Figure()
    fig_np = nprop_system.plot_as_frame(ts,fig_np)

    fig_np.update_layout({'title': 'Non-Proportional Damping'})
    fig_np.show()


def main_simple_figure(speed_up=False, springs=False):
    """Build and show the side-by-side Mode 1 / Mode 2 animation."""
    fig = plot_simplest_oscillator(speed_up=speed_up, springs=springs)
    fig.show()


def _displacement_eigvec(osc, k):
    """
    Return the phase-normalised complex displacement eigenvector (2D)
    for mode k of the given SimplestOscillator.
    Mirrors the rotation in mode_ics but returns the raw complex phi
    (before extracting real part) so it can be plotted in the complex plane.
    """
    z1, z2, z3, D = osc.zeta1, osc.zeta2, osc.zeta3, osc.Delta
    lam1 = -z2 + 1j * np.sqrt(max(D  - z2**2, 0.0))
    lam2 = -z1 + 1j * np.sqrt(max(1.0 - z1**2, 0.0))
    if k == 1:
        psi = np.array([0, np.conj(lam1) / D, 0, 1], dtype=complex)
    elif k == 2:
        if z3 == 0:
            psi = np.array([1, 0, lam2, 0], dtype=complex)
        else:
            Gamma = (z1 - z2) / z3
            alpha = (1 - D) * np.conj(lam2) / (2 * z3)
            psi   = np.array([np.conj(lam2) * (alpha + Gamma),
                              np.conj(lam2), alpha + Gamma, 1], dtype=complex)
    phi = psi[:2]
    # j   = int(np.argmax(np.abs(phi)))
    j   = len(phi) - 1  # use the last component (y) for rotation
    phi = np.exp(-1j * np.angle(phi[j])) * phi
    return phi / np.linalg.norm(phi)


def export_plotly_animation(preset_idx=1, output_path='animation.gif',
                            fps=15, nt=200, tf=50.0, springs=False):
    """
    Export one preset's Plotly animation to GIF or MP4/MOV.

    This renders the exact same Plotly figure you see interactively — no
    matplotlib maintenance required.

    Requires on your local machine (not the sandbox):
        pip install kaleido
        pip install imageio[ffmpeg]   # for .mp4/.mov output

    Parameters
    ----------
    preset_idx : int   — 0 = Ex 1, 1 = Ex 2 (default), 2 = Ex 3
    output_path : str  — e.g. 'animation.gif' or 'animation.mp4'
    fps         : int  — frames per second in the output video
    nt          : int  — time steps; fewer = faster export (try 60 for a quick GIF)
    tf          : float — simulation time
    springs     : bool — passed through to plot_simplest_oscillator
    """
    import io as _io
    import copy
    import os

    try:
        import imageio.v3 as iio
        _iio_version = 3
    except ImportError:
        import imageio as iio
        _iio_version = 2

    try:
        import plotly.io as pio
    except ImportError:
        raise ImportError('plotly is required')

    # Verify kaleido is available before spending time building the figure
    try:
        pio.to_image(go.Figure(), format='png', width=10, height=10)
    except Exception as e:
        if 'kaleido' in str(e).lower() or 'orca' in str(e).lower():
            raise RuntimeError(
                'kaleido is not installed. Run:  pip install kaleido\n'
                'Then retry export_plotly_animation().'
            ) from e
        raise

    print('Building Plotly figure …')
    fig = plot_simplest_oscillator(tf=tf, nt=nt, springs=springs)

    # Index all frames by name
    frame_map = {f.name: f for f in fig.frames}

    # Base layout with animation controls stripped (not needed for static renders)
    import json
    base_layout_dict = json.loads(fig.layout.to_json())
    base_layout_dict.pop('updatemenus', None)
    base_layout_dict.pop('sliders', None)

    images = []
    print(f'Rendering {nt} frames (this takes a minute) …')
    for k in range(nt):
        frame = frame_map[f'ex{preset_idx}_t{k}']

        # Start from base trace dicts, then apply per-frame updates
        snap_data_dicts = [json.loads(t.to_json()) for t in fig.data]
        for ti, td in zip(frame.traces, frame.data):
            frame_td_dict = json.loads(td.to_json())
            snap_data_dicts[ti].update(frame_td_dict)

        # Apply per-frame layout (annotations etc.)
        snap_layout_dict = copy.deepcopy(base_layout_dict)
        if frame.layout:
            frame_layout_dict = json.loads(frame.layout.to_json())
            snap_layout_dict.update(frame_layout_dict)

        snap = go.Figure(data=snap_data_dicts, layout=snap_layout_dict)

        img_bytes = pio.to_image(snap, format='png',
                                 width=base_layout_dict.get('width', 820),
                                 height=base_layout_dict.get('height', 840))
        if _iio_version == 3:
            images.append(iio.imread(_io.BytesIO(img_bytes)))
        else:
            images.append(iio.imread(_io.BytesIO(img_bytes), format='png'))

        if (k + 1) % 20 == 0:
            print(f'  {k + 1}/{nt} frames done')

    ext = os.path.splitext(output_path)[1].lower()
    print(f'Writing {output_path} …')
    if ext == '.gif':
        if _iio_version == 3:
            iio.imwrite(output_path, images, plugin='pillow',
                        duration=int(1000 / fps), loop=0)
        else:
            iio.mimsave(output_path, images, fps=fps)
    else:  # .mp4 / .mov
        if _iio_version == 3:
            iio.imwrite(output_path, images, plugin='ffmpeg', fps=fps)
        else:
            writer = iio.get_writer(output_path, fps=fps)
            for img in images:
                writer.append_data(img)
            writer.close()

    print(f'Saved to {output_path}')


def plot_example_trajectory(preset_idx=0, tf=50.0, nt=600, output_path=None,
                            show=False, figsize=(8, 4.2), ic_label='      ',
                            legend=True, halo=False, point=False,
                            ticks=False, ticklabels=False):
    """
    Plot the configuration-space (x-y) trajectory of Mode 1 and Mode 2 side by
    side for one example preset, as a static matplotlib figure.

    The trajectory is the path the mass traces in the (x, y) displacement plane
    under the pure-mode initial conditions returned by SimplestOscillator.mode_ics.
    A straight line through the origin ⇒ monophase mode; a curved (decaying loop)
    path ⇒ non-monophase mode.  Line color runs from dark (τ=0) to light (τ=tf)
    to show the decay over time.

    Parameters
    ----------
    preset_idx  : int    0 = Ex 1 (Δ=1), 1 = Ex 2 (Δ=1.5), 2 = Ex 3 (ζ₃=0)
    tf          : float   Final dimensionless time τ
    nt          : int     Number of time samples along the trajectory
    output_path : str or None   If given, save a PNG here (e.g. 'ex1_traj.png')
    show        : bool    If True, call plt.show()
    figsize     : tuple   Figure size in inches

    Returns
    -------
    fig, axes : matplotlib Figure and the two Axes (Mode 1, Mode 2)
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.patheffects as pe
    import matplotlib.colors as mcolors
    from matplotlib.collections import LineCollection
    from matplotlib.legend_handler import HandlerPatch

    # Legend handler that draws the IC entry as an arrow (not a plain line).
    def _legend_arrow(legend, orig_handle, xdescent, ydescent,
                      width, height, fontsize):
        return mpatches.FancyArrow(0, 0.5 * height, width, 0,
                                   length_includes_head=True,
                                   head_width=0.5 * height, head_length=0.4 * width,
                                   color='black')

    p   = PRESETS[preset_idx]
    osc = SimplestOscillator(p['z1'], p['z2'], p['z3'], p['D'])
    ts  = np.linspace(0, tf, nt)

    # Solve each mode under its pure-mode initial conditions.
    traj = []
    for k in (1, 2):
        u0, up0 = osc.mode_ics(k)
        sol = SimplestOscillator(p['z1'], p['z2'], p['z3'], p['D'],
                                 list(u0), list(up0)).solve(ts)
        traj.append((np.real(sol[0]), np.real(sol[1])))   # (xs, ys)

    # Unified, square, origin-centered axis limits across both modes.
    all_x = np.concatenate([t[0] for t in traj])
    all_y = np.concatenate([t[1] for t in traj])
    # lim   = 1.05 * max(np.abs(all_x).max(), np.abs(all_y).max(), 1e-3)
    xlim = 2.3
    ylo = -2.8
    yhi = 1.8

    # Time-gradient colormap. For a single solid color instead, set
    #     traj_cmap = None;  traj_color = '#1f77b4'
    # and pass color=traj_color (drop cmap/array) to the LineCollection
    traj_cmap  = mcolors.LinearSegmentedColormap.from_list(
        'traj', plt.get_cmap('viridis')(np.linspace(0.0, 1.0, 256)))

    fig, axes = plt.subplots(2, 1, figsize=figsize)
    for ax, (xs, ys) in zip(axes, traj):
        pts  = np.column_stack([xs, ys]).reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        lc   = LineCollection(segs, cmap=traj_cmap, array=ts[:-1],
                              linewidths=2.0, zorder=2)
        ax.add_collection(lc)
        # Vector from the origin to the starting point (the initial condition).
        ann = ax.annotate('', xy=(xs[0], ys[0]), xytext=(0, 0),
                        arrowprops=dict(arrowstyle='-|>', color='black', lw=1),
                        zorder=6)
        if halo:
            # Thin white halo to keep the arrow legible over the trajectory
            ann.arrow_patch.set_path_effects([pe.withStroke(linewidth=3.5, foreground='white')])
        if point:
            if halo:
                ax.plot(xs[0], ys[0], 'o', color='black', ms=2, zorder=7,
                        path_effects=[pe.withStroke(linewidth=2.5, foreground='white')])
            else:
                ax.plot(xs[0], ys[0], 'o', color='black', ms=2, zorder=7)
        ax.axhline(0, color='gray', lw=0.6, zorder=0)
        ax.axvline(0, color='gray', lw=0.6, zorder=0)
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(ylo, yhi)
        ax.set_aspect('equal')
        if ticks:
            ax.set_xticks([-2, -1, 0, 1, 2])
            ax.set_yticks([-2, -1, 0, 1, 2])
        else:
            ax.set_xticks([])
            ax.set_yticks([])
        if ticklabels:
            ax.set_xticklabels(['', '', '', '', ''])
            ax.set_yticklabels(['', '', '', '', ''])
        # Small coordinate frame (x, y directions) in the lower-left corner
        ox, oy, fl = -1.75, -2.25, 0.8
        for dx, dy in [(fl, 0.0), (0.0, fl)]:
            ax.annotate('', xy=(ox + dx, oy + dy), xytext=(ox, oy),
                        arrowprops=dict(arrowstyle='-|>', color='black', lw=1.2,
                                        shrinkA=0, shrinkB=0),
                        zorder=8)
        if legend:
            # Legend entry for the IC vector (blank label for tikz annotation).
            ic_proxy = mpatches.FancyArrow(0, 0, 1, 0, color='black')
            ax.legend([ic_proxy], [ic_label],
                      handler_map={mpatches.FancyArrow: HandlerPatch(patch_func=_legend_arrow)},
                      loc='lower right', frameon=True, handlelength=2.0)

    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches='tight')
        print(f'Saved {output_path}')
    if show:
        plt.show()
    return fig, axes


if __name__ == "__main__":
    main_simple_figure(speed_up=True)
    names = ['ex1_classical', 'ex2_monophase', 'ex3_non-monophase']
    for i,exname in enumerate(names):
        plot_example_trajectory(i, output_path=f'{exname}_traj.png',legend=False)
        export_plotly_animation(preset_idx=i, output_path=f'{exname}.gif', fps=15, nt=200)
"""
Minimal partially coherent (Hopkins/Abbe) scalar imaging model for periodic masks.
Units: lengths in nm; frequencies in 1/nm. Scalar, thin-mask, no defocus unless specified.

Periodic cell of size (Lx, Ly) with Nx x Ny mask pixels. Mask m real in [-1,1] (PSM) or [0,1] (binary).
Diffraction orders f_k = k/L, k integer. Pupil passes |f| <= NA/lambda.
Source: set of point sources s_j (in direction-cosine units, normalized so |s|<=sigma*NA) with weights w_j.
Abbe: I(x) = sum_j w_j |sum_k P(f_k + s_j/lambda) M(f_k) e^{2 pi i f_k x}|^2
"""
import numpy as np

class Litho1D:
    def __init__(self, L, N, wavelength, NA, source_pts, source_w, defocus=0.0):
        self.L, self.N, self.lam, self.NA = L, N, wavelength, NA
        self.k = np.arange(-(N//2), N - N//2)            # order indices
        self.f = self.k / L                              # spatial frequencies (1/nm)
        self.src = np.asarray(source_pts, float)         # source positions in sigma units * NA (direction cosines)
        self.w = np.asarray(source_w, float) / np.sum(source_w)
        self.x = np.arange(N) * L / N
        self.defocus = defocus
        # DFT matrix: M = F m  (orders from pixels), with pixel integration (sinc) ignored -> point samples
        n = np.arange(N)
        self.F = np.exp(-2j*np.pi*np.outer(self.k, n)/N)/N   # M_k = (1/N) sum_n m_n e^{-2pi i k n/N}
        # Pupil for each source point: pass if |f_k + s_j/lam| <= NA/lam
        self.P = np.zeros((len(self.src), N), complex)
        for j, s in enumerate(self.src):
            g = self.f + s/self.lam                      # shifted frequency (direction cosine / lam)
            inside = np.abs(g)*self.lam <= self.NA
            ph = np.zeros(N)
            if defocus != 0.0:
                # scalar defocus phase: exp(i 2pi/lam * dz * (sqrt(1-(g lam)^2)-1))
                gl = np.clip(g*self.lam, -1, 1)
                ph = 2*np.pi/self.lam*defocus*(np.sqrt(1-gl**2)-1)
            self.P[j] = inside*np.exp(1j*ph)
        # Image sampling matrix: E_{x,k} = e^{2 pi i f_k x}
        self.E = np.exp(2j*np.pi*np.outer(self.x, self.f))
        # Coherent system matrices: a_j[x,:] = E[x,:]*P[j,:]  -> field_j(x) = a_j[x,:] @ M
        self.A = np.stack([self.E*self.P[j][None,:] for j in range(len(self.src))])  # (S, N, N)
        # Pixel-domain quadratic forms Q[x] (N x N real symmetric PSD): I(x) = m^T Q[x] m
        self.Q = self._pixel_forms()

    def orders(self, m):
        return self.F @ m

    def image(self, m):
        M = self.orders(m)
        fields = self.A @ M                              # (S, Nx)
        return np.einsum('j,jx->x', self.w, np.abs(fields)**2).real

    def _pixel_forms(self):
        # I(x) = sum_j w_j |b_j(x)^H m|^2 with b_j(x) = (a_j[x,:] @ F)^*  -> Q[x] = sum_j w_j Re(c c^H) (real m)
        S, N, _ = self.A.shape
        C = np.einsum('jxk,kn->jxn', self.A, self.F)     # (S, Nx, Npix)  : field_j(x) = C[j,x,:] @ m
        Q = np.zeros((N, N, N))
        for x in range(N):
            Qx = np.zeros((N, N), complex)
            for j in range(S):
                c = C[j, x, :]
                Qx += self.w[j]*np.outer(c.conj(), c)
            Q[x] = Qx.real                               # m real -> m^T Re(Qx) m
            Q[x] = 0.5*(Q[x]+Q[x].T)
        return Q

    def image_from_forms(self, m):
        return np.einsum('n,xnm,m->x', m, self.Q, m)

if __name__ == "__main__":
    # sanity: dense lines, coherent, check cutoff pitch = lam/NA (binary mask) 
    lam, NA = 193.0, 1.35
    for pitch in [100, 120, 140, 145, 160, 200]:
        lt = Litho1D(L=pitch, N=32, wavelength=lam, NA=NA, source_pts=[0.0], source_w=[1.0])
        m = (np.arange(32) < 16).astype(float)           # 50% duty binary
        I = lt.image(m)
        I2 = lt.image_from_forms(m)
        print(f"pitch={pitch:4d}  contrast={(I.max()-I.min())/(I.max()+I.min()+1e-12):.3f}  form-consistency={np.abs(I-I2).max():.2e}  cutoff lam/NA={lam/NA:.1f}")

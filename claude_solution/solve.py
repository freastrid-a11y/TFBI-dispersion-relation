"""Independent root finder for Eq. 72: argument-principle count of unstable roots
in a box of the upper half omega-plane + grid-minimum seeding + complex Newton."""
import os, sys, argparse, numpy as np, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import indep as I
from scipy.special import wofz

def Dvec(w, kx, ky, p):
    """Vectorised Eq. 72, valid for Im(omega_i) > 0 (Z-tilde = Z)."""
    w = np.asarray(w, complex); k2 = kx*kx+ky*ky
    wi = w + 1j*p['nu_i']; ka = ky*p['ai0']
    s1 = np.sqrt(2*(k2*p['VTn2']+1j*ka)); d2 = kx**2*p['VTx2']+ky**2*p['VTy2']+1j*ka; s2=np.sqrt(2*d2)
    z2 = (wi-ky*p['U0'])/s2
    Iion = k2*p['VTn2']/(p['betaT']*d2)*(1+z2*I.Z(z2))/(1+1j*p['nu_i']/s1*I.Z(wi/s1))
    kae = ky*p['ae0']; A = k2*p['VTe2']-2j*kae; Bd = A-p['ae0']**2/p['VTe2']
    Dw = kx*p['Vex']+ky*p['Vey']-w; a1 = 5j*kae/A
    q4 = -(3j*p['nu0']*Dw)*I.Om_e**2/(4*p['nu0']**2*Bd)
    r = q4[None]**0.25*np.array([1,1j,-1,-1j])[:,None]
    q = np.take_along_axis(r, np.argmax(np.minimum(r.real,r.imag),0)[None],0)[0]
    Iel = A/Bd*(1+2*a1+(q4+a1*q*q)/q*I.Z(q)-1j*(q4-a1*q*q)/q*I.Z(1j*q))
    return Iion+Iel+k2*p['lD2']

def winding(f, corners, n=4000):
    pts=[]
    for a,b in zip(corners, corners[1:]+corners[:1]):
        pts.append(a+(b-a)*np.linspace(0,1,n,endpoint=False))
    z=np.concatenate(pts); ph=np.angle(f(z))
    d=np.diff(np.concatenate([ph,ph[:1]])); 
    if np.max(np.abs((d+np.pi)%(2*np.pi)-np.pi))>2.0:  # under-resolved
        return winding(f, corners, n*4) if n<256000 else np.nan
    return int(round(np.sum((d+np.pi)%(2*np.pi)-np.pi)/(2*np.pi)))

def newton(f, w, it=60):
    for _ in range(it):
        h=1e-6*max(abs(w),1.0); fw=f(np.array([w]))[0]
        dfw=(f(np.array([w+h]))[0]-f(np.array([w-h]))[0])/(2*h)
        if dfw==0 or not np.isfinite(dfw): return None
        step=fw/dfw; w=w-step
        if abs(step)<1e-10*max(abs(w),1): return w
    return None

def roots_k(kx, ky, p, gmin=1.0):
    k=np.hypot(kx,ky); R=1.5*k*2100+20*p['nu_i']; G=k*2100+10*p['nu_i']
    f=lambda w: Dvec(w,kx,ky,p)
    box=[complex(-R,gmin),complex(R,gmin),complex(R,G),complex(-R,G)]
    N=winding(f,box)
    wr=np.linspace(-R,R,241); wg=np.linspace(gmin,G,121); W=wr[None,:]+1j*wg[:,None]
    A=np.abs(f(W.ravel())).reshape(W.shape)
    loc=(A[1:-1,1:-1]<=np.minimum.reduce([A[:-2,1:-1],A[2:,1:-1],A[1:-1,:-2],A[1:-1,2:]]))
    seeds=W[1:-1,1:-1][loc]
    found=[]
    for s in seeds:
        r=newton(f,s)
        if r is None or not (gmin<r.imag<G and abs(r.real)<R): continue
        if abs(f(np.array([r]))[0])>1e-8*max(1,np.abs(f(np.array([s]))[0])): continue
        if all(abs(r-x)>1e-6*abs(r) for x in found): found.append(r)
    return N, sorted(found,key=lambda w:-w.imag)

def fine_seeds(kx, ky, p, gmin=1.0):
    """Denser seed grid hugging the real axis, for weakly growing roots near
    the gamma = 0 boundary that the coarse grid in roots_k can step over."""
    k = np.hypot(kx, ky); R = 1.5*k*2100 + 20*p['nu_i']; G = 0.1*(k*2100 + 10*p['nu_i'])
    wr = np.linspace(-R, R, 801); wg = np.geomspace(gmin, G, 60)
    W = wr[None, :] + 1j*wg[:, None]
    A = np.abs(Dvec(W.ravel(), kx, ky, p)).reshape(W.shape)
    loc = (A[1:-1, 1:-1] <= np.minimum.reduce([A[:-2, 1:-1], A[2:, 1:-1], A[1:-1, :-2], A[1:-1, 2:]]))
    return list(W[1:-1, 1:-1][loc])


def repair(om, cnt, nf, ks, p, passes=20):
    """Where the argument principle counted more unstable roots than were found,
    re-seed Newton from neighbours, from the mirror point -k
    (omega(-k) = -conj(omega(k))), and finally from a denser seed grid.
    Repeats until nothing changes, so fixes propagate along a missed edge."""
    n = len(ks)
    for _ in range(passes):
        changed = False
        for iy, ix in np.argwhere(cnt > nf):
            kx, ky = ks[ix], ks[iy]
            f = lambda w: Dvec(w, kx, ky, p)
            seeds = [om[iy+a, ix+b] for a in (-1, 0, 1) for b in (-1, 0, 1)
                     if 0 <= iy+a < n and 0 <= ix+b < n]
            seeds.append(-np.conj(om[n-1-iy, n-1-ix]))
            rs = [r for s in seeds if np.isfinite(s)
                  for r in [newton(f, s)] if r is not None and r.imag > 1.0]
            if not rs and not np.isfinite(om[iy, ix]):
                rs = [r for s in fine_seeds(kx, ky, p)
                      for r in [newton(f, s)] if r is not None and r.imag > 1.0]
            if rs:
                best = max(rs, key=lambda w: w.imag)
                if not np.isfinite(om[iy, ix]) or best.imag > om[iy, ix].imag + 1e-9*abs(best):
                    om[iy, ix] = best; changed = True
                nf[iy, ix] = max(nf[iy, ix], len({complex(np.round(r, 3)) for r in rs}))
        if not changed:
            break
    return om


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Independent solve of Eq. 72 on a k grid")
    ap.add_argument("--alt", type=int, required=True, choices=[104, 110],
                    help="effective altitude in km (104 -> run 5, 110 -> run 4)")
    ap.add_argument("--points", type=int, default=100)
    ap.add_argument("--kmax", type=float, default=100.0)
    ap.add_argument("--E-scale", type=float, default=1.0, help="multiply E0 (2.0 = double E)")
    ap.add_argument("--ae0-sign", type=float, default=1.0, choices=[1.0, -1.0])
    ap.add_argument("--out", default=None, help="output .npz (default: out/indep_<alt>.npz)")
    a = ap.parse_args()
    run = {110: 4, 104: 5}[a.alt]; n = a.points
    p = I.params(run, E_scale=a.E_scale, ae0_sign=a.ae0_sign)
    ks = np.linspace(-a.kmax, a.kmax, n)
    om = np.full((n, n), np.nan+0j); cnt = np.full((n, n), -1.0); nf = np.zeros((n, n), int)
    t = time.time()
    with np.errstate(all="ignore"):
        for iy, ky in enumerate(ks):
            for ix, kx in enumerate(ks):
                if kx == 0 and ky == 0: continue
                N, r = roots_k(kx, ky, p); cnt[iy, ix] = N; nf[iy, ix] = len(r)
                if r: om[iy, ix] = r[0]
            print(f"row {iy+1}/{n}  ({time.time()-t:.0f} s)", flush=True)
        om = repair(om, cnt, nf, ks, p)
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", f"indep_{a.alt}.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, omega=om, count=cnt, nfound=nf, k=ks, alt=a.alt)
    print("saved", out)
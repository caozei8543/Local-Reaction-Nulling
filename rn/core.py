"""Equations (2)-(10), (S4)-(S8), (S12)-(S21). SI units; complex linear S."""
import math
import numpy as np
EPS0 = 8.8541878128e-12
F0 = 5.8e9
DESIGNS = ("ND", "DLO", "RN")
CONDITIONS = ("common", "unequal")
PARAMETERS = ("l1", "l2", "d1", "d2", "ws", "xf", "lp")
BOUNDS = np.array([[.8,2.2],[.8,2.2],[.6,1.8],[.6,1.8],[1.5,1.9],[4.5,7.0],[10.5,11.5]])

def polar(mag, deg):
    return complex(mag * np.exp(1j*np.deg2rad(deg)))

def angle(z):
    return float(np.angle(z, deg=True))

def pair(z):
    return [float(np.real(z)), float(np.imag(z))]

def relative_error(value, reference):
    if abs(reference) == 0: raise ValueError("Zero reference; relative error undefined")
    return abs(value-reference)/abs(reference)

def phase_difference(a,b):
    if abs(a)==0 or abs(b)==0: raise ValueError("Phase of zero undefined")
    return abs(angle(a/b))

def db(z):
    with np.errstate(divide='ignore'):
        return 20*np.log10(np.abs(z))

def linear(db_value):
    return 10**(db_value/20)

def material_eta(eps_real, loss_tangent):
    return eps_real*(1-1j*loss_tangent)-1

def reaction(e1,e2,weights,g=1):
    """Complex bilinear dot product (NOT vdot/conjugation), V/m fields and m^3 weights."""
    e1,e2,w = np.asarray(e1,complex),np.asarray(e2,complex),np.asarray(weights,float)
    if e1.ndim!=2 or e1.shape[1]!=3 or e1.shape!=e2.shape or w.shape!=(len(e1),):
        raise ValueError("Require E1/E2 N×3 and N volume weights")
    g=np.broadcast_to(np.asarray(g,complex),w.shape)
    if not all(np.isfinite(v).all() for v in (e1,e2,w,g)) or np.any(w<=0):
        raise ValueError("Nonfinite field/profile or nonpositive volume weight")
    return complex(np.sum(w*g*np.sum(e1*e2,axis=1)))

def derivative_from_reaction(r, frequency_hz=F0, phasor="rms", p1_w=1., p2_w=1.):
    """S5: convert to unit incident power, zero-phase, 50-ohm reference waves."""
    if phasor not in ("rms","peak") or min(frequency_hz,p1_w,p2_w)<=0: raise ValueError("Invalid phasor or power/frequency")
    return -1j*2*np.pi*frequency_hz*EPS0*r/( (2 if phasor=="rms" else 4)*math.sqrt(p1_w*p2_w))

def central_difference(plus,minus,delta):
    if not np.isfinite(delta) or delta<=0: raise ValueError("Positive finite step required")
    return (plus-minus)/(2*delta)

def differential_ratio(dl,dr):
    if abs(dl+dr)==0: raise ValueError("Zero common-mode slope; ratio undefined")
    return abs(dl-dr)/abs(dl+dr)

def delta_linear(dl,dr,eta_l,eta_r):
    return dl*eta_l + dr*eta_r

def common_differential(eta_l,eta_r):
    return (eta_l+eta_r)/2,(eta_l-eta_r)/2

def worst_case_first_order(dl,dr,r):
    if r<0: raise ValueError("Radius must be nonnegative")
    return r*(abs(dl)+abs(dr))

def remainder(s_loaded,s0,dl,dr,eta_l,eta_r):
    return s_loaded-s0-delta_linear(dl,dr,eta_l,eta_r)

def cancellation(a,b):
    ma,mb=abs(a),abs(b)
    if not ma or not mb: raise ValueError("Both contributions must be nonzero")
    phi=180-phase_difference(a,b)
    imbalance=(ma-mb)**2
    phase_term=4*ma*mb*np.sin(np.deg2rad(phi)/2)**2
    total=imbalance+phase_term
    return {"sum_re_im":pair(a+b),"sum_magnitude":abs(a+b),"sum_phase_deg":angle(a+b),"a_magnitude":ma,"b_magnitude":mb,"phase_separation_deg":180-phi,"phase_deficit_deg":phi,"cres":abs(a+b)/(ma+mb),"imbalance_fraction":float(imbalance/total) if total else None}

def total_efficiency(eta_rad,s11,s21):
    """Port 1 incident, other port matched; assumes eta_rad referenced to accepted (non-port) power."""
    return eta_rad*(1-abs(s11)**2-abs(s21)**2)

def feasible(p,s11,s21,s22,eta_rad):
    p=np.asarray(p,float)
    return bool(p.shape==(7,) and np.isfinite(p).all() and np.all(p>=BOUNDS[:,0]) and np.all(p<=BOUNDS[:,1]) and abs(s21)<=linear(-20) and max(abs(s11),abs(s22))<=linear(-15) and .8<=eta_rad<=1)

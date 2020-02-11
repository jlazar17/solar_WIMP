import numpy as np

class Hist():

    def __init__(self, ch, m):

        self.ch   = ch
        self.m    = m

        self.hist     = np.load("/data/user/jlazar/solar_WIMP/data/e_d_theta_hist/ch%d_m%d_f1.000000_f_tot_e_d_theta_hist.npy" % (ch, m))
        self.pmf      = self.hist / np.sum(self.hist)
        self.thetapmf = np.sum(self.hist, axis=0)
        self.epmf     = np.sum(self.hist, axis=1)[::-1]
        self.thetacmf = np.concatenate([[0], np.cumsum(self.thetapmf)])
        self.ecmf     = np.concatenate([[0], np.cumsum(self.epmf)])

    def get_e_cdf(self, val):
        return cdf(val, e_bins, self.epmf)

    def get_theta_cdf(self, val):
        return cdf(val, theta_bins, self.thetapmf)

    def get_e_pmf(self, val):
        return pmf(val, e_bins, self.epmf)


def cmf(pmf,axis=None):
    if axis is None:
        return np.concatenate([[0], np.cumsum(pmf)])
    else:
        shape = np.shape(pmf)
        z_shape = np.copy(shape)
        z_shape[axis] = 1
        z = np.zeros(z_shape)
        return np.concatenate([z, np.cumsum(pmf, axis=axis)], axis=axis)
def pmf(cmf, axis=-1):
    return np.diff(cmf, axis=axis)
def interp(x, x0, x1, y0, y1):
    r = np.empty(shape=np.shape(x))
    mask_0 = x == x0
    mask_1 = x == x1
    mask_2 = ~np.logical_or(mask_0, mask_1)
    r[mask_0] = y0[mask_0]
    r[mask_1] = y1[mask_1]
    r[mask_2] = y0[mask_2]+(x[mask_2]-x0[mask_2])*(y1[mask_2]-y0[mask_2])/(x1[mask_2]-x0[mask_2])
    return r


def cmf_interp_helper(cmf_0, edges_0, y_val, i_0):
    x0s = []
    # Determine where we are with respect to cmf_0
    while i_0 < len(cmf_0)-1 and cmf_0[i_0] < y_val:
        i_0 += 1
    # cmf[i_0] is greater than or equal to y_val
    if cmf_0[i_0] == y_val:
        # Check for duplicate values
        i_0_upper = i_0
        while i_0_upper < len(cmf_0)-1 and cmf_0[i_0_upper] == y_val:
            i_0_upper += 1
        i_0_upper -= 1
        if i_0 == i_0_upper:
            # We don't have duplicate values
            x0s.append(edges_0[i_0])
        else:
            # We have duplicate values
            x0s.append(edges_0[i_0])
            x0s.append(edges_0[i_0_upper])
    else:
        # Use cmf_0[i_0 - 1] for interpolation
        x0 = interp(y_val, cmf_0[i_0-1], cmf_0[i_0], edges_0[i_0-1], edges_0[i_0])
        x0s.append(x0)
    return x0s, i_0

def compute_interp_cmf(cmf_0, cmf_1, edges_0, edges_1, x):
    y_vals = np.unique(np.concatenate([cmf_0, cmf_1]))
    i_0 = 0
    i_1 = 0
    x_pairs = []
    return_y = []
    for y_val in y_vals:
        x0s, i_0 = cmf_interp_helper(cmf_0, edges_0, y_val, i_0)
        x1s, i_1 = cmf_interp_helper(cmf_1, edges_1, y_val, i_1)
        pairs = [(x0, x1) for x0 in x0s for x1 in x1s]
        x_pairs.extend(pairs)
        return_y.extend([y_val]*len(pairs))
    return_x = [x0*(1.0-x) + x1*x for x0, x1 in x_pairs]
    return (np.array(x) for x in zip(*sorted(zip(return_y, return_x), key=lambda x: x[1])))



def rebin_cmf(cmf_0, edges_0, edges_1):
    r = np.empty(shape=len(edges_1))
    mask_less = edges_1 < np.amin(edges_0)
    mask_more = edges_1 >= np.amax(edges_0)
    mask = ~np.logical_or(mask_less, mask_more)
    r[mask_less] = 0.0
    r[mask_more] = 1.0
    
    x_pos = np.digitize(edges_1[mask], bins=edges_0) - 1
    y_vals = interp(edges_1[mask],
           edges_0[x_pos], edges_0[x_pos+1],
           cmf_0[x_pos], cmf_0[x_pos])
    r[mask] = y_vals
    return r

def interp_e_d_theta(ch, m0, m1, x, reverse_order=False):
    
    if m0 > m1:
        tmp = m1
        m1  = m0
        m0  = tmp
        
    h0   = Hist(ch, m0)
    h1   = Hist(ch, m1)
    n0   = np.sum(h0.hist)
    n1   = np.sum(h1.hist)
    norm = (n1-n0)*x+n0
    
    if reverse_order:
        n           = 360
        m           = 60
        pmf0        = np.sum(h0.pmf, axis=1)[::-1]
        pmf1        = np.sum(h1.pmf, axis=1)[::-1]
        edges       = np.logspace(0.5, 6.5, m+1)
        edges_prime = np.linspace(0, 180, n+1)
        we0         = h0.pmf
        we1         = h1.pmf
        
    else:
        n           = 60
        m           = 360
        pmf0        = np.sum(h0.pmf, axis=0)
        pmf1        = np.sum(h1.pmf, axis=0)
        edges       = np.linspace(0, 180, m+1)
        edges_prime = np.logspace(6.5, 0.5, n+1)
        we0         = h0.pmf.T
        we1         = h1.pmf.T
        
    int_cmf, int_edges = compute_interp_cmf(cmf(pmf0), cmf(pmf1), edges, edges, x)
    int_cmf = rebin_cmf(int_cmf, int_edges, edges)
    int_pmf = pmf(int_cmf)
#    plt.plot(edges[:-1], int_pmf)
#    plt.semilogx()
#    plt.show()
    
    
    
    
    pmfs0 = []
    pmfs1 = []
    for arr0, arr1 in zip(we0, we1):
        s0 = np.sum(arr0)
        s1 = np.sum(arr1)
        if s0 !=0:
            pmf0 = arr0 / np.sum(arr0)
            pmfs0.append(pmf0)
        else:
            pmfs0.append(None)
        if s1 !=0:
            pmf1 = arr1 / np.sum(arr1)
            pmfs1.append(pmf1)
        else:
            pmfs1.append(None)
    
    interp_pmfs = []
    for pmf0_, pmf1_ in zip(pmfs0, pmfs1):
        if pmf0_ is None and  pmf1_ is None:
            interp_pmfs.append(np.zeros(n))
        elif pmf0_ is None:
            interp_pmfs.append(pmf1_*x)
        elif pmf1_ is None:
            interp_pmfs.append(pmf0_*(1-x))
        else:
            cmf_, edges_ = compute_interp_cmf(cmf(pmf0_), cmf(pmf1_), edges_prime, edges_prime, x)
            cmf_ = rebin_cmf(cmf_, edges_, edges_prime)
            p = pmf(cmf_)
            p = np.where(p<0, 0,p )
            interp_pmfs.append(p)

    h = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            if not reverse_order:
                h[i,j] = int_pmf[i]*interp_pmfs[i][j]
            else:
                h[i,j] = int_pmf[i]*interp_pmfs[-i][j]

            
    if reverse_order:
        h =h[::-1]
    else:
        h = h.T
    h = np.where(np.isnan(h), 0, h)
#    plt.imshow(
#           np.log10(5*h)[18:55, :81],
##            np.log10(5*h),
#           extent=[0,20, 1,4.7],
#           aspect="auto",
#           vmin=-4,
#           vmax=0,
#           )
#    plt.show()
    return h*norm
        

#if __name__==__main__:
    

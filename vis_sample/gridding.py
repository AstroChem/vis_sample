import numpy as np

# assumes m = 6
def spheroid(eta=0, alpha=1):
    # alpha can be [0, 0.5, 1, 1.5, 2]
    spheroid_dict_1 = {
        0 : spheroid_00_1,
        0.5 : spheroid_05_1,
        1 : spheroid_10_1,
        1.5 : spheroid_15_1,
        2 : spheroid_20_1,
    }

    spheroid_dict_2 = {
        0 : spheroid_00_2,
        0.5 : spheroid_05_2,
        1 : spheroid_10_2,
        1.5 : spheroid_15_2,
        2 : spheroid_20_2,
    }

    etalim = 0.75  # specific for m = 6

    eta = abs(eta)

    if eta <= etalim:
        nn = eta**2 - etalim**2
        return spheroid_dict_1[alpha](nn)

    elif eta <= 1.00000000001:
        nn = eta**2 - 1.0
        return spheroid_dict_2[alpha](nn)

    else:
        print "The spheroid is only defined for alpha = 0.0, 0.5, 1.0, 1.5, and 2.0"

    return 0.

def corrfun_single(eta, alpha):
    return spheroid(eta, alpha)

def corrfun(etas, alpha):
    return [corrfun_single(eta, alpha) for eta in etas]

def gcffun_single(eta, alpha):
    return ((abs(1 - eta**2))**alpha)*spheroid(eta, alpha)

def gcffun(etas, alpha):
    return [gcffun_single(eta, alpha) for eta in etas]


# Apply the correction function to the image.
def apply_corrfun(img, alpha, mu_RA, mu_DEC, corr="", return_corr=False):
    if (corr==""):
        ndec, nra, nvel = img.data.shape
        # The size of one half-of the image.
        # sometimes ra and dec will be symmetric about 0, othertimes they won't
        # so this is a more robust way to determine image half-size
        del_ra = abs(img.ra[1] - img.ra[0])
        del_dec = abs(img.dec[1] - img.dec[0])
        maxra = del_ra * nra/2
        maxdec = del_dec * ndec/2

        # If the image will be later offset via a phase shift, then this means that
        # the corrfunction will need to be applied *as if the image were already
        # offset.*

        corr = np.zeros([ndec, nra])
        for i in range(nra):
            for j in range(ndec):
                eta_x = (img.ra[i] + mu_RA)/maxra
                # not sure why -del_dec needs to be applied - has to do with pixel shift I think
                eta_y = (img.dec[j] - del_dec + mu_DEC)/maxdec
                if (abs(eta_x) <= 1.0 and abs(eta_y) <= 1.0):
                    corr[j, i] = 1.0/(corrfun_single(eta_x, alpha) * corrfun_single(eta_y, alpha))

        for k in range(nvel):
            img.data[:,:,k] = img.data[:,:,k]*corr

        if (return_corr==True):
            return corr
    else:
        for k in range(nvel):
            img.data[:,:,k] = img.data[:,:,k]*corr


def spheroid_00_1(nn):
    p_0 = [0.3303194, -0.6324887, 0.6256387, -0.3019847, 0.05613913]
    p_1 = [0.2535284, 0.9077644, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))

def spheroid_05_1(nn):
    p_0 = [0.27657, -0.5829747, 0.6302307, -0.3342119, 0.06843713]
    p_1 = [0.22914, 0.8626056, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))

def spheroid_10_1(nn):
    p_0 = [0.2312756, -0.5335581, 0.627866, -0.3644705, 0.08203343]
    p_1 = [0.2078043, 0.8212018, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))

def spheroid_15_1(nn):
    p_0 = [0.1934013, -0.485747, 0.6197133, -0.3922489, 0.09675562]
    p_1 = [0.1890848, 0.7831755, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))

def spheroid_20_1(nn):
    p_0 = [0.1618978, -0.4405326, 0.6069622, -0.4172349, 0.1124069]
    p_1 = [0.1726085, 0.7481828, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))

def spheroid_00_2(nn):
    p_0 = [0.07747182, -0.1109391, 0.06888533, -0.01616105, 0.0008531865]
    p_1 = [0.3858544, 1.10127, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))

def spheroid_05_2(nn):
    p_0 = [0.07094106, -0.1170228, 0.08595213, -0.02558954, 0.00206076]
    p_1 = [0.3337648, 1.025431, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))

def spheroid_10_2(nn):
    p_0 = [0.06412774, -0.1201436, 0.1021332, -0.03697768, 0.004028559]
    p_1 = [0.2918724, 0.9599102, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))

def spheroid_15_2(nn):
    p_0 = [0.0574421, -0.1207733, 0.1168451, -0.04994202, 0.006887946]
    p_1 = [0.2575337, 0.9025276, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))

def spheroid_20_2(nn):
    p_0 = [0.05112822, -0.1194208, 0.1297386, -0.06404749, 0.01071895]
    p_1 = [0.2289667, 0.851747, 1.0]
    return (np.polyval(p_0,nn)/np.polyval(p_1,nn))





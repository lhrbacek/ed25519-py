curve25519_P = 2 ** 255 - 19
curve25519_A = 486662

# from https://gist.github.com/nickovs/cc3c22d15f239a2640c185035c06f8a3
def _point_add(point_n, point_m, point_diff):
    """Given the projection of two points and their difference, return their sum"""
    (xn, zn) = point_n
    (xm, zm) = point_m
    (x_diff, z_diff) = point_diff
    x = (z_diff << 2) * (xm * xn - zm * zn) ** 2
    z = (x_diff << 2) * (xm * zn - zm * xn) ** 2
    return x % curve25519_P, z % curve25519_P


# from https://gist.github.com/nickovs/cc3c22d15f239a2640c185035c06f8a3
def _point_double(point_n):
    """Double a point provided in projective coordinates"""
    (xn, zn) = point_n
    xn2 = xn ** 2
    zn2 = zn ** 2
    x = (xn2 - zn2) ** 2
    xzn = xn * zn
    z = 4 * xzn * (xn2 + curve25519_A * xzn + zn2)
    return x % curve25519_P, z % curve25519_P


# from https://gist.github.com/nickovs/cc3c22d15f239a2640c185035c06f8a3
def _const_time_swap(a, b, swap):
    """Swap two values in constant time"""
    index = int(swap) * 2
    temp = (a, b, b, a)
    return temp[index:index+2]


# from https://gist.github.com/nickovs/cc3c22d15f239a2640c185035c06f8a3
def _raw_curve25519(base, n):
    """Raise the point base to the power n"""
    zero = (1, 0)
    one = (base, 1)
    mP, m1P = zero, one

    for i in reversed(range(256)):
        bit = bool(n & (1 << i))
        mP, m1P = _const_time_swap(mP, m1P, bit)

        if i >= 247:
            print("x1: ", hex(mP[0]))
            print("z1: ", hex(mP[1]))
            print("x2: ", hex(m1P[0]))
            print("z2: ", hex(m1P[1]))
            print("--------------------")

        mP, m1P = _point_double(mP), _point_add(mP, m1P, one)
        mP, m1P = _const_time_swap(mP, m1P, bit)

    x, z = mP
    x1, z1 = m1P
    #inv_z = pow(z, curve25519_P - 2, curve25519_P)
    #return (x * inv_z) % curve25519_P
    return x, z, x1, z1


# from https://gist.github.com/nickovs/cc3c22d15f239a2640c185035c06f8a3
# base - x coordinate of point as bytearray/bytes in big endian (msb on left)
# n - scalar
#
# Usage:
#    point = bytearray(bytes.fromhex('7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed'))
#    g = scamult_once_curve25519(point, 0)
#    print(hex(g))
#
def scamult_once_curve25519(base, n):
    """Raise the point base to the power n, return x of mP without conoversion back to affine"""
    Px = int.from_bytes(base, byteorder='big')
    zero = (1, 0)
    one = (Px, 1)
    mP, m1P = zero, one

    for i in reversed(range(1)):
        bit = bool(n & (1 << i))
        mP, m1P = _const_time_swap(mP, m1P, bit)
        mP, m1P = _point_double(mP), _point_add(mP, m1P, one)
        mP, m1P = _const_time_swap(mP, m1P, bit)

    x, z = mP
    x1, z1 = m1P
    #inv_z = pow(z, curve25519_P - 2, curve25519_P)
    #return (x * inv_z) % curve25519_P
    return x

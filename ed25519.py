from Crypto.PublicKey import ECC
from Crypto.Hash import SHA512
from Crypto.Signature import eddsa
from cryptography.hazmat.primitives.asymmetric import ed25519
from cryptography.hazmat.primitives import serialization

from curve25519 import _raw_curve25519

class Ed25519:
    p: int = 2 ** 255 - 19
    p_bytes: bytes = bytes.fromhex('7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed')
    a_bytes: bytes = bytes.fromhex('7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffec')
    d_bytes: bytes = bytes.fromhex('52036cee2b6ffe738cc740797779e89800700a4d4141d8ab75eb4dca135978a3')
    G: ECC.EccPoint = ECC.EccPoint(int.from_bytes(bytes.fromhex('216936D3CD6E53FEC0A4E231FDD6DC5C692CC7609525A7B2C9562D608F25D51A'), byteorder="big"),
                                   int.from_bytes(bytes.fromhex('6666666666666666666666666666666666666666666666666666666666666658'), byteorder="big"),
                                   curve='Ed25519')
    G_ma_u: int = 9
    G_ma_v: int = 14781619447589544791020593568409986887264606134616475288964881837755586237401
    L_bytes: bytes = bytes.fromhex('1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed')
    scaling_factor_pos: int = 16416487832636737118837039172820900612695230415163812779824790760673067034857  # TODO https://math.stackexchange.com/questions/1392277/point-conversion-between-twisted-edwards-and-montgomery-curves
    scaling_factor_neg: int = -16416487832636737118837039172820900612695230415163812779824790760673067034857 # TODO https://math.stackexchange.com/questions/1392277/point-conversion-between-twisted-edwards-and-montgomery-curves


# Taken and modified from RFC 8032
def point_compress(P: ECC.EccPoint, p: int):
    return int.to_bytes((int(P.y) % p) | (((int(P.x) % p) & 1) << 255), length=32, byteorder='little')


# Taken and modified from RFC 8032
def modp_inv(x):
    return pow(x, Ed25519.p - 2, Ed25519.p) # x^(p-2) % p


# Square root of -1 from RFC 8032
modp_sqrt_m1 = pow(2, (Ed25519.p - 1) // 4, Ed25519.p)


# Compute corresponding x-coordinate, with low bit corresponding to sign, or return None on failure
# Taken and modified from RFC 8032
def recover_x(y, sign):
    p = Ed25519.p
    if y >= p:
        return None
    x2 = (y*y-1) * modp_inv(int.from_bytes(Ed25519.d_bytes, byteorder='big')*y*y+1)
    if x2 == 0:
        if sign:
            return None
        else:
            return 0

    # Compute square root of x2
    x = pow(x2, (p+3) // 8, p)
    if (x*x - x2) % p != 0:
        x = x * modp_sqrt_m1 % p
    if (x*x - x2) % p != 0:
        return None

    if (x & 1) != sign:
        x = p - x
    return x


# Taken and modified from RFC 8032
def point_decompress(s):
    if len(s) != 32:
        raise Exception("Invalid input length for decompression")
    y = int.from_bytes(s, "little")
    sign = y >> 255
    y &= (1 << 255) - 1

    x = recover_x(y, sign)
    if x is None:
        return None
    else:
        return ECC.EccPoint(x, y, curve='Ed25519')


def point_conversion_ea_ma(P: ECC.EccPoint):
    mu = ((1 + int(P.y)) * pow(1 - int(P.y), Ed25519.p - 2, Ed25519.p)) % Ed25519.p
    mv = ((1 + int(P.y)) * pow(((1 - int(P.y)) * int(P.x)), Ed25519.p - 2, Ed25519.p)) % Ed25519.p
    mv = mv*pow(Ed25519.scaling_factor_pos, Ed25519.p - 2, Ed25519.p) % Ed25519.p
    return mu, mv


# just one inversion
def point_conversion_mp_ea(U, V, W):
    Vs = (V * Ed25519.scaling_factor_pos) % Ed25519.p

    U_plus_W = (U + W) % Ed25519.p
    
    # Montgomery trick
    T = (Vs * U_plus_W) % Ed25519.p
    R = pow(T, Ed25519.p - 2, Ed25519.p) # T^-1

    x = (U * R * U_plus_W) % Ed25519.p # U / Vs = U * R*(U+W)
    y = ((U - W) * R * Vs) % Ed25519.p # (U-W) / (U+W) = (U-W) * R*Vs

    return x, y


def point_conversion_mp_ma(U, V, W):
    inv_W = pow(W, Ed25519.p - 2, Ed25519.p) 
    u = (U * inv_W) % Ed25519.p
    v = (V * inv_W) % Ed25519.p
    return u, v


def point_conversion_ma_ea(u, v):
    x = (u * pow(v, Ed25519.p - 2, Ed25519.p)) % Ed25519.p
    y = ((u - 1) * pow(u + 1, Ed25519.p - 2, Ed25519.p)) % Ed25519.p
    x = x*pow(Ed25519.scaling_factor_pos, Ed25519.p - 2, Ed25519.p) % Ed25519.p
    return x, y


# According to https://eprint.iacr.org/2017/212.pdf Algorithm 5
# which is according to Okeya-Sakurai y-coordinate recovery algorithm 1
def y_recovery(x, y, X1, Z1, X2, Z2): # P=(x,y,1), [k]P=(X1:Z1), [k+1]P=(X2:Z2)
    t1 = (x * Z1)  % Ed25519.p
    t2 = (X1 + t1) % Ed25519.p
    t3 = (X1 - t1) % Ed25519.p
    t3 = (t3 * t3) % Ed25519.p
    t3 = (t3 * X2) % Ed25519.p

    t1 = (2*486662 * Z1) % Ed25519.p
    t2 = (t2 + t1)       % Ed25519.p
    t4 = (x * X1)        % Ed25519.p
    t4 = (t4 + Z1)       % Ed25519.p
    t2 = (t2 * t4)       % Ed25519.p

    t1 = (t1 * Z1)    % Ed25519.p
    t2 = (t2 - t1)    % Ed25519.p
    t2 = (t2 * Z2)    % Ed25519.p
    recY1 = (t2 - t3) % Ed25519.p
    t1 = (2*1 * y)    % Ed25519.p

    t1 = (t1 * Z1)    % Ed25519.p
    t1 = (t1 * Z2)    % Ed25519.p
    recX1 = (t1 * X1) % Ed25519.p
    recZ1 = (t1 * Z1) % Ed25519.p

    return recX1, recY1, recZ1


# b=256
# 3 hashes, 1 scamult, 1 coordinates conversion (1 inversion), 2 or 3 modulos, arithmetics in the end
def my_eddsa_sign(B: ECC.EccPoint, priv_pub_key: bytearray, M: bytes):
    #print("========== Python Sign() ==========")

    priv_k = priv_pub_key[:32]
    A_comp = priv_pub_key[32:]

    L = int.from_bytes(Ed25519.L_bytes, byteorder='big')

    h = SHA512.new(data=priv_k).digest()                                    # 1) h = H(priv_k)
    s = int.from_bytes(h[0:32], byteorder='little')                         # s = h[0:32]
                                                                            # 1.5) pruning according to RFC8032
    
    #print("s: ", hex(s))

    s &= (1 << 254) - 8                                                     # clear the lowest three bits of the first octet
    s |= (1 << 254)                                                         # set the second highest bit of the last octet (the highest bit is already clear I suppose)

    #print("s cleared: ", hex(s))

    r = SHA512.new(data=(h[32:64] + M)).digest()                            # 2) r = H(k[32:64] || M)
    #print("r: ", r.hex())
    r = int.from_bytes(r, byteorder='little')
    r = r % L                                                               # for efficiency according to RFC8032 (and now it has to be here because of mont scamult)
    #print("r: ", hex(r))
    #r = 43690
    #print("r_mod_l: ", bytearray(int.to_bytes(r, length=32, byteorder="little")).hex())

                                                                            # 3) R = [r]B
    R = r * B                                                               # just for comparison with library

    Bmu = Ed25519.G_ma_u                                                    # precomputed generator G in Montgomery affine
    Bmv = Ed25519.G_ma_v
    rBmU, rBmW, r1BmU, r1BmW = _raw_curve25519(Bmu, r)                      # Curve25519 scamult, returns [r]B, [r+1]B, in Montgomery projective
    rBmU, rBmV, rBmW = y_recovery(Bmu, Bmv, rBmU, rBmW, r1BmU, r1BmW)       # recovery of y (V) coordinate of [r]B
    Rx, Ry = point_conversion_mp_ea(rBmU, rBmV, rBmW)                       # conversion of [r]B to Twisted Edwards affine

    R = ECC.EccPoint(Rx, Ry, curve='Ed25519')                               
    
    #print("R edwards (pycryptodome) scamult:    ", int(R.x), int(R.y))
    #print("R montgomery (my conversion) scamult:", Rx, Ry)
    print('x_ea', bytearray(int.to_bytes(Rx, length=32, byteorder="little")).hex())
    print('y_ea', bytearray(int.to_bytes(Ry, length=32, byteorder="little")).hex())
    #print('Rx_ea', bytearray(int.to_bytes(int(R.x) % Ed25519.p, length=32, byteorder="little")).hex())
    #print('Ry_ea', bytearray(int.to_bytes(int(R.y) % Ed25519.p, length=32, byteorder="little")).hex())


    R_comp = point_compress(R, Ed25519.p)                                   # 3.5) encoded R' = r*B
    #print('R_comp:', R_comp.hex())
    
    k = SHA512.new(data=(R_comp + A_comp + M)).digest()                     # 4) k = H(R'||A'||M)
    #print("k: ", k.hex())
    k = int.from_bytes(k, byteorder='little') % L                           # modulo for efficiency according to RFC8032
    
    #print("k_mod_l: ", hex(k))
    k_mul_s = (k * s) % L
    #print("k_mul_s: ", hex(k_mul_s))
    S = (r + k_mul_s) % L
    S = (r + k * s) % L                                                     # 5) S = (r + H(R'||A'||M)*s) mod l
    
    return (R_comp + (int.to_bytes(S, length=32, byteorder='little')))      # 6) return (R, S) concatenated


def my_eddsa_verify(signature: bytes, A_comp: bytes, B: ECC.EccPoint, M: bytes):
    R_comp = signature[0:32]                                     # 1) R' = signature[0:32]
    R = point_decompress(R_comp)                                 # R = decode R' only to check if point encoding is valid
    if R is None:
        raise Exception('Error when decompressing in verify')
    
    S = int.from_bytes(signature[32:64], byteorder='little')     # S = signature[32:64]
    if S >= int.from_bytes(Ed25519.L_bytes, byteorder='big'):
        raise Exception('Error in verify, S >= L')

    A = point_decompress(A_comp)                                 # A = decode A'
    if A is None:
        raise Exception('Error when decompressing in verify')
    
    k = SHA512.new(data=(R_comp + A_comp + M)).digest()          # 2) k = H(R'||A'||M)
    k = int.from_bytes(k, byteorder='little')

    SB = S * B                                                   # SB = [S]B
    kA = k * A                                                   # kA = [k]A
    return SB == (R + kA)                                        # 3) [S]B == R + [k]A



if __name__ == "__main__":
    testvec1_priv_k = bytes.fromhex('9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60')
    cryptography_private_key = ed25519.Ed25519PrivateKey.from_private_bytes(testvec1_priv_k)
    cryptodome_key = eddsa.import_private_key(testvec1_priv_k)
    #cryptography_private_key = ed25519.Ed25519PrivateKey.generate()
    #cryptodome_key = eddsa.import_private_key(cryptography_private_key.private_bytes(serialization.Encoding.Raw, serialization.PrivateFormat.Raw, serialization.NoEncryption()))

    message = b'this is mesidz'
    message = b''
    message = bytes.fromhex('aa11ccdd')
    #message = b''

    print('msg:', message.hex())

    A_comp =  point_compress(cryptodome_key.pointQ, Ed25519.p)


    print('===MY WHOLE===')
    my_priv_key = bytearray(cryptodome_key.seed)
    #print("my_priv_key:", bytes(my_priv_key).hex())
    #print("A_comp:", A_comp.hex())
    h = SHA512.new(data=my_priv_key).digest()
    #print("hash:", h.hex())
    my_priv_pub_key = my_priv_key.copy()
    my_priv_pub_key.extend(bytearray(A_comp))
    #print("my_priv_pub_key:", bytes(my_priv_pub_key).hex())
    my_signature = my_eddsa_sign(Ed25519.G, my_priv_pub_key, message)
    print(my_eddsa_verify(my_signature, A_comp, Ed25519.G, message))


    print('===CRYPTODOME WHOLE===')
    signer = eddsa.new(cryptodome_key, 'rfc8032')
    signature = signer.sign(message)
    verifier = eddsa.new(cryptodome_key, 'rfc8032')
    try:
        verifier.verify(message, signature)
        print("The message is authentic.")
    except ValueError:
        print("The message is not authentic.")


    print('===CRYPTODOME SIGN MY VERIFY===')
    print(my_eddsa_verify(signature, A_comp, Ed25519.G, message))


    print('===MY SIGN CRYPTODOME VERIFY===')
    try:
        verifier.verify(message, my_signature)
        print("The message is authentic.")
    except ValueError:
        print("The message is not authentic.")


    print('===CRYPTOGRAPHY WHOLE===')
    cryptography_signature = cryptography_private_key.sign(message)
    public_key = cryptography_private_key.public_key()
    public_key.verify(cryptography_signature, message)
    print('OK if no exception')


    print("-----------------------")
    print('My signature:           ', my_signature.hex())
    print('Pycryptodome signature: ', signature.hex())
    print('Cryptography signature: ', cryptography_signature.hex())


    print("-----------------------")
    #print('scaling_factor_pos:', bytearray(int.to_bytes(Ed25519.scaling_factor_pos % Ed25519.p, length=32, byteorder="little")).hex())
    #print('scaling_factor_neg:', bytearray(int.to_bytes(Ed25519.scaling_factor_neg % Ed25519.p, length=32, byteorder="little")).hex())
    #print(hex(Ed25519.scaling_factor_pos % Ed25519.p))
    #print(hex(Ed25519.scaling_factor_neg % Ed25519.p))
    #hash_in = bytes.fromhex('0102030000000000000000000000000000000000000000000000000000ccbbaa')
    #print(SHA512.new(data=(hash_in)).digest().hex())





import random

 R.<X> = Integers(0)[]

def convolution(f,g,N):
    #This could also be done by list comprehension for better performance
    return f*g %(X^N -1)

def generate_key(dpos,dneg,N):
    # We first randomly place dpos 1s and then dneg -1s,
    # the rest are 0 
    poly = [0]*N
    while dpos > 0:
        i = randrange(N)
        if poly[i] == 0:
            poly[i] = 1
            dpos -= 1
    while dneg >0:
        i = randrange(N)
        if poly[i] == 0:
            poly[i] = -1
            dneg -= 1
    return R(poly)

def poly_mod(f,q):
    # This can be done by list comprehension for better performance
    Rq.<X> = Integers(q)[]
    coefflist = f.list()
    return Rq(coefflist)

def center_lift(f,p):
    # Define the center lift of a polynomial
    coefflist = R(f).list()
    return R([((coefflist[i] + p//2) %p) - p//2 for i in range(len(coefflist))])

def inverse_prime(f,j,N):
    #find the inverse of the polynomialm f modulo a prime j
    # We could also do this using the extended euclidian algorithm 
    # since // is implemented modulo primes
    Rj.<X> = Integers(j)[]
    RR = Rj.quotient(X^(N) -1)
    return R(lift(RR(f.list())^(-1)))

def inverse_prime_pow(f,t,q,N):
    # find the inverse of the polynomialm f modulo a prime power q of t
    # The idea here is that we first find an inverse g mod t. 
    # Then we know that g(x)*f(x) = 1 - tr(x). 
    # This means that (1 + tr(x))g(x)*f(x) = 1 - t^2r^2(x)
    # We iterate this process until we get g'(x)f(x) = 1 - qr^m(x), 
    # and g'(x) is the inverse of f(x) mod q
    
    assert q.is_power_of(t)
    g = inverse_prime(f,t,N)
    assert poly_mod(convolution(f,g,N),t) == 1
    while True:
        r = center_lift(convolution(f,g,N),q)
        if r == 1:
            return g
        g = center_lift(convolution(g,2-r,N),q)

def inverse(f,q,N):
    # find the inverse of f modulo q where q is composite
    # We first find the inverse modulo the prime powers composing q
    # Then we use the CRT to find the inverse mod the product of those prime powers, aka q
    
    # prime factorization of q
    fact = factor(q)
    factors = list(fact)
    
    # inverses mod prime powers
    inverse_dict = {}
    for l in factors:
        inverse_dict[l[0]] = inverse_prime_pow(f,l[0],(l[0]^l[1]),N)
        
    # CRT to find the inverse mod q
    inverse = 0
    for i in factors:
        Ni = q // (i[0]^i[1])
        inverse += inverse_dict[i[0]]*Ni*xgcd(Ni,i[0]^i[1])[1]
    return inverse

def NTRU_setup(N, p, q, d):

    # check that we don't have to create impossible polynomials
    assert 2*d + 1 <= N

    # We could check that the decryption will always happen accurately by writing:
    # assert q > (6*d+1)*p
    # But in practice the decryption is accurate almost all the time even if q < (6*d+1)*p.
    # This is because it's highly unlikely that all the coefficients will line up to attain that bound 
    
    #generate key f invertible mod p and mod q of degree N-1 with coefficients 1, 0, -1 
    while True:
        try:
            f = generate_key(d+1,d,N)
            Fp = inverse(f,p,N)
            Fq = inverse(f,q,N)
            break
        except: pass
    
    assert poly_mod(convolution(f,Fp,N),p) == 1
    assert poly_mod(convolution(f,Fq,N),q) == 1

    #generate key g with coefficients 1,0,-1
    g = generate_key(d,d,N)
    
    #create the public key
    h = poly_mod(convolution(Fq,g,N),q)

    return f, Fp, h

def NTRU_encrypt(m,p,q,N,h):

    #generate a random perturbation
    r = generate_key(d,d,N)
    
    #return the encrypted message
    return poly_mod(p*convolution(h,r,N) + m,q)

def NTRU_decrypt(f,Fp,e,p,q,N):

    a = R(center_lift(convolution(f,e,N),q))
    return center_lift(convolution(a,Fp,N),p)

#change parameters to play with it as you see fit :)
# the current parameters are considered "safe" fi we assume that the key generation is truly random
p = 3
q = 4096
R.<X> = Integers(0)[]
N = 821
d = 271
f_Fp_h = NTRU_setup(N,p,q,d)
f = f_Fp_h[0]
Fp = f_Fp_h[1]
h = f_Fp_h[2]
#print(h)
m = generate_key(d,d,N)
print(m)
e = NTRU_encrypt(m,p,q,N,h)
#print(e)
NTRU_decrypt(f,Fp,e,p,q,N)

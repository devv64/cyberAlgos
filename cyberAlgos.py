import time
import math
import hashlib
import secrets
from random import randrange

def gcdExt(a, b):
    global x, y
    if a == 0:
        x = 0
        y = 1
        return b
    gcd = gcdExt(b%a, a)
    x1 = x
    y1 = y
    x = y1 - (b//a) * x1
    y = x1
    return gcd

def isNotPrime(a):
    for i in range(2,a):
        if (a%i) == 0:
            return True
    return False

def modInv(a,b):
    gcdExt(b,a)
    # if isNotPrime(a):
    #     return "a is not prime"
    if (0>=b or b>=a):
        return "b is less than 0 or b is greater than a"
    else:
        return (x%a+a)%a

naiveTime = 0
efficientTime = 0

def naiveApproach(a, b, n):
    start_time = time.time()
    x = 1
    for i in range(b):
        x = (x*a) % n
    print("Naive approach takes %s seconds" % (time.time() - start_time))
    print("It returns: ", x)
    return x

def efficientApproach(a, b, n):
    start_time = time.time()
    x = 1
    a = a%n
    if(a==0):
        return 0
    while b:
        if((b&1)==1):
            x = (x*a) % n
        b = b>>1
        a = (a*a) % n
    efficientTime = time.time() - start_time
    print("Efficient approach takes %s seconds" % (time.time() - start_time))
    print("It returns: ", x)
    return x

def main(a,b,n):
    naive_time = time.time()
    naiveApproach(a,b,n)
    naiveTime = time.time() - naive_time
    efficient_time = time.time()
    efficientApproach(a,b,n)
    efficientTime = time.time() - efficient_time
    print("The difference in computation time is: ")
    print(naiveTime - efficientTime)

def BSGS(g, b, p):
    n = int(math.sqrt(p) + 1)
    res = 1
    for i in range(n):
        res = (res * g) % p
    store = [0] * p
    point = res
    for i in range(1, n+1):
        if(store[point] == 0):
            store[point] = i
            point = (point * res) % p
    print(point)
    point = b
    for i in range(n+1):
        if(store[point]>0):
            ans = store[point] * n - i
            if (ans < p):
                return ans
        point = (point * g) % p
    return -1

def sieve(n, prime) :
    p = 2
    while( p * p <= n ):
        if (prime[p] == True) :
            for i in range(p * 2, n, p) :
                prime[i] = False       
        p += 1
                  
def SophieGermainNumbers(n) :
    nums = 0
    prime = [True]*(2 * n + 1)
    sieve(2 * n + 1, prime)
    for i in range(2, n + 1) :
        if(len(prime) > (2 * i + 1)):
            if (prime[i] and prime[2 * i + 1]):
                nums = i
    return nums

def smallestG(q, p, ZpS):
    for g in ZpS:
        if efficientApproach(g, q, p) == 1: 
            return g

def calcXY(g, a, b, p):
    x = int(pow(g,a,p)) 
    y = int(pow(g,b,p))
    return x,y

def main():
    n = randrange(1,1000)
    P = SophieGermainNumbers(n)
    Q = (SophieGermainNumbers(n) * 2 + 1)
    print(SophieGermainNumbers(n))
    ZpS = list(map(lambda x: x, range(2,P)))
    ZqS = list(map(lambda x: x, range(2,Q)))
    g = smallestG(P, Q, ZpS)
    print(g)
    a = randrange(1,1000)
    b = randrange(1,1000)
    x = 0
    y = 0
    xy = calcXY(g, a, b, P)
    x = xy[0]
    y = xy[1]
    xb = pow(x, b) % P
    ya = pow(y, a) % P
    if xb == ya:
        print("Verified that x^b mod P = y^a mod P!")

def isPrime(n):
    if (n <= 1):
        return False
    for i in range(2, int(math.sqrt(n))+1):
        if (n % i == 0):
            return False
    return True

def genPQ():
    PQ = 0
    while not isPrime(PQ):
        PQ = randrange(pow(2,31), pow(2,32))
    return PQ

P = genPQ()
Q = genPQ()

N = P * Q

PHIN = (P-1) * (Q-1)

def computeGCD(x, y):
    while(y):
       x, y = y, x % y
    return abs(x)

def rande(PHIN):
    bool = True
    while bool:
        e = randrange(1,PHIN)
        if computeGCD(e, PHIN) == 1:
            bool = False
            return e

e = rande(PHIN)
d = modInv(PHIN, e)

def bitsof(bt, nbits):
    neededbytes = (nbits+7)//8
    i = int.from_bytes(bt[:neededbytes], 'big')
    if nbits % 8:
        i >>= 8 - nbits % 8
    return i

M = "important message"
def sender(M):
    digest = hashlib.sha256(M.encode('utf-8')).digest()
    digest_60 = bitsof(digest, 60)
    DS = pow(digest_60,d,N)
    receiver(M, DS)

def receiver(M, DS):
    RecievedM = M
    DSprime = DS
    digest = hashlib.sha256(RecievedM.encode('utf-8')).digest()
    digest_60 = bitsof(digest, 60)
    print(digest_60)
    print(DSprime)
    print(e)
    print(N)
    DScheck = pow(DSprime,e,N)
    if digest_60 == DScheck:
        print("M' is signed by sender")
    else: print("signature DS' is invalid.")

digestM = hashlib.sha256(M.encode('utf-8')).digest()

sk=[[],[]]
pk=[[],[]]
ds = []
d = hashlib.sha256(M.encode('utf-8')).hexdigest()
binD = bin(int(d, 16))[2:]

def senderkey():
    for i in range(2):
        for j in range(256):
            sk[i].append(secrets.token_bytes(32))

def publickey():
    for i in range(2):
        for j in range(256):
            pk[i].append(hashlib.sha256(str(sk[i][j]).encode('utf-8')).hexdigest())

def dsGen():
        for i in range(len(binD)):
            ds.append(sk[int(binD[i])][i])

def receiver():
    for i in range(len(binD)):
        if hashlib.sha256(str(ds[i]).encode('utf-8')).hexdigest() != pk[int(binD[i])][i]:
            print("not verified")
            return

    print("DS is verified")

def server():
    k = str(randrange(pow(2,31), pow(2,32)))
    n = 10
    hash = hashlib.sha256(k.encode('utf-8')).hexdigest()
    hashes = [hash]

    for i in range(n):
        hash = hashlib.sha256(hash.encode('utf-8')).hexdigest()
        hashes.append(hash)
    
    client(list(reversed(hashes[:-1])), hashes[-1])
    pass

def client(hashes, serverkey):
    clientkey = hashes[0]
    if serverkey == hashlib.sha256(clientkey.encode('utf-8')).hexdigest():
        print('authenticated')
        serverkey = clientkey
    else:
        print('not authenticated')
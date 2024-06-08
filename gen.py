import numpy as np
import math
import bisect
import time
import sympy
import random
import os
from functools import lru_cache

class PrimeCalc:
    def __init__(self):
        self.primes = [2,3]
    def load_data(self,lim=-1):
        PATH = "primes.npy"
        if(not os.path.exists(PATH)): return
        arr = np.load(PATH)
        inx = np.searchsorted(arr,lim,side="right")
        self.primes = list(arr[:inx])
    
    def calc_next_prime(self):
        v = self.primes[-1]+2
        while(True):
            lim = int(math.sqrt(v+1))
            good = True
            for p in self.primes:
                if(p > lim): break
                if(v % p ==0): good = False; break
            if(good): break
            v += 2            
        self.primes.append(v)
        return v
    def get_prime_after(self,n):
        n += 1 + n%2
        while(not self.is_prime(n)): n += 2
        return n

    def is_prime(self,n:int) -> bool:
        if(n < 2): return False
        if(n <= self.primes[-1]):
            index = bisect.bisect_left(self.primes, n)
            return (0 <= index < len(self.primes)) and self.primes[index] == n
        lim = int(math.sqrt(n)+1)
        for p in self.primes:
            if(p > lim):  return True
            if(n % p == 0): return False
        while(self.calc_next_prime()<=lim):
            if(n % self.primes[-1] == 0): return False
        return True
    
    def get_prime_number(self,n) -> int: # return the n prime number
        if(n < len(self.primes)): return self.primes[n]
        l = n - len(self.primes)+1
        for i in range(l): self.calc_next_prime()
        return self.primes[-1]
        
    def load_primes_until(self,n) -> list: # return all the primes that are smaller then n
        if(self.primes[-1]>n):
            index = bisect.bisect_left(self.primes, n+1)
            return self.primes[:index]
        while(True):
            v = self.calc_next_prime()
            #print(f"\r{v} ({round(100*v/n,2)}%)   ",end="")
            if(v > n): return self.primes[:-1]
    def count_primes_until(self,n) -> int:
        if(n > self.primes[-1]): return len(self.load_primes_until(n))
        return bisect.bisect_left(self.primes,n+1)
    def get_prime_index(self,p):
        loc = bisect.bisect_left(self.primes,p)
        if(self.primes[loc] != p): return -1
        return loc

    def get_factors_primes(self,n): #return prime factors without duplication (exp for 4 return [2] for 12 return [2,3])
        return sympy.primefactors(n)
        if(n < 2): return []
        if(self.is_prime(n)): return [n]
        l,lim = [], int(math.sqrt(n))
        self.load_primes_until(lim)
        for p in self.primes:
            if(p > lim):  break
            if(n % p == 0): 
                l.append(p)
        r = []
        for f in l:
            d = n//f
            if(d != f and self.is_prime(d)): r.append(d)
        r.reverse()
        l += r
        
        return l
    
    def get_factors(self,n):
        if(True):
            f,l = sympy.factorint(n),[]
            f = dict(f)
            if(len(f)==1): return [list(f.keys())[0]**(list(f.values())[0]-1)]
            for i,v in f.items():
                l.append(i**v)
            return l

        f = self.get_factors_primes(n)
        for i,v in enumerate(f):
            m = v
            while(True):
                t = m*v
                if(n % t != 0 or n == t): f[i]=m; break
                m = t
        return f
    
    def phi(self,n):
        s,fac = n,self.get_factors_primes(n)
        for f in fac: s *= (1-1/f)
        return int(s)

class TimePredicter:
    def __init__(self) -> None:
        self.start = time.perf_counter()
    def predict_left(self,a,b) -> float: #predict the time left if done a from b
        return ((time.perf_counter()-self.start)/a)*(b-a)

class FactorsSaver:
    def __init__(self, primeCalc:PrimeCalc):
        self.factors = [{},{}]
        self.calc = primeCalc
    
    def add_to_dict(self,dic: dict, key, amount=1):
        res = dic.get(key)
        if (res is None): dic[key] = amount
        else: dic[key] += amount

    def get_factors(self,n) -> dict:
        if(n < len(self.factors)): return self.factors[n]

        factors = {}
        primes_devs = self.calc.load_primes_until(math.floor(math.sqrt(n)))

        for prime in primes_devs:
            if (n % prime != 0): continue
            self.add_to_dict(factors,prime)
        
        if(len(factors) == 0): factors = {n:1}
        else:
            facs_prod = np.prod(np.array(list(factors.keys())))
            dev_factors = self.get_factors(n // facs_prod)
            for dev_factor, amount in dev_factors.items():
                self.add_to_dict(factors,dev_factor, amount)

        if(n == len(self.factors)): self.factors.append(factors)
        return factors

    def get_prime_factors(self, n):
        return list(self.get_factors(n).keys())

    def get_full_factors(self, n):
        factors = self.get_factors(n)
        res = []
        for fac, amount in factors.items():
            res.append(fac**amount)
        return res


def get_digits(n) -> list:
    l = []
    while(n > 0):
        n,m = divmod(n,10)
        l.append(m)
    l.reverse()
    return l
def digts_to_num(l:list)-> int:
    s = 0
    for v in l:
        s *= 10
        s += v
    return s

def is_same_digits(a,b):
    return is_same_digits_l(get_digits(a),get_digits(b))
def is_same_digits_l(l1,l2):
    s = 0
    for v in l1: s ^= v
    for v in l2: s ^= v
    if(s != 0): return False
    l1.sort()
    l2.sort()
    return l1 == l2

def binary_search(arr,n) -> int: #return index of element if found else return -1
    index = bisect.bisect_left(arr, n)
    return -1 if (index==len(arr) or arr[index] != n) else index

def phi(n,PR:PrimeCalc=None):
    if(PR is None): PR = PrimeCalc()
    s = n
    PR.load_primes_until(n)
    for p in PR.primes:
        if(n % p ==0): s *= (1-1/p)
    return int(s)

def written_as_two(arr,n) -> list: # return in num can be written as two integers from array
    for v in arr:
        i = binary_search(arr,n-v)
        if(i != -1): return [v,arr[i]]
    return None
        
def get_permuation(arr:list,n) -> list:
    if(len(arr) == 0): return []
    s = math.factorial(len(arr)-1)
    d,m = divmod(n,s)
    v = arr.pop(d)
    l = get_permuation(arr,m)
    arr.insert(d,v)
    return [v] + l

def is_list_equal(l1,l2) -> bool:
    if(len(l1) != len(l2)): return False
    for i in range(len(l1)):
        if(l1[i] != l2[i]): return False
    return True
def is_palindrom(l,d=0):
    if(len(l)-(d<<1)<2): return True
    if(l[d] != l[-(d+1)]): return False
    return is_palindrom(l,d+1)

def is_2_nums_equal_target(arr, target,a=0): #if arr is sorted return if there pair at arr the equal to target
    for i in range(a,len(arr)):
        inx = bisect.bisect_left(arr,target-arr[i],i+1,len(arr))
        if(inx < len(arr) and arr[inx]+arr[i]==target): return True
    return False


def sort2(list,met):
    def swap(a,b): list[a],list[b] = list[b],list[a]
    def prap(l1,l2):
        val,loc = list[l2],l1
        for i in range(l1,l2):
            if(met(list[i],val)):
                swap(i,loc)
                loc += 1
        swap(loc,l2)
        return loc
    
    def inside_sort(a,b):
        if(b-a<2):
            if(b-a==1 and list[b]<list[a]): swap(a,b)
            return
        q = prap(a,b)
        inside_sort(a,q-1)
        inside_sort(q+1,b)
    inside_sort(0,len(list)-1)
    return list

def remove_duplicate_lists(arr,a=0): #from list of sorted lists remove duplicates  
    unique_lists = set(tuple(lst) for lst in arr)
    return [list(tpl) for tpl in unique_lists]

    res = []
    table = {}
    temp = []
    for v in arr:
        if(len(v)<=a): 
            res.append(v)
            continue
        elif(table.get(v[a]) is None):
            inx = len(table)
            table.__setitem__(v[a],inx)
            temp.append([])
        else: inx = table[v[a]]
        temp[inx].append(v) 
    
    for i,t in enumerate(temp):
        temp[i] = remove_duplicate_lists(t,a+1)
        res.extend(temp[i])

    return res

def cond_binary_search(a,b,test):
    if(b-a<=1): return a
    q = (a+b)//2
    if(test(q)): return cond_binary_search(q,b,test)
    return cond_binary_search(a,q,test)


class Frac:
    def __init__(self,a,b=1):
        self.a = a
        self.b = b
        self.minimaize()
    
    def minimaize(self):
        self.val = self.a/self.b
        c = math.gcd(abs(self.a),abs(self.b))
        if(self.a < 0 and self.b < 0): c *= -1
        self.a //= c
        self.b //= c
        
    def __add__(self,o):
        if(type(o) == Frac):
            self.a = self.a*o.b+o.a*self.b
            self.b *= o.b
        else:
            self.a += o
        self.minimaize()
        return self
    def __sub__(self,o):
        return self.__add__(o*-1)

    def __mul__(self, o):
        if(type(o) == Frac):
            self.a *= o.a
            self.b *= o.b
        else: self.a *= o
        self.minimaize()
        return self
    def __truediv__(self, o):
        if(type(o) == Frac):
            self.a *= o.b
            self.b *= o.a
        else: self.b *= o
        self.minimaize()
        return self

    def multiply(self,x):
        return x*self.val
    
    def __gt__(self, o) -> bool:
        if(type(o) == Frac): return self.val > o.val
        return self.val > o    
    def __ge__(self,o) -> bool:
        return self.__gt__(o) or self.__eq__(o)

    def __eq__(self, __value) -> bool:
        if(type(__value) == Frac): return  __value.a == self.a and __value.b == self.b
        return self.val == __value

    def __str__(self) -> str:
        return f"{self.a}/{self.b}"

class Matrix:
    def __init__(self, arr) -> None:
        self.arr = arr

    def get_clone(self):
        arr = []
        for l in self.arr:
            arr.append([])
            for r in l: arr[-1].append(r)
        return Matrix(arr)
    
    def __mod__(self, o):
        res = []
        for i in range(len(self.arr)):
            res.append([])
            for g in range(len(self.arr[i])):
                res[-1].append(self.arr[i][g] % o)
        return Matrix(res)
        

    def __add__(self, o):
        result = []
        for i in range(len(self.arr)):
            result.append([])
            for g in range(len(self.arr[i])):
                result[-1].append(self.arr[i][g] + o.arr[i][g])
        return Matrix(result)

    def __mul__(self, o):
        if not isinstance(o, Matrix):
            raise ValueError("Multiplication is only supported between two matrices.")
        if len(self.arr[0]) != len(o.arr):
            raise ValueError("Number of columns in the first matrix must match the number of rows in the second matrix.")

        result = []
        for i in range(len(self.arr)):
            row = []
            for j in range(len(o.arr[0])):
                sum = 0
                for k in range(len(o.arr)):
                    sum += self.arr[i][k] * o.arr[k][j]
                row.append(sum)
            result.append(row)
        return Matrix(result)

    def __pow__(self, num: int, mod=None):
        if (type(num) == int):
            bits = []
            while (num > 0):
                bits.append(num & 1)
                num >>= 1
            res = Matrix([[1 if g == i else 0 for g in range(len(self.arr))] for i in range(len(self.arr))])
            adder = self.get_clone()
            for i,b in enumerate(bits):
                if (b): res *= adder
                if (i != len(bits) - 1): adder *= adder
                if (mod is not None):
                    res %= mod
                    adder %= mod
            return res
        else:
            raise RuntimeError(f"can pow matrix with object from type {type(num)}")


    def __str__(self):
        result = ""
        for row in self.arr:
            result += " ".join(map(str, row)) + "\n"
        return result

if(__name__ == "__main__"):
    mat = Matrix([[1,21],[1,0]])
    mat = mat.__pow__(8,3)
    print(mat)
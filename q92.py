import gen
import numpy as np
import math
import time
import bisect
import sympy as sp
import sys
from matplotlib import pyplot as plt

sys.setrecursionlimit(100000)
from functools import cmp_to_key

PR = gen.PrimeCalc()

def read_data():
    with open("data.txt","r") as f:
        txt = f.read()
        f.close()
    return txt
def write_data(txt):
    with open("data.txt","w") as f:
        f.write(txt)
        f.close()


if(False): #q92

    
    def getV(n):
        s = 0
        for d in gen.get_digits(n): s += d*d
        return s
    def c89(n):
        if(mem[n] == 0): mem[n] = 1 if c89(getV(n)) else -1
        return mem[n]==1
        
    pred = gen.TimePredicter()
    c,S = 0,1
    L = 10**7
    mem = np.zeros(L+1)
    mem[1] = -1
    mem[89] = 1

    for i in range(1,L):
        good = c89(i)
        if(good): c += 1; mem[i] = 1
        else: mem[i] = -1
        S += 1
        print(f"\r{i}\\{L} ({round(100*i/L,2)}%) pred {pred.predict_left(i,L)}   ",end="")
    print(c)

if(False): #q100
    v1,v2 = 1,4
    L = 10**12
    while(True):
        v1,v2 = v2,6*v2-v1-2
        if(v2 >= L): break
    print(v2)

if(False): #q79
    exm = set(read_data().split("\n"))
    def does_pass(password:str):
        def sear(p:str,c:str):
            i = p.find(c)
            return len(p) if i==-1 else i

        for e in exm:
            loc = [sear(password,v) for v in e]
            if(loc[0] > loc[1] or loc[1] > loc[2]): return False
        return True
    def does_done(password:str):
        for e in exm:
            for v in e:
                if(password.find(v)==-1): return False
        return True
    
    opts = [""]
    while(True):
        vals = []
        for opt in opts:
            for i in range(10):
                r = opt + str(i)
                if(does_pass(r)): vals.append(r)
        opts = vals
        for o in opts:
            if(does_done(o)):
                print(o)
                break

if(False): #q700
    a = 1504170715041707
    b = 4503599627370517
    n = 0
    vals = []
    inx = []
    

    if(0):
        for i in range(1,1000000):
            n = (a*i) % b
            if(len(vals)==0 or n<vals[-1]): vals.append(n); inx.append(i)
    else:
        vals = [b]
        inx = [0]
        i = 1
        while(vals[-1] != 1):
            n = (a*i)%b
            vals.append(n)
            inx.append(i)
            i *= vals[-2]//vals[-1]
            i += inx[-1]-inx[-2]
    
    for i,v in enumerate(vals):
        print(f"{inx[i]} -> {v} [{((b if i==0 else vals[i-1])//v)}] <{b//v}> [{(inx[i]+i-1)%3}]")
    for i in range(1,len(vals)):
        print(f"{inx[i]} - {inx[i-1]} = {inx[i]-inx[i-1]}")
        
    print(sum(vals)-vals[0])
    

if(False): #q816
    s,L,M = 290797,2000000,50515093
    T=M**2
    GRID_PARTS = 400
    points = np.zeros((L,2),dtype=np.int64)
    for i in range(L):
        points[i][0] = s
        s *= s
        s %= M
        points[i][1] = s
        s *= s
        s %= M
    
    #loc = np.zeros((L,2))
    grids = []
    for y in range(GRID_PARTS):
        grids.append([])
        for x in range(GRID_PARTS): grids[y].append([])

    for i,p in enumerate(points):
        grids[(p[0]*GRID_PARTS)//M][((p[1]*GRID_PARTS))//M].append(p)
        #loc[i][0] = (p[0]*GRID_PARTS)//M
        #loc[i][1] = (p[1]*GRID_PARTS)//M
    
    LOWEST = len(grids[0][0])
    for g1 in grids:
        for g2 in g1:
            LOWEST = min(LOWEST,len(g2))
    print(f"smallest grid is {LOWEST}")
    
    def get_min_dis2(p1,p2):
        dis = T
        for v1 in p1:
            for v2 in p2:
                d = (v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])
                if(dis == -1 or d < dis): dis = d
        return dis
    def get_min_dis(p1):
        dis = T
        for i,v1 in enumerate(p1):
            for g in range(i+1,len(p1)):
                v2 = p1[g]
                d = (v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])
                if(dis == -1 or d < dis): dis = d
        return dis


    t1 = time.perf_counter()
    dis = T
    CL = gen.TimePredicter()
    for y in range(GRID_PARTS):
        for x in range(GRID_PARTS):
            #print(f"[{len(grids[y][x])}]",end="")
            dis = min(dis,get_min_dis(grids[y][x]))
            if(y != GRID_PARTS-1): dis = min(dis,get_min_dis2(grids[y][x],grids[y+1][x]))
            if(x == GRID_PARTS-1): continue
            dis = min(dis,get_min_dis2(grids[y][x],grids[y][x+1]))
            if(y != GRID_PARTS-1): dis = min(dis,get_min_dis2(grids[y][x],grids[y+1][x+1]))
            if(y != 0): dis = min(dis,get_min_dis2(grids[y][x],grids[y-1][x+1]))
            print(f"\r{y*GRID_PARTS+x}\\{GRID_PARTS**2} time left {CL.predict_left(y*GRID_PARTS+x+1,GRID_PARTS**2)}",end="")
        #print()

    print(f"\ntook {time.perf_counter()-t1}")
    print(f"\n\nres: {math.sqrt(dis)}")

if(False): #q?
    primes = []
    for l in range(1,8):
        v = 10**l
        for i in range(v):
            num = v+i
            if(PR.is_prime(num)): primes.append(num)

    res = []
    for p in primes:
        s = p*p
        digits = gen.get_digits(s)
        digits.reverse()
        r = gen.digts_to_num(digits)
        if(s == r): continue
        r_sqrt = math.sqrt(r)
        if(r_sqrt.is_integer() and PR.is_prime(int(r_sqrt)) and (not s in res)):
            res.append(s)
            res.append(r)
    res.sort()
    print(res)
    print(len(res))

if(False): #q71
    data = read_data()
    grid = []
    for line in data.split("\n"):
        grid.append([int(v) for v in line.split(",")])
    points_sums = []
    for i in range(len(grid)*2-1): points_sums.append([])
    print(points_sums)
    for y in range(len(grid)):
        for x in range(len(grid)):
            points_sums[y+x].append(grid[y][x])
    for g in range(len(points_sums)-2,-1,-1):
        points = points_sums[g]
        prev = points_sums[g+1]
        for i,p in enumerate(points):
            bigger = len(points) > len(prev)
            if(not bigger):
                points[i] += min(prev[i],prev[i+1])
                continue
            if(i==0): points[i] += prev[0]
            elif(i==len(prev)): points[i] += prev[i-1]
            else: points[i] += min(prev[i],prev[i-1])
    print(points_sums[0][0])


if(False): #q102
    class Line:
        def __init__(self, p1,p2):
            x1,y1 = p1
            x2,y2 = p2
            self.var = x1 == x2
            if self.var:
                self.B = x1
                self.X1 = min(y1, y2)
                self.X2 = max(y1, y2)
                return
            self.M = (y2 - y1) / (x2 - x1)
            self.B = y1 - self.M * x1
            self.X1 = min(x1, x2)
            self.X2 = max(x1, x2)
        
        def is_in_range(self,x):
            return (self.B==x) if self.var else (self.X1 < x < self.X2)

        def get_y(self, x):
            return self.M * x + self.B

    def originIn(p1,p2,p3):
        lines = [Line(p1,p2),Line(p2,p3),Line(p1,p3)]
        down,up = 0,0
        for line in lines:
            if(not line.is_in_range(0)): continue
            if(line.get_y(0) > 0): up += 1
            else: down += 1
        return down == 1 and up == 1

    data = read_data()
    c = 0
    for line in data.split("\n"):
        p = [int(v) for v in line.split(",")]
        orin = originIn((p[0],p[1]),(p[2],p[3]),(p[4],p[5]))
        if(orin): c += 1
    print(c)

if(False): #q104
    L=1000000
    M1=10**9
    a,b = 1,1
    dn = [1,2,3,4,5,6,7,8,9]
    for i in range(L):
        a,b = b,a+b
        d1 = gen.get_digits(b%M1)
        print(f"\r{i}\\{L}   ",end="")
        if(not gen.is_same_digits_l(d1,dn)): continue
        if(not gen.is_same_digits_l(gen.get_digits(b)[:9],dn)): continue
        print(f"\n\nat {i}")
        break

if(False): #q91
    V = 51
    L = V**2
    c = 0
    for p1i in range(1,L):
        y1,x1 = divmod(p1i,V)
        d1 = y1*y1+x1*x1
        for p2i in range(p1i+1,L):
            y2,x2 = divmod(p2i,V)
            d2 = y2*y2+x2*x2
            d3 = (y1-y2)*(y1-y2)+(x1-x2)*(x1-x2)
            good = False
            if(d1 > d2):
                if(d1 > d3): good = d2+d3==d1
                else: good = d1+d2==d3
            else:
                if(d2 > d3): good = d1+d3==d2
                else: good = d1+d2==d3
            if(good): c += 1
        print(f"\r{p1i}\\{L}  ",end="")

    print(f"\n\n{c}")

if(False): #q82
    data = read_data()
    grid = []
    for line in data.split("\n"):
        grid.append([int(v) for v in line.split(",")])
    T = 10**9

    for x in range(len(grid)-2,-1,-1):
        for y in range(0,len(grid)):
            val = T
            if(y != 0): val = grid[y-1][x]+grid[y][x]
            step = 0
            for y2 in range(y,len(grid)):
                step += grid[y2][x]
                if(step > val): break
                val = min(step+grid[y2][x+1],val)
            grid[y][x] = val
    #
    # for y in grid:  print(y)
    
    print(min(np.array(grid)[:,0]))

if(False): #q87
    p1,p2,p3 = 2,2,2
    v1,v2,v3 = 4,8,16
    L,s = 50*(10**6),0
    nums = []
    while(v3 < L):
        p2,v2 = 2,8
        while(v3 + v2 < L):
            p1,v1 = 2,4
            s1 = v3+v2
            while(s1 + v1 < L):
                #print(f"{p1} + {p2} + {p3} [{v1} + {v2} + {v3}] = {v1+v2+v3}")
                nums.append(s1+v1)
                p1 = PR.get_prime_after(p1)
                v1 = p1*p1
            p2 = PR.get_prime_after(p2)
            v2 = p2*p2*p2
        p3 = PR.get_prime_after(p3)
        v3 = p3**4
        #print(f"\r v1 is {v1} with {p1} ",end="   ")
    print("\n\n" + str(s))
    print(len(set(nums)))

if(False): #q105
    def isRule1(arr,si=None):
        if(len(arr)==1): return True
        if(not si is None):
            v1 = arr[si]
            for i,v2 in enumerate(arr):
                if(v1==v2 and i!=si): 
                    print(f"{v1} = {v2}")
                    return False
        tes = np.zeros(len(arr)-1)
        for i1 in range(len(arr)):
            for i2 in range(i1+1,len(arr)):
                i3 = 0
                for i,v in enumerate(arr):
                    if(i != i1 and i != i2):
                        tes[i3] = v
                        i3 += 1
                tes[i3] = arr[i1]+arr[i2]
                if(not isRule1(tes,i3)): return False
        return True
    def can_add(arr,target): #if there substring of arr the equal target
        if(len(arr)==1): return arr[0]==target
        arr = arr[:bisect.bisect_right(arr,target)]
        sub = np.array(arr[1:])
        for i,v in enumerate(arr):
            if(v==target): return True
            if(can_add(sub,target-v)): return True
            if(i != len(sub)): sub[i]=v
        return False
    
    def isRule1B(arr,min=0):
        sm = sum(arr)
        if(sm % 2 == 0):
            arr2 = np.multiply(arr,2)
            if(can_add(arr2,sm)): 
                #print(arr)
                return False
        for i in range(len(arr)):
            if(arr[i] < min): continue
            arr2 = arr[:i] + arr[i+1:]
            if(not isRule1B(arr2,arr[i])): return False
        return True
    def isRule2(arr):
        s1,s2 = arr[0],0
        for i in range(1,(len(arr)+1)//2):
            s1 += arr[i]
            s2 += arr[-i]
            if(s2 >= s1): return False
        return True

    def isGood(arr):
        return isRule2(arr) and isRule1B(arr)

    if(False): #q105
        data = read_data()
        sets = []
        for line in data.split("\n"): sets.append([int(v) for v in line.split(",")])
        sets = np.array(sets,dtype=object)
        tt = 0
        for i,set in enumerate(sets):
            set.sort()
            print(f"\rset {i}\\{len(sets)} is {set} len {len(set)}",end="                                       ")
            if(not isRule2(set)): continue
            if(not isRule1B(set)): continue
            #print(f"{set}")
            tt += sum(set)
        print(f"\n\n{tt}")

    if(False): #q103
        S,LIM = 7,70
        
        base = []
        def get_arrs(l=S,min=1,ms=235):
            base = []
            for i in range(min,LIM): base.append([i])
            #print(f"base is {base}")
            if(l == 1): return base
            res = []
            for b in base:
                L = ms-b[0]
                for v in get_arrs(l-1,b[0]+1):
                    if(sum(v)>L): continue
                    res.append(b+v)
            return res
        
        SV = 19
        arrs = get_arrs(6,SV+1)
        print(f"build {len(arrs)} arrs")
        tps = []
        for arr in arrs:
            arr = [SV] + arr
            tps.append((arr,sum(arr)))
        tps.sort(key=lambda x:x[1])
        arrs = [v[0] for v in tps]
        #arrs = arrs[int(len(arrs)*0.5):]
        for i,arr in enumerate(arrs):
            print(f"\r{i}\\{len(arrs)} ({round(100*i/len(arrs),2)}%)",end="     ")
            if(not isGood(arr)): continue
            print(f"best {S} is {arr} ")
            break

if(False): #q95
        print(sp.divisors(28))
        L = 10**6
        chain_ids = np.zeros(L,dtype=np.int)
        if(False):
            sums = np.zeros(L)
            for i in range(1,L):
                #devs = PR.get_factors(i)
                #for d in devs: sums[i-1] += sums[d-1] + d
                sums[i-1] = sum(sp.divisors(i))-i
                print(f"\rfactor {i}\\{L}   ",end="")
            np.save("arr.npy",sums)
        else: sums = np.array(np.load("arr.npy"),dtype=np.uint)
        
        print(sums)
        L = 1000000
        sums = np.insert(sums,0,0)
        sums = sums[:L]

        id = 0
        for i in range(L):
            print(f"\r{i}\\{L}",end="     ")
            if(chain_ids[i] != 0): continue
            id += 1
            chain_ids[i] = id
            inxs = [i,sums[i]]
            while(inxs[-1]<L and chain_ids[inxs[-1]] == 0): 
                chain_ids[inxs[-1]] = id
                inxs.append(sums[inxs[-1]])
            mark = -1 if inxs[-1]>=L else min(chain_ids[inxs[-1]],id)
            for inx in inxs: 
                if(inx >= L): continue
                chain_ids[inx] = mark
        
        chins = []
        print(f"max is {max(chain_ids)+1}")
        for i in range(max(chain_ids)+1): chins.append([])
        for i,s in enumerate(sums):
            if(chain_ids[i] == -1): continue
            chins[chain_ids[i]].append(i)
        chins2 = []
        for c in chins:
            if(len(c)==0 or c[0]==0): continue
            c:list
            j1,j2 = sums[c[0]],sums[sums[c[0]]]
            while(j1 != j2):
                j1,j2 = sums[j1],sums[sums[j2]]
            
            arr = []
            while(True):
                arr.append(j1)
                j1 = sums[j1]
                if(j1 == j2): break
            chins2.append(arr); continue
    
        chins2.sort(key=lambda x:len(x))
        for c in chins2: print(c)
        print(f"\nres is {min(chins2[-1])}")

        
if(False): #54
    def get_rank(cards):
        kinds = np.zeros(4)
        for card in cards: kinds[card[1]] += 1
        row = True
        for i in range(1,len(cards)):
            if(cards[i-1][0]+1!=cards[i][0]): row = False; break
        if(row and max(kinds)==5): return (0,cards[-1][0]) if (cards[-1][0]==14) else (1,cards[-1][0])#return royal flash or stait falsh
        comps = [cards[i][0]==cards[i+1][0] for i in range(len(cards)-1)]
        comps_count = sum(comps)
        if(comps_count==3 and (cards[0][0]==cards[-2][0] or cards[1][0]==cards[-1][0])): return (2,cards[-2][0] if (cards[0][0]==cards[-2][0]) else cards[-1][0]) #four of a kind
        if((cards[0][0]==cards[2][0] and cards[3][0] == cards[4][0]) or (cards[0][0]==cards[1][0] and cards[2][0]==cards[4][0])): return 3,(cards[-1][0] if cards[-1][0]==cards[-3][0] else cards[-3][0]) #full house
        if(max(kinds)>=4):
            last_odd = cards[-1][1] != cards[-2][1] and cards[-1][1] != cards[-3][1]
            return 4,(cards[-2][0] if last_odd else cards[-1][0]) #flush
        if(row): return 5,cards[0][-1] #row
        if(cards[0][0]==cards[2][0] or cards[2][0]==cards[4][0]): return 6,(cards[2][0] if cards[0][0]==cards[2][0] else cards[4][0]) #tree of a kind
        if(comps_count==2): return 7,(cards[4][0] if comps[-1] else cards[3][0]) #two pairts
        if(comps_count==1): 
            inx = comps.index(True)
            return 8,cards[inx+1][0] #pair
        return 9,cards[-1][0]

    data = read_data()
    wins = 0
    kind_dict = {"S":0,"C":1,"D":2,"H":3}
    num_dict = {"T":10,"J":11,"Q":12,"K":13,"A":14}
    for i in range(1,10): num_dict[str(i)] = i
    for line in data.split("\n"):
        marks = line.split(" ")
        cards = [(num_dict[m[0]],kind_dict[m[1]]) for m in marks]
        cards1 = cards[:5]
        cards2 = cards[5:]
        cards1.sort(key=lambda x:x[0])
        cards2.sort(key=lambda x:x[0])
        r1,l1 = get_rank(cards1)
        r2,l2 = get_rank(cards2)
        if(r1 == r2 and l1 != l2):
            if(l1 > l2): r1 -= 1
            else: r2 -= 1
        elif(r1 == r2):
            w1 = True
            for i in range(len(cards1)-1,-1,-1):
                if(cards1[i][0]==cards2[i][0]): continue
                w1 = cards1[i][0]>cards2[i][0]
                break
            if(w1): r1 -= 1
            else: r2 -= 1
        print(f"{cards1 if r1<r2 else cards2} win {cards2 if r1<r2 else cards1} with {min(r1,r2)} comare to {max(r1,r2)} won [{1 if r1<r2 else 2}]")
        if(r1<r2): wins += 1
    print(wins)
        
if(False): #85
    def get_rect(w,h):
        return ((w*w+w)//2)*((h*h+h)//2)
    G = 2*(10**6)
    best = 1
    dis = G
    for h in range(1,100):
        for w in range(h,100):
            val = get_rect(h,w)
            if(abs(G-val) < dis): dis,best = abs(G-val),w*h
    print(f"area {best} with dis {dis}")

if(False): #112
    def isB(n):
        down,up = False,False
        n,m = divmod(n,10)
        while(n > 0):
            n,c = divmod(n,10)
            if(c == m): continue
            if(c > m): up = True
            else: down = True
            if(up and down): return True
            m = c
        return False

    print(isB(66420))

    bon = 0
    for i in range(1,10000000):
        if(isB(i)): bon += 1
        v = i-bon
        if(bon == 0 or v == 0): continue
        p = bon/(v*(1+bon/v))
        if(p==0.99):
            print(f"at {i} with {bon}")
            break
    print(f"ratio is {p} with {bon}")

if(False): #q74
    def get_digits_fac(n):  return sum([math.factorial(d) for d in gen.get_digits(n)])
    if(0):
        nums = np.zeros(2177281,dtype=np.uint)
        for i in range(len(nums)): nums[i] = get_digits_fac(i)
        np.save("arr.npy",nums)
    else: nums = np.array(np.load("arr.npy"),dtype=np.uint)
    

    hits = np.zeros(len(nums))
    L = min(10**6,len(nums))

    mx = 0
    res = 0
    for i in range(1,L):
        ln,num = 0,i
        while(True):
            if(hits[num]==i): break
            hits[num] = i
            num = nums[num]
            ln += 1
        mx = max(mx,ln)
        if(ln == 60): res += 1
        print(f"\r{i}\\{L} with biggest {mx}",end="     ")
    print(f"\n\nmx is {mx}")
    print(f"res is {res}")
        
if(False): #q75-1
    L = 15*(10**3)+1
    kinds = np.zeros(L,dtype=np.uint)
    a,v1=1,1
    while(a<L):
        b,v2 = a+1,(a+1)*(a+1)
        while(a + b < L):
            sm = v1 + v2
            b += 1
            v2 = b*b
            c = math.sqrt(sm)
            if(not c.is_integer()): continue
            inx = int(a+b+c-1)
            if(inx > len(kinds)): break
            kinds[inx] += 1
        a += 1
        v1 = a*a
    for i,v in enumerate(kinds[:200]): print(f"[{i}] have {v}")
    print(sum([v==1 for v in kinds]))
if(False): #q75-2
    L = 15*(10**5)+1
    trip = []
    for s in range(1,L):
        ss = s*s
        if(ss >= L): break
        for t in range(1,s):
            if((s % 2 == 1 and  t % 2 == 1) or (math.gcd(s,t) != 1)): continue
            tt = t*t
            b = ss-tt
            if(b < 0): break
            a,c = 2*s*t,ss+tt
            if(a+b+c > L): continue
            if(a>b): trip.append((b,a,c))
            else: trip.append((a,b,c))

    TL = len(trip)
    SM = trip[:TL]
    for i in range(TL):
        t,m = trip[i],2
        while(True):
            v = (t[0]*m,t[1]*m,t[2]*m)
            if(sum(v)<=L): trip.append(v)
            elif(sum(v)>L): break
            m += 1
    #print(trip)
    print(len(trip))
        
    con = np.zeros(L)
    for t in trip: con[sum(t)] += 1
    for i,c in enumerate(con[:200]): print(f"[{i}] is [{c}]")
    print(sum([c==1 for c in con]))
    

if(False): #q88

    class TList(list):
        def __hash__(self) -> int:
            return sum(self).__hash__()

    def get_splits(arr,a=0):
        prev = -1
        arr.sort()
        res = [([1]+arr)]
        for i in range(a,len(arr)):
            p = arr[i]
            if(prev==p): continue
            prev = p
            if(p==1 or PR.is_prime(p)): continue
            temp = arr[:i] + arr[i+1:]
            splits = get_mult_res(p,2)
            for sp in splits: res.append(temp+sp)
        return res
   
    def get_mult_res(num,size=-1,sm=-1,se=True):
        lm = int(math.sqrt(num))+1
        prods = []
        for i in range(2,lm):
            if(num % i == 0): prods.append([i,num//i])

        full = []
        itr = 2
        while(True):
            if(len(prods) == 0): break
            if(itr==size): break
            itr += 1
            new = set() if se else []
            for prod in prods:
                for spl in get_splits(prod):
                    if(sm!=-1 and sum(spl)>sm): continue
                    if(spl[0]==1): 
                        if(len(full)!=0 and full[-1]==spl): continue
                        full.append(spl)
                    else:
                        #if(len(new)!=0): continue
                        if(se): new.add(TList(spl))
                        else: new.append(spl)
            prods = new


        for i,f in enumerate(full):
            amount = size-len(f)
            f:list
            for g in range(amount): f.insert(0,1)
            if(sm == -1 or sum(f)<=sm):
                if(se):  prods.add(TList(f))
                else: prods.append(f)
        
        if(size==-1):
            for p in prods: p.pop(0)

        return prods

    for m in get_mult_res(36): print(m)
    print("\n---\n")
    for m in get_mult_res(36,se=False): print(m)
    
    L = 1000000
    TEM = np.zeros(12001)
    adds = 2
    for i in range(L):
        print(f"\rrun {i}\\{L} fill {adds}\\{len(TEM)} ({round(100*adds/len(TEM),2)}%)",end="    ")
        spl = list(get_mult_res(i,sm=i+1))
        sums = [sum(s) for s in spl]
        for g,s in enumerate(sums):
            k = i+len(spl[g])-s
            if(k >= len(TEM)): continue
            if(TEM[k]==0): TEM[k] = i; adds += 1
        if(adds == len(TEM)): break

    for i,v in enumerate(TEM):
        print(f"k[{i}] = [{v}]")

    print(f"res is {sum(set(TEM))}")
    exit()

    res = []
    for k in range(2,12000):
        print(f"\rat {k}",end="    ")
        found = False
        for i in range(4,10000):
            opts = get_mult_res(i,k,i)
            opts = [o for o in opts if sum(o)==i]
            if(len(opts) == 0): continue
            res.append(i)
            found = True
            #print(f"for k={k} have [{i}] sum {sum(opts[0])}")
            break
        if(not found): print(f"error didn't found for {k}")
    print(f"\nfinal result is {sum(set(res))}")

    L = 12
    minN = np.zeros(L)
    for i in range(2,L): pass
        
if(False): #q83
    data = read_data()
    C = 10**6
    finl = []
    grid = []
    for line in data.split("\n"):
        grid.append([int(v) for v in line.split(",")])
        finl.append([int(v) for v in line.split(",")])

    def get_at(y,x):
        if(x < 0 or y < 0 or x >= len(grid) or y >= len(grid)): return C
        return grid[y][x]

    def get_min(x1,y1,x2,y2):
        xmin,xmax,ymin,ymax = min(x1,x2),max(x1,x2),min(y1,y2),max(y1,y2)
        
    
    def calc(first):
        for s in range(2*len(grid)-3,-1,-1):
            x = min(len(grid)-1,s)
            while(x >= 0):
                y = s-x
                if(y >= len(grid)): break
                #if(x==len(grid)-1): val = grid[y+1][x]
                #elif(y==len(grid)-1): val = grid[y][x+1]
                #else: val = min(grid[y+1][x],grid[y][x+1])
                if(first): vals = [get_at(y+1,x),get_at(y,x+1)]
                if(not first):
                    vals = []
                    for i in range(1,2):
                        px,py = i,0
                        while(px >= 0):
                            if(py != 0): vals.extend([get_at(y+py,x),get_at(y-py,x)])
                            if(px != 0): vals.extend([get_at(y,x+px),get_at(y,x-px)])
                            px -= 1; py += 1

                val = min(vals)
                grid[y][x] = finl[y][x] + val
                x -= 1
    calc(True)
    prev = grid[0][0]
    while(True):
        print(prev)
        calc(False)
        if(grid[0][0]==prev): break
        prev = grid[0][0]
    print(grid[0][0])

if(False): #q77
    mem = [[],[]]
    lens = [0,0]

    def compar(item1, item2):
        for a,b in zip(item1,item2):
            if(a==b): continue
            return a < b
        return len(item1)<len(item2)


    def ways_to_write(n):
        if(n < len(mem)): return mem[n]
        res = []
        if(PR.is_prime(n)): res.append([n])
        for i in range(1,n//2+1):
            v1,v2 = ways_to_write(i), ways_to_write(n-i)
            for nv1 in v1:
                for nv2 in v2:
                    res.append(nv1+nv2)
        
        if(len(res)==0): return []
        
        for val in res: val.sort()
        
        if(0):
            gen.sort2(res,compar)
            res2 = [res[0]]
            for i in range(1,len(res)):
                if(res2[-1]==res[i]): continue
                res2.append(res[i])
        else: res2 = gen.remove_duplicate_lists(res)

        return res2

    sm = 0
    for i in range(2,1000):
        w = ways_to_write(i)
        sm += len(w)
        if(len(w) > 5000):
            print(f"\n\nfound {i} with {len(w)}")
            break
        mem.append(w)
        lens.append(len(w))
        print(f"\rat [{i}] max ways is {max(lens)}  ",end="  ")
    
    print(f"\nsm is {sm}")


if(False): #q101
    def get_len(pol):
        for i in range(len(pol)-1,-1,-1): 
            if(pol[i] != 0): return i+1
        return 0

    def multiply_pols(p1,p2):
        l1,l2 = get_len(p1), get_len(p2)
        res = np.zeros(l1+l2-1)
        for i1 in range(l1):
            v1 = p1[i1]
            for i2 in range(l2):
                res[i1+i2] += v1*p2[i2]
        return res
    def calculate(pol,x):
        s = 0
        for i,p in enumerate(pol): s += p*(x**i)
        return s

    def generate_polynom(points):
        pol = np.zeros(len(points))
        for i,p in enumerate(points):
            sub_pol = np.zeros(len(points))
            sub_pol[0] = 1
            for g,p2 in enumerate(points):
                if(i==g): continue
                sub_pol = multiply_pols(sub_pol,[-g,1])
            sub_pol /=  calculate(sub_pol,i)
            sub_pol *= p
            pol += sub_pol
        return pol

    UP = []
    for i in range(11): UP.append(1 if i % 2 == 0 else -1)
    #UP = [0,0,0,1]
    UP = np.array(UP)
    VALS = [calculate(UP,i) for i in range(1,len(UP)+1)]
    print(VALS)
    S = 0
    for i in range(1,len(UP)):
        pol = generate_polynom(VALS[:i])
        print(f"pol degree [{i}] generate [{calculate(pol,i)}] with vals {VALS[:i]}")
        S += calculate(pol,i)
    print(S)

if(False): #q108
    def get_ways(n):
        c,l = 0,n*2
        for i in range(n+1,l+1):
            if(i*n % (i-n) == 0): c += 1
        return c
    def get_ways2(n):
        pass
        
    for i in range(1,10):
        v1,v2 = get_ways(i),get_ways2(i)
        if(v1 != v2):
            print(f"error at {i} with {v1} and {v2}")
            break
    exit()

    L = 10000000000000
    mx = 1
    for i in range(0,L,4324320):
        print(f"\r{i}\\{L} with max [{mx}]",end="     ")
        w = get_ways(i)
        if(w > mx): print(f"good at {i} with {w}")
        mx = max(w,mx)
        if(w > 1000000):
            print(f"\n\nfound at {i} with {w}")
            break


if(False): #q107
    data = read_data()
    linesN = data.split("\n")
    lines = []
    nodes = []
    for i in range(len(linesN)): nodes.append([])
    print(len(lines))
    for i,line in enumerate(linesN):
        vals = line.split(",")
        for g,v in enumerate(vals):
            if(v=='-' or g < i): continue
            lines.append([i,g,int(v),0,0])
            nodes[i].append(lines[-1])
            nodes[g].append(lines[-1])
            #nodes[-1].append((min(i,g),max(i,g),int(v),0,0))
    
    print(nodes)

    def get_sum(all=False):
        return sum([line[2] for line in lines if (all or line[-1])])

    START_SUM = get_sum(True)

    #mark only the shortes line
    for node in nodes:
        mi,mv = 0,node[0][2]
        for i,v in enumerate(node):
            if(v[2] < mv): mi,mv = i,v[2]
        node[mi][-1] = 1
    
    MARKS = np.zeros(len(nodes))
    global MARK
    MARK = 0

    def get_nab(i,mark):
        if(MARKS[i] == mark): return []
        MARKS[i] = mark
        nabs = [i]
        for n in nodes[i]:
            if(not n[-1]): continue
            nabs.extend(get_nab(n[0] if n[1]==i else n[1],mark))
        return nabs

    def get_blocks(): 
        global MARK
        blocks = []
        MARK += 1
        for i in range(len(nodes)):
            nab = get_nab(i,MARK)
            if(len(nab) == 0): continue
            blocks.append(nab)
        return blocks
    
    def connect_blocks(blocks):
        targ = 0
        for i,b in enumerate(blocks):
            if(len(b) < len(blocks[targ])): targ = i
        
        mt,mti,mtl = 0,0,10000
        for t in blocks[targ]:
            for i,n in enumerate(nodes[t]):
                if(n[-1]): continue
                if(n[2] < mtl): mt,mti,mtl = t,i,n[2]
        nodes[mt][mti][-1] = 1

        

    bl = get_blocks()
    while(len(bl) > 1):
        print(f"len is {len(bl)}")
        connect_blocks(bl)
        bl = get_blocks()

    print(f"blocks are {bl} with sum of [{get_sum()}] diffrent [{START_SUM-get_sum()}]")

if(False): #q114 & 115
    M = 3
    TEMP = np.zeros(200)
    def get_options(s):
        if(s < M): return 1
        if(s < len(TEMP) and TEMP[s] != 0): return TEMP[s]
        su = get_options(s-1)
        for i in range(M,s+1):
            su += get_options(s-i-1)
        if(s < len(TEMP)): TEMP[s] = su
        return su

    t1 = time.perf_counter()
    get_options(50)
    t2 = time.perf_counter()
    print(f"took {t2-t1}")

    for i in range(200):
        r = get_options(i)
        if(r > 1000000):
            print(f"found at {i}")
            break

if(False): #q116 & 117
    M = [2,3,4]
    TEMP = np.zeros(200)
    def get_options(s):
        if(s < M[0]): return 1
        if(s < len(TEMP) and TEMP[s] != 0): return TEMP[s]
        
        su = get_options(s-1)
        for m in M:
            if(s < m): break
            su += get_options(s-m)

        if(s < len(TEMP)): TEMP[s] = su
        return su
    
    print(get_options(50))

if(False): #q119
    for i in range(10,0):
        sm = sum(gen.get_digits(i))
        if(sm==1): continue
        if(math.log(i,sm).is_integer()):
            print(f"found {i} with {sm}^{int(math.log(i,sm))}")

    nums = []
    for a in range(2,100):
        for b in range(2,100):
            n = a**b
            if(sum(gen.get_digits(n))==a):
                nums.append(n)
                print(f"found {n} to be {a}^{b}")
    
    print(f"found {nums}\nlength is [{len(nums)}]")
    nums.sort()
    print(nums[1])
    print(nums[29])
            

if(False): #q120
    def r_max(a):
        M1,M2,AS = a-1,a+1,a*a
        V1,V2 = 1,1
        MARK,MX = -1,2
        for i in range(0,10000):
            R = (V1+V2)%AS
            if(MARK == -1 and R != 2): MARK = R
            elif(R == MARK): return MX
            MX = max(R,MX)
            #print(f"for n={i} mod {R}")
            V1 *= M1; V2 *= M2
            
    
    sm = 0
    for i in range(3,1001):
        sm += r_max(i)
    print(f"sm is {sm}")

if(False): #q123
    mx = 0
    L = 10**10
    for i in range(2000,100000):
        print(f"\rat {i} with mx [{mx}]",end="    ")
        p = PR.get_prime_number(i-1)
        M = ((p-1)**i + (p+1)**i) % (p*p)
        mx = max(M,mx)
        if(M > L):
            print(f"found at {i} with p {p}")
            break

if(False): #124
    rads = []
    def rad(n):
        devs = PR.get_factors_primes(n)
        return np.prod(devs)
    for i in range(1,100001):
        print(f"\rbuild {i}",end="   ")
        rads.append((i,rad(i)))
    rads.sort(key=lambda x:x[1])
    print(rads[9999][0])

if(False): #q125
    print(np.sum(np.arange(1,6,dtype=np.uint64)*np.arange(1,6,dtype=np.uint64)))
    LIM = 10**8
    def get_squares_sums():
        nums = []
        a,b,s=1,1,0
        while(True):
            AS = a*a
            if(AS >= LIM): break
            b,s = a,AS
            while(s < LIM):
                #print(f"{np.arange(a,b+1)} = {s}")
                if(s != AS): nums.append(s)
                b += 1
                s += b*b
            a += 1
        return nums



    nums = get_squares_sums()
    print(len(nums))
    nums = [n for n in nums if gen.is_palindrom(gen.get_digits(n))]
    print(nums)
    print(len(nums))
    print(sum(set(nums)))

if(False): #q122
    def fix_best(n):
        s = int(math.log2(n))-1
        while(n > 0):
            if(n&1): s += 1
            n>>=1
        return s
    def brute_force(nums,target,size):
        if(size == len(nums)): return size
        #if(size == len(nums)+1): return len(nums) if gen.is_2_nums_equal_target(nums,target,a=len(nums)-1) else size
        best = size
        for i,n1 in enumerate(nums):
            for g in range(i,len(nums)):
                if(g != len(nums)-1): continue
                n2 = nums[g]
                if((n1+n2 <= nums[-1]) or (n1+n2>target) or ((n1+n2) in nums)): continue
                if(n1+n2==target): return len(nums)
                best = min(best,brute_force(nums+[n1+n2],target,size))
        return best
    
    print(brute_force([1],15,7))
    sm = 0
    for i in range(1,201):
        print(f"\r{i}\\{200} next fix is [{fix_best(i)}] ",end="")
        fix = fix_best(i)
        sm += brute_force([1],i,fix)
    print(f"\n\nsm is [{sm}]")

if(False): #q98
    EL = 10000
    arr = np.zeros(EL)
    for i in range(10000000):
        arr[(i*i)%EL] += 1
    ender = [i for i in range(len(arr)) if arr[i]>0]
    #print(arr)
    #print(ender)
    #print(len(ender))

    def word_to_code(word):
        arr = np.zeros(26)
        for c in word: arr[ord(c)-65] = 1
        n,m = 0,1
        for v in arr:
            if(v): n += m
            m *= 2
        return n
    def is_words_pair(w1:str,w2:str):
        c1,c2 = [c for c in w1],[c for c in w2]
        c1.sort()
        c2.sort()
        return c1 == c2
    
    def get_map_val(word,mp):
        s = 0
        for i in range(len(word)):
            s *= 10
            s += mp[ord(word[i])]
        return math.sqrt(s)
    def is_valid_map(word,mp): return get_map_val(word,mp).is_integer()
    def get_all_maps(word):
        mp = np.zeros(91,dtype=np.uint8)
        maps = []
        inxs = []
        for i in range(len(word)-1,-1,-1):
            if(ord(word[i]) in inxs): continue
            inxs.append(ord(word[i]))
        

        if(False and len(inxs) >= math.log10(EL)): #use enders
            startV = int(math.log10(EL))
            for i in range(len(ender)):
                print(f"\ruse enders {i}\\{len(ender)} ({round(100*i/len(ender),2)}%)",end="     ")
                v = ender[i]
                for i in range(startV):
                    mp[inxs[i]] = v % 10
                    v = v//10
                while(True):
                    #if(is_valid_map(word,mp)): nums.append(int(get_map_val(word,mp)))
                    if(is_valid_map(word,mp)): maps.append(mp.copy())
                    v = startV
                    while(v < len(inxs)):
                        mp[inxs[v]] += 1
                        if(mp[inxs[v]] < 10): break
                        mp[inxs[v]] = 0
                        v += 1
                    if(v == len(inxs)): break
        else:
            g,L = 1,10**len(word)
            while(True):
                print(f"\r{g}\\{L} ({round(100*g/L,2)}%)",end="     ")
                g += 1
                #if(is_valid_map(word,mp)): nums.append(int(get_map_val(word,mp)))
                if(is_valid_map(word,mp)): maps.append(mp.copy())
                v = 0
                while(v < len(inxs)):
                    mp[inxs[v]] += 1
                    if(mp[inxs[v]] < 10): break
                    mp[inxs[v]] = 0
                    v += 1
                if(v == len(inxs)): break
        return maps
    
    if(False):
        t1 = time.perf_counter()
        MPS = get_all_maps('ABCD')
        t2 = time.perf_counter()
        print(f"found maps {MPS}\nwith len {len(MPS)}")
        print(f"took [{t2-t1}]")
        exit()



    data = read_data()
    words = data.split(",")
    words = [w[1:-1] for w in words]
    max_len = max([len(w) for w in words])
    words_by_letter = []
    for i in range(max_len+1): words_by_letter.append([])
    for w in words: words_by_letter[len(w)].append((word_to_code(w),w))
    angs = []
    for arr in words_by_letter:
        arr.sort(key=lambda x:x[0])
        lastM = 0
        for i in range(1,len(arr)):
            if(arr[i-1][0] == arr[i][0]):
                for g in range(lastM,i):
                    if(not is_words_pair(arr[g][1],arr[i][1])): continue
                    #print(f"anagram {arr[g][1]} <-> {arr[i][1]}")
                    angs.append((arr[g][1],arr[i][1]))
            else: lastM = i
    print(angs[-1])
    angs[-1] = (angs[-1][1],angs[-1][0])

    best = 0
    for g,ang in enumerate(angs):
        if(g != len(angs)-1): continue
        print(f"\rat {g}\\{len(angs)} with [{ang[0]}] best is [{best}]            ",end="")
        #if(ang[0] != "RACE" and ang[1] != "RACE"): continue
        maps1 = get_all_maps(ang[0])
        for i in range(len(maps1)-1,-1,-1):
            if(is_valid_map(ang[1],maps1[i])):
                val = get_map_val(ang[1],maps1[i])
                if(val > best): best = val

    print(f"\n\nbest is {best**2}")
        
if(False): #q205
    def get_values_change_p(probs,dices):
        if(dices == 1): return probs
        return np.convolve(get_values_change_p(probs,dices-1),probs)
    
    def get_values_change(values,dices):
        pobs = get_values_change_p(np.full(values,1/values,dtype=np.float64),dices)
        pobs = np.concatenate((np.zeros(dices),pobs))
        return pobs
    
    def print_pobs(pobs):
        for i,p in enumerate(pobs):
            print(f"for [{i}] have [{p}]")

    probs1 = get_values_change(4,9)
    probs2 = get_values_change(6,6)
    sm = 0
    for i1,p1 in enumerate(probs1):
        sm += p1*np.sum(probs2[:i1])
    print(f"sm is {sm}")
    print(f"sm round to 7 is {round(sm,7)}")


if(False): #q127
    """sm = 0
    def rad(n):
        nm = 1
        #for v in FC.get_prime_factors(n): nm *= v
        for v in PR.get_factors_primes(n): nm *= v
        return nm
        #return np.prod(PR.get_factors_primes(n))
    t1 = time.perf_counter()
    for c in range(3,2000):
        #facs = PR.get_factors_primes(c)
        if(PR.is_prime(c)): continue
        sbrad = rad(c)
        if(sbrad >= c): continue
        for a in range(1,c//2):
            b = c-a
            if(math.gcd(c,a) != 1 or math.gcd(c,b) != 1 or math.gcd(b,a) != 1): continue
            if(rad(a*b)*sbrad >= c): continue
            print(f"found {(a,b,c)}")
            sm += c
    print(f"\nsm is [{sm}]")
    print(f"took {time.perf_counter()-t1}")"""

    t1 = time.perf_counter()
    L = 120000
    factors = []
    factors_mult = []
    for i in range(L): 
        factors.append(sp.primefactors(i))
        mv = 1
        for f in factors[-1]: mv *= f
        factors_mult.append(mv)

    sm = 0
    for c in range(L):
        print(f"\rat {c}\\{L} ",end="")
        if(factors_mult[c] >= c): continue
        for a in range(1,c//2):
            if(math.gcd(c,a) != 1): continue
            if(factors_mult[c]*factors_mult[a]*factors_mult[c-a] >= c): continue
            #print(f"found {(a,c-a,c)}")
            sm += c
    print(f"\n\nsm is [{sm}]")
    print(f"took [{time.perf_counter()-t1}]")

if(False): #q121
    #N = 5
    #vals = np.zeros(N)
    #for i in range(2,2+N): vals[i] = 1/i

    def play(a,b,left):
        if(left == 0): return 1
        if(b == a): return 0
        p = 1/a
        return p*play(a+1,b,left-1)+(1-p)*play(a+1,b,left)

    num = play(2,17,8)
    print(num)
    print(1/num)
    print(f"would pay {math.floor(1/num)}")

if(False): #q187
    L = 10**8    
    LS = math.sqrt(L)
    PR.load_data(L)
    primes = np.array(PR.load_primes_until(L))
    sm = 0
    for i,p in enumerate(primes):
        print(f"\rat {i}\\{len(primes)}",end="    ")
        if(p > LS): break
        inx = np.searchsorted(primes,L/p,side="left")
        sm += inx-i
    print(f"\nsm is {sm}")

if(False): #q131
    L = 10**6
    cubs = []
    for i in range(L): cubs.append(i*i*i)
    cubs = np.array(cubs)
    def is_cube(n):
        inx = np.searchsorted(cubs,n)
        return 0 <= inx < len(cubs) and cubs[inx]==n

    CS,CI = 0,0
    PR.load_data(L)
    primes = PR.load_primes_until(L)
    sm = 0
    for i2,p in enumerate(primes):
        for i in range(CS,min(len(cubs),CS+i2-CI)):
            n = cubs[i]
            #if(n > p*p): break
            if(is_cube(n+p)):
                print(f"found p={p} and n={n} test {round(math.pow(n,1/3))} cs is {CS} and limit is {CS+i2-CI}")
                sm += 1
                CS,CI = i+1,i2
    print(f"found {sm}")
            

if(False): #q204
    def real_force(mx,lim):
        sm = 1
        for i in range(2,lim):
            if(max(PR.get_factors_primes(i)) <= mx): 
                print(f"aprove [{i}]")
                sm += 1
        return sm
            
    def brute_force(primes,lim,a=0):
        #if(len(primes)-a == 0 or primes[a] >= lim): return 1
        if(len(primes)-a == 1): return math.ceil(math.log(lim,primes[a]))
        v = 1
        am = 0
        while(v < lim):
            am += brute_force(primes,lim/v,a+1)
            #print(f"v is {v} and add {brute_force(primes,lim/v,a+1)}")
            v *= primes[a]
        return am
    primes = PR.load_primes_until(100)
    print(primes)
    print(brute_force(primes,10**9+1))
    #print(real_force(primes[-1],27))
         
if(False): #q203
    def get_pasc(prev,l=1):
        arr = [1]
        for i in range(1,len(prev)):  arr.append(prev[i-1]+prev[i])
        arr.append(1)
        if(l<=1): return arr
        return get_pasc(arr,l-1) + arr
   
    nums = set(get_pasc([1],50))
    mx = max(nums)
    print(f"max is {mx}")
    primes = PR.load_primes_until(math.sqrt(mx))
    squares = []
    for p in primes: squares.append(p*p)

    sm = 0
    for num in nums:
        print(f"\rat {num}\\{mx}                           ",end="")
        good = True
        for sq in squares:
            if(num % sq == 0): good = False; break
        if(good): sm += num
    print(f"sm is {sm}")

if(False): #q132
    VALSR = [0]
    def R(n): 
        v = 0
        for i in range(n):
            v *= 10
            v += 1
        return v
   
    
    facs = []
    inxs = {}
    par = [[]]
    N = 10**9

    primes = PR.load_primes_until(200000)
    def get_pasabo_factors(n):
        facs = []
        for p in primes:
            if(p > n): break
            if(n % p == 0): facs.append(p)
        return facs

    for i in range(1,100000):
        if(N % i != 0): continue
        #facs.append(PR.get_factors_primes(R(i)))
        facs.append(get_pasabo_factors(R(i)))
        print(f"R({i}) factors are {facs[-1]}")
        mambers = []
        for fac in facs[-1]:
            if(not inxs.get(fac) is None): continue
            inxs[fac] = []
            mambers.append(fac)
        par.append(mambers)
    
    devs = []
    for i,p in enumerate(par):
        #if(i == 0 or N % i != 0): continue
        devs.extend(p)
    print(par)
    devs.sort()
    print(f"devs are {devs} have {len(devs)}")
    print(f"max dev is {max(devs)}")
    print(f"full sum dev is {sum(devs)}")
    print(f"sum part dev is {sum(devs[:min(40,len(devs))])}")


            
    
    def f(n): return n*3-1
    i = 2
    per = set(facs[f(1)])
    while(f(i)<len(facs)):
        per.intersection_update(set(facs[f(i)]))
        i += 1

if(False): #q134 (fail)
    PR.load_data(1000004)
    primes = PR.load_primes_until(1000003)
    print(primes[-5:])

    def get_pasbo_factors(n,prms):
        facs = []
        for p in prms:
            if(p > n): break
            if(n % p == 0): facs.append(p)
        return facs

    def S(i):
        p1,p2 = primes[i],primes[i+1]
        #print(f"p1={p1}  p2={p2}")
        M = 10**(len(str(p1)))
        #M = 10**(math.floor(math.log10(p1))+1)
        #brute fource
        if(0):
            v = p2
            while(v%M != p1): v += p2
            return v
        else:
            return ((pow(int(p2),-1,int(M))*p1)%M)*p2

    
    print(S(4))
    print(pow(7,-1,4))

    sm = 0
    for i in range(2,len(primes)-1):
        print(f"\r{i}\\{len(primes)}   ",end="")
        sm += S(i)
        #if(i > 2000): break
    print(f"\n\nsm is {sm}")


if(False): #q135
    def solve_quad(a,b,c):
        root_in = b*b-4*a*c
        if(root_in < 0): return (-1,-1)
        ab,root = -b/(2*a),math.sqrt(b*b-4*a*c)/(2*a)
        return (ab-root,ab+root) if root>0 else (ab+root,ab-root)

    LIM = 10**6

    def calculate_sum(a,d):
        if(a-d-d <= 0 or d<=0): return 0
        return a*a-(a-d)*(a-d)-(a-d-d)*(a-d-d)
    def get_sums(A):
        z1,z2 = solve_quad(-5,A*6,-A*A)
        t1,t2 = solve_quad(-5,A*6,-A*A-LIM)
        z1,z2,t1,t2 = math.floor(z1),math.ceil(z2),math.ceil(t1),math.floor(t2)
        if(t1 == -1): t1,t2 = z2,z2
        #print(f"from {z1}-{t1} and {z2}-{t2}")
        sums = []
        for i in range(z1,t1): sums.append(calculate_sum(A,i))
        for i in range(t2,z2): sums.append(calculate_sum(A,i))
        #if(27 in sums): print(f"at {A}")
        return sums
    
    print(get_sums(3))
    #exit()

    L = 10**6
    vals = []
    for i in range(1,L):
        print(f"\r{i}\\{L}",end="    ")
        vals.extend(get_sums(i))
    #print(vals)
    print(f"have {len(vals)} vals")
    count = np.zeros(L)
    for v in vals: 
        if(v >= L or v < 0): continue
        count[v] += 1
        
    sm = 0
    for i,c in enumerate(count):
        if(c == 10):
            print(f"first is {i}")
            sm += 1
    print(f"\nsm is {sm}")

if(False): #q197
    N = 30.403243784
    M = math.pow(10,-9)
    def f(x,n=1):
        v = math.floor(math.pow(2,N-x*x))*M
        if(n == 1): return v
        return f(v,n-1)
    v = -1
    vals = []
    for i in range(1000):
        print(f"u[{i}] = {v + f(v)}")
        v = f(v)
        vals.append(v)
    
    vals2 = []
    for i in range(1,len(vals)):
        vals2.append(vals[i-1]+vals[i])
    vals = vals2

    plt.plot(np.arange(1,len(vals)+1),np.array(vals))
    plt.grid()
    plt.show()

if(False): #q188

    M = 10**8
    n = 9
    for i in range(100000):
        if(n > M and ((n % M) % 3)== 0): print(f"found at {i}"); print(((n))); break
        n *= 3
    print("didnt found")
    exit()

    def tat(a,b):
        if(b == 1): return a
        return a**tat(a,b-1)
    def tatM(a,b,M):
        if(b == 1): return a
        p = tat(a,b-1)%M
        return a**p

    sys.set_int_max_str_digits(500000)
    M = math.ceil(math.log(10**8,3))
    print(M)
    for i in range(1,6):
        v1 = str(tat(2,i))
        if(len(v1) > 8): v1 = v1[-8:]
        print(v1)
        v2 = str(tatM(2,i,M))
        if(len(v2) > 8): v2 = v2[-8:]
        print(v2)

if(False):
    L = 10**7
    def get_d(n):
        return len(sp.divisors(n))
    sm = 0
    pd = get_d(1)
    for i in range(2,L):
        cd = get_d(i)
        if(pd == cd): sm += 1
        #if(pd == cd): print(f"found {i-1} - {i}")
        pd = cd
        print(f"\r{i}\\{L} ({round(100*i/L,2)}%) with {sm}    ",end="")
    print(f"\n\nfound {sm}")

if(False): #q162
    def fac(n,a=0):
        if(n <= a): return 1
        return n*fac(n-1,a)
    #def fix_put(size,n):
    #    return fac(size,size-n)
    def brute_force(n):
        L = 16**n
        sm = 0
        for i in range(L):
            hx = hex(i)
            good = hx.find('a') != -1 and (hx.find('0',2) != -1 ) and hx.find('1') != -1
            if(good): sm += 1
            #if(good): print(f"{hx} {'+' if good else ''}")
        return sm


    def solve_nf(n):
        return (n-1)*(n-1)*(n-2)*(16**(n-3))
    
    def solve_n(n:int) -> int:
        fl = 15*(16**(n-1))
        wl = 13**n #without a b and c
        wl += 2*(14*(15**(n-1))) + 15**n #without a or without b or without c
        wl -= 2*(14**n) + 13*(14**(n-1))
        return fl-wl        
        fs:int = 2*(n-1)*(n-2)
        fs *= pow(16,n-3)
        if(n == 3): ss = 0
        else:
            ss:int = (n-1)*(n-2)*(n-3)
            ss *= 15
            ss *= pow(16,n-4)
        return fs + ss
    
    sm = sum([solve_n(i) for i in range(3,17)])
    print(sm)
    #print(solve_nf(4))
    #print(hex(solve_nf(4)).upper())
    #print(brute_force(5))

    print(hex(sm).upper())

if(False): #q164
    blocks = []
    for i in range(0,1000):
        block = (i//100,(i//10)%10,i%10)
        if(sum(block) <= 9): blocks.append(block)
    nexts = []
    for block in blocks:
        conts = []
        for i,b2 in enumerate(blocks):
            if(b2[0] == block[1] and b2[1] == block[2]): conts.append(i)
        nexts.append(conts)
    
    BREAK_INXS2 = [4,8,12,16]
    BREAK_INXS = []
    BREAK_COUNTS = []

    def count(inx,amount):
        #print(f"count {inx} and {amount}")
        if(amount in BREAK_INXS):
            return BREAK_COUNTS[BREAK_INXS.index(amount)][inx]

        nx = nexts[inx]
        if(amount == 1): return len(nx)
        sm = 0
        for i in nx: sm += count(i,amount-1)
        return sm
    
    for brk in BREAK_INXS2:
        print(f"load brk {brk}")
        vals = []
        for i in range(len(blocks)):
            vals.append(count(i,brk))
        BREAK_COUNTS.append(vals)
        BREAK_INXS.append(brk)
        

    S = 1
    print(blocks[S])
    for v in nexts[S]: print(f"[{blocks[v]}]")
    
    sm = 0
    for i,block in enumerate(blocks):
        if(block[0] == 0): continue
        sm += count(i,17)
    print(f"sm is {sm}")

if(False): #q165
    global SRAND
    SRAND = 290797
    def next_point():
        global SRAND
        #vals = [27,44,12,32,46,53,17,62,46,70,22,40]
        #if(SRAND > 20): SRAND = 0
        #v = vals[SRAND];SRAND += 1;return v
        SRAND = SRAND*SRAND % 50515093
        return SRAND % 500
    

    lines = []
    slops = []
    bvals = []
    STR = []

    def true_point(i1,i2):
        if(slops[i1] == slops[i2]): return None
        x = (bvals[i2]-bvals[i1])/(slops[i1]-slops[i2])
        if(x <= gen.Frac(max(lines[i1][0][0],lines[i2][0][0])) or x >= gen.Frac(min(lines[i1][1][0],lines[i2][1][0]))): return None
        return (x,slops[i1]*x+bvals[i1])

    def true_point_str(st,li):
        line = lines[li]
        stx = STR[st][0][0]
        if(line[0][0] >= stx or line[1][0] <= stx): return None
        ly = slops[li]*stx + bvals[li]
        sty1 = min(STR[st][0][1],STR[st][1][1])
        sty2 = max(STR[st][0][1],STR[st][1][1])
        if(ly <= gen.Frac(sty1) or ly >= gen.Frac(sty2)): return None
        return (stx,ly)

    for i in range(5000):
        line = [[next_point(),next_point()],[next_point(),next_point()]]
        if(line[0][0] > line[1][0]): line = [line[1],line[0]]
        if(line[0][0] == line[1][0]):
            STR.append(line)
            continue
        #slops.append((line[0][1]-line[1][1])/(line[0][0]-line[1][0]))
        #bvals.append(line[0][1]-slops[-1]*line[0][0])
        slops.append(gen.Frac(line[0][1]-line[1][1],line[0][0]-line[1][0]))
        bvals.append((slops[-1]*line[0][0]*-1)+line[0][1])
        lines.append(line)
    print([str(b) for b in bvals[:5]])
    points = []

    for i in range(len(STR)):
        for g in range(len(lines)):
            p = true_point_str(i,g)
            if(p is not None): points.append(p)

    for i in range(len(lines)):
        for g in range(i+1,len(lines)):
            p = true_point(i,g)
            if(p is not None): points.append(p)
        print(f"\r{i}\\{len(lines)} ({round(100*i/len(lines),2)}%)   ",end="")

    for i in range(len(points)): points[i] = (round(points[i][0],10),round(points[i][1],10))

    #print(points)    
    print(f"\n{len(points)}")
    points = list(set(points))
    print(f"\n{len(points)}")
    points.sort(key=lambda x:x[0]*x[0]+x[1]*x[1])
    #print(points)
    print(f"min points x: {min([((points[i][0]-points[i-1][0])**2+(points[i][1]-points[i-1][1])**2) for i in range(1,len(points))])}")
    #points.sort(key=lambda x:x[1])
    #print(f"min points y: {min([points[i][1]-points[i-1][1] for i in range(1,len(points))])}")
    #print(lines[0])
    #print(lines[1])
    #print(lines[2])

if(False): #q216
    def t(n): return 2*(n*n)-1


    def brute_force(l):
        sm = 0
        for i in range(2,l):
            if(PR.is_prime(t(i))): sm += 1
        return sm
    def brute1(l):
        prs = []
        for i in range(8,math.ceil(math.sqrt(2)*l),8):
            if(PR.is_prime(i-1)): prs.append(i-1)
            if(PR.is_prime(i+1)): prs.append(i+1)
        return prs
        sm = 0
        for i in range(2,l):
            n = t(i)
            lm = math.sqrt(n)
            good = True
            for p in prs:
                if(p > lm): break
                if(n % p == 0): good = False; break
            if(good): sm += 1
        return sm

            
            
        
    facs = []
    for i in range(2,100):
        n = t(i)
        facs.extend(sp.divisors(n))
        print(f"t({i}) = {n} -> {'' if PR.is_prime(n) else sp.divisors(n)[1:-1]}")
    facs = list(set(facs))
    facs.sort()
    print(facs[-1])
    #print(facs)

    ts = set(brute1(1000))
    ts.difference_update(set(facs))
    ts = list(ts)
    ts.sort()
    print(ts)
    exit()
    
    t1 = time.perf_counter()
    print(brute_force(10000))
    t2 = time.perf_counter()
    print(brute1(10000))
    t3 = time.perf_counter()
    print(f"took [{round(t2-t1,3)}] and [{round(t3-t2,3)}]")

if(False): #q235
    def calc(r):
        sm = 0
        for i in range(1,5001):
            sm += (900-3*i)*math.pow(r,i-1)
        return sm
    
    T = -6*(10**11)
    jump = 0.25
    val = 1.05
    for i in range(1000):
        res = calc(val)
        if(res > T): val += jump
        else: val -= jump
        jump /= 2
    print(calc(val))
    print(val)

if(False): #q243
    L = 15499/94744
    G = 10**9
    MN = 1
    for i in range(0,G,910):
        if(i == 0): continue
        v = PR.phi(i)/(i-1)
        if(v < MN): print(f"at {i}")
        MN = min(v,MN)
        print(f"\r{i}\\{G} with min {round(MN,10)} r: [{round(MN/L,10)}]   ",end="")
        if(v < L):
            print(f"\nfound at {i}")
            break

if(False): #q76B
    BEL = [[0],[0,1]]
    def get_ways2(n):return get_ways(n,n)
    def get_ways(n,mx):
        if(n < len(BEL)): return BEL[n][min(n,mx)]
        ways = [0]
        sm = 0
        for i in range(1,n):
            sm += get_ways(n-i,i)
            ways.append(sm)
        ways.append(sm+1)
        BEL.append(ways)
        return sm
    
    
    for i in range(2,10001):
        print(f"w({i}) = {get_ways2(i)}")

    #print(BEL)

if(False): #q800
    def find_for_p(L,p):
        #targ = L*math.log(L,p)
        targ = math.log(L,p)*L
        x = math.ceil(targ)
        
        
        """while(x + p*math.log(x,p) > targ and (not PR.is_prime(x))): x -= 1
        inx = PR.get_prime_index(x)
        if(inx != -1):
            while(x + p*math.log(x,p) > targ):
                x = PR.primes[inx - 1]
                inx -= 1
                if(inx == -1): return 0"""
        
        def test(v): 
            return v + p*math.log(v,p) <= targ
        x = gen.cond_binary_search(0,x,test)

        #print(x)
        #x = PR.get_prime_number(x)
        
        #mults = PR.load_primes_until(x)
        #if(len(mults)==0): return 0 
        #size = len(mults)
        size = PR.count_primes_until(x)
        #ize = len(mu)
        if(size == 0): return 0
        if(p <= x): size -= 1
        #print(f"{mults} - [{size}] with p = {p}")
        return size
    
    #print(find_for_p(800,19))
    #exit()
    L = 800800
    p = 2
    sm = 0
    LL = math.ceil(math.log2(L)*L)
    PR.load_data(LL)
    while(p < LL):
        sm += find_for_p(L,p)
        p = PR.get_prime_after(p)
        print(f"\r{p}\\{L}   ",end="")
    sm /= 2
    print(sm)

if(False): #q684

    MOD = 10**9+7
    print(MOD)

    def s(n):
        nines = (n-1)//9
        return pow(10,nines,MOD)*(n-nines*9+1) - 1
        

        M = 10**nines
        N1 = M-1
        return (M*(n-nines*9)+N1)
        
        #if(n < 10): return 10 + (n-1)
        digs_len = (n-1)//9 + 1
        digs = np.full(digs_len,9)
        s = 9*digs_len
        inx = 0
        while(s-8>=n):
            digs[inx] = 1
            s -= 8
        digs[inx] -= s-n
        #print(f"for n={n} sum digs are [{sum(digs)}] digts are {digs} and num is {gen.digts_to_num(digs)}")
        return gen.digts_to_num(digs)
    
    PREV_S = [0]

    def ls(n):
        return sum([s(i) for i in range(1,n+1)])%MOD
        if(n < len(PREV_S)): return PREV_S[n]
        #prev = ls(n-1)
        prev = ls(len(PREV_S)-1)
        for i in range(len(PREV_S),n+1):
            prev += s(i)
            prev %= MOD
            PREV_S.append(prev)
        return prev

    def ls2(n):
        dev,md = divmod(n,9)
        primes = 6*pow(10,dev,MOD)-6-9*dev
        for i in range(n-md+1,n+1): primes += s(i)
        return primes % MOD


    print(ls2(1118))
    print(ls(1118))


    #print(s(18))
    #print([s(i) for i in range(1,21)])
    #print(sum([s(i) for i in range(1,21)]))

    a,b = 1,1
    sm = 0
    for i in range(2,91):
        print(f"\rat {i}\\{91}   ",end="")
        a,b = b,a+b
        print(f"f_{i} = {a}")
        sm += ls2(a)
    #print(sm)
    print(sm % MOD)

if(False): #q686
    
    def get_min_jump(num,targ,goal1):
        dev,goal = 1,goal1
        p_goal = 2; p_val = 1
        v1,v2 = num,num
        for i in range(1,1000):
            v1 <<= 2
            #v2 <<= 2; v2 += 1
            if(v1 > goal): dev  *= 10; goal *= 10
            if(i > p_goal): p_val += 1; p_goal <<= 2
            d1 = v1//dev
            if(d1 > targ): continue
            if(d1+p_val < targ): continue
            return i


    def brute_force1(l,n):
        size = len(str(l))
        start = math.floor(math.log2(size))
        goal1 = 10**(size)
        dev,goal = 1,10**(size)
        num = 1
        c = 0
        for i in range(1,10000000):
            num <<= 1
            #num_size = math.floor(math.log10(num))+1
            #val = num // (10**(num_size-size))
            ext = num > goal
            if(ext): 
                dev  *= 10; goal *= 10
            val = num//dev
            if(val == l): c += 1
            if(c == n): return i
            if(ext): i += get_min_jump(val,l,goal1)

    t1 = time.perf_counter()
    print(brute_force1(123,678))
    print(f"took [{round(time.perf_counter()-t1,2)}]")

if(False): #q853
    def get_phi(n):
        a,b = 1,1
        while(b < n):
            a,b = b,a+b
        a,b = b%n,(a+b)%n
        t1,t2 = a,b

        for i in range(1000000):
            a,b = b,(a+b)%n
            if(a == t1 and b == t2): return i+1
            #print(a)
    #print(f"res is {get_phi(13*19*5)}")

    print(get_phi(5*19*5))
    exit()

    primes = PR.load_primes_until(100000)
    allowed_primes = []
    for p in primes:
        ph = get_phi(p)
        print(f"f({p}) = {ph}")
        if(ph <= 120): allowed_primes.append(p)
    
    print(allowed_primes)
    print(len(allowed_primes))


if(False):
    pass
    LIM = 100
    RES = np.zeros(LIM)
    nums = []
    sm = 0
    for i in range(2,LIM):
        if(RES[i] != 0): continue
        RES[i] = i
        exe = [i]
        for n in nums:
            v = i*n
            if(v >= LIM): break
            RES[v] = i
            exe.append(v)
        nums.extend(exe)
        nums.sort()
    
    for i,r in enumerate(RES):
        print(f"f({i}) = {r}")

if(False): #q128
    lim = 1000
    TOP = np.zeros(lim)
    TOP[1] = 1
    v,tinx = 2,1
    while(v < lim):
        TOP[v] = 1
        print(v)
        v += 6*tinx
        tinx += 1
    
    def is_top(n): return TOP[n]==1


if(False): #q86
    def is_int(a:float): return abs(abs(a%1-0.5)-0.5)<0.00001
    def d(a,b): return math.sqrt(a*a+b*b)
    def get_shortets(a,b,c):
        x = a*b/(b+c)
        return d(x,b) + d(a-x,c)
    
    L = 10**6
    M = 0
    sm = 0
    while(sm < L):
        print(f"\rm is [{M}] sm is [{sm}\\{L}]",end="")
        M += 1
        for b in range(1,M+1):
            for c in range(1,b+1):
                if(is_int(get_shortets(M,b,c))): sm += 1
    print(sm)
    print(M)

if(False): #q173
    L = 1000001
    sm = 0
    start = 1

    a = 1
    while(a < L):
        #print(a)
        print(f"\r{a}\\{L}   ",end="")
        cs = start + (0 if (start%2==a%2) else 1)
        if(a*a-cs*cs < L): sm += (a-cs)//2
        else:
            sa = a*a
            while(sa-start*start > L): start += 1
            #print("<>")
            a -= 1
        a += 1


    print(sm)
        
if(False): #q174
    LIM = 1000000
    def resiveL(n):
        arr = []
        limit = LIM + n*n
        num = n + 2
        while(num * num <= limit):
            arr.append(num*num - n*n)
            num += 2
        return arr

    lists = []
    for i in range(1,LIM):
        res = resiveL(i)
        if(len(res) == 0): break
        lists.append(res)

    print(len(lists))

    appr = np.zeros(LIM+1)
    for arr in lists:
        for v in arr:
            appr[v] += 1
    
    def N(n):
        sm = 0
        for v in appr:
            if(v == n): sm += 1
        return sm


    print(sum([N(i) for i in range(1,11)]))

if(False): #q188
    MOD = 10**8
    exper = 1777
    ml = exper
    for i in range(1854):
        ml = pow(exper,ml,MOD)
    print(ml)

if(False):
    text = read_data()
    data = []
    for line in text.split("\n"):
        data.append(([int(v) for v in list(line[:line.index(" ")])],int(line[line.index(';')+1])))
    data.sort(key = lambda d: d[1])
    print(data)

    number_size = len(str(data[0][0]))
    
    def is_possible(options): #for some part of number return if possible
        for example in data:
            minC,maxC = 0,0
            for i,digit in enumerate(example):
                if(options == -1): maxC += 1;
                elif(num[i] == digit): 
                    minC += 1; maxC += 1;
            if(example[1] < minC or example[1] > maxC): return False
        return True

    options = np.zeros((number_size,10))

if(False): #q103
    
    def CanAdd(arr,n):
        if(len(arr) == 0): return True
        if(len(arr) == 1): return arr[-1] < n        
        if(len(arr) % 2 == 0):
            total_sum = sum(arr)
            first_sum = sum(arr[:len(arr)//2 + 1])
            if(2*first_sum <= total_sum + n): return False
        else:
            sum1 = sum(arr[:(len(arr)//2 + 1)])
            sum2 = sum(arr[(len(arr)//2 + 2):])
            if(sum1 <= sum2 + n): return False
        return pass_rule1(arr + [n])
    
    def get_all_sums(arr, n):
        if(n == 1): return arr
        sums = []
        for i,v in enumerate(arr):
            part = get_all_sums(arr[i+1:],n-1)
            for p in part:
                sums.append(p+v)
        return sums

    def pass_rule1(arr):
        for i in range(2,len(arr)-1):
            allSums = get_all_sums(arr,i)
            if(len(allSums) != len(set(allSums))): return False
        return True
    
    MAX = 50
    BIG_INT = 10000
    def extandsSet(arr:list, size, limit=BIG_INT):
        if(size == 0): return [v for v in arr]
    
        maxVal = MAX if len(arr) < 2 else min(MAX,arr[0] + arr[1])
        maxVal = min(maxVal, limit - sum(arr))

        best, bestSum = None,limit
        for i in range((0 if len(arr) == 0 else arr[-1]) + 1,maxVal):
            if(not CanAdd(arr,i)): continue
            arr.append(i)
            extend = extandsSet(arr,size - 1,bestSum)
            arr.pop(-1)
            if (extend is None):
                continue
            if (size == 1): return extend
            if (bestSum > sum(extend)):
                best,bestSum = extend, sum(extend)
        return best
    
    print(CanAdd([2,3],4))
    t1 = time.time()
    res = extandsSet([],7)
    print(res)
    print(f"took {time.time() - t1}")
    for i in res:
        print(i,end="")

if(False):
    def isProg(n):
        lim = math.sqrt(n)
        for v in range(n-1,int(lim),-1):
            rat = (v * v) / n
            if(rat.is_integer() and rat < v): return True
        return False
    
    def isProg2(n):
        lim = int(math.pow(n,2/3))+1
        for d in range(2,lim):
            q,r = divmod(n,d)
            if(r == 0): continue
            if(d*r == q*q): return True
        return False
    
    def isProgForRemainder(n,r):
        mn = (n-r)/r
        d = round(math.pow(r*mn*mn,1/3),9)
        return d.is_integer()

    def isProg3(n):
        if(n < 100000): return isProg2(n)
        C = 729
        C = max(min(C,int(math.sqrt(n))-2),0)
        
        for i in range(1,C+1):
            if(isProgForRemainder(n,i)): return True
        
        C += 1
        lim = int(math.pow((n-C)*(n-C)/C,1/3)) + 1

        for d in range(2,lim):
            q,r = divmod(n,d)
            if(r == 0): continue
            if(d*r == q*q): return True
        return False
        

    for i in range(2,10000):
        num = i
        if(isProg2(num) != isProg3(num)):
            print(f"found at {num} with p2: {isProg2(num)}")

    M = 10**5
    MLIM = math.sqrt(M)
    t1 = time.time()
    res = sum([(n*n if isProg3(n*n) else 0) for n in range(2,int(MLIM))])
    print(time.time() - t1)
    print(res)

if(False): #q191
    arrs = read_data().split(" ")
    
    print(arrs)

    print(sum([(0 if name.find("L") != -1 else 1) for name in arrs]))
    print([(0 if name.find("L") != -1 else 1) for name in arrs])

    def brute_possible_a(size,rs):
        arrs = np.zeros((2**size,size))
        for i in range(len(arrs)):
            v = i
            for g in range(size):
                arrs[i][g] = v & 1
                v //= 2
        sm = 0
        for arr in arrs:
            tmp = sum(arr[:rs])
            for i in range(rs,size):
                if(tmp == rs): sm += 1; break
                tmp -= arr[i-rs]
                tmp += arr[i]
                if(tmp == rs): sm += 1; break
        return sm
        
    def get_possibles_a(size,rs):
        if(size == rs): return 1
        if(size < rs): return 0
        sm = 0
        for i in range(size-rs+1):
            v1 = 2**(size-rs-i)
            v2 = 1
            if(i != 0):
                v2 = 2**(i-1) - get_possibles_a(i-1,rs)
            sm += v1*v2
        return sm
    
    def get_without_l(days,rs):
        return (2**days) - get_possibles_a(days,rs)
    
    DAYS = 30
    RS = 3
    sm = get_without_l(DAYS,RS)
    
    for i in range(DAYS):
        v1 = get_without_l(i,RS)
        v2 = get_without_l(DAYS - 1 - i,RS)
        sm += v1*v2
    print(sm)

if(False): #q231
    FC = gen.FactorsSaver(PR)

    def get_factors_sum(num):
        facs = FC.get_factors(num)
        #facs = sp.factorint(num)
        facs:dict
        sm = 0
        for key,val in facs.items():
            sm += key*val
        return sm

    def solve(n,k):
        for i in range(n):
            if(i % 100 == 0): print(f"\rfac {i}\\{n}  ",end="")
            FC.get_factors(i)

        prods = np.arange(n-k + 1,n+1)
        sm1 = sum([get_factors_sum(v) for v in prods])
        devs = np.arange(2,k+1)
        sm2 = sum([get_factors_sum(v) for v in devs])
        #print(prods)
        #print(devs)
        return sm1 - sm2
    
    a,b = 20000000,15000000
    #tot = math.comb(a,b)
    t1 = time.perf_counter()
    r1 = solve(a,b)
    t2 = time.perf_counter()
    r2 = 2
    print(f"got {r1} and {r2} at [{t2-t1}] and {time.perf_counter() - t2}")

if(False): #q250
    L = 250250
    M = 5
    GM = 10**9

    elementsMod = []
    for i in range(1,L+1):
        elementsMod.append(pow(i,i,M))
    counts = np.zeros(M,dtype=np.uint)
    for ele in elementsMod:
        counts[ele] += 1
    
    counts = np.full(5,400)

    #create factorial data
    print(max(counts[1:]))
    factorial_vals = [1]
    for i in range(1,max(counts[1:])+1):
        factorial_vals.append(factorial_vals[-1]*i)

    def comb(n,k) -> int:
        return factorial_vals[n] // (factorial_vals[k]*factorial_vals[n-k])
    
    def rool(arr:list,amount:int) -> list[int]:
        res = []
        for i in range(len(arr)):
            res.append(arr[(i-amount) % len(arr)])
        return res
    
    def mul_array(arr, t) -> list[int]:
        return [int((v*t)) for v in arr]
    def add_array(arr, addArr) -> list[int]:
        return [int(arr[i] + addArr[i]) for i in range(len(arr))]

    def get_ways_arr(value,amount,mod):
        ways = list(np.zeros(mod,np.uint64))
        for i in range(amount + 1):
            #print(f"\rway {i}\\{amount}  ",end="")
            ways[value*i % mod] += comb(amount,i) 
            #ways[value*i % mod] %= GM
        return ways
   
    def get_ways_to_rep(amounts, mod):
        rep = list(np.zeros(mod,np.uint64))
        rep[0] = 1
        for i,a in enumerate(amounts):
            print(f"\r{i} \\ {len(amounts)} ",end="")
            result = list(np.zeros(mod,np.uint64))
            ways = get_ways_arr(i,a,mod)
            for g,w in enumerate(ways):
                result = add_array(result,mul_array(rool(rep,g),w))
            print(result)
            for i,v in enumerate(result):
                rep[i] = v % GM
        return rep
    
    def get_ways_to_rep2(amounts, mod):
        pol = np.poly1d(1)
        for i,a in enumerate(amounts):
            arr = np.zeros(i+1)
            arr[0] = 1
            arr[i] = a
            pol *= np.poly1d(arr)
        print(np.array(pol,dtype=np.uint64))

    get_ways_to_rep2(counts,250)
    exit()

    #print(get_ways_arr(0,counts[0],250))
    zeorsAmount = counts[0]
    #counts[0] = 0
    zlways = get_ways_to_rep(counts,M)
    print(zlways[0])
    
    #add options of zeros
    #zeros = pow(2,int(zeorsAmount),GM)
    #print(f"zeros {zeros}")
    #print(f"res {(zlways[0]*zeros) % GM}")

if(False): #q265
    def get_possible_bits(size):
        if (size == 1): return [[0],[1]]
        prev = get_possible_bits(size - 1)
        arr = [[0] + p for p in prev]
        arr.extend([[1] + p for p in prev])
        return arr
    
    def get_opts(n):
        size = 2**n
        arrs = []
        temp = [[0]*n]
        
        #for the first start
        print(size - 2*n - 1)
        options = get_possible_bits(size - 2*n - 1)
        arr = np.zeros(size,dtype=np.uint8)
        arr[:n] = 0
        arr[n:2*n] = 1
        arr[-1] = 1
        for opt in options:
            arr[2*n:-1] = opt
            arrs.append(list(arr))

        for i in range(n+2,size-n):
            print(f"{i}\\{size}")
            arr = np.zeros(size,dtype=np.uint8)
            arr[:n] = 0
            arr[n] = 1
            arr[i:i+n] = 1
            arr[-1] = 1
            opts = get_possible_bits(size - 2*n - 2)
            for opt in opts:
                arr[n+1:i] = opt[:(i-n-1)]
                arr[i+n:-1] = opt[(i-n-1):]
                arrs.append(list(arr))
        
        return arrs
    
    print(get_opts(5))

if(False): #q348
    LIM = 10**7
    squares = []
    i = 2
    while(True):
        res = i*i
        if(res > LIM): break
        squares.append(res)
        i += 1
    cubes = []
    i = 2
    while(True):
        res = i*i*i
        if(res > LIM): break
        cubes.append(res)
        i += 1

    nums = list(set(squares).union(set(cubes)))
    nums.sort()
    res = []
    for num in nums:
        if(not gen.is_palindrom(gen.get_digits(num))): continue
        res.append(num)
        if(len(res) == 5): break
    print(res)
    print(sum(res))

if (False): #q365
    prms = PR.load_primes_until(5000)
    prms = [v for v in prms if v > 1000]
    muls = []
    for a in range(len(prms)):
        for b in range(a,len(prms)):
            for c in range(b,len(prms)):
                muls.append(prms[a]*prms[b]*prms[c])
    print(len(muls))
    
if (False): #q345
    txt = read_data()
    lines = txt.split("\n")
    colms = [line.split(" ") for line in lines]
    matrix = [[int(v) for v in line if v != ''] for line in colms]
    
    ords = np.arange(0,len(matrix))
    def get_sum():
        return sum([matrix[i][ords[i]] for i in len(matrix)])

if (False): #q346
    LIM = 10**12
    bases = int(math.sqrt(4*LIM-3)-1)
    print(bases)
    nums = []
    unc = []
    for b in range(2,bases):
        ls,a = [1],1
        while (True):
            a *= b
            nx = ls[-1]+a
            if (nx >= LIM): break
            if (a > b): unc.append(nx)
            ls.append(nx)
        nums.append(ls)
    #print(nums)
    #print(unc)
    print(f"sm {sum(set(unc))+1}")

if (False): #q347
    def M(p,q,N):
        if (p*q < 0):
            print(f"neg {p} - {q} - {p*q}")
            input("a")
        if (p*q > N): return 0
        a = int(math.log(N/q,p))
        mul = (p**a)*q
        mx = mul
        while(a > 1):
            mul //= p
            a -= 1
            while (True):
                nex = mul*q
                if (nex <= N): mul = nex
                else: break
            mx = max(mx,mul)
        return mx
    
    L = 10**7
    sm = 0
    PR.load_data(L)
    primes = PR.load_primes_until(L)
    primes = [int(p) for p in primes]
    for i,p1 in enumerate(primes):
        print(f"\r{i}/{len(primes)}   ",end="")
        for g in range(i+1, len(primes)):
            p2 = primes[g]
            if (p1*p2 > L): break
            sm += M(p1,p2,L)
    print(f"sm is {sm}")

if (False): #q348
    XR = 30000
    YR = 2000
    nums = []
    for x in range(2,XR):
        for y in range(2,YR):
            #if (x*x+y*y*y == 19225):print(f"found {x} - {y}")
            nums.append(x*x+y*y*y)
    nums.sort()
    print(len(nums))
    sc,p = 1,0
    fd = []
    for i,n in enumerate(nums):
        if (not gen.is_palindrom(gen.get_digits(n))): continue
        if (n == p): 
            sc += 1
            if (sc == 4 and nums[i+1] != p):
                print(f"found {n} {nums[i-5:i+5]}")
                fd.append(n)
                if (len(fd) == 5): break
        else: 
            sc,p = 1,n
    print(f"sm: {sum(fd)}")

if (False): #q686
    def get_n_first_digits(num, digits):
        if (num < 10**digits): return num
        num_len = int(math.log10(num))+1
        dever = 10**(num_len - digits)
        return num // dever

    def brute(L,n):
        i,c,SL,num = 0,0,str(L),1
        inx = []
        while (c < n):
            print(f"\rc is {c}\\{n}  ",end="")
            num *= 2
            i += 1
            if (get_n_first_digits(num,len(SL)) == L):
                c += 1
                inx.append(i)

        print([inx[i]-inx[i-1] for i in range(1,len(inx))])
        return i

    def spes(n):
        c = 90
        vals = []
        num = 2**c
        adders = [196,289,485]
        muls = [2**a for a in adders]
        for i in range(n-1):
            print(f"\r{i}\\{n}   ",end="")
            for g,a in enumerate(adders):
                temp = num*muls[g]
                if (get_n_first_digits(temp,3) == 123):
                    num = temp
                    vals.append(c)
                    c += a
                    break
        print(vals)
        difs = [vals[i]-vals[i-1] for i in range(1,len(vals))]
        st = set(difs)
        for s in st:
            print(f"[{s}] -> {[i for i in range(len(difs)) if difs[i] == s]}")
        at = [i for i in range(len(difs)) if difs[i] == 485]
        print([(2**vals[i]) for i in at if i < 10])

        return c

    #print(brute(123,678910)) 
    print(spes(910))

if (False): #q479
    def solve_k(k):
        s1 = math.sqrt(4*(k**11)+27*(k**10)-18*(k**8)-k**6+4*(k**3))
        s1 *= 3*math.sqrt(3)
        s1 += 2*(k**6) + 27*(k**5) - 9*(k**3)
        s1 = math.pow(s1,1/3)

        r1 = s1/(3*math.pow(2,1/3)*k)
        r1 -= math.pow(2,1/3)*(3*k-k**4)/(3*s1*k)
        r1 += k/3

        r2 = -complex(1,-math.sqrt(3))/(6*math.pow(2,1/3)*k)*s1
        r2 += (3*k-k**4)*complex(1,math.sqrt(3))/(3*math.pow(2,2/3)*k*s1)
        r2 += k/3

        r3 = -complex(1,math.sqrt(3))/(6*math.pow(2,1/3)*k)*s1
        r3 += (3*k-k**4)*complex(1,-math.sqrt(3))/(3*math.pow(2,2/3)*k*s1)
        r3 += k/3

        return r1, r2, r3
    
    def sum_mod(n,q,m):
        if (n == 1): return q % m
        if (n % 2 != 0): return (sum_mod(n-1,q,m) + pow(q,n,m)) % m
        return (pow(q,n//2,m) + 1) * sum_mod(n//2,q,m)


    
    M = 1000000007
    def S(n):
        sm = 0
        for k in range(1,n+1):
            print(f"\r{k}\\{n}   ({round(100*k/n,2)}%)  ",end="")
            #a,b,c = solve_k(k)
            #num = round(((a+b)*(b+c)*(a+c)).real)
            num = (k+1)*(1-k)
            sm += sum_mod(n,num,M)
            #sm += (num*(pow(num,n)-1)//(num-1)) % M
            #for p in range(1,n+1):   sm += num**p
        return sm
    res = S(10**6)
    print(f"\n\nres is {res % M}")
    
    #print(res % 1000000007)

if (False): #q523
    def count_sort(L):
        a,c = 1,0
        while (a < len(L)):
            if (L[a] < L[a-1]):
                t = L[a]
                for i in range(a,0,-1): L[i] = L[i-1]
                L[0] = t 
                a = 1
                c += 1
            else: a += 1
        return c
    def guess_force(n,LIM=10):
        sm = 0
        for i in range(LIM):
            nums = np.arange(n)
            np.random.shuffle(nums)
            sm += count_sort(nums)
        return sm / LIM
    def brute_force(n):
        def get_prems(nums):
            if (len(nums) == 1): return [nums]
            pv = get_prems(nums[1:])
            res = []
            for p in pv:
                for i in range(len(p)+1):
                    res.append(p[:i] + [nums[0]] + p[i:])
            return res
        arrs = get_prems([i for i in range(1,n+1)])
        #arrs = [a for a in arrs if a[2] == 4]
        #for a in arrs: print(f"{a} -> {count_sort(a)}")
        sm = sum([count_sort(arr) for arr in arrs])
        return sm / len(arrs)
    def test1(n):
        if (n == 1): return 0
        return (n*test1(n-1) + sum([i for i in range(n)])) / n
    def test2(n):
        if (n <= 0): return 0
        if (n == 4): return 3.25
        sm = 0
        for i in range(0,n):
            v = test2(i-1)
            for g in range(i,n):
                for h in range(g+1): v += h/(h+1)
            sm += v / n
        return sm
            
        
    print(brute_force(4)) 
    
    for i in range(1,9):
        print(f"{i} -> {brute_force(i)} | {test2(i)}")
    
if (False): #q849
    def brute_force(n):
        coms = []
        for i in range(n):
            for g in range(i+1,n):
                v1,v2,v3 = tuple([np.zeros(n,dtype=np.int32) for h in range(3)])
                v1[i] = 2
                v2[i] = 1; v2[g] = 1
                v3[g] = 2
                coms.append([v1,v2,v3])
        def get_outs(arr):
            if (len(arr) == 1): return set([tuple(a) for a in arr[0]])
            prev = get_outs(arr[1:])
            res = set([])
            for p in prev:
                for a in arr[0]:
                    res.add(tuple(np.add(p,a)))
            return res
        print(f"coms len is {len(coms)}")
        pos = get_outs(coms)
        return len(pos)
    print(brute_force(7))

if (False): #q752
    def brute_force(x):
        a,b,n = 1,1,1
        while (n < x*x):
            if (a % x == 1 and b % x == 0): return n
            n += 1
            a,b = a+7*b, a+b
        return 0
    def g2(x):
        a,b,n,l = 1,1,1,x*x
        while (n < l):
            if (a % x == 1 and b % x == 0): return n
            n += 1
            a,b = a+7*b, a+b
            a %= x
            b %= x
        return 0
    def g3(x):
        if (not sp.isprime(x)): return g2(x)
        if (x == 7): return 7
        def get(n):
            A = gen.Matrix([[1,1],[1,0]])
            B = gen.Matrix([[2,1],[1,0]])
            C = gen.Matrix([[2,1],[6,0]])
            C = C.__pow__(n,x)
            A *= C
            B *= C
            return A.arr[0][1],B.arr[1][1]
        L=x*x-1
        for i in range(4*x,0,-1):
            num = L//i
            #if (i == 1): num = L-1
            a,b = get(num)
            if (a % x == 1 and b % x == 0): return num
        return 0

    def test2(x):
        if (x % 5 != 0 and x % 7 != 0 and (not sp.isprime(x))): return 0
        a,b,n = 1,1,1
        while (n < x*x):
            if (a % x == 1 and b % x == 0): return n
            n += 1
            a,b = a+7*b, a+b
        return 0
    def brute_sum(n):
        sm = 0
        for i in range(2,n+1):
            sm += brute_force(i)
        return sm
    def sum2(n):
        primes = PR.load_primes_until(n)
        pv = []
        for i,p in enumerate(primes):
            print(f"\rprime {i}\\{len(primes)}   ",end="")
            v = g3(p)
            if (v != 0): pv.append((p,v))
        print("get primes")
        def get_values(prev):
            res = []
            for p1 in pv:
                for p2 in prev:
                    inx = p1[0]*p2[0]
                    if (inx > n): break
                    res.append((inx,math.gcd(p1[0],p2[0])*p1[1]*p2[1]//math.gcd(p1[1],p2[1])))
            res = list(set(res))
            res.sort(key=lambda x: x[0])
            return res
        
        total = 0
        last = pv
        while (len(last) > 0):
            total += sum([p[1] for p in last])
            last = get_values(last)
        return total
    
    print(sum([g2(i) for i in range(100)]))

    print(sum(PR.get_factors_primes(70)) == 12)
    for i in range(2,100):
        if (brute_force(i) != g3(i)): print(f"dif at {i}") 
        print(f"{i} -> {brute_force(i)} | {g3(i)} | {i*i/brute_force(i) if brute_force(i) and False != 0 else ''}")
    print(sum2(10**3))
    #print([i for i in range(2,200) if brute_force(i) != 0])
    #print(sum([g3(i) for i in range(1000)]))

if (False): #q634
    def sub_find(n,bc):
        #print(f"got {n/bc}")
        return max(int(math.sqrt(n/bc)-1),0)
    def brute(n,bc):
        i = 2
        while (i*i*bc <= n): i += 1
        return i - 2
    def get_mutches(n):
        v1,v2,i1,i2 = [],[],2,2
        while (i1*i1 < n): v1.append(i1*i1); i1 += 1
        while(i2*i2*i2 < n): v2.append(i2*i2*i2); i2 += 1
        return len(set(v1).intersection(set(v2)))
        
    
    M,sm,i = 3*(10**6),0,2
    print(get_mutches(M))
    while(True):
        c = i*i*i
        if (c >= M): break
        #sm += sub_find(M,c)
        #sm += brute(M,c)
        print(f"for b={i} b^3={c} there are {sub_find(M,c)} below {M}")
        i += 1
    print(sm)

if (False): #q216
    def brute_force(n):
        sm = 0
        for i in range(2,n):
            print(f"\r{i}\\{n}   ({round(100*i/n,2)}%)",end="")
            num = 2*i*i-1
            #print(f"{i} -> {num} | \t{sp.isprime(num)}")
            if (sp.isprime(num)): sm += 1
        return sm
    
    print(brute_force(5*10**7))

if (False): #q282
    def A(m,n):
        if (m == 0): return n+1
        if (n == 0): return A(m-1,1)
        return A(m-1,A(m,n-1))
    def cont(a,b,n,m=None):
        if (n == 1): return pow(a,b,m)
        val = a
        for i in range(b-1):
            val = cont(a,val,n-1,m)
            if (not m is None): val %= m
        return val

    print(cont(3,3,2,5))

if (True): #q401
    def s1(n,a=1):
        return sum([i**a for i in range(1,n+1) if n % i == 0])
    
    def SIG1(n):
        return sum(s1(i,2) for i in range(n+1))
    
    def SIG2_s(n, i):
        return (n//i)*(i*i)
    def SIG2(n,l = None):
        if (l is None): l = n
        return sum([SIG2_s(n,i) for i in range(1,l+1)])
    
    def SIG3(n):
        def ss(k):
            return (k-1)*k*(2*k-1)//6
        x1 = int(math.sqrt(n))*5
        sm = SIG2(n,x1)
        #x1,sm = 1,0
        while(x1 <= n):
            print(f"\r{x1}\\{n} ({round(100*x1/n,2)}%)   ",end="")
            dv = n//x1
            x2 = math.ceil((n+1)/dv)
            #print(f"dv = {dv}, i = {i}, x1 = {x1}, x2 = {x2}")
            #input("a")
            sm += (ss(x2)-ss(x1))*dv
            x1 = x2
        
        return sm
            

    
    N = 10
    for i in range(1,10):  print(f"s({i}) = {s1(i,2)}")
    #for i in range(1,10): print(f"SIG({i}) = {SIG2(i)} | SIG2({i}) = {SIG3(i)}")
    print(f"\n{SIG3(10**12)}")    
    #print(SIG2_s(10**15,int(math.sqrt(10**15))*5))

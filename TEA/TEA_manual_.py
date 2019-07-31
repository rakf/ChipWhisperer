import numpy as np
traces   = np.load(r'C:\Users\user\chipwhisperer\projects\TEA_fixed_key_random_plaintext_data\traces\2019.08.01-02.23.59_traces.npy')
pt       = np.load(r'C:\Users\user\chipwhisperer\projects\TEA_fixed_key_random_plaintext_data\traces\2019.08.01-02.23.59_textin.npy')
knownkey = np.load(r'C:\Users\user\chipwhisperer\projects\TEA_fixed_key_random_plaintext_data\traces\2019.08.01-02.23.59_keylist.npy')
#print knownkey
numtraces = np.shape(traces)[0]-1
numpoint = np.shape(traces)[1]

# Use less than the maximum traces by setting numtraces to something smaller
numtraces = 1000
bestguess = 0
   
def intermediate_1(v, k):
    v1 = ( v[4] << 24
         | v[5] << 16
         | v[6] << 8
         | v[7])
    k0 = ( k[0] << 24
         | k[1] << 16
         | k[2] << 8
         | k[3])
    return ((v1 << 4) + k0) & 0xffffffff

def intermediate_2(v, k):
    v1 = ( v[4] << 24
         | v[5] << 16
         | v[6] << 8
         | v[7])
    k1 = ( k[4] << 24
         | k[5] << 16
         | k[6] << 8
         | k[7])
    return ((v1 >> 5) + k1) & 0xffffffff
    
def intermediate_3(v, k):
    v0 = ( v[0] << 24
         | v[1] << 16
         | v[2] << 8
         | v[3])
    v1 = ( v[4] << 24
         | v[5] << 16
         | v[6] << 8
         | v[7])
    k0 = ( k[0] << 24
         | k[1] << 16
         | k[2] << 8
         | k[3])
    k1 = ( k[4] << 24
         | k[5] << 16
         | k[6] << 8
         | k[7])
    int1 = ((v1 << 4) + k0)  & 0xffffffff
    int2 = (v1 + 0x9e3779b9) & 0xffffffff
    int3 = ((v1 >> 5) + k1)  & 0xffffffff
    return v0 + (int1 ^ int2 ^ int3)
HW = [bin(n).count("1") for n in range(0,256)]
cpaoutput = [0]*256
maxcpa = [0]*256
for kguess in range(0, 256):
    #print "Key = %02x: "%(kguess),
    
    # Build our key guess 
    k = [0, 0, 0, 0, kguess, 0, 0, 0]

    #Initialize arrays and variables to zero
    sumnum = np.zeros(numpoint)
    sumden1 = np.zeros(numpoint)
    sumden2 = np.zeros(numpoint)

    # Calculate hypothesis = Hamming weight of intermediate values
    hyp = np.zeros(numtraces)
    for tnum in range(0, numtraces):
        int = (intermediate_1(pt[tnum], k) >> 0) & 0xff
        hyp[tnum] = HW[int]

    #Means of hypothesis and traces
    meanh = np.mean(hyp, dtype=np.float64)
    meant = np.mean(traces, axis=0, dtype=np.float64)

    #For each trace, calculate correlation
    for tnum in range(0, numtraces):
        hdiff = (hyp[tnum] - meanh)
        tdiff = traces[tnum,:] - meant

        sumnum = sumnum + (hdiff*tdiff)
        sumden1 = sumden1 + hdiff*hdiff 
        sumden2 = sumden2 + tdiff*tdiff

    cpaoutput[kguess] = sumnum / np.sqrt( sumden1 * sumden2 )
        
    maxcpa[kguess] = max(abs(cpaoutput[kguess]))
    #print maxcpa[kguess]

bestguess = np.argmax(maxcpa)
cparefs = np.argsort(maxcpa)[::-1]

print "Best Key Guess: %02x"%(bestguess)
#print [hex(c) for c in cparefs]

#!/usr/bin/env python

#####################################
# Edited by Yen-Lung Chen 2023.04.06
#####################################



from serial import Serial, EIGHTBITS, PARITY_NONE, STOPBITS_ONE
from sys import argv
from struct import unpack, pack

assert len(argv) == 2
s = Serial(
    port=argv[1],
    baudrate=115200,
    bytesize=EIGHTBITS,
    parity=PARITY_NONE,
    stopbits=STOPBITS_ONE,
    xonxoff=False,
    rtscts=False
)

fp_pat  = open('./random_pattern.bin','rb') 
fp_gold = open('./pattern_ans.txt','r')
fp_out  = open('output_result.bin','wb')
assert fp_pat and fp_gold and fp_out

pat = fp_pat.read()
lines = fp_gold.readlines()
assert len(pat) % 96 == 0       #hap 32 bytes, read 32 bytes, read_BQ 32 bytes
set_num = (len(pat) / 96)
print ("Open files successfully!")
print ("Total pattern num: " + str(set_num))
print ('============')

err_count = 0
for i in range(0, len(pat), 96):
    print ('Pattern ' + str(i/96) +'\n')
    s.write(pat[i:i+32])       # Send in haplotype seq
    s.write(pat[i+32:i+64])    # Send in read seq
    s.write(pat[i+64:i+96])    # Send in read base quality seq
    
    out = s.read(5)           # Receive output result
    fp_out.write(out)
    
    count   = 0
    vals    = unpack("{}B".format(5), out)
    score   = 0
    for val in vals:
        if (count > 1):
            score = (score << 8) | val
        count = count + 1    
    
    print ('*Golden*')
    print ('Score: \t' + str(lines[int(i/96)]))
    
    print ('*Output*')
    print ('Score: \t' + str(score))
    
    if (int(lines[int(i/)]) != score):
        err_count = err_count + 1
        print('error')
    else:
        print('correct!')
    print ('============')

if (err_count == 0):
    print('\nYou have passed all patterns!\n')
elif (err_count == 1):
    print('\nThere is 1 error...\n')
else:
    print('\nThere are ' + str(err_count) + ' errors...\n')

fp_pat.close()
fp_gold.close()
fp_out.close()

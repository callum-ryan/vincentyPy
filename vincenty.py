import numpy as np
import sys

c1 = sys.argv[1].split(',')
c2 = sys.argv[2].split(',')

def to_rad(lat_long):
	for x, y in enumerate(lat_long):
		lat_long[x] = (np.float64(y) * np.pi) / 180
	return lat_long

c1 = to_rad(c1)
c2 = to_rad(c2)

a = 6378.137
f = 1/298.257223563
b = 6356.7523142

L = c1[1] - c2[1]
u1 = np.arctan((1-f)*np.tan(c1[0]))
u2 = np.arctan((1-f)*np.tan(c2[0]))

lmd = L
lmdd = 0
count = 0

while abs(lmd-lmdd) > 1e-12 and count <= 100:
	sinSqsig = (np.cos(u2)*np.sin(lmd))*(np.cos(u2)*np.sin(lmd)) + (np.cos(u1)*np.sin(u2) - np.sin(u1)*np.cos(u2)*np.cos(lmd))*(np.cos(u1)*np.sin(u2) - np.sin(u1)*np.cos(u2)*np.cos(lmd))
	sinsig = np.sqrt(sinSqsig)
	cossig = np.sin(u1)*np.sin(u2) + np.cos(u1)*np.cos(u2)*np.cos(lmd)
	sig = np.arctan2(sinsig,cossig)
	sina = np.cos(u1)*np.cos(u2)*np.sin(lmd) / sinsig
	cosSqa = 1 - sina * sina
	cos2sigM = cossig - 2*np.sin(u1)*np.sin(u2)/cosSqa
	C = f/16*cosSqa*(4+f*(4-3*cosSqa))
	lmdd = lmd
	lmd = L + (1-C) * f * sina * (sig + C*sinsig*(cos2sigM+C*cossig*(-1+2*cos2sigM*cos2sigM)))
	count += 1

uSq = cosSqa * (a*a - b*b) / (b*b)
A = 1 + uSq/16384. * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
B = uSq / 1024. * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
dsig = B * sinsig * (cos2sigM + B/4. * (cossig * (-1 + 2 * cos2sigM*cos2sigM) - B / 6. * cos2sigM * (-3 + 4 * sinsig * sinsig) * (-3 + 4 * cos2sigM * cos2sigM)))

s = b * A * (sig-dsig)
a1 = np.arctan2(np.cos(u2)*np.sin(lmd), np.cos(u1)*np.sin(u2)-np.sin(u1)*np.cos(u2)*np.cos(lmd))
a2 = np.arctan2(np.cos(u1)*np.sin(lmd), -np.sin(u1)*np.cos(u2)+np.cos(u1)*np.sin(u2)*np.cos(lmd))

print('distance        : '+'{0:.2f}'.format(s)+' km')
print('initial bearing : '+'{0:.2f}'.format(np.degrees(a1))+' °')
print('final bearing   : '+'{0:.2f}'.format(np.degrees(a2))+' °')
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

λ = L
λd = 0
count = 0

while abs(λ-λd) > 1e-12 and count <= 100:
	sinSqσ = (np.cos(u2)*np.sin(λ))*(np.cos(u2)*np.sin(λ)) + (np.cos(u1)*np.sin(u2) - np.sin(u1)*np.cos(u2)*np.cos(λ))*(np.cos(u1)*np.sin(u2) - np.sin(u1)*np.cos(u2)*np.cos(λ))
	sinσ = np.sqrt(sinSqσ)
	cosσ = np.sin(u1)*np.sin(u2) + np.cos(u1)*np.cos(u2)*np.cos(λ)
	σ = np.arctan2(sinσ,cosσ)
	sinα = np.cos(u1)*np.cos(u2)*np.sin(λ) / sinσ
	cosSqα = 1 - sinα * sinα
	cos2σM = cosσ - 2*np.sin(u1)*np.sin(u2)/cosSqα
	C = f/16*cosSqα*(4+f*(4-3*cosSqα))
	λd = λ
	λ = L + (1-C) * f * sinα * (σ + C*sinσ*(cos2σM+C*cosσ*(-1+2*cos2σM*cos2σM)))
	count += 1

uSq = cosSqα * (a*a - b*b) / (b*b)
A = 1 + uSq/16384. * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
B = uSq / 1024. * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
Δσ = B * sinσ * (cos2σM + B/4. * (cosσ * (-1 + 2 * cos2σM*cos2σM) - B / 6. * cos2σM * (-3 + 4 * sinσ * sinσ) * (-3 + 4 * cos2σM * cos2σM)))

s = b * A * (σ-Δσ)
a1 = np.arctan2(np.cos(u2)*np.sin(λ), np.cos(u1)*np.sin(u2)-np.sin(u1)*np.cos(u2)*np.cos(λ))
a2 = np.arctan2(np.cos(u1)*np.sin(λ), -np.sin(u1)*np.cos(u2)+np.cos(u1)*np.sin(u2)*np.cos(λ))

print('distance        : '+'{0:.2f}'.format(s)+' km')
print('initial bearing : '+'{0:.2f}'.format(np.degrees(a1))+' °')
print('final bearing   : '+'{0:.2f}'.format(np.degrees(a2))+' °')
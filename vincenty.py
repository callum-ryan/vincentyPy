import numpy as np

c1 = [44.0528021, 8.22575059999997]
c2 = [51.457347, -0.196769]

for x,y in enumerate(c1):
	c1[x] = (y * np.pi) / 180

for x,y in enumerate(c2):
	c2[x] = (y * np.pi) / 180

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
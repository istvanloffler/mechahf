import numpy as np
import sympy as sym

#Adatok
R=0.454;
L=0.6*10**(-4);
k_m=0.0116;
J_r=4.42*10**(-6);
omega_n=356.047;
i_n=3.96;
N=27;
J_t=10**(-5);
B_r=2*10**(-6);
U_be=6;
dt = 0.05


#1 e, feladat
t = sym.Symbol('t')
s = sym.Symbol('s')

P11 = L*(J_r*N**2+J_t)
P12 = R*(J_r*N**2+J_t)+L*B_r*N**2
P13 = R*B_r*N**2+k_m**2*N**2
nevezo = P11*s**2+P12*s+P13
szamlalo = k_m*N
gsz = 1/s #gerjesztés

W_1 = (k_m*N)/(P11*s**2+P12*s+P13)*gsz
print("WUom(s) = ",W_1)

W_1_invL = sym.integrals.inverse_laplace_transform(W_1, s, t)
print("omegat és Ube közötti átmeneti fgv = ",W_1_invL)
#sym.plot(W_1_invL,(t,0,0.2), ylabel='szögsebesség[rad/s]', xlabel='idő[s]')

#1 f, feladat
s12 = sym.solve(nevezo,s)
print("Gyökök (S1 ; s2) = ",s12)

W_1s = szamlalo/((s-s12[0])*(s-s12[1]))
T1 = -1/s12[1]
T2 = -1/s12[0]
print("WUom = ",W_1s)
print("Időállandók (T1 ; T2) = ",[T1,T2])


#2 b, feladat
W_UOm = W_1/gsz
W_Ui = (B_r*N**2+(J_r*N**2+J_t)*s)/nevezo
W_MOm = (-R-L*s)/nevezo
W_Mi = (k_m*N)/nevezo

lim_W_UOm = sym.limit(W_UOm,s,0)
lim_W_Ui = sym.limit(W_Ui,s,0)
lim_W_MOm = sym.limit(W_MOm,s,0)
lim_W_Mi = sym.limit(W_Mi,s,0)

print("limes WUomega = ",lim_W_UOm)
print("limes WUI = ",lim_W_Ui)
print("limes WMomega = ",lim_W_MOm)
print("limes WMI = ",lim_W_Mi)

#2 c, fleadat
W_1s_invLalg = sym.integrals.inverse_laplace_transform(W_1s,s,t)
print("Wuomega = ",W_1s_invLalg)

#2 d, feladat
W_Omt = W_1*U_be
W_Omt = sym.integrals.inverse_laplace_transform(W_Omt,s,t)
#sym.plot(W_Omt,(t,0,0.2), ylabel='szögsebesség[rad/s]', xlabel='idő[s]')


#2 e, feladat
Mmax = abs((i_n-U_be*lim_W_Ui)/lim_W_Mi)
print("Mmax = ",Mmax)


#3 a, feladat
A11 = -R*dt/L-1
A12 = -dt*k_m/L
A21 = dt*k_m/J_r
A22 = 1
B1 = dt/L
B2 = 0
C1 = 0
C2 = 1
D = 0

print("A11a = ",A11)
print("A12a = ",A12)
print("A21a = ",A21)
print("A22a = ",A22)
print("B1a = ",B1)


#3 b, feladat
Ab = sym.Matrix(2,2,[-R/L,-k_m/L,k_m/J_r,0])
I = sym.Matrix(2,2,[1,0,0,1])
M1 = I-Ab*dt
Abd = M1.inv()
Bb = sym.Matrix(2,1,[(1/L),0])
Bbd = Abd*Bb*dt
Cb = sym.Matrix([0,1])
Cbd = Abd*Cb
Db = sym.Matrix([0,0])
Dbd = sym.Matrix([Cb[0]*Bbd[0],Cb[1]*Bbd[1]])

print("Ad Matrix = ",Abd)
print("Bd Matrix = ",Bbd)
print("Cd Matrix = ",Cbd)
print("Dd Matrix = ",Dbd)

#3 c, feladat
n = 5
omega1 = [0]*6
i1 = [0]*6
for k in range(n):
    i1[k+1] = A11*i1[k]+A12*omega1[k] + B1*U_be
    omega1[k+1] = A21*i1[k] + A22*omega1[k]

omega2 = [0]*6
i2 = [0]*6
for k in range(n):
    i2[k+1] = Abd[0,0]*i2[k]+Abd[0,1]*omega2[k] + Bbd[0]*U_be
    omega2[k+1] = Abd[1,0]*i2[k] + Abd[1,1]*omega2[k] + Bbd[1]*U_be

print("Forward Euler egységugrásra adott válasz =", omega1)
print("Backward Euler egységugrásra adott válasz =", omega2)


#4 szorgalmi feladat
P11_2 = L*(J_r*N**2+J_t)
P12_2 = R*(J_r*N**2+J_t)+L*10*B_r*N**2
P13_2 = R*10*B_r*N**2+k_m**2*N**2
nevezo_2 = P11_2*s**2+P12_2*s+P13_2

W_Ui2 = (10*B_r*N**2+(J_r*N**2+J_t)*s)/nevezo_2
W_Mi2 = (k_m*N)/nevezo_2

lim_W_Ui2 = sym.limit(W_Ui2,s,0)
lim_W_Mi2 = sym.limit(W_Mi2,s,0)

Mmax2 = abs((i_n-U_be*lim_W_Ui2)/lim_W_Mi2)
print("Mt B_r*10 estén = ",Mmax2)


P11_3 = L*(J_r*N**2+J_t)
P12_3 = R*(J_r*N**2+J_t)+L*B_r/10*N**2
P13_3 = R*B_r/10*N**2+k_m**2*N**2
nevezo_3 = P11_3*s**2+P12_3*s+P13_3

W_Ui3 = (B_r/10*N**2+(J_r*N**2+J_t)*s)/nevezo_3
W_Mi3 = (k_m*N)/nevezo_3

lim_W_Ui3 = sym.limit(W_Ui3,s,0)
lim_W_Mi3 = sym.limit(W_Mi3,s,0)

Mmax3 = abs((i_n-U_be*lim_W_Ui3)/lim_W_Mi3)
print("Mt B_r*/10 estén = ",Mmax3)

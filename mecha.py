import numpy as np
import sympy as sym

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


#1.5-os feladat

t = sym.Symbol('t')
s = sym.Symbol('s')

P11 = L*(J_r*N**2+J_t)
P12 = R*(J_r*N**2+J_t)+L*B_r*N**2
P13 = R*B_r*N**2+k_m**2*N**2
nevezo = P11*s**2+P12*s+P13
szamlalo = k_m*N
print(nevezo)

gsz = 1/s #gerjesztés

W_1 = (k_m*N)/(P11*s**2+P12*s+P13)*gsz
print(W_1)

W_2 = sym.integrals.inverse_laplace_transform(W_1, s, t)
print(W_2)


#1.6-os feladat

s12 = sym.solve(nevezo,s)
print(s12)

W_1s = szamlalo/((s-s12[0])*(s-s12[1]))
T1 = -1/s12[1]
T2 = -1/s12[0]
print(W_1s)
print([T1,T2])


#2.2-es feladat
W_UOm = W_1/gsz
W_Ui = (B_r*N**2+(J_r*N**2+J_t)*s)/nevezo
W_MOm = (-R-L*s)/nevezo
W_Mi = (k_m*N)/nevezo

lim_W_UOm = sym.limit(W_UOm,s,0)
lim_W_Ui = sym.limit(W_Ui,s,0)
lim_W_MOm = sym.limit(W_MOm,s,0)
lim_W_Mi = sym.limit(W_Mi,s,0)

print(lim_W_UOm)
print(lim_W_Ui)
print(lim_W_MOm)
print(lim_W_Mi)









#W_1 = 783/(2500*s*((7326498479328773*s^2)/37778931862957161709568 +(1691909079830335*s)/1152921504606846976 + 7116132150715785/72057594037927936))
x1=7116132150715785/72057594037927936*2500
x2=(1691909079830335)/1152921504606846976*2500
x3=2500*((7326498479328773)/37778931862957161709568)


#print(783/1000,"/ (",x1/1000," + ",x2/1000,"s + ",x3/1000,"s^2)")
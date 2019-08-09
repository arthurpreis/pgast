from sympy import *
import numpy
from einsteinpy.symbolic import MetricTensor, ChristoffelSymbols, RiemannCurvatureTensor, RicciTensor

def Rm_independent_element(el1, el2):
    #Checa se o termo é independente
    a1, b1, c1, d1 = el1
    a2, b2, c2, d2 = el2
    if (a1 == a2) and (b1 == b2) and (c1 == d2) and (d1 == c2):
        return False
    if (a1 == b2) and (b1 == a2) and (c1 == c2) and (d1 == d2):
        return False
    if (a1 == c2) and (b1 == d2) and (c1 == a2) and (d1 == b2):
        return False
    if (a1 == b2) and (b1 == a2) and (c1 == d2) and (d1 == c2):
        return False
    else:
        return True

def Ric_independent_element(el1, el2):
    #Checa se o termo é independente
    a1, b1 = el1
    a2, b2 = el2
    if (a1 == b2) and (b1 == a2):
        return False
    else:
        return True

def print_tex(T, t='riem'):
    #printa em formato LaTeX
    index, Tel = T
    s = ''
    if t == 'riem':
        (a,b,c,d) = index
        s = 'R_{' + a + b + c + d + '} &= '
    elif t == 'ricc':
        (a,b) = index
        s = 'R_{' + a + b + '} &= '
    elif t == 'eins':
        (a,b) = index
        s = 'G_{' + a + b + '} &= '
    elif t=='scal':
        s = 'R &= '
    print(s + latex(Tel) + '\\\\')

#Define as coordenadas
t, r, theta, phi = symbols('t r theta phi')
nu = Function('nu', real = True)(t,r)
lamb = Function('lambda', real = True)(t,r)

#Define a métrica
metric = [[0 for i in range(4)] for i in range(4)]
metric[0][0] = exp(nu)
metric[1][1] = - exp(lamb)
metric[2][2] = - r**2
metric[3][3] = - r**2 * sin(theta)**2

g = MetricTensor(metric, (t, r, theta, phi))
gup = g.change_config('uu')
#Calcula o tensor de Riemann, tensor Ricci, escalar de Ricci
#e tensor de einstein
rm = RiemannCurvatureTensor.from_metric(g)
rm = rm.change_config('llll') #Deixa covariante

ric = RicciTensor.from_metric(g)

#---------------Escrevendo os termos não nulos independentes de Riemann-----
non_zero = []
for i, el in enumerate(rm.tensor()):
    if el != 0: #Salva só os termos não nulos
        e = simplify(el)
        a = str(int(i/64))
        b = str(int((i%64)/16))
        c = str(int((i%16)/4))
        d = str(int(i%4))
        R_index = (a,b,c,d)
        non_zero.append((R_index,e))

unique = []
for el in non_zero:
    #Salva só os termos independetes
    elindex, Rel = el
    flag = True
    for uniqel in unique:
        uniqindex, Run = uniqel
        if not Rm_independent_element(elindex, uniqindex):
            flag = False
    if flag:
        unique.append(el)

for el in unique:
    print_tex(el, 'riem')


#---------------Escrevendo os termos não nulos independentes de Ricci-----
non_zero = []
for i, el in enumerate(ric.tensor()):
    if el != 0: #Salva só os termos não nulos
        e = simplify(el)
        a = str(int((i%16)/4))
        b = str(int(i%4))
        ric_index = (a,b)
        non_zero.append((ric_index,e))

unique = []
for el in non_zero:
    #Salva só os termos independetes
    elindex, Rel = el
    flag = True
    for uniqel in unique:
        uniqindex, Run = uniqel
        if not Ric_independent_element(elindex, uniqindex):
            flag = False
    if flag:
        unique.append(el)

for el in unique:
    print_tex(el, 'ricc')


#---------------Escrevendo o escalar de Ricci-----
s = 0
for i in range(4):
    for j in range(4):
        s += gup[i*4 + j]*ric[i*4+j]
ricci_scalar = simplify(s)

print_tex((ricci_scalar,ricci_scalar),'scal')

#---------------Escrevendo o tensor de einstein---------
G = [0 for i in range(16)]
uniq = []
for i in range(16):
    G[i] = simplify(ric[i] - (1/2)*g[i]*ricci_scalar)
    if G[i] != 0:
        a = str(int((i%16)/4))
        b = str(int(i%4))
        flag = True
        for u in unique:
            if not Ric_independent_element((a,b), u):
                flag = False
        if flag:
            unique.append((a,b))
            print_tex(((a,b),G[i]),t='eins')



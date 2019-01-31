import numpy as np
import pylab as pl
from scipy.optimize import curve_fit as cf

def fill_between(x, y1, y2=0, ax=None, **kwargs):
    """
    Plot filled region between `y1` and `y2`.
    
    This function works exactly the same as matplotlib's fill_between, except
    that it also plots a proxy artist (specifically, a rectangle of 0 size)
    so that it can be added it appears on a legend.
    """
    ax = ax if ax is not None else pl.gca()
    ax.fill_between(x, y1, y2, **kwargs)
    p = pl.Rectangle((0, 0), 0, 0, **kwargs)
    ax.add_patch(p)
    return p


def f(x, a, b):
	return a/x + b

#define all the catalogues
cat_3C273 = open('3C273/3C273_sensitivity.txt')
cat_3C279 = open('3C279/3C279_sensitivity.txt')
cat_3C84 = open('3C84/3C84_sensitivity.txt')
cat_uranus = open('Uranus/uncorrected/Uranus_sensitivity.txt')
cat_CB68 = open('CB68/CB68_sensitivity.txt')

#the 725.0/571.0 conversion is to apply the POL-2 FCF = 725 and get rid of the SCUBA-2 FCF = 537 (erronusly written in the code as 571)

#import 3C273 data
T_850_3C273 = []
e_T_850_3C273 = []
NEFD_850_3C273 = []
e_NEFD_850_3C273 = []
t_elp_3C273 = []
t_exp_3C273 = []
e_t_exp_3C273 = []
for line in cat_3C273.readlines()[2:]:
    tmp = line.split()
    T_850_3C273.append(float(tmp[0]))
    e_T_850_3C273.append(float(tmp[1]))
    NEFD_850_3C273.append(725.0/571.0 * float(tmp[2]))
    e_NEFD_850_3C273.append(float(tmp[5]))
    t_elp_3C273.append(float(tmp[3]))
    t_exp_3C273.append(float(tmp[4]))
    e_t_exp_3C273.append(float(tmp[6]))

#import 3C279 data
T_850_3C279 = []
e_T_850_3C279 = []
NEFD_850_3C279 = []
e_NEFD_850_3C279 = []
t_elp_3C279 = []
t_exp_3C279 = []
e_t_exp_3C279 = []
for line in cat_3C279.readlines()[2:]:
    tmp = line.split()
    T_850_3C279.append(float(tmp[0]))
    e_T_850_3C279.append(float(tmp[1]))
    NEFD_850_3C279.append(725.0/571.0 * float(tmp[2]))
    e_NEFD_850_3C279.append(float(tmp[5]))
    t_elp_3C279.append(float(tmp[3]))
    t_exp_3C279.append(float(tmp[4]))
    e_t_exp_3C279.append(float(tmp[6]))

#import 3C84 data
T_850_3C84 = []
e_T_850_3C84 = []
NEFD_850_3C84 = []
e_NEFD_850_3C84 = []
t_elp_3C84 = []
t_exp_3C84 = []
e_t_exp_3C84 = []
for line in cat_3C84.readlines()[2:]:
    tmp = line.split()
    T_850_3C84.append(float(tmp[0]))
    e_T_850_3C84.append(float(tmp[1]))
    NEFD_850_3C84.append(725.0/571.0 * float(tmp[2]))
    e_NEFD_850_3C84.append(float(tmp[5]))
    t_elp_3C84.append(float(tmp[3]))
    t_exp_3C84.append(float(tmp[4]))
    e_t_exp_3C84.append(float(tmp[6]))

#import uranus data
T_850_uranus = []
e_T_850_uranus = []
NEFD_850_uranus = []
e_NEFD_850_uranus = []
t_elp_uranus = []
t_exp_uranus = []
e_t_exp_uranus = []
for line in cat_uranus.readlines()[2:]:
    tmp = line.split()
    T_850_uranus.append(float(tmp[0]))
    e_T_850_uranus.append(float(tmp[1]))
    NEFD_850_uranus.append(725.0/571.0 * float(tmp[2]))
    e_NEFD_850_uranus.append(float(tmp[5]))
    t_elp_uranus.append(float(tmp[3]))
    t_exp_uranus.append(float(tmp[4]))
    e_t_exp_uranus.append(float(tmp[6]))

#import CB68 data
T_850_CB68 = []
e_T_850_CB68 = []
NEFD_850_CB68 = []
e_NEFD_850_CB68 = []
t_elp_CB68 = []
t_exp_CB68 = []
e_t_exp_CB68 = []
for line in cat_CB68.readlines()[2:]:
    tmp = line.split()
    T_850_CB68.append(float(tmp[0]))
    e_T_850_CB68.append(float(tmp[1]))
    NEFD_850_CB68.append(725.0/571.0 * float(tmp[2]))
    e_NEFD_850_CB68.append(float(tmp[5]))
    t_elp_CB68.append(float(tmp[3]))
    t_exp_CB68.append(float(tmp[4]))
    e_t_exp_CB68.append(float(tmp[6]))

cat_3C273.close()
cat_3C279.close()
cat_3C84.close()
cat_uranus.close()
cat_CB68.close()

#create "all" arrays
T_850_all = np.concatenate([T_850_3C273,T_850_3C279,T_850_3C84,T_850_uranus,T_850_CB68])
e_T_850_all = np.concatenate([e_T_850_3C273,e_T_850_3C279,e_T_850_3C84,e_T_850_uranus,e_T_850_CB68])
NEFD_850_all = np.concatenate([NEFD_850_3C273,NEFD_850_3C279,NEFD_850_3C84,NEFD_850_uranus,NEFD_850_CB68])
e_NEFD_850_all = np.concatenate([e_NEFD_850_3C273,e_NEFD_850_3C279,e_NEFD_850_3C84,e_NEFD_850_uranus,e_NEFD_850_CB68])
t_elp_all = np.concatenate([t_elp_3C273,t_elp_3C279,t_elp_3C84,t_elp_uranus,t_elp_CB68])
t_exp_all = np.concatenate([t_exp_3C273,t_exp_3C279,t_exp_3C84,t_exp_uranus,t_exp_CB68])
e_t_exp_all = np.concatenate([e_t_exp_3C273,e_t_exp_3C279,e_t_exp_3C84,e_t_exp_uranus,e_t_exp_CB68])

#find relation between exposure time and elapsed time
const = np.mean(t_exp_all/t_elp_all)
e_const = np.std(t_exp_all/t_elp_all)

print const, e_const

#curve fit
params, cov = cf(f, T_850_all, NEFD_850_all)
a,b = params[0],params[1]
e_a,e_b = cov[0,0]**0.5,cov[1,1]**0.5

print a*1000,b*1000,e_a*1000,e_b*1000

#$\frac{str(int(round(a*1000,0))) + '\pm' + str(int(round(e_a*1000,0)))}{T{850}} + (str(int(round(b*1000,0))) + '\pm' + str(int(round(e_b*1000,0)))) [mJy\sqrt{s}]$

#define Transmission space
Tspace = np.linspace(0.4,0.9,500)

#NEFD plots
pl.figure()
pl.xlabel(r'$T_{850}$ (WVM$_{225 GHz}$)')
pl.ylabel(r'NEFD (Jy $\sqrt{s}$)')
pl.title("Variance map method (PI)")

pl.plot(T_850_3C273,NEFD_850_3C273,'bo',label='3C273')
pl.errorbar(T_850_3C273,NEFD_850_3C273,yerr=e_NEFD_850_3C273,xerr=e_T_850_3C273,color='b',ls='.')

pl.plot(T_850_3C279,NEFD_850_3C279,'go',label='3C279')
pl.errorbar(T_850_3C279,NEFD_850_3C279,yerr=e_NEFD_850_3C279,xerr=e_T_850_3C279,color='g',ls='.')

pl.plot(T_850_3C84,NEFD_850_3C84,'mo',label='3C84')
pl.errorbar(T_850_3C84,NEFD_850_3C84,yerr=e_NEFD_850_3C84,xerr=e_T_850_3C84,color='m',ls='.')

pl.plot(T_850_uranus,NEFD_850_uranus,'ro',label='Uranus')
pl.errorbar(T_850_uranus,NEFD_850_uranus,yerr=e_NEFD_850_uranus,xerr=e_T_850_uranus,color='r',ls='.')

pl.plot(T_850_CB68,NEFD_850_CB68,'co',label='CB68')
pl.errorbar(T_850_CB68,NEFD_850_CB68,yerr=e_NEFD_850_CB68,xerr=e_T_850_CB68,color='c',ls='.')

pl.plot(Tspace,f(Tspace,a,b),'k-')#,label="NEFD"+r"=$\frac{226\pm21}{T_{850}} + (-58\pm31)[mJy\sqrt{s}]$")
#fill_between(Tspace,f(Tspace,a+e_a,b),f(Tspace,a-e_a,b),facecolor='black',alpha=0.2)
fill_between(Tspace,f(Tspace,a,b+e_b),f(Tspace,a,b-e_b),facecolor='black',alpha=0.2)

pl.legend(fontsize=10,ncol=1,frameon=True,numpoints=1,scatterpoints=1,loc=1)

pl.savefig('figures/all_NEFD_850_var.png',bbox_inches='tight')
pl.savefig('figures/all_NEFD_850_var.pdf',bbox_inches='tight')
pl.show()
pl.close()

#actual elapsed vs predicted elapsed time plot
pl.figure()
pl.xlabel(r'Actual elapsed time (s)')
pl.ylabel(r'Predicted elapsed time (s)')
pl.title("Variance map method (PI)")

pl.plot(t_elp_3C273,t_exp_3C273/const,'b.',label='3C273')
pl.errorbar(t_elp_3C273,t_exp_3C273/const,yerr=e_t_exp_3C273/const,color='b',ls='.')

pl.plot(t_elp_3C279,t_exp_3C279/const,'g.',label='3C279')
pl.errorbar(t_elp_3C279,t_exp_3C279/const,yerr=e_t_exp_3C279/const,color='g',ls='.')

pl.plot(t_elp_3C84,t_exp_3C84/const,'m.',label='3C84')
pl.errorbar(t_elp_3C84,t_exp_3C84/const,yerr=e_t_exp_3C84/const,color='m',ls='.')

pl.plot(t_elp_uranus,t_exp_uranus/const,'r.',label='Uranus')
pl.errorbar(t_elp_uranus,t_exp_uranus/const,yerr=e_t_exp_uranus/const,color='r',ls='.')

pl.plot(t_elp_CB68,t_exp_CB68/const,'c.',label='CB68')
pl.errorbar(t_elp_CB68,t_exp_CB68/const,yerr=e_t_exp_CB68/const,color='c',ls='.')

pl.plot(range(500,2500),range(500,2500),'k--',label='1:1 ratio')

pl.xlim(400,2600)
pl.ylim(400,2600)

pl.legend(fontsize=12,ncol=1,frameon=True,numpoints=1,scatterpoints=1,loc=2)

pl.savefig('figures/t_elp_ratio_var.png',bbox_inches='tight')
pl.savefig('figures/t_elp_ratio_var.pdf',bbox_inches='tight')
pl.show()
pl.close()


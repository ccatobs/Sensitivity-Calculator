import numpy as np
import matplotlib.pyplot as plt
import sys

# Theta [deg.]
# Phi   [deg.]
#Abs(E   )[dBV/m ]
#Abs(Horiz)[dBV/m ]
# Phase(Horiz)[deg.]
#Abs(Verti)[dBV/m ]
# Phase(Verti)[deg.]
#Ax.Ratio[dB    ]

labels = ['Abs(E)[dBV/m]',
          'Abs(Horiz)[dBV/m]',
          'Phase(Horiz)[deg.]',
          'Abs(Verti)[dBV/m]',
          'Phase(Verti)[deg.]',
          'Ax.Ratio[dB]',
          ]

half_angle = 13.4


def power2(db):
    return 10.0**(db/10.0)


def power(db):
    # return 10.0**(db/10.0)
    return db**2


d = np.genfromtxt('data/tolTEC_staircase_singleHorn_280GHz.txt', skip_header=2)

print(d.shape)
# print d[:,0][:721]
# exit()

n = 721
d = d.reshape(-1, n, 8)

phis = d[:, 0, 1]
# print phis[6]
# print d[0][:,0]
# print phis
# exit()

# ebeam = horn_resp_right(th,0,f)
# hbeam = horn_resp_right(th,np.pi/2.,f)
# # tot = np.sum((ebeam+hbeam)*np.sin(th)*th)
# tot = np.trapz((ebeam+hbeam)*np.sin(th),th)
# th_cutoff = np.where(th<np.radians(half_angle))
# ebeam = horn_resp_right(th[th_cutoff],0,f)
# hbeam = horn_resp_right(th[th_cutoff],np.pi/2.,f)

# # beam = np.sum((ebeam+hbeam)*np.sin(th[th_cutoff])*th[th_cutoff])
# beam = np.trapz((ebeam+hbeam)*np.sin(th[th_cutoff]),th[th_cutoff])#*th[th_cutoff])
# # print beam, tot, beam/tot
# col.append(beam/tot)

# print phis[6]
# exit()
th = np.radians(d[0, :, 0])
n_th = len(th)

# ns = n_th/2
# ne = ns+ns
# print ne
# print ns
# print n_th
# exit()
tot = np.trapz(power2(d[0, :n, 3]) /
               np.max(power2(d[0, :n, 3]))*np.sin(th[:n]), th[:n])

th_cutoff = np.where(np.abs(th[:n]) < np.radians(half_angle))[0]

beam = np.trapz(power2(d[0, :n, 3][th_cutoff])/np.max(power2(d[0, :n, 3]
                [th_cutoff]))*np.sin(th[:n][th_cutoff]), th[:n][th_cutoff])

spill_eff = beam/tot
# print beam

# print power(d[6,:n,3][th_cutoff])/np.max(power(d[6,:n,3][th_cutoff]))

toltec_280 = np.genfromtxt('data/beam_280.txt')

tot = np.trapz(power(toltec_280[:, 1])/np.max(power(toltec_280[:, 1]))
               * np.sin(toltec_280[:, 0]), toltec_280[:, 0])
th_cutoff = np.where(toltec_280[:, 0] < np.radians(half_angle))[0]
beam = np.trapz(power(toltec_280[th_cutoff, 1])/np.max(power(toltec_280[th_cutoff, 1]))
                * np.sin(toltec_280[th_cutoff, 0]), toltec_280[th_cutoff, 0])
spill_eff_tol = beam/tot

plt.plot(np.degrees(th[:n]), power2(d[0, :n, 3])/np.max(power2(d[0, :n, 3])),
         label='Doug phi = 0 beam\nspill_eff = %.2f' % spill_eff, linewidth=2)
plt.plot(np.degrees(toltec_280[:, 0]), power(toltec_280[:, 1])/np.max(power(
    toltec_280[:, 1])), label='Sara beam 280GHz\nspill_eff = %.2f' % spill_eff_tol, linewidth=2)

# plt.plot(th[:n][th_cutoff],power(d[6,:n,3])[th_cutoff]/np.max(power(d[6,:n,3][th_cutoff])))
plt.axvline(x=half_angle, color='k', linewidth=2, label='Lyot stop angle')
fac = float(sys.argv[1])/1.2
ind = np.where(np.abs(half_angle*fac-np.degrees(th[:n])) == np.min(
    np.abs(half_angle*fac-np.degrees(th[:n]))))[0]

x = 10*np.log10((power2(d[0, :n, 3])/np.max(power2(d[0, :n, 3])))[ind])

tot = np.trapz(power2(d[0, :n, 3]) /
               np.max(power2(d[0, :n, 3]))*np.sin(th[:n]), th[:n])

th_cutoff = np.where(np.abs(th[:n]) < np.radians(half_angle*fac))[0]

beam = np.trapz(power2(d[0, :n, 3][th_cutoff])/np.max(power2(d[0, :n, 3]
                [th_cutoff]))*np.sin(th[:n][th_cutoff]), th[:n][th_cutoff])

spill_eff = beam/tot
print("beam,beam_db,spill", (power2(
    d[0, :n, 3])/np.max(power2(d[0, :n, 3])))[ind], '%.2f' % (x[0]), '%.2f' % spill_eff)
# exit()

plt.legend(loc=0)
# plt.xlim(0,90)
plt.xlim(0, 180)
plt.yscale('log')
plt.xlabel('angle [deg]')
plt.ylabel('normalized beam')
# plt.savefig('figs/doug_beam_spill_200625.png')
plt.show()
plt.clf()
exit()


TE = power(toltec_280[:, 1])/np.max(power(toltec_280[:, 1]))
TH = power(toltec_280[:, 2])/np.max(power(toltec_280[:, 2]))
DE = power2(d[0, :n, 3])/np.max(power2(d[0, :n, 3]))
DH = power2(d[2, :n, 3])/np.max(power2(d[2, :n, 3]))


plt.subplot(211)
plt.plot(np.degrees(toltec_280[:, 0]), power(
    toltec_280[:, 1])/np.max(power(toltec_280[:, 1])), label='Sara E^2', linewidth=2)
plt.plot(np.degrees(toltec_280[:, 0]), power(
    toltec_280[:, 2])/np.max(power(toltec_280[:, 2])), label='Sara H^2', linewidth=2)


plt.plot(np.degrees(th[:n]), power2(d[0, :n, 3]) /
         np.max(power2(d[0, :n, 3])), label='Doug E^2', linewidth=2)
plt.plot(np.degrees(th[:n]), power2(d[2, :n, 3]) /
         np.max(power2(d[2, :n, 3])), label='Doug H^2', linewidth=2)

plt.legend(loc=0)
plt.axvline(x=13.4, color='k')
plt.yscale('log')
plt.ylim(1e-1, 1)
plt.xlim(0, 20)
plt.subplot(212)

plt.plot(np.degrees(toltec_280[:, 0]), (TE-TH)/(TE+TH),
         label='Sara (E^2-H^2)/(E^2+H^2)', linewidth=2)
plt.plot(np.degrees(th[:n]), (DE-DH)/(DE+DH),
         label='Doug (E^2-H^2)/(E^2+H^2)', linewidth=2)
plt.legend(loc=0)
plt.axhline(y=0, color='k', alpha=0.8)
plt.axvline(x=13.4, color='k')
# plt.yscale('log')
plt.ylim(-0.25, 0.1)
plt.xlim(0, 20)
plt.xlabel('angle [deg]')
plt.savefig('figs/beam_symmetry_200625.png')

exit()


# plt.plot(np.degrees(th[:n]),(d[6,:n,2]),label='2',linewidth=2)
# plt.plot(np.degrees(th[:n]),(d[6,:n,3]),label='3',linewidth=2)
# plt.plot(np.degrees(th[:n]),(d[6,:n,4]),label='4',linewidth=2)
# plt.plot(np.degrees(th[:n]),(d[6,:n,5]),label='5',linewidth=2)
# plt.plot(np.degrees(th[:n]),(d[6,:n,6]),label='6',linewidth=2)
# plt.plot(np.degrees(th[:n]),(d[6,:n,7]),label='7',linewidth=2)

# plt.legend()
# plt.show()


# exit()

# print th_cutoff
# print beam
# print tot
# print (beam/tot)

# exit()
for j in [6, 12]:
    for i in [2]:
        plt.plot(d[j, :, 0], power(d[j, :, i])/np.max(power(d[j, :, i])),
                 label='%s %i' % (labels[i-2], phis[j]), linewidth=2)

plt.xlim(-120, 120)
# plt.ylim(-5,35)
plt.legend()
plt.show()
